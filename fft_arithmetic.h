#ifndef FFT_ARITHMETIC_H
#define FFT_ARITHMETIC_H

#include "complex_arithmetic.h"

class ComputeFft
{
public:
    ComputeFft(const ComplexVector& vec)
    {
        init(vec, true);
    }

    ComputeFft(const std::vector<double>& realVec)
    {
        ComplexVector vec(realVec.size());

        for (unsigned i = 0; i < realVec.size(); i++)
        {
            ComplexNumber z(realVec[i], 0.0);
            if (std::isnan(realVec[i]))
                throw std::invalid_argument(
                        "ComputeFft: incompatible real data w/NaN");
            vec.setElement(i, z);
        }

        init(vec, false);
    }

private:
    void init(const ComplexVector& vec, bool checkIsNan)
    {
        if (checkIsNan)
        {
            for (unsigned i = 0; i < vec.size(); i++)
            {
                ComplexNumber z;
                vec.getElement(i, z);
                if (std::isnan(z.re()) || std::isnan(z.im()))
                    throw std::invalid_argument(
                            "ComputeFft: incompatible complex data w/NaN");
            }
        }

        _input = vec;
    }

    // Cooley-Tukey FFT: radix-2 decimation-in-time (DIT)
    // The count of samples is 2^nsampBits.
    ComplexVector getPow2Fft(
            const ComplexVector& input,
            unsigned i0,
            unsigned nsampBits,
            unsigned s,
            bool doIfft)
    {
        if (_input.size() == 0)
            return ComplexVector(0);

        double nsampBits0 = std::log2(input.size());
        double di;
        if (modf(nsampBits0, &di) != 0.0)
            return ComplexVector(0);

        unsigned N = std::pow(2, nsampBits);
        ComplexVector vec(N);

        if (nsampBits == 0)
        {
            ComplexNumber z;
            input.getElement(i0, z);
            vec.setElement(0, z);
        
            return vec;
        }

        ComplexVector vec1 = getPow2Fft(input, i0, nsampBits - 1, 2 * s, doIfft);
        ComplexVector vec2 = getPow2Fft(input, i0 + s, nsampBits - 1, 2 * s, doIfft);

        double alpha = -2.0 * M_PI / N;
        if (doIfft)
            alpha = -alpha;
        for (unsigned i = 0; i < N/2; i++)
        {
            ComplexNumber p;
            ComplexNumber q;
            vec1.getElement(i, p);
            vec2.getElement(i, q);
            q = q * ComplexNumber(std::cos(alpha * i), std::sin(alpha * i));
            vec.setElement(i, p + q);
            vec.setElement(i + N/2, p - q);
        }

        if (vec.size() == input.size() && doIfft)
            vec = vec / vec.size();

        return vec;
    }

    // Bluestein's FFT algorithm:
    // Using the convolution theorem and the Cooley-Tukey FFT
    ComplexVector getSelectFft(bool doIfft)
    {
        if (_input.size() == 0)
            return ComplexVector(0);

        unsigned nsamp0 = _input.size();
        double nsampBits0 = std::log2(nsamp0);
        double di;
        if (modf(nsampBits0, &di) == 0.0)
        {
            const unsigned nsampBits = std::log2(_input.size());
            return getPow2Fft(_input, 0, nsampBits0, 1, doIfft);
        }

        unsigned nsamp1 = 2*nsamp0 - 1;
        unsigned nsampBits2 = std::ceil(std::log2(nsamp1));
        unsigned nsamp2 = std::pow(2, nsampBits2);

        ComplexNumber z0(0.0, 0.0);
        ComplexVector a(nsamp2, z0);
        ComplexVector b(nsamp2, z0);

        ComplexNumber z;
        ComplexNumber zz;
        double alpha = -M_PI/nsamp0;
        if (doIfft)
            alpha = -alpha;
        ComplexVector chirp(nsamp0);
        for (unsigned i = 0; i < nsamp2; i++)
        {
            if (i < nsamp0)
            {
                _input.getElement(i, z);
                double x = std::cos(alpha*i*i);
                double y = std::sin(alpha*i*i);
                chirp.setElement(i, x, y);
                a.setElement(i, z * ComplexNumber(x, y));
                b.setElement(i, ComplexNumber(x, -y));
            }
            else
            {
                unsigned j = nsamp2 - i;
                if (j < nsamp0)
                {
                    b.getElement(j, z);
                    b.setElement(i, z);
                }
            }
        }

        ComplexVector fa
                = getPow2Fft(a, 0, nsampBits2, 1, false);
        ComplexVector fb
                = getPow2Fft(b, 0, nsampBits2, 1, false);

        ComplexVector ab
                = getPow2Fft(fa*fb, 0, nsampBits2, 1, true);

        ComplexVector vec(nsamp0);
        for (unsigned i = 0; i < nsamp0; i++)
        {
            ab.getElement(i, z);
            chirp.getElement(i, zz);
            vec.setElement(i, z * zz);
        }

        if (doIfft)
            vec = vec / _input.size();

        return vec;
    }

public:
    ComplexVector getInput() { return _input; }

    ComplexVector getFft()
    {
        if (_fft.size() == 0)
        {
            _fft = getSelectFft(false);
        }

        return _fft;
    }

    ComplexVector getIfft()
    {
        if (_ifft.size() == 0)
        {
            _ifft = getSelectFft(true);
        }

        return _ifft;
    }

private:
    ComplexVector   _input;
    ComplexVector   _fft;
    ComplexVector   _ifft;
};

#endif
