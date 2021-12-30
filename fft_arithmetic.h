#ifndef FFT_ARITHMETIC_H
#define FFT_ARITHMETIC_H

#include "complex_arithmetic.h"

class ComputeFft
{
public:
    ComputeFft(const ComplexVector& vec)
    {
        if (vec.size() < 2)
            throw std::runtime_error("ComputeFft: incompatible data");

        for (unsigned i = 0; i < vec.size(); i++)
        {
            ComplexNumber z;
            vec.getElement(i, z);
            if (std::isnan(z.re()) || std::isnan(z.im()))
                throw std::runtime_error("ComputeFft: incompatible data w/NaN");
        }

        _nsampBits = 0;
        double nsampBits = std::log2(vec.size());
        double di;
        if (modf(nsampBits, &di) == 0.0)
            _nsampBits = nsampBits;

        _input = vec;
    }

    // The count of samples is 2^nsampBits.
    ComplexVector getPow2Fft(
            unsigned i0,
            unsigned nsampBits,
            unsigned s,
            bool doIfft)
    {
        unsigned N = std::pow(2, nsampBits);
        ComplexVector vec(N);

        if (nsampBits == 0)
        {
            ComplexNumber z;
            _input.getElement(i0, z);
            vec.setElement(0, z);
        
            return vec;
        }

        ComplexVector vec1 = getPow2Fft(i0, nsampBits - 1, 2 * s, doIfft);
        ComplexVector vec2 = getPow2Fft(i0 + s, nsampBits - 1, 2 * s, doIfft);

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

        return vec;
    }

    ComplexVector getInput() { return _input; }

    ComplexVector getFft()
    {
        // TODO: Still has to be pow-2 count samples
        if (_nsampBits == 0)
            return ComplexVector(0);

        if (_fft.size() == 0)
        {
            _fft = getPow2Fft(0, _nsampBits, 1, false);
        }

        return _fft;
    }

    ComplexVector getIfft()
    {
        // TODO: Still has to be pow-2 count samples
        if (_nsampBits == 0)
            return ComplexVector(0);

        if (_ifft.size() == 0)
        {
            _ifft = getPow2Fft(0, _nsampBits, 1, true);
        }

        return _ifft;
    }

private:
    ComplexVector _input;
    unsigned _nsampBits;
    ComplexVector _fft;
    ComplexVector _ifft;
};

#endif
