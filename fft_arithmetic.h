#ifndef FFT_ARITHMETIC_H
#define FFT_ARITHMETIC_H

#include "complex_arithmetic.h"

class ComputeFft
{
public:
    ComputeFft(const PairVector& pvec)
    {
        init(pvec, true);
    }

    ComputeFft(const std::vector<double>& realVec)
    {
        PairVector pvec(realVec.size());

        for (unsigned i = 0; i < pvec.size(); i++)
        {
            if (std::isnan(realVec[i]))
            {
                throw std::invalid_argument(
                        "ComputeFft: incompatible real data w/NaN");
            }
            pvec[i] = std::make_pair(realVec[i], 0.0);
        }

        init(pvec, false);
    }

    ComputeFft(const ComplexVector& vec)
    {
        PairVector pvec(vec.size());

        for (unsigned i = 0; i < vec.size(); i++)
        {
            ComplexNumber z;
            vec.getElement(i, z);
            if (std::isnan(z.re()) || std::isnan(z.im()))
            {
                throw std::invalid_argument(
                        "ComputeFft: incompatible complex data w/NaN");
            }
            pvec[i] = std::make_pair(z.re(), z.im());
        }

        init(pvec, false);
    }

private:
    void init(const PairVector& pvec, bool checkIsNan)
    {
        if (checkIsNan)
        {
            for (unsigned i = 0; i < pvec.size(); i++)
            {
                if (std::isnan(pvec[i].first) || std::isnan(pvec[i].second))
                {
                    throw std::invalid_argument(
                            "ComputeFft: incompatible complex data w/NaN");
                }
            }
        }

        _input = pvec;
    }

    // Cooley-Tukey FFT: radix-2 decimation-in-time (DIT)
    // The count of samples is 2^nsampBits.
    PairVector getPow2Fft(
            const PairVector& input,
            unsigned i0,
            unsigned N,
            unsigned s,
            bool doIfft)
    {
        if (N == 1)
        {
            return { input[i0] };
        }

        PairVector pvec(N);

        N /= 2;

        PairVector pvec1 = getPow2Fft(input, i0, N, 2 * s, doIfft);
        PairVector pvec2 = getPow2Fft(input, i0 + s, N, 2 * s, doIfft);

        double alpha = -M_PI / N;
        if (doIfft)
            alpha = -alpha;
        for (unsigned i = 0; i < N; i++)
        {
            double arg = alpha * i;
            double c = std::cos(arg);
            double s = std::sin(arg);
            DoublePair q;
            q.first = pvec2[i].first * c - pvec2[i].second * s;
            q.second = pvec2[i].first * s + pvec2[i].second * c;
            pvec[i].first = pvec1[i].first + q.first;
            pvec[i].second = pvec1[i].second + q.second;
            pvec[i + N].first = pvec1[i].first - q.first;
            pvec[i + N].second = pvec1[i].second - q.second;
            if (pvec.size() == input.size() && doIfft)
            {
                pvec[i].first /= input.size();
                pvec[i].second /= input.size();
                pvec[i + N].first /= input.size();
                pvec[i + N].second /= input.size();
            }
        }

        return pvec;
    }

    // Bluestein's FFT algorithm:
    // Using the convolution theorem and the Cooley-Tukey FFT
    PairVector getSelectFft(bool doIfft)
    {
        if (_input.size() == 0)
            return PairVector(0);

        unsigned nsamp0 = _input.size();
        double nsampBits0 = std::log2(nsamp0);
        double di;
        if (modf(nsampBits0, &di) == 0.0)
        {
            const unsigned nsampBits = std::log2(_input.size());
            return getPow2Fft(_input, 0, nsamp0, 1, doIfft);
        }

        unsigned nsamp1 = 2*nsamp0 - 1;
        unsigned nsampBits2 = std::ceil(std::log2(nsamp1));
        unsigned nsamp2 = std::pow(2, nsampBits2);

        PairVector a(nsamp2);
        PairVector b(nsamp2);

        double alpha = -M_PI/nsamp0;
        if (doIfft)
            alpha = -alpha;
        PairVector chirp(nsamp0);
        for (unsigned i = 0; i < nsamp2; i++)
        {
            if (i < nsamp0)
            {
                double arg = alpha*i*i;
                double x = std::cos(arg);
                double y = std::sin(arg);
                chirp[i] = std::make_pair(x, y);
                a[i].first = _input[i].first * x - _input[i].second * y;
                a[i].second = _input[i].first * y + _input[i].second * x;
                b[i].first = x;
                b[i].second = -y;
            }
            else
            {
                a[i].first = 0.0;
                a[i].second = 0.0;
                unsigned j = nsamp2 - i;
                if (j < nsamp0)
                {
                    b[i] = b[j];
                }
                else
                {
                    b[i].first = 0.0;
                    b[i].second = 0.0;
                }
            }
        }

        a = getPow2Fft(a, 0, nsamp2, 1, false);
        b = getPow2Fft(b, 0, nsamp2, 1, false);

        PairVector c(nsamp2);
        for (unsigned i = 0; i < nsamp2; i++)
        {
            c[i].first = a[i].first * b[i].first - a[i].second * b[i].second;
            c[i].second = a[i].first * b[i].second + a[i].second * b[i].first;
        }
        a.clear();
        b.clear();

        c = getPow2Fft(c, 0, nsamp2, 1, true);

        PairVector pvec(nsamp0);
        for (unsigned i = 0; i < nsamp0; i++)
        {
            pvec[i].first = c[i].first * chirp[i].first - c[i].second * chirp[i].second;
            pvec[i].second = c[i].first * chirp[i].second + c[i].second * chirp[i].first;
            if (doIfft)
            {
                pvec[i].first /= _input.size();
                pvec[i].second /= _input.size();
            }
        }

        return pvec;
    }

public:
    PairVector getInput()
    {
        return _input;
    }

    ComplexVector getInputComplexVector()
    {
        if (_input.size() > 0)
        {
            if (_inputComplexVector.size() == 0)
            {
                _inputComplexVector = ComplexVector(_input);
            }
        }

        return _inputComplexVector;
    }

    PairVector getFft()
    {
        if (_fft.size() == 0)
        {
            _fft = getSelectFft(false);
        }

        return _fft;
    }

    ComplexVector getFftComplexVector()
    {
        if (_fft.size() == 0)
        {
            _fft = getSelectFft(false);
        }

        if (_fftComplexVector.size() == 0)
        {
            _fftComplexVector = ComplexVector(_fft);
        }

        return _fftComplexVector;
    }

    PairVector getIfft()
    {
        if (_ifft.size() == 0)
        {
            _ifft = getSelectFft(true);
        }

        return _ifft;
    }

    ComplexVector getIfftComplexVector()
    {
        if (_ifft.size() == 0)
        {
            _ifft = getSelectFft(true);
        }

        if (_ifftComplexVector.size() == 0)
        {
            _ifftComplexVector = ComplexVector(_ifft);
        }

        return _ifftComplexVector;
    }

private:
    PairVector      _input;
    ComplexVector   _inputComplexVector;

    PairVector      _fft;
    ComplexVector   _fftComplexVector;

    PairVector      _ifft;
    ComplexVector   _ifftComplexVector;
};

#endif
