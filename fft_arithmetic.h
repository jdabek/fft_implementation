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
    void getPow2Fft(
            const PairVector& input,
            PairVector& output,
            const PairVector& cosSin,
            unsigned j0,
            unsigned i0,
            unsigned N,
            unsigned s,
            bool doIfft)
    {
        if (N == 1)
        {
            output[j0] = input[i0];
            return;
        }

        N /= 2;

        getPow2Fft(input, output, cosSin, j0, i0, N, 2 * s, doIfft);
        getPow2Fft(input, output, cosSin, j0 + N, i0 + s, N, 2 * s, doIfft);

        unsigned jump = input.size()/(2*N);
        for (unsigned i = 0; i < N; i++)
        {
            double c = cosSin[i*jump].first;
            double s = cosSin[i*jump].second;
            if (doIfft)
                s = -s;
            DoublePair q;
            q.first = output[j0+N+i].first * c - output[j0+N+i].second * s;
            q.second = output[j0+N+i].first * s + output[j0+N+i].second * c;
            output[j0+N+i].first = output[j0+i].first - q.first;
            output[j0+N+i].second = output[j0+i].second - q.second;
            output[j0+i].first = output[j0+i].first + q.first;
            output[j0+i].second = output[j0+i].second + q.second;
        }
    }

    // Bluestein's FFT algorithm:
    // Using the convolution theorem and the Cooley-Tukey FFT
    void getSelectFft(bool doIfft)
    {
        if (_input.size() == 0)
            return;

        unsigned nsamp0 = _input.size();
        double nsampBits0 = std::log2(nsamp0);
        double di;
        if (modf(nsampBits0, &di) == 0.0)
        {
            PairVector cosSin(nsamp0/2);
            double alpha0 = -2.0*M_PI/nsamp0;
            for (unsigned i = 0; i < nsamp0/2; i++)
            {
                cosSin[i].first = std::cos(i*alpha0);
                cosSin[i].second = std::sin(i*alpha0);
            }

            const unsigned nsampBits = std::log2(_input.size());
            if (!doIfft)
            {
                _fft.resize(nsamp0);
                getPow2Fft(_input, _fft, cosSin, 0, 0, nsamp0, 1, doIfft);
            }
            else
            {
                _ifft.resize(nsamp0);
                getPow2Fft(_input, _ifft, cosSin, 0, 0, nsamp0, 1, doIfft);
                for (unsigned i = 0; i < nsamp0; i++)
                {
                    _ifft[i].first /= nsamp0;
                    _ifft[i].second /= nsamp0;
                }
            }
            return;
        }

        unsigned nsamp1 = 2*nsamp0 - 1;
        unsigned nsampBits2 = std::ceil(std::log2(nsamp1));
        unsigned nsamp2 = std::pow(2, nsampBits2);

        PairVector cosSin(nsamp2/2);
        double alpha2 = -2.0*M_PI/nsamp2;
        for (unsigned i = 0; i < nsamp2/2; i++)
        {
            cosSin[i].first = std::cos(i*alpha2);
            cosSin[i].second = std::sin(i*alpha2);
        }

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

        PairVector aout(nsamp2);
        getPow2Fft(a, aout, cosSin, 0, 0, nsamp2, 1, false);
        a.clear();
        PairVector bout(nsamp2);
        getPow2Fft(b, bout, cosSin, 0, 0, nsamp2, 1, false);
        b.clear();

        PairVector c(nsamp2);
        for (unsigned i = 0; i < nsamp2; i++)
        {
            c[i].first = aout[i].first * bout[i].first - aout[i].second * bout[i].second;
            c[i].second = aout[i].first * bout[i].second + aout[i].second * bout[i].first;
        }
        aout.clear();
        bout.clear();

        PairVector cout(nsamp2);
        getPow2Fft(c, cout, cosSin, 0, 0, nsamp2, 1, true);
        c.clear();

        PairVector pvec(nsamp0);
        for (unsigned i = 0; i < nsamp0; i++)
        {
            pvec[i].first = cout[i].first * chirp[i].first - cout[i].second * chirp[i].second;
            pvec[i].second = cout[i].first * chirp[i].second + cout[i].second * chirp[i].first;
            
            pvec[i].first /= nsamp2;
            pvec[i].second /= nsamp2;
        }

        if (!doIfft)
        {
            std::swap(_fft, pvec);
        }
        else
        {
            std::swap(_ifft, pvec);
            for (unsigned i = 0; i < nsamp0; i++)
            {
                _ifft[i].first /= nsamp0;
                _ifft[i].second /= nsamp0;
            }
        }
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
            getSelectFft(false);
        }

        return _fft;
    }

    ComplexVector getFftComplexVector()
    {
        if (_fft.size() == 0)
        {
            getSelectFft(false);
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
            getSelectFft(true);
        }

        return _ifft;
    }

    ComplexVector getIfftComplexVector()
    {
        if (_ifft.size() == 0)
        {
            getSelectFft(true);
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
