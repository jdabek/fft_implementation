#include <iostream>
#include <chrono>
#include "complex_arithmetic.h"
#include "fft_arithmetic.h"

#define ENSURE(x) { if (!(x)) { std::cout << __FUNCTION__ << " FAILED: line " << __LINE__ << std::endl; return 1; } }
#define SUCCESS { std::cout << __FUNCTION__ << " PASSED" << std::endl; return 0; }

class FftImplementationUnitTests
{
public:
    static int testComplexArithmetic()
    {
        ComplexNumber z(3.0, 4.0);
        ENSURE(z.re() == 3 && z.im() == 4);
        ENSURE(z.conj().re() == 3 && z.conj().im() == -4);
        ComplexNumber znan(std::nan("0"), std::nan("0"));
        ENSURE(std::isnan(znan.conj().re()) && std::isnan(znan.conj().im()));
        ENSURE(z.sqr() == 25);
        ENSURE(std::isnan(znan.sqr()));
        ENSURE(z.abs() == 5);
        ENSURE(std::isnan(znan.abs()));
        ENSURE(z.phi() == std::atan2(4.0, 3.0));
        ENSURE(std::isnan(znan.phi()));

        ComplexNumber z2(5.0, 0.5);
        ENSURE((-z2).re() == -5 && (-z2).im() == -0.5);
        ENSURE(std::isnan((-znan).re()) && std::isnan((-znan).im()));
        ENSURE((z + z2).re() == 8 && (z + z2).im() == 4.5);
        ENSURE(std::isnan((z + znan).re()) && std::isnan((z + znan).im()));
        ENSURE(std::isnan((znan + z).re()) && std::isnan((znan + z).im()));
        ENSURE(std::isnan((znan + znan).re()) && std::isnan((znan + znan).im()));
        ENSURE((z - z2).re() == -2 && (z - z2).im() == 3.5);
        ENSURE(std::isnan((z - znan).re()) && std::isnan((z - znan).im()));
        ENSURE(std::isnan((znan - z).re()) && std::isnan((znan - z).im()));
        ENSURE(std::isnan((znan - znan).re()) && std::isnan((znan - znan).im()));
        ENSURE((z * z2).re() == 13 && (z * z2).im() == 21.5);
        ENSURE(std::isnan((z * znan).re()) && std::isnan((z * znan).im()));
        ENSURE(std::isnan((znan * z).re()) && std::isnan((znan * z).im()));
        ENSURE(std::isnan((znan * znan).re()) && std::isnan((znan * znan).im()));
        ENSURE((z / z2).re() == 17.0/25.25 && (z / z2).im() == 18.5/25.25);
        ENSURE(std::isnan((z / znan).re()) && std::isnan((z / znan).im()));
        ENSURE(std::isnan((znan / z).re()) && std::isnan((znan / z).im()));
        ENSURE(std::isnan((znan / znan).re()) && std::isnan((znan / znan).im()));
        ComplexNumber z0(0.0, 0.0);
        ENSURE(std::isnan((z / z0).re()) && std::isnan((z / z0).im()));

        ComplexNumber z3(1.0, 2.0);
        ComplexNumber z4(2.0, 1.0);
        ENSURE((z3^2).re() == -3 && (z4^2).im() == 4);
        ENSURE(std::isnan((znan^2).re()) && std::isnan((znan^2).im()));

        ENSURE(ComplexNumber(3.0, 0.0) == ComplexNumber(3.0, 0.0));
        ENSURE(!(ComplexNumber(3.0, 0.0) == ComplexNumber(0.0, 0.0)));
        ENSURE(!(ComplexNumber(3.0, 0.0) != ComplexNumber(3.0, 0.0)));
        ENSURE(ComplexNumber(3.0, 0.0) != ComplexNumber(0.0, 0.0));
        ENSURE(ComplexNumber(3.0, 0.0) == 3);
        ENSURE(3 == ComplexNumber(3.0, 0.0));
        ENSURE(!(ComplexNumber(3.0, 0.0) == 4));
        ENSURE(!(4 == ComplexNumber(3.0, 0.0)));

        ComplexNumber z5(3.0, 6.0);
        double d = 3.0;

        ENSURE(z5 + d == ComplexNumber(6.0, 6.0));
        ENSURE(d + z5 == ComplexNumber(6.0, 6.0));
        ENSURE(z5 - d == ComplexNumber(0.0, 6.0));
        ENSURE(d - z5 == ComplexNumber(0.0, -6.0));
        ENSURE(z5 * d == ComplexNumber(9.0, 18.0));
        ENSURE(d * z5 == ComplexNumber(9.0, 18.0));
        ENSURE(z5 / d == ComplexNumber(1.0, 2.0));
        ENSURE(d / z5 == ComplexNumber(0.2, -0.4));

        SUCCESS;
    }

    static int testComplexVector()
    {
        ComplexVector vec(3);
        ENSURE(vec.size() == 3u);
        try
        {
            vec[5];
            ENSURE(false);
        }
        catch (...)
        {
            ENSURE(true);
        }

        ENSURE(ComplexVector(3, 1.0) == ComplexNumber(1.0, 0.0));

        ComplexNumber z;
        ENSURE(vec.getElement(0, z));
        ENSURE(std::isnan(z.re()) && std::isnan(z.im()));
        ENSURE(vec.getElement(1, z));
        ENSURE(std::isnan(z.re()) && std::isnan(z.im()));
        ENSURE(vec.getElement(2, z));
        ENSURE(std::isnan(z.re()) && std::isnan(z.im()));
        ENSURE(!vec.getElement(5, z));

        PairVector vecPair;
        vecPair.emplace_back(1.0, 5.0);
        vecPair.emplace_back(2.0, 4.0);
        vec = ComplexVector(vecPair);
        ENSURE(vec[0] == ComplexNumber(1.0, 5.0));
        ENSURE(vec[1] == ComplexNumber(2.0, 4.0));

        PairVector vecPair2;
        vecPair2 = vec.getPairVector();
        ENSURE(vecPair[0] == std::make_pair(1.0, 5.0));
        ENSURE(vecPair[1] == std::make_pair(2.0, 4.0));

        vec = ComplexVector(3, ComplexNumber(1.0, 2.0));
        ENSURE(vec[0] == ComplexNumber(1.0, 2.0));
        ENSURE(vec[1] == ComplexNumber(1.0, 2.0));
        ENSURE(vec[2] == ComplexNumber(1.0, 2.0));

        ENSURE(vec.setElement(0, 1.0, -5.0));
        ENSURE(vec.setElement(1, 2.0, 0.0));
        ENSURE(vec.setElement(2, 3.0, 5.0));
        ENSURE(vec.setElement(2, ComplexNumber(3.0, 5.0)));
        ENSURE(!vec.setElement(3, 4.0, 10.0));

        ENSURE(vec.getElement(0, z));
        ENSURE(z.re() == 1 && z.im() == -5);
        ENSURE(vec.getElement(1, z));
        ENSURE(z.re() == 2 && z.im() == 0);
        ENSURE(vec.getElement(2, z));
        ENSURE(z.re() == 3 && z.im() == 5);

        ComplexVector vec2(3);
        ENSURE(vec2.setElement(0, ComplexNumber(1.0, -5.0)));
        ENSURE(vec2.setElement(1, ComplexNumber(2.0, 0.0)));
        ENSURE(vec2.setElement(2, ComplexNumber(3.0, 5.0)));
        ENSURE(!vec2.setElement(3, ComplexNumber(4.0, 10.0)));

        std::vector<double> onev;
        onev = vec.re();
        ENSURE(onev.size() == 3u);
        ENSURE(onev[0] == 1 && onev[1] == 2 && onev[2] == 3);
        onev = vec.im();
        ENSURE(onev.size() == 3u);
        ENSURE(onev[0] == -5 && onev[1] == 0 && onev[2] == 5);

        ComplexVector vec3 = vec.conj();
        ENSURE(vec3.getElement(0, z));
        ENSURE(z.re() == 1 && z.im() == 5);
        ENSURE(vec3.getElement(1, z));
        ENSURE(z.re() == 2 && z.im() == 0);
        ENSURE(vec3.getElement(2, z));
        ENSURE(z.re() == 3 && z.im() == -5);

        ComplexVector vecnan(3);
        vecnan.setElement(1, 0.0, 2.0);
        ComplexVector vecnanconj = vecnan.conj();
        ENSURE(vecnanconj.getElement(0, z));
        ENSURE(std::isnan(z.re()) && std::isnan(z.im()));
        ENSURE(vecnanconj.getElement(1, z));
        ENSURE(z.re() == 0 && z.im() == -2);
        ENSURE(vecnanconj.getElement(2, z));
        ENSURE(std::isnan(z.re()) && std::isnan(z.im()));

        onev = vec.sqr();
        ENSURE(onev[0] == 26 && onev[1] == 4 && onev[2] == 34);
        
        onev = vecnan.sqr();
        ENSURE(std::isnan(onev[0]) && onev[1] == 4 && std::isnan(onev[2]));
        
        onev = vec.abs();
        ENSURE(onev[0] == std::sqrt(26.0));
        ENSURE(onev[1] == 2);
        ENSURE(onev[2] == std::sqrt(34));

        onev = vecnan.abs();
        ENSURE(std::isnan(onev[0]) && onev[1] == 2 && std::isnan(onev[2]));

        onev = vec.phi();
        ENSURE(std::atan2(-5.0, 1.0));
        ENSURE(M_PI/2.0);
        ENSURE(std::atan2(5.0, 3.0));

        onev = vecnan.phi();
        ENSURE(std::isnan(onev[0]) && onev[1] == M_PI/2.0 && std::isnan(onev[2]));

        ComplexVector vec4 = vec + vec;
        ENSURE((-vec4).getElement(0, z));
        ENSURE(z.re() == -2 && z.im() == 10);
        ENSURE((-vec4).getElement(1, z));
        ENSURE(z.re() == -4 && z.im() == 0);
        ENSURE((-vec4).getElement(2, z));
        ENSURE(z.re() == -6 && z.im() == -10);
        
        ComplexVector vecnanplus = vecnan + vec;
        ENSURE(vecnanplus.getElement(0, z));
        ENSURE(std::isnan((-z).re()) && std::isnan((-z).im()));
        ENSURE(vecnanplus.getElement(1, z));
        ENSURE((-z).re() == -2 && (-z).im() == -2);
        ENSURE(vecnanplus.getElement(2, z));
        ENSURE(std::isnan((-z).re()) && std::isnan((-z).im()));
        
        ComplexVector plusvecnan = vec + vecnan;
        ENSURE(plusvecnan.getElement(0, z));
        ENSURE(std::isnan(z.re()) && std::isnan(z.im()));
        ENSURE(plusvecnan.getElement(1, z));
        ENSURE(z.re() == 2 && z.im() == 2);
        ENSURE(plusvecnan.getElement(2, z));
        ENSURE(std::isnan(z.re()) && std::isnan(z.im()));

        ComplexVector vec5 = vec - vec.conj();
        ENSURE(vec5.getElement(0, z));
        ENSURE(z.re() == 0 && z.im() == -10);
        ENSURE(vec5.getElement(1, z));
        ENSURE(z.re() == 0 && z.im() == 0);
        ENSURE(vec5.getElement(2, z));
        ENSURE(z.re() == 0 && z.im() == 10);
        
        ComplexVector vecnanminus = vecnan - vec;
        ENSURE(vecnanminus.getElement(0, z));
        ENSURE(std::isnan(z.re()) && std::isnan(z.im()));
        ENSURE(vecnanminus.getElement(1, z));
        ENSURE(z.re() == -2 && z.im() == 2);
        ENSURE(vecnanminus.getElement(2, z));
        ENSURE(std::isnan(z.re()) && std::isnan(z.im()));
        
        ComplexVector minusvecnan = vec - vecnan;
        ENSURE(minusvecnan.getElement(0, z));
        ENSURE(std::isnan(z.re()) && std::isnan(z.im()));
        ENSURE(minusvecnan.getElement(1, z));
        ENSURE(z.re() == 2 && z.im() == -2);
        ENSURE(minusvecnan.getElement(2, z));
        ENSURE(std::isnan(z.re()) && std::isnan(z.im()));

        ComplexVector vec6 = vec * vec;
        ENSURE(vec6.getElement(0, z));
        ENSURE(z.re() == -24 && z.im() == -10);
        ENSURE(vec6.getElement(1, z));
        ENSURE(z.re() == 4 && z.im() == 0);
        ENSURE(vec6.getElement(2, z));
        ENSURE(z.re() == -16 && z.im() == 30);
        
        ComplexVector vecnantimes = vecnan * vec;
        ENSURE(vecnantimes.getElement(0, z));
        ENSURE(std::isnan(z.re()) && std::isnan(z.im()));
        ENSURE(vecnantimes.getElement(1, z));
        ENSURE(z.re() == 0 && z.im() == 4);
        ENSURE(vecnantimes.getElement(2, z));
        ENSURE(std::isnan(z.re()) && std::isnan(z.im()));
        
        ComplexVector timesvecnan = vec * vecnan;
        ENSURE(timesvecnan.getElement(0, z));
        ENSURE(std::isnan(z.re()) && std::isnan(z.im()));
        ENSURE(timesvecnan.getElement(1, z));
        ENSURE(z.re() == 0 && z.im() == 4);
        ENSURE(timesvecnan.getElement(2, z));
        ENSURE(std::isnan(z.re()) && std::isnan(z.im()));

        ComplexVector vec7 = vec / vec;
        ENSURE(vec7.getElement(0, z));
        ENSURE(z.re() == 1 && z.im() == 0);
        ENSURE(vec7.getElement(1, z));
        ENSURE(z.re() == 1 && z.im() == 0);
        ENSURE(vec7.getElement(2, z));
        ENSURE(z.re() == 1 && z.im() == 0);
        
        ComplexVector vecnanper = vecnan / vec;
        ENSURE(vecnanper.getElement(0, z));
        ENSURE(std::isnan(z.re()) && std::isnan(z.im()));
        ENSURE(vecnanper.getElement(1, z));
        ENSURE(z.re() == 0 && z.im() == 1);
        ENSURE(vecnanper.getElement(2, z));
        ENSURE(std::isnan(z.re()) && std::isnan(z.im()));
        
        ComplexVector pervecnan = vec / vecnan;
        ENSURE(pervecnan.getElement(0, z));
        ENSURE(std::isnan(z.re()) && std::isnan(z.im()));
        ENSURE(pervecnan.getElement(1, z));
        ENSURE(z.re() == 0 && z.im() == -1);
        ENSURE(pervecnan.getElement(2, z));
        ENSURE(std::isnan(z.re()) && std::isnan(z.im()));

        ComplexVector vec8 = ComplexNumber(2.0, 4.0) * vec;
        ENSURE(vec8.getElement(0, z));
        ENSURE(z.re() == 22 && z.im() == -6);
        ENSURE(vec8.getElement(1, z));
        ENSURE(z.re() == 4 && z.im() == 8);
        ENSURE(vec8.getElement(2, z));
        ENSURE(z.re() == -14 && z.im() == 22);

        vecnantimes = vecnan * ComplexNumber(0.0, 2.0);
        ENSURE(vecnantimes.getElement(0, z));
        ENSURE(std::isnan(z.re()) && std::isnan(z.im()));
        ENSURE(vecnantimes.getElement(1, z));
        ENSURE(z.re() == -4 && z.im() == 0);
        ENSURE(vecnantimes.getElement(2, z));
        ENSURE(std::isnan(z.re()) && std::isnan(z.im()));
        
        timesvecnan = ComplexNumber(0.0, 2.0) * vecnan;
        ENSURE(timesvecnan.getElement(0, z));
        ENSURE(std::isnan(z.re()) && std::isnan(z.im()));
        ENSURE(timesvecnan.getElement(1, z));
        ENSURE(z.re() == -4 && z.im() == 0);
        ENSURE(timesvecnan.getElement(2, z));
        ENSURE(std::isnan(z.re()) && std::isnan(z.im()));

        ComplexVector vec9 = vec / ComplexNumber(0.0, 1.0);
        ENSURE(vec9.getElement(0, z));
        ENSURE(z.re() == -5 && z.im() == -1);
        ENSURE(vec9.getElement(1, z));
        ENSURE(z.re() == 0 && z.im() == -2);
        ENSURE(vec9.getElement(2, z));
        ENSURE(z.re() == 5 && z.im() == -3);

        vec9.setElement(0, 1.0, 0.0);
        vec9.setElement(1, 0.0, 1.0);
        vec9.setElement(2, 3.0, 4.0);
        vec9 = ComplexNumber(0.0, 5.0) / vec9;
        ENSURE(vec9.getElement(0, z));
        ENSURE(z.re() == 0 && z.im() == 5);
        ENSURE(vec9.getElement(1, z));
        ENSURE(z.re() == 5 && z.im() == 0);
        ENSURE(vec9.getElement(2, z));
        ENSURE(z.re() == 4.0/5.0 && z.im() == 3.0/5.0);

        vecnanper = vecnan / ComplexNumber(0.0, 2.0);
        ENSURE(vecnanper.getElement(0, z));
        ENSURE(std::isnan(z.re()) && std::isnan(z.im()));
        ENSURE(vecnanper.getElement(1, z));
        ENSURE(z.re() == 1 && z.im() == 0);
        ENSURE(vecnanper.getElement(2, z));
        ENSURE(std::isnan(z.re()) && std::isnan(z.im()));
        
        pervecnan = ComplexNumber(0.0, 2.0) / vecnan;
        ENSURE(pervecnan.getElement(0, z));
        ENSURE(std::isnan(z.re()) && std::isnan(z.im()));
        ENSURE(pervecnan.getElement(1, z));
        ENSURE(z.re() == 1 && z.im() == 0);
        ENSURE(pervecnan.getElement(2, z));
        ENSURE(std::isnan(z.re()) && std::isnan(z.im()));

        ComplexVector vec10(2);
        vec10.setElement(0, -0.5, 0.0);
        vec10.setElement(1, 0.0, -0.5);

        ComplexVector vec11 = vec10^3;
        ENSURE(vec11.getElement(0, z));
        ENSURE(z.re() == -0.125 && std::abs(z.im()) < 1e-10);
        ENSURE(vec11.getElement(1, z));
        ENSURE(std::abs(z.re()) < 1e-10 && z.im() == 0.125);

        ComplexVector vec12 = vecnan^3;
        ENSURE(vec12.getElement(0, z));
        ENSURE(std::isnan(z.re()) && std::isnan(z.im()));
        ENSURE(vec12.getElement(1, z));
        ENSURE(std::abs(z.re()) < 1e-10 && z.im() == -8);
        ENSURE(vec12.getElement(2, z));
        ENSURE(std::isnan(z.re()) && std::isnan(z.im()));

        ENSURE((vecnan + vec11).size() == 0u);
        ENSURE((vecnan - vec11).size() == 0u);
        ENSURE((vecnan * vec11).size() == 0u);
        ENSURE((vecnan / vec11).size() == 0u);

        ENSURE((vec11 + vecnan).size() == 0u);
        ENSURE((vec11 - vecnan).size() == 0u);
        ENSURE((vec11 * vecnan).size() == 0u);
        ENSURE((vec11 / vecnan).size() == 0u);

        ENSURE((vecnan + vec11).size() == 0u);
        ENSURE((vecnan - vec11).size() == 0u);
        ENSURE((vecnan * vec11).size() == 0u);
        ENSURE((vecnan / vec11).size() == 0u);

        try
        {
            vec11[5];
            ENSURE(false);
        }
        catch (std::invalid_argument e)
        {
            std::cout << e.what() << std::endl;
            ENSURE(true);
        }

        try
        {
            const ComplexVector vec13(3);
            vec13[5];
            ENSURE(false);
        }
        catch (std::invalid_argument e)
        {
            std::cout << e.what() << std::endl;
            ENSURE(true);
        }

        ComplexVector vec13a(2, ComplexNumber(1.0, 0.0));
        ComplexVector vec13b(2, ComplexNumber(1.0, 0.0));

        ENSURE(vec13a == vec13b);
        ENSURE(!(vec13a != vec13b));
        ENSURE(vec13a == ComplexNumber(1.0, 0.0));
        ENSURE(!(vec13a != ComplexNumber(1.0, 0.0)));
        ENSURE(ComplexNumber(1.0, 0.0) == vec13a);
        ENSURE(!(ComplexNumber(1.0, 0.0) != vec13a));
        ENSURE(vec13a == 1);
        ENSURE(!(vec13a != 1));
        ENSURE(1 == vec13a);
        ENSURE(!(1 != vec13a));
        vec13b.setElement(1, 0.0, 0.0);
        ENSURE(!(vec13a == vec13b));
        ENSURE(vec13a != vec13b);
        ENSURE(!(vec13b == ComplexNumber(1.0, 0.0)));
        ENSURE(vec13b != ComplexNumber(1.0, 0.0));
        ENSURE(!(ComplexNumber(1.0, 0.0) == vec13b));
        ENSURE(ComplexNumber(1.0, 0.0) != vec13b);
        ENSURE(!(vec13b == 1));
        ENSURE(vec13b != 1);
        ENSURE(!(1 == vec13b));
        ENSURE(1 != vec13b);

        double d = 3.0;

        ComplexVector vec14(2);
        vec14.setElement(0, 3.0, 6.0);
        vec14.setElement(1, 6.0, 3.0);
        ComplexVector vec15 = vec14;

        vec15.setElement(0, 6.0, 6.0);
        vec15.setElement(1, 9.0, 3.0);
        ENSURE(vec14 + d == vec15);
        ENSURE(d + vec14 == vec15);

        vec15.setElement(0, 0.0, 6.0);
        vec15.setElement(1, 3.0, 3.0);
        ENSURE(vec14 - d == vec15);
        vec15.setElement(0, 0.0, -6.0);
        vec15.setElement(1, -3.0, -3.0);
        ENSURE(d - vec14 == vec15);

        vec15.setElement(0, 9.0, 18.0);
        vec15.setElement(1, 18.0, 9.0);
        ENSURE(vec14 * d == vec15);
        ENSURE(d * vec14 == vec15);

        vec15.setElement(0, 1.0, 2.0);
        vec15.setElement(1, 2.0, 1.0);
        ENSURE(vec14 / d == vec15);
        vec15.setElement(0, 0.2, -0.4);
        vec15.setElement(1, 0.4, -0.2);
        ENSURE(d / vec14 == vec15);

        SUCCESS;
    }

    static int testFftImplementation()
    {
        ComplexNumber z;
   
        ComplexVector data(16);
        int f = 4;
        for (unsigned i = 0; i < data.size(); i++)
        {
            data.setElement(
                    i,
                    std::cos(f * 2.0 * M_PI * i / data.size()),
                    std::sin(f * 2.0 * M_PI * i / data.size()));
        }

        ComputeFft cfft(data);

        ComplexVector fft = cfft.getFftComplexVector();
        ComplexVector ifft = cfft.getIfftComplexVector();

        for (unsigned i = 0; i < data.size(); i++)
        {
            fft.getElement(i, z);
            if (i == f)
            {
                ENSURE((z - ComplexNumber(1.0*data.size(), 0.0)).abs() < 1e-10);
            }
            else
            {
                ENSURE(z.abs() < 1e-10);
            }
            ifft.getElement(i, z);
            if (i == data.size() - f)
            {
                ENSURE((z - ComplexNumber(1.0, 0.0)).abs() < 1e-10);
            }
            else
            {
                ENSURE(z.abs() < 1e-10);
            }
        }

        f = 3;
        data = data * ComplexNumber(0.0, 0.0);
        data.setElement(f, 0.0, 1.0);

        cfft = ComputeFft(data);

        fft = cfft.getFftComplexVector();
        ifft = cfft.getIfftComplexVector();

        for (unsigned i = 0; i < data.size(); i++)
        {
            ComplexNumber z;
            fft.getElement(i, z);
            ENSURE((z - ComplexNumber(
                    std::sin(f * 2.0 * M_PI * i / data.size()),
                    std::cos(f * 2.0 * M_PI * i / data.size()))).abs() < 1e-10);
            ifft.getElement(i, z);
            ENSURE((z*ComplexNumber(data.size(), 0.0) - ComplexNumber(
                    -std::sin(f * 2.0 * M_PI * i / data.size()),
                    std::cos(f * 2.0 * M_PI * i / data.size()))).abs() < 1e-10);
        }

        std::vector<double> realVec(16);
        for (unsigned i = 0; i < realVec.size(); i++)
        {
            realVec[i] = std::sin(f * 2.0 * M_PI * i / realVec.size());
        }

        cfft = ComputeFft(realVec);
        fft = cfft.getFftComplexVector();
        ifft = cfft.getIfftComplexVector();

        for (unsigned i = 0; i < realVec.size(); i++)
        {
            if (i == f)
            {
                fft.getElement(i, z);
                ENSURE((z - ComplexNumber(0.0, -0.5*realVec.size())).abs() < 10e-6);
                ifft.getElement(i, z);
                ENSURE((z - ComplexNumber(0.0, 0.5)).abs() < 10e-6);
            }
            else if (i == realVec.size() - f)
            {
                fft.getElement(i, z);
                ENSURE((z - ComplexNumber(0.0, 0.5*realVec.size())).abs() < 10e-6);
                ifft.getElement(i, z);
                ENSURE((z - ComplexNumber(0.0, -0.5)).abs() < 10e-6);
            }
            else
            {
                fft.getElement(i, z);
                ENSURE(z.abs() < 10e-6);
                ifft.getElement(i, z);
                ENSURE(z.abs() < 10e-6);
            }
        }

        realVec = std::vector<double>(16, 0.0);
        realVec[f] = 1.0;

        cfft = ComputeFft(realVec);
        fft = cfft.getFftComplexVector();
        ifft = cfft.getIfftComplexVector();

        for (unsigned i = 0; i < realVec.size(); i++)
        {
            double c = std::cos(f * 2.0 * M_PI * i / realVec.size());
            double s = std::sin(f * 2.0 * M_PI * i / realVec.size());
            fft.getElement(i, z);
            ENSURE((z - ComplexNumber(c, -s)).abs() < 10e-6);
            ifft.getElement(i, z);
            ENSURE((z*ComplexNumber(realVec.size(), 0.0)
                    - ComplexNumber(c, s)).abs() < 10e-6);
        }

        z = ComplexNumber(1.0, 2.0);
        ComplexVector data2 = ComplexVector(1, z);
        cfft = ComputeFft(data2);
        ENSURE(cfft.getFftComplexVector()[0] == z && cfft.getIfftComplexVector()[0] == z);

        data2 = ComplexVector(0);
        cfft = ComputeFft(data2);
        ENSURE(cfft.getInput().size() == 0);
        ENSURE(cfft.getFft().size() == 0);
        ENSURE(cfft.getIfft().size() == 0);

        try
        {
            std::vector<double> realVec = { 1.0, std::nan("0") };
            cfft = ComputeFft(realVec);
            ENSURE(false);
        }
        catch (std::invalid_argument e)
        {
            std::cout << e.what() << std::endl;
            ENSURE(true);
        }
        
        try
        {
            data.setElement(1,
                    ComplexNumber(std::nan("0"), std::nan("0")));
            cfft = ComputeFft(data);
            ENSURE(false);
        }
        catch (std::invalid_argument e)
        {
            std::cout << e.what() << std::endl;
            ENSURE(true);
        }

        try
        {
            data = ComplexVector(1);
            cfft = ComputeFft(data);
            ENSURE(false);
        }
        catch (std::invalid_argument e)
        {
            std::cout << e.what() << std::endl;
            ENSURE(true);
        }

        for (unsigned k = 0; k < 4; k++)
        {
            data = ComplexVector(127+k);
            for (unsigned i = 0; i < data.size(); i++)
            {
                data.setElement(i, ComplexNumber(0.0, 0.0));
            }
            data.setElement(f, ComplexNumber(0.0, 1.0));

            cfft = ComputeFft(data);

            fft = cfft.getFftComplexVector();
            for (unsigned i = 0; i < data.size(); i++)
            {
                fft.getElement(i, z);
                ENSURE((z - ComplexNumber(
                        std::sin(f * 2.0 * M_PI * i / data.size()),
                        std::cos(f * 2.0 * M_PI * i / data.size()))).abs() < 1e-10);
            }

            data = ComplexVector(127+k);
            for (unsigned i = 0; i < data.size(); i++)
            {
                data.setElement(i, ComplexNumber(
                        std::cos(f * 2.0 * M_PI * i / data.size()),
                        std::sin(f * 2.0 * M_PI * i / data.size())));
            }

            cfft = ComputeFft(data);

            cfft = ComputeFft(cfft.getFft());

            ifft = cfft.getIfftComplexVector();
            for (unsigned i = 0; i < data.size(); i++)
            {
                ifft.getElement(i, z);
                ENSURE((z - ComplexNumber(
                        std::cos(f * 2.0 * M_PI * i / data.size()),
                        std::sin(f * 2.0 * M_PI * i / data.size()))).abs() < 1e-10);
            }
        }

        SUCCESS;
    }

    static int testFftPerformance(unsigned nbits)
    {
        unsigned pow2 = std::pow(2, nbits);

        ComplexVector datapow2(pow2);
        int f = 50;
        for (unsigned i = 0; i < datapow2.size(); i++)
        {
            datapow2.setElement(
                    i,
                    std::cos(f * 2.0 * M_PI * i / datapow2.size()),
                    std::sin(f * 2.0 * M_PI * i / datapow2.size()));
        }

        ComputeFft cfftpow2(datapow2);

        auto tpow2fft1 = std::chrono::high_resolution_clock::now();
        cfftpow2.getFft();
        auto tpow2fft2 = std::chrono::high_resolution_clock::now();

        auto tpow2ifft1 = std::chrono::high_resolution_clock::now();
        cfftpow2.getIfft();
        auto tpow2ifft2 = std::chrono::high_resolution_clock::now();

        ComplexVector datapow2p1(pow2+1);
        for (unsigned i = 0; i < datapow2p1.size(); i++)
        {
            datapow2p1.setElement(
                    i,
                    std::cos(f * 2.0 * M_PI * i / datapow2p1.size()),
                    std::sin(f * 2.0 * M_PI * i / datapow2p1.size()));
        }

        ComputeFft cfftpow2p1(datapow2p1);

        auto tpow2p1fft1 = std::chrono::high_resolution_clock::now();
        cfftpow2p1.getFft();
        auto tpow2p1fft2 = std::chrono::high_resolution_clock::now();

        auto tpow2p1ifft1 = std::chrono::high_resolution_clock::now();
        cfftpow2p1.getIfft();
        auto tpow2p1ifft2 = std::chrono::high_resolution_clock::now();

        auto uspow2fft = std::chrono::duration_cast<std::chrono::microseconds>(
                tpow2fft2 - tpow2fft1);
        auto uspow2ifft = std::chrono::duration_cast<std::chrono::microseconds>(
                tpow2ifft2 - tpow2ifft1);

        auto uspow2p1fft = std::chrono::duration_cast<std::chrono::microseconds>(
                tpow2p1fft2 - tpow2p1fft1);
        auto uspow2p1ifft = std::chrono::duration_cast<std::chrono::microseconds>(
                tpow2p1ifft2 - tpow2p1ifft1);

        std::cout << "Duration fft (microseconds) pow2 (" << pow2 << "): "
                  << uspow2fft.count() << std::endl;
        std::cout << "Duration ifft (microseconds) pow2 (" << pow2 << "): "
                  << uspow2ifft.count() << std::endl;

        std::cout << "Duration fft (microseconds) pow2+1 (" << pow2+1 << "): "
                  << uspow2p1fft.count() << std::endl;
        std::cout << "Duration ifft (microseconds) pow2+1 (" << pow2+1 << "): "
                  << uspow2p1ifft.count() << std::endl;
        
        int scalerThreshold = 20;
        ENSURE(uspow2p1fft.count() < scalerThreshold*uspow2fft.count());
        ENSURE(uspow2p1ifft.count() < scalerThreshold*uspow2ifft.count());

        SUCCESS;
    }

};

int main(void)
{
    int nfail = 0;

    nfail += FftImplementationUnitTests::testComplexArithmetic();
    nfail += FftImplementationUnitTests::testComplexVector();

    nfail += FftImplementationUnitTests::testFftImplementation();
    nfail += FftImplementationUnitTests::testFftPerformance(2);
    nfail += FftImplementationUnitTests::testFftPerformance(4);
    nfail += FftImplementationUnitTests::testFftPerformance(8);
    nfail += FftImplementationUnitTests::testFftPerformance(16);

    if (!nfail)
        std::cout << "All tests PASSED." << std::endl;
    else
        std::cout << "There were " << nfail << " FAILED tests." << std::endl;

    return 0;
}

