#include <iostream>
#include "complex_arithmetic.h"
#include "fft_arithmetic.h"

#define ENSURE(x) { if (!(x)) { std::cout << __FUNCTION__ << " FAILED: line " << __LINE__ << std::endl; return; } }
#define SUCCESS { std::cout << __FUNCTION__ << " PASSED" << std::endl; }

class FftImplementationUnitTests
{
public:
    static void testComplexArithmetic()
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

        SUCCESS;
    }

    static void testComplexVector()
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

        ComplexNumber z;
        ENSURE(vec.getElement(0, z));
        ENSURE(std::isnan(z.re()) && std::isnan(z.im()));
        ENSURE(vec.getElement(1, z));
        ENSURE(std::isnan(z.re()) && std::isnan(z.im()));
        ENSURE(vec.getElement(2, z));
        ENSURE(std::isnan(z.re()) && std::isnan(z.im()));
        ENSURE(!vec.getElement(5, z));

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

        SUCCESS;
    }

    static void testFft2PowImplementation()
    {
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

        ComplexVector fft = cfft.getFft();
        ComplexVector ifft = cfft.getIfft();

        for (unsigned i = 0; i < data.size(); i++)
        {
            ComplexNumber z;
            fft.getElement(i, z);
            if (i == f)
            {
                ENSURE((z - ComplexNumber(16.0, 0.0)).abs() < 1e-10);
            }
            else
            {
                ENSURE(z.abs() < 1e-10);
            }
            ifft.getElement(i, z);
            if (i == 16-f)
            {
                ENSURE((z - ComplexNumber(16.0, 0.0)).abs() < 1e-10);
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

        fft = cfft.getFft();
        ifft = cfft.getIfft();

        for (unsigned i = 0; i < data.size(); i++)
        {
            ComplexNumber z;
            fft.getElement(i, z);
            ENSURE((z - ComplexNumber(
                    std::sin(f * 2.0 * M_PI * i / data.size()),
                    std::cos(f * 2.0 * M_PI * i / data.size()))).abs() < 1e-10);
            ifft.getElement(i, z);
            ENSURE((z - ComplexNumber(
                    -std::sin(f * 2.0 * M_PI * i / data.size()),
                    std::cos(f * 2.0 * M_PI * i / data.size()))).abs() < 1e-10);
        }

        data = ComplexVector(3);
        data.setElement(0, 1.0, 0.0);
        data.setElement(1, 2.0, 1.0);
        data.setElement(2, 3.0, 2.0);
        
        cfft = ComputeFft(data);

        ENSURE(cfft.getFft().size() == 0);
        ENSURE(cfft.getIfft().size() == 0);

        SUCCESS;
    }
};

int main(void)
{
    FftImplementationUnitTests::testComplexArithmetic();
    FftImplementationUnitTests::testComplexVector();

    FftImplementationUnitTests::testFft2PowImplementation();

    return 0;
}

