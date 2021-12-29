#include <iostream>
#include "fft_implementation.h"

#define ENSURE(x) { if (!(x)) { std::cout << __FUNCTION__ << " failed on line " << __LINE__ << std::endl; return; } }
#define SUCCESS { std::cout << __FUNCTION__ << " PASSED" << std::endl; }

class FftImplementationUnitTests
{
public:
    static void testComplexArithmetics()
    {
        ComplexNumber z(3.0, 4.0);
        ENSURE(z.re() == 3 && z.im() == 4);
        ENSURE(z.conj().re() == 3 && z.conj().im() == -4);
        ENSURE(z.sqr() == 25);
        ENSURE(z.abs() == 5);
        ENSURE(z.phi() == std::atan2(4.0, 3.0));

        ComplexNumber z2(5.0, 0.5);
        ENSURE((-z2).re() == -5 && (-z2).im() == -0.5);
        ENSURE((z + z2).re() == 8 && (z + z2).im() == 4.5);
        ENSURE((z - z2).re() == -2 && (z - z2).im() == 3.5);
        ENSURE((z*z2).re() == 13 && (z*z2).im() == 21.5);
        ENSURE((z/z2).re() == 17.0/25.25 && (z/z2).im() == 18.5/25.25);
        
        ComplexNumber z3(1.0, 2.0);
        ComplexNumber z4(2.0, 1.0);
        ENSURE((z3^2).re() == -3 && (z4^2).im() == 4);

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
        ENSURE(vec.getElement(1, z));
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

        onev = vec.sqr();
        ENSURE(onev[0] == 26 && onev[1] == 4 && onev[2] == 34);
        
        onev = vec.abs();
        ENSURE(onev[0] == std::sqrt(26.0));
        ENSURE(onev[1] == 2);
        ENSURE(onev[2] == std::sqrt(34));

        onev = vec.phi();
        ENSURE(std::atan2(-5.0, 1.0));
        ENSURE(M_PI/2.0);
        ENSURE(std::atan2(5.0, 3.0));

        ComplexVector vec4 = vec + vec;
        ENSURE((-vec4).getElement(0, z));
        ENSURE(z.re() == -2 && z.im() == 10);
        ENSURE((-vec4).getElement(1, z));
        ENSURE(z.re() == -4 && z.im() == 0);
        ENSURE((-vec4).getElement(2, z));
        ENSURE(z.re() == -6 && z.im() == -10);
        
        ENSURE(vec4.getElement(0, z));
        ENSURE(z.re() == 2 && z.im() == -10);
        ENSURE(vec4.getElement(1, z));
        ENSURE(z.re() == 4 && z.im() == 0);
        ENSURE(vec4.getElement(2, z));
        ENSURE(z.re() == 6 && z.im() == 10);
        
        ComplexVector vec5 = vec - vec.conj();
        ENSURE(vec5.getElement(0, z));
        ENSURE(z.re() == 0 && z.im() == -10);
        ENSURE(vec5.getElement(1, z));
        ENSURE(z.re() == 0 && z.im() == 0);
        ENSURE(vec5.getElement(2, z));
        ENSURE(z.re() == 0 && z.im() == 10);
        
        ComplexVector vec6 = vec * vec;
        ENSURE(vec6.getElement(0, z));
        ENSURE(z.re() == -24 && z.im() == -10);
        ENSURE(vec6.getElement(1, z));
        ENSURE(z.re() == 4 && z.im() == 0);
        ENSURE(vec6.getElement(2, z));
        ENSURE(z.re() == -16 && z.im() == 30);
        
        ComplexVector vec7 = vec / vec;
        ENSURE(vec7.getElement(0, z));
        ENSURE(z.re() == 1 && z.im() == 0);
        ENSURE(vec7.getElement(1, z));
        ENSURE(z.re() == 1 && z.im() == 0);
        ENSURE(vec7.getElement(2, z));
        ENSURE(z.re() == 1 && z.im() == 0);
        
        ComplexVector vec8 = ComplexNumber(2.0, 4.0) * vec;
        ENSURE(vec8.getElement(0, z));
        ENSURE(z.re() == 22 && z.im() == -6);
        ENSURE(vec8.getElement(1, z));
        ENSURE(z.re() == 4 && z.im() == 8);
        ENSURE(vec8.getElement(2, z));
        ENSURE(z.re() == -14 && z.im() == 22);

        ComplexVector vec9 = vec / ComplexNumber(0.0, 1.0);
        ENSURE(vec9.getElement(0, z));
        ENSURE(z.re() == -5 && z.im() == -1);
        ENSURE(vec9.getElement(1, z));
        ENSURE(z.re() == 0 && z.im() == -2);
        ENSURE(vec9.getElement(2, z));
        ENSURE(z.re() == 5 && z.im() == -3);

        ComplexVector vec10(2);
        vec10.setElement(0, -0.5, 0.0);
        vec10.setElement(1, 0.0, -0.5);

        ComplexVector vec11 = vec10^3;
        ENSURE(vec11.getElement(0, z));
        ENSURE(z.re() == -0.125 && std::abs(z.im()) < 1e-10);
        ENSURE(vec11.getElement(1, z));
        ENSURE(std::abs(z.re()) < 1e-10 && z.im() == 0.125);

        ENSURE((vec+vec11).size() == 0u);
        ENSURE((vec-vec11).size() == 0u);
        ENSURE((vec*vec11).size() == 0u);
        ENSURE((vec/vec11).size() == 0u);

        std::cout << vec;

        SUCCESS;
    }
};

int main(void)
{
    FftImplementationUnitTests::testComplexArithmetics();
    FftImplementationUnitTests::testComplexVector();

    return 0;
}

