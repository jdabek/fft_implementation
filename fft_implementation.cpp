#include <iostream>
#include <fstream>
#include "fft_arithmetic.h" 

int main(int argc, char* argv[])
{
    unsigned nsamp = 127;
    ComplexVector data(nsamp);
    std::vector<std::pair<double, double>>
            freqPhase = { std::make_pair(5.7, 2.1),
                          std::make_pair(29.2, 1.5) };

    for (unsigned i = 0; i < nsamp; i++)
    {
        double x = 0.0;
        double y = 0.0;
        for (const auto& fp : freqPhase)
        {
            x += std::cos(fp.second + fp.first * 2.0 * M_PI * i / nsamp);
            y += std::sin(fp.second + fp.first * 2.0 * M_PI * i / nsamp);
        }
        data.setElement(i, x, y);
    }

    ComputeFft cfft(data);

    ComplexVector fft = cfft.getFft();

    std::ofstream ofs;
    ofs.open("fft.txt");
    ofs << fft;
    ofs.close();

    return 0;
}

