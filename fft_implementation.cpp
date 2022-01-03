#include <iostream>
#include <fstream>
#include <limits>
#include <iomanip>
#include "fft_arithmetic.h" 

PairVector readNumericData(const std::string fin)
{
    PairVector vals;

    std::ifstream fins(fin);
    if (!fins.good())
    {
        std::cout << "File does not exist: " << fin << std::endl;
        return vals;
    }

    unsigned nsamp = 0;
    std::string line;
    double a, b;
    for (unsigned nparse = 2; nparse > 0; nparse--)
    {
        if (nsamp > 0)
            break;
        while (std::getline(fins, line))
        {
            std::istringstream iss(line);
            if (nparse == 2)
            {
                if (!(iss >> a >> b))
                    break;
            }
            else
            {
                if (!(iss >> a))
                    break;
                b = 0.0;
            }
            nsamp++;
            vals.emplace_back(a, b);
        }
        fins.clear();
        fins.seekg(0);
    }

    return vals; 
}

int main(int argc, char* argv[])
{
    if (argc != 4)
    {
        std::cout << "Please give 3 arguments: [i]fft input_file.txt output_file.txt" << std::endl;
        std::cout << "[i]fft: fft for foward transformation, ifft for backward" << std::endl;
        std::cout << "input_file.txt : real or complex list of input values" << std::endl;
        std::cout << "output_file.txt : write list of complex output values" << std::endl;
        std::cout << "List type: one (real) or two (complex) floats per line" << std::endl;
        return 0;
    }

    std::string transform = argv[1];
    bool doIfft = false;
    if (transform == "ifft")
        doIfft = true;
    if (!doIfft && transform != "fft")
    {
        std::cout << "Unknown transformation: " << transform << std::endl;
        return 1;
    }
    std::string fin = argv[2];
    std::string fout = argv[3];

    PairVector vals = readNumericData(fin);

    if (vals.empty())
    {
        std::cout << "Could not parse file: " << fin << std::endl;
        return 1;
    }

    unsigned nsamp = vals.size();

    ComputeFft cfft(vals);

    PairVector trans;
    if (!doIfft)
        trans = cfft.getFft();
    else
        trans = cfft.getIfft();

    std::ofstream ofs;
    ofs.open(fout);

    using dbl = std::numeric_limits<double>;
    ofs << std::setprecision(dbl::max_digits10);

    for (DoublePair pair : trans)
        ofs << pair.first << " " << pair.second << std::endl;
    ofs.close();

    return 0;
}

