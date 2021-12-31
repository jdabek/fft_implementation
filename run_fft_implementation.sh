#!/bin/bash

# First, compile and run unit tests
g++ fft_implementation_unit_tests.cpp -o fft_implementation_unit_tests.o
./fft_implementation_unit_tests.o

# Then, compile and run the fft implementation
g++ fft_implementation.cpp -o fft_implementation.o
# Note: This script saves the fft output to "fft.txt"
./fft_implementation.o

# Finally, plot the data in "fft.txt"
# Note: This python script saves the fft output as an image "fft.png"
python3 plot_ffts.py

