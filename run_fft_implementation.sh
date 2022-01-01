#!/bin/bash

function run_impl()
{
    # First, compile and run unit tests
    g++ fft_implementation_unit_tests.cpp -o fft_implementation_unit_tests.o
    ./fft_implementation_unit_tests.o

    # Then, compile and run the fft implementation
    g++ fft_implementation.cpp -o fft_implementation.o

    # Generate some data
    echo "1.0 0.0" >data.txt
    echo "0.0 1.0" >>data.txt
    echo "-1.0 0.0" >>data.txt
    echo "0.0 -1.0" >>data.txt
    echo "1.0 0.0" >>data.txt

    echo "Input data:"
    cat data.txt

    # Note: This script now reads from data.txt and outputs either fft.txt or ifft.txt
    ./fft_implementation.o fft data.txt fft.txt
    ./fft_implementation.o ifft data.txt ifft.txt

    echo "Output fft:"
    cat fft.txt
    echo "Output ifft:"
    cat ifft.txt

    # Validate the implementation
    # Note: This python script uses text files to operate on the transformations
    python3 validate_implementation.py 1 10
    python3 validate_implementation.py 120 130
}

# Note: the tests save log.txt
run_impl 2>&1 | tee log.txt

if [[ ! -z $(grep "FAILED" log.txt) ]]; then
    echo "Some tests FAILED"
else
    echo "All tests PASSED"
fi
