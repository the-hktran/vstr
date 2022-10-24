#!/bin/bash

g++ -fPIC --shared -std=c++11 -Wall -O3 -fopenmp vhci_jf_functions.cpp vhci_jf_functions_wrapper.cpp \
    -I${CONDA_PREFIX}/include \
    -I${CONDA_PREFIX}/include/python3.10 \
    -I${CONDA_PREFIX}/include/eigen3 \
    -I/home/henry/include/spectra-0.6.2/include \
    -o vhci_jf_functions.so
