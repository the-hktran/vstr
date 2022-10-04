#!/bin/bash

g++ -fPIC --shared -std=c++11 -Wall -O3 vhci_functions.cpp vhci_functions_wrapper.cpp \
    -I${CONDA_PREFIX}/include \
    -I${CONDA_PREFIX}/include/python3.10 \
    -I${CONDA_PREFIX}/include/eigen3 \
    -o vhci_functions.so
