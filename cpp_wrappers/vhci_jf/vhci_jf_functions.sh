#!/bin/bash

g++ -fPIC --shared -std=c++11 -Wall -O3 -fopenmp vhci_jf_functions.cpp vhci_jf_functions_wrapper.cpp \
    -I${CONDA_PREFIX}/include \
    -I${CONDA_PREFIX}/include/python3.10 \
    -I${CONDA_PREFIX}/include/eigen3 \
    -I/insomnia001/home/hkt2112/include/spectra-0.6.2/include/ \
    -I/insomnia001/home/hkt2112/include/boost_1_66_0/ \
    -o vhci_jf_functions.so
    #-I/burg/home/hkt2112/include/spectra-0.6.2/include \
    #-I/burg/home/hkt2112/include/libtorch/include \
    #-I/burg/home/hkt2112/include/libtorch/include/torch/csrc/api/include \
