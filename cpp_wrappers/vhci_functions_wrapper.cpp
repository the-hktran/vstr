#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include "vhci_functions.hpp"

PYBIND11_MODULE(vhci_functions, m)
{
    m.doc() = "Module for C++ implementations for VHCI";
    m.def("BasisConnectionsCPP", BasisConnectionsCPP, "Generates the connections between basis sets given the anharmonic potential");
    m.def("HamVCPP", HamVCPP, "Generates the vibrational Hamiltonian, dense.");
}
