#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include "vhci_functions.hpp"
#include <vector>

PYBIND11_MODULE(vhci_functions, m)
{
    m.doc() = "Module for C++ implementations for VHCI";
    pybind11::class_<VDeriv>(m, "VDerivCPP")
        .def(pybind11::init<double, std::vector<int>, bool>())
        .def_readwrite("QIndices", &VDeriv::QIndices)
        .def_readwrite("W", &VDeriv::W)
        .def_readwrite("Order", &VDeriv::Order)
        .def_readwrite("QUnique", &VDeriv::QUnique)
        .def_readwrite("QPowers", &VDeriv::QPowers);
    m.def("FormBasisConnectionsCPP", FormBasisConnectionsCPP, "Generates the connections between basis sets given the anharmonic potential");
    m.def("HamVCPP", HamVCPP, "Generates the vibrational Hamiltonian, dense.");
}
