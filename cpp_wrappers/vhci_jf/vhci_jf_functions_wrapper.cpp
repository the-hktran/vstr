#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
//#include <pybind11/spectra.h>
#include <pybind11/stl.h>
#include "VCI_headers.h"
#include <vector>

PYBIND11_MODULE(vhci_jf_functions, m)
{
    m.doc() = "Module for C++ implementations for VHCI based on JF's code.";
    pybind11::class_<FConst>(m, "FConst")
        .def(pybind11::init<double, std::vector<int>, bool>())
        .def_readwrite("QIndices", &FConst::QIndices)
        .def_readwrite("fc", &FConst::fc)
        .def_readwrite("fcpow", &FConst::fcpow)
        .def_readwrite("Order", &FConst::Order)
        .def_readwrite("QUnique", &FConst::QUnique)
        .def_readwrite("QPowers", &FConst::QPowers);
    pybind11::class_<HOFunc>(m, "HOFunc")
        .def_readwrite("Freq", &HOFunc::Freq)
        .def_readwrite("Quanta", &HOFunc::Quanta);
    pybind11::class_<WaveFunction>(m, "WaveFunction")
        .def(pybind11::init<std::vector<int>, std::vector<double>>())
        .def_readwrite("M", &WaveFunction::M)
        .def_readwrite("Modes", &WaveFunction::Modes);

    m.def("DenseDiagonalizeCPP", DenseDiagonalizeCPP, "Forms and diagonalizes dense vibrational Hamiltonian");
    m.def("SparseDiagonalizeCPP", SparseDiagonalizeCPP, "Forms and diagonalizes sparse virbational Hamiltonian");
    m.def("AddStatesHB", AddStatesHB, "Screens for states above the HB threshold");
    m.def("HeatBath_Sort_FC", HeatBath_Sort_FC, "Sorts the force constants from highest to lowest");
    m.def("DoPT2", DoPT2, "Runs PT2 corrections");
    m.def("DoSPT2", DoSPT2, "Runs SPT2 corrections");
}