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
    pybind11::class_<ptr_wrapper<double>>(m, "ptr_double");

    //m.def("DenseDiagonalizeCPP", DenseDiagonalizeCPP, "Forms and diagonalizes dense vibrational Hamiltonian");
    //m.def("SparseDiagonalizeCPP", SparseDiagonalizeCPP, "Forms and diagonalizes sparse virbational Hamiltonian");
    m.def("GenerateHamV", GenerateHamV, "Forms dense vibrational Hamiltonian");
    m.def("GenerateHam0V", GenerateHam0V, "Forms dense zeroth order vibrational Hamiltonian");
    m.def("GenerateSparseHamV", GenerateSparseHamV, "Forms sparse vibrational Hamiltonian");
    m.def("GenerateSparseHamVOD", GenerateSparseHamVOD, "Forms sparse vibrational Hamiltonian");
    m.def("GenerateSparseHamAnharmV", GenerateSparseHamAnharmV, "Forms sparse vibrational Anharmonic Hamiltonian");
    m.def("GenerateHamAnharmV", GenerateHamAnharmV, "Forms vibrational Anharmonic Hamiltonian");
    m.def("AddStatesHB", AddStatesHB, "Screens for states above the HB threshold");
    m.def("HeatBath_Sort_FC", HeatBath_Sort_FC, "Sorts the force constants from highest to lowest");
    m.def("DoPT2", DoPT2, "Runs PT2 corrections");
    m.def("DoSPT2", DoSPT2, "Runs SPT2 corrections");
    m.def("GetVEffSLOW1CPP", GetVEffSLOW1CPP, "Generates the effective potential for each mode.");
    m.def("MakeCTensorsCPP", MakeCTensorsCPP, "Generates all C Tensors for each anharmonic term.");
    m.def("MakeCTensorCPP", MakeCTensorCPP, "Generates the C Tensor given the list of modes.");
    m.def("GetVEffSLOW2CPP", GetVEffSLOW2CPP, "Generates the effective potential for each mode using memory efficient implementation.");
    m.def("GetVEffCPP", GetVEffCPP, "Generates the effective potential for each mode using mode independent memory efficient implementation.");
    m.def("ContractedAnharmonicPotential", ContractedAnharmonicPotential, "Creates Y matrices.");
    m.def("ContractedHOTerms", ContractedHOTerms, "Creates X matrices.");
    m.def("VCIHamFromVSCF", VCIHamFromVSCF, "Generates H in the modal basis.");
    m.def("VCISparseHamFromVSCF", VCISparseHamFromVSCF, "Generates H in the modal basis.");
    m.def("AddStatesHBWithMax", AddStatesHBWithMax, "Screens for states above the HB threshold with a maximum on Quanta per mode.");
    m.def("AddStatesHBFromVSCF", AddStatesHBFromVSCF, "Screens for states above the HB with exact matrix elements in VSCF.");
    m.def("DoPT2FromVSCF", DoPT2FromVSCF, "Runs PT2 corrections in the modal basis.");
    m.def("DoSPT2FromVSCF", DoSPT2FromVSCF, "Runs stochastic PT2 corrections in the modal basis.");
    m.def("ProdU", ProdU, "Multiplies a list of 2 x 2 matrices.");
    m.def("SetUij", SetUij, "Makes a rotation metrix.");
    m.def("SetUs", SetUs, "Makes a list of rotation matrices");
    m.def("ContractFCCPP", [](pybind11::array_t<double> buffer3, pybind11::array_t<double> buffer4, pybind11::array_t<double> buffer5, pybind11::array_t<double> buffer6, Eigen::MatrixXd U, int N)
        {
            pybind11::buffer_info info3 = buffer3.request(); 
            pybind11::buffer_info info4 = buffer4.request(); 
            pybind11::buffer_info info5 = buffer5.request(); 
            pybind11::buffer_info info6 = buffer6.request(); 

            return ContractFCCPP(static_cast<double*>(info3.ptr),static_cast<double*>(info4.ptr),static_cast<double*>(info5.ptr),static_cast<double*>(info6.ptr), U, N);
        });
    m.def("SpectralFrequencyPrune", SpectralFrequencyPrune, "Prunes basis based on how close w is to Hnn-E0");
}
