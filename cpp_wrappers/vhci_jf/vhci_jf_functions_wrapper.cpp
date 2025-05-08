#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
//#include <pybind11/spectra.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
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
        .def_readwrite("Modes", &WaveFunction::Modes)
        .def(pybind11::self == pybind11::self)
        .def(pybind11::self != pybind11::self);
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
    m.def("AddStatesHBStoreCoupling", AddStatesHBStoreCoupling, "Add states based on HB but keeps coupling elements too");
    m.def("SpectralFrequencyPrune", SpectralFrequencyPrune, "Prunes basis based on how close w is to Hnn-E0");
    m.def("SpectralFrequencyPruneFromVSCF", SpectralFrequencyPruneFromVSCF, "Prunes basis based on how close w is to Hnn-E0 using VSCF modals");
    m.def("DoSpectralPT2", DoSpectralPT2, "Runs spectral PT2 corrections");
    m.def("VCISparseHamNMode", VCISparseHamNMode, "Generates H using n-Mode potential.");
    m.def("VCISparseHamNModeArray", [](std::vector<WaveFunction> BasisSet1, std::vector<WaveFunction> BasisSet2, std::vector<double> Frequencies, double V0, pybind11::array_t<double> buffer2, pybind11::array_t<double> buffer3, pybind11::array_t<double> buffer4, pybind11::array_t<double> buffer5, pybind11::array_t<double> buffer6, bool DiagonalBlock, int MaxNMode, int MaxQ)
            {
                pybind11::buffer_info info2 = buffer2.request();
                pybind11::buffer_info info3 = buffer3.request();
                pybind11::buffer_info info4 = buffer4.request();
                pybind11::buffer_info info5 = buffer5.request();
                pybind11::buffer_info info6 = buffer6.request();
                return VCISparseHamNModeArray(BasisSet1, BasisSet2, Frequencies, V0, static_cast<double*>(info2.ptr), static_cast<double*>(info3.ptr), static_cast<double*>(info4.ptr), static_cast<double*>(info5.ptr), static_cast<double*>(info6.ptr), DiagonalBlock, MaxNMode, MaxQ);
            });
    m.def("VCISparseHamNModeFromOM", VCISparseHamNModeFromOM, "Generates H using n-Mode potential in one mode eigenbasis.");
    m.def("VCISparseHamNModeFromOMArray", [](std::vector<WaveFunction> BasisSet1, std::vector<WaveFunction> BasisSet2, std::vector<double> Frequencies, double V0, std::vector<Eigen::VectorXd> OneMode_Eig, pybind11::array_t<double> buffer3, pybind11::array_t<double> buffer4, pybind11::array_t<double> buffer5, pybind11::array_t<double> buffer6, bool DiagonalBlock, int MaxNMode, int MaxQ)
            {
                pybind11::buffer_info info3 = buffer3.request();
                pybind11::buffer_info info4 = buffer4.request();
                pybind11::buffer_info info5 = buffer5.request();
                pybind11::buffer_info info6 = buffer6.request();
                return VCISparseHamNModeFromOMArray(BasisSet1, BasisSet2, Frequencies, V0, OneMode_Eig, static_cast<double*>(info3.ptr), static_cast<double*>(info4.ptr), static_cast<double*>(info5.ptr), static_cast<double*>(info6.ptr), DiagonalBlock, MaxNMode, MaxQ);
            });
    m.def("ConnectedStatesCIPSI", ConnectedStatesCIPSI, "Finds all connected configurations given an n-mode potential to a space of configurations.");
    m.def("AddStatesCIPSI", AddStatesCIPSI, "Selects configurations based on the CIPSI criterion.");
    m.def("AddStatesHB2Mode", AddStatesHB2Mode, "Selects configurations based on 2-mode potential sorting.");
    m.def("AddStatesHB2ModeArray", [](std::vector<WaveFunction> BasisSet1, pybind11::array_t<double> buffer2, pybind11::array_t<int> buffer3, Eigen::VectorXd C, double eps, bool ExactSingles, int NModes, int MaxQ)
            {
                pybind11::buffer_info info2 = buffer2.request();
                pybind11::buffer_info info3 = buffer3.request();
                return AddStatesHB2ModeArray(BasisSet1, static_cast<double*>(info2.ptr), static_cast<int*>(info3.ptr), C, eps, ExactSingles, NModes, MaxQ);
            });
    m.def("DoSpectralPT2NMode", DoSpectralPT2NMode, "Runs spectral PT2 corrections for nMode potential");
    m.def("DownFoldPT2NModeFromOMArray", [](Eigen::MatrixXd C, std::vector<WaveFunction> BasisSet, double V0, std::vector<Eigen::VectorXd> OneModeEig, pybind11::array_t<double> buffer2, pybind11::array_t<double> buffer3, pybind11::array_t<double> buffer4, pybind11::array_t<double> buffer5, pybind11::array_t<int> buffer6, int MaxNMode, int MaxQ, double PT2_Eps, double w)
        {
            pybind11::buffer_info info2 = buffer2.request();
            pybind11::buffer_info info3 = buffer3.request();
            pybind11::buffer_info info4 = buffer4.request();
            pybind11::buffer_info info5 = buffer5.request();
            pybind11::buffer_info info6 = buffer6.request();
            return DownFoldPT2NModeFromOMArray(C, BasisSet, V0, OneModeEig, static_cast<double*>(info2.ptr), static_cast<double*>(info3.ptr), static_cast<double*>(info4.ptr), static_cast<double*>(info5.ptr), static_cast<int*>(info6.ptr), MaxNMode, MaxQ, PT2_Eps, w);
        });
    m.def("VCISparseHamDiagonalNModeFromOM", VCISparseHamDiagonalNModeFromOM, "Generates H diagonal elements using n-Mode potential in one mode eigenbasis");
    m.def("VCISparseT", VCISparseT, "Generates kinetic energy in HO basis.");
    //m.def("VCISparseHamTCI", VCISparseHamTCI, "Generates H using TCI potential.");
}
