#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include "vscf_functions.hpp"
#include "vhci_jf/VCI_headers.h"
#include <vector>

PYBIND11_MODULE(vhci_functions, m)
{
    m.doc() = "Module for C++ implementations for VSCF";
    m.def("GetVEffCPP", GetVEffCPP, "Generates the effective potential for each mode.");
}
