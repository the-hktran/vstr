# ***V***ibrational ***STR***ucture

`vstr` is a Python-based package of vibrational structure methods based on advances in electronic structure theory.

## Dependencies

### Python Dependencies
* NumPy
* SciPy
* Cython
* Pybind11
* Eigen
* Matplotlib
* Numdifftools
* PySCF and its dependencies
* h5py
<!--* maxvolpy
* TnTorch-->
* xfacpy

### C++ Dependencies
* Eigen
* Boost
* Spectra

## Features

* Vibrational Self-Consistent Field (**VSCF**)
* Vibrational Heat-Bath Configuration Interaction (**VHCI**) in harmonic oscillator product bases and VSCF modal bases using both Taylor Series PES and n-Mode PES
* Epstein-Nesbet Pertubation Theory of Second Order on top of VHCI solutions done deterministically, stochastically, and semi-stochastically (**VHCI+PT2**)
* Optimized Coordinates for VSCF (**OC-VSCF**)
* Numerical PES derivatives computed at various levels of theory up to sextic through the use of PySCF
* IR spectral intensities and frequencies determined from VHCI using a dipole surface expansion or an n-Mode dipole surface

## Author

Written by Henry K. Tran (hkt2112@columbia.edu) based on code written by Jonathan H. Fetherolf at https://github.com/berkelbach-group/VHCI. 

## Publications

* Semistochastic VHCI+PT2 and VSCF references: H. K. Tran and T. C. Berkelbach, "Vibrational heat-bath configuration interaction with semistochastic perturbation theory using harmonic oscillator or VSCF modals", J. Chem. Phys. 159, 194101 (2023) https://doi.org/10.1063/5.0172702
* VHCI and VHCI+PT2: J. H. Fetherolf and T. C. Berkelbach, "Vibrational heat-bath configuration interaction", J. Chem. Phys. 154, 074104 (2021) https://doi.org/10.1063/5.0035454 
