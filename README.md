# ***V***ibrational ***STR***ucture

`vstr` is a Python-based package of vibrational structure methods based on advances in electronic structure theory.

## Dependencies

### Python Dependencies
* NumPy
* SciPy
* Pybind11
* Eigen
* Matplotlib
* PySCF and its dependencies

### C++ Dependencies
* Eigen
* Boost
* Spectra

## Features

* Vibrational Self-Consistent Field (**VSCF**)
* Vibrational Heat-Bath Configuration Interaction (**VHCI**) in harmonic oscillator product bases and VSCF modal bases
* Epstein-Nesbet Pertubation Theory of Second Order on top of VHCI solutions done deterministically, stochastically, and semi-stochastically (**VHCI+PT2**)
* Optimized Coordinates for VSCF (**OC-VSCF**)
* Numerical PES derivatives computed at various levels of theory up to sextic through the use of PySCF
* IR spectral intensities and frequencies determined from VHCI using a dipole surface expansion

## Author

Written by Henry K. Tran (hkt2112@columbia.edu) based on code written by Jonathan H. Fetherolf at https://github.com/berkelbach-group/VHCI. 

## Publications

* VHCI and VHCI+PT2: Jonathan H. Fetherolf and Timothy C. Berkelbach, "Vibrational heat-bath configuration interaction", J. Chem. Phys. 154, 074104 (2021) https://doi.org/10.1063/5.0035454 
