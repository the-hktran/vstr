import sys
import numpy as np
import scipy
import numdifftools as nd
import h5py
from vstr.utils import init_funcs, constants
from vstr.ff.force_field import ScaleFC_me
from vstr.cpp_wrappers.vhci_jf.vhci_jf_functions import VCISparseHamNMode
from vstr.nmode.mol import NormalModes
from vstr.spectra.dipole import GetDipole
from vstr.nmode.mol import Molecule, get_qmat_ho, get_tmat_ho
from vstr.utils.perf_utils import TIMER
from pyscf import gto, scf, cc
import h5py

#import tntorch as tn
#import torch

A2B = 1.88973
AU2CM = 219474.63 
AMU2AU = 1822.888486209

class TCIMolecule(Molecule):
    '''
    Atomic units are used throughout. 
    '''

    def __init__(self, potential_cart, natoms, mass, **kwargs):

        self.potential_cart = potential_cart
        self.dipole_cart = None
        self.natoms = natoms
        self.Nm = 3 * natoms - 6
        self.LowFrequencyCutoff = None
        self.mass = mass
        self.onemode_coeff = []
        self.onemode_eig = []

        self.ngridpts = 12
        self.tt_method = 'tci'
        self.rank = 121
        self.tci_tol = 1e-6
        self.Order = 2
        self.OrderPlus = None
        self.calc_integrals = True
        self.calc_dipole = False
        self.IntsFile = "./ints.h5"
        self.ReadGeom = False
        self.doGeomOpt = True
        self.doShiftPotential = True

        self.__dict__.update(kwargs)

        self.Timer = TIMER(3)
        self.TimerNames = ["Opt + NM", "Tensor Train", "ContractTT"]

    def __str__(self):
        mol_str = ""
        mol_str += "_______________________\n"
        mol_str += "|                     |\n"
        mol_str += "|      Molecule       |\n"
        mol_str += "|_____________________|\n\n"
        mol_str += "Number of Atoms       : %d\n" % self.natoms
        mol_str += "Normal Mode Frequences:\n"
        for w in self.Frequencies:
            mol_str += "    %.6f\n" % w
        mol_str += "Number of Grid Points : %d\n" % self.ngridpts
        mol_str += "Dipole                : %s\n" % self.calc_dipole
        mol_str += "Intergral Path        : %s\n" % self.IntsFile
        mol_str += "Minimum Energy        : %.6f cm-1\n" % self.V0
        mol_str += "Rank                  : %d\n" % self.core_tensors[0].shape[1]
        mol_str += "TCI Tolerance         : %.6f\n" % self.tci_tol
        mol_str += "Minimum Geometry (A)  :\n"
        for i in range(self.x0.shape[0]):
            mol_str += "    %.6f    %.6f    %.6f\n" % (self.x0[i, 0], self.x0[i, 1], self.x0[i, 2])
        return mol_str

    def _potential(self, x):
        """
        Calculate the potential with 3N-dimensional argument.
        """
        return self.potential_cart(x.reshape((self.natoms,3)))
    
    def _dipole(self, x):
        return self.dipole_cart(x.reshape((self.natoms, 3)))

    def get_heg(self, ngridpts, optimized=False, ngridpts0=None):
        """
        Calculate the HEG grid points and functions by diagonalizing the Q
        matrix in a basis of harmonic oscillator eigenfunctions.

        Parameters:
        ngridpts (int or list of int): The number of harmonic oscillator
            eigenfunctions (per mode) to use when diagonalizing Q.
        """
        if isinstance(ngridpts, (int, np.int)):
            ngridpts = [ngridpts] * self.nm.nmodes
        gridpts = list()
        coeff = list()
        for i in range(self.nm.nmodes):
            nho = ngridpts[i]
            omega = self.nm.freqs[i]
            qmat = get_qmat_ho(omega, nho)
            q, u = scipy.linalg.eigh(qmat)

            if optimized:
                q0, u0 = q, u 

                tmat = get_tmat_ho(omega, nho)
                vgrid = np.array([self.nm.potential_1mode(i,qi) for qi in q0]) # this should be vectorized
                vmat = np.dot(u0, np.dot(np.diag(vgrid), u0.T))
                e, c = scipy.linalg.eigh(tmat+vmat)
                c = c[:,:ngridpts[i]]
                qmat = np.dot(c.T, np.dot(qmat, c))
                q, u = scipy.linalg.eigh(qmat)

            gridpts.append(q)
            coeff.append(u)

        self.gridpts, self.dvr_coeff = gridpts, coeff
        return self.gridpts, self.dvr_coeff 

    def CalcTT(self, tt_method = 'tci', rank = 121, tci_tol = 1e-6):
        gridpts, dvr_coeff = self.get_heg(self.ngridpts)
        self.Timer.start(1)
        if tt_method.upper() == 'TCI':
            cores = self.nm.do_tci(gridpts, maxit = rank, tol = tci_tol)
            print("Tensor Ranks")
            for core in cores:
                print(core.shape)
        else:
            cores = self.nm.do_tt(gridpts)
        self.Timer.stop(1)
        #self.core_tensors = [np.einsum('irj,nr,mr->ijnm', core.detach().numpy(), dvr_c, dvr_c, optimize=True) for core, dvr_c in zip(cores, dvr_coeff)]
        self.Timer.start(2)
        self.core_tensors = [np.einsum('irj,nr,mr->ijnm', core, dvr_c, dvr_c, optimize=True) for core, dvr_c in zip(cores, dvr_coeff)]
        self.Timer.stop(2)

    def ContractTT(self):
        VTT = self.cores[0]
        for i in range(1, len(self.cores)):
            VTT = np.tensordot(VTT, self.cores[i], axes=(-1, 0))
        return VTT[0, ..., 0]

    def ExplicitVTensor(self):
        gridpts_mesh = np.meshgrid(*self.gridpts, indexing='ij')
        for i in range(len(gridpts_mesh)):
            gridpts_mesh[i] = gridpts_mesh[i].flatten()
        return self.nm.potential_nm_vec(*gridpts_mesh).reshape([self.ngridpts]*self.nm.nmodes)
        #gridpts_mesh_torch = [torch.tensor(gridpts_mesh[i], dtype=torch.float64) for i in range(len(gridpts_mesh))]
        #del gridpts_mesh
        #return self.nm.potential_nm_torch(*gridpts_mesh_torch).reshape([self.ngridpts]*self.nm.nmodes).detach().numpy()

    def SaveCoreTensors(self, IntsFile = None):
        if IntsFile is None:
            IntsFile = self.IntsFile

        with h5py.File(IntsFile, "a") as f:
            if "x0" in f:
                del f["x0"]
            f.create_dataset("x0", data = self.nm.x0)
            if "V0" in f:
                del f["V0"]
            f.create_dataset("V0", data = self.nm.V0)
            if "mu0" in f:
                del f["mu0"]
            f.create_dataset("mu0", data = self.nm.mu0)
            if "nm_coeff" in f:
                del f["nm_coeff"]
            f.create_dataset("nm_coeff", data = self.nm.nm_coeff)
            if "freq" in f:
                del f["freq"]
            f.create_dataset("freq", data = self.Frequencies)

            if "core_tensors" in f:
                del f["core_tensors"]
            for i, core in enumerate(self.core_tensors):
                f.create_dataset("core_tensors/%d" % i, data = self.core_tensors[i])
    
    def ReadCoreTensors(self, IntsFile = None):
        if IntsFile is None:
            IntsFile = self.IntsFile

        with h5py.File(IntsFile, "r") as f:
            self.core_tensors = [f["core_tensors/%d" % i][()] for i in range(len(self.Nm))]

    def kernel(self, x0 = None):
        self.Timer.start(0)
        self.CalcNM(x0 = x0)
        self.Timer.stop(0)
        self.CalcTT(tt_method = self.tt_method, rank = self.rank, tci_tol = self.tci_tol)
        print(self)

        self.Timer.report(self.TimerNames)

if __name__ == "__main__":
    from vstr.examples.potentials.h2o_n.h2o_n_pot import calc_h2o_n_pot
    from pyscf import gto, scf

    mass_h = 1836.152697 # in au
    mass_o = 32810.46286
    natoms = 3
    mass = [mass_h, mass_h, mass_o]
    '''
    natoms = 6
    mass = [mass_h, mass_h, mass_h, mass_h, mass_o, mass_o]
    '''
    '''
    def potential_cart_h2o(coords):
        if np.array(coords).ndim == 3:
            return calc_h2o_pot(coords, len(coords))
        else:
            return calc_h2o_pot([coords], 1)[0]
    '''

    def potential_cart_h2o(coords):
        return calc_h2o_n_pot(coords)[0]

    mol = TCIMolecule(potential_cart_h2o, natoms, mass, ngridpts = 5)

    x0 = np.array([
        [0, -0.757, 0.587],
        [0,  0.757, 0.587],
        [0,  0    , 0    ]]) 
    '''
    x0 = np.array([
      [-0.50113604,  -0.00000008,  -0.09491969],
      [-1.79430900,   0.00000593,  -0.89798097],
      [1.78424679,   0.75999552,   0.49046049],
      [1.78424142,  -0.76001421,   0.49044679],
      [-1.46030412,  -0.00000130,  -0.00111733],
      [1.46515384,  -0.00000381,   0.001]])
    x0 *= A2B
    '''

    mol.kernel(x0=x0)

    from vstr.vhci.vhci import TCIVHCI
    print("3")
    mCI = TCIVHCI(mol, MaxTotalQuanta = 3)
    mCI.kernel(doVHCI = False)
    
    print("4")
    mCI = TCIVHCI(mol, MaxTotalQuanta = 4)
    mCI.kernel(doVHCI = False)

    print("5")
    mCI = TCIVHCI(mol, MaxTotalQuanta = 5)
    mCI.kernel(doVHCI = False)

    print("6")
    mCI = TCIVHCI(mol, MaxTotalQuanta = 6)
    mCI.kernel(doVHCI = False)

    print("7")
    mCI = TCIVHCI(mol, MaxTotalQuanta = 7)
    mCI.kernel(doVHCI = False)


 
