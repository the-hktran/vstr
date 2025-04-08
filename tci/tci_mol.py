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
from vstr.mf.lo import NMBoys
from vstr.utils.perf_utils import TIMER
from pyscf import gto, scf, cc
import h5py

#import tntorch as tn
#import torch

A2B = 1.88973
AU2CM = 219474.63 
AMU2AU = 1822.888486209

def traveling_salesman_brute_force(distance_matrix):
    import itertools
    """
    Solves the Traveling Salesman Problem (TSP) using brute force.

    Args:
        distance_matrix: A 2D list representing the distances between cities.
                       distance_matrix[i][j] is the distance from city i to city j.

    Returns:
        A tuple containing:
            - The shortest path (a list of city indices).
            - The total distance of the shortest path.
            Returns None if the input is invalid (e.g., non-square matrix).
    """

    num_cities = len(distance_matrix)

    # Input validation: Check if it's a square matrix
    if any(len(row) != num_cities for row in distance_matrix):
        return None  # Or raise an exception: raise ValueError("Distance matrix must be square.")


    if num_cities <= 1: # trivial cases
      return [0], 0 if num_cities==1 else float('inf')

    min_distance = float('inf')
    shortest_path = None

    # Try all possible permutations of cities (except starting city which we fix at 0)
    for path_permutation in itertools.permutations(range(1, num_cities)):
        path = [0] + list(path_permutation) + [0]  # Start and end at city 0
        current_distance = 0

        for i in range(num_cities):
            current_distance += distance_matrix[path[i]][path[i + 1]]

        if current_distance < min_distance:
            min_distance = current_distance
            shortest_path = path

    return shortest_path[:-1], min_distance

def fiedler_vector(I):
    D = np.diag(np.sum(I, axis=0))
    L = D - I
    w, v = np.linalg.eigh(L)
    return v[:, 1]

class TCIMolecule(Molecule):
    '''
    Atomic units are used throughout. 
    '''

    def __init__(self, potential_cart, natoms, mass, loc_method = None, order_method = None, **kwargs):

        self.potential_cart = potential_cart
        self.dipole_cart = None
        self.natoms = natoms
        self.Nm = 3 * natoms - 6
        self.LowFrequencyCutoff = None
        self.NonFrzCoords = None
        self.mass = mass
        self.onemode_coeff = []
        self.onemode_eig = []

        self.ngridpts = 12
        self.tt_method = 'tci'
        self.rank = 121
        self.tci_tol = 1e-6
        self.loc_method = loc_method
        self.order_method = order_method
        self.Order = 2
        self.OrderPlus = None
        self.calc_integrals = True
        self.calc_dipole = False
        self.IntsFile = "./ints.h5"
        self.ReadGeom = False
        self.doGeomOpt = True
        self.doShiftPotential = True
        self.ReadTensors = False
        self.isLinear = False

        self.__dict__.update(kwargs)

        if self.isLinear == True:
            self.Nm = 3 * natoms - 5

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
        self.cores = cores

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

    def PlotPotential1D(self, i):
        import matplotlib.pyplot as plt
        gridpts, dvr_coeff = self.get_heg(self.ngridpts)
        Vgrid = np.zeros(self.ngridpts)
        Vtens = np.zeros(self.ngridpts)
        for ni, qi in enumerate(gridpts[i]):
            q = np.zeros(self.nm.nmodes)
            q[i] = qi
            for j in range(self.nm.nmodes):
                if j != i:
                    q[j] = gridpts[j][self.ngridpts // 2]
            Vgrid[ni] = self.nm.potential_nm(q)

            nq = np.full((len(self.cores)), self.ngridpts // 2, dtype=int)
            nq[i] = ni
            G = self.cores[0][:, nq[0], :]
            for k in range(1, len(self.cores)):
                G = G @ self.cores[k][:, nq[k], :]
            Vtens[ni] = G[0, 0]
            print(ni, Vtens[ni])

        plt.plot(gridpts[i], Vgrid, label='exact')
        plt.plot(gridpts[i], Vtens, label='tensor')
        plt.legend()
        plt.savefig('1d_potential_%d.png' % i)
        plt.clf()

    def PlotPotential2D(self, i, j):
        import matplotlib.pyplot as plt
        gridpts, dvr_coeff = self.get_heg(self.ngridpts)
        Vgrid = np.zeros((self.ngridpts, self.ngridpts))
        Vtens = np.zeros((self.ngridpts, self.ngridpts))
        for ni, qi in enumerate(gridpts[i]):
            for nj, qj in enumerate(gridpts[j]):
                q = np.zeros(self.nm.nmodes)
                q[i] = qi
                q[j] = qj
                for k in range(self.nm.nmodes):
                    if k != i and k != j:
                        q[k] = gridpts[k][self.ngridpts // 2]
                Vgrid[ni, nj] = self.nm.potential_nm(q)

                nq = np.full((len(self.cores)), self.ngridpts // 2, dtype=int)
                nq[i] = ni
                nq[j] = nj
                G = self.cores[0][:, nq[0], :]
                for k in range(1, len(self.cores)):
                    G = G @ self.cores[k][:, nq[k], :]
                Vtens[ni, nj] = G[0, 0]

        # plot Vgrid with gridpts
        plt.pcolormesh(gridpts[i], gridpts[j], Vgrid, cmap='RdBu')
        plt.xlabel('q%d' % i)
        plt.ylabel('q%d' % j)
        plt.colorbar()
        plt.savefig('potential_%d_%d.png' % (i, j))
        plt.clf()

        # plot Vtens with gridpts
        plt.pcolormesh(gridpts[i], gridpts[j], Vtens, cmap='RdBu')
        plt.xlabel('q%d' % i)
        plt.ylabel('q%d' % j)
        plt.colorbar()
        plt.savefig('tens_%d_%d.png' % (i, j))
        plt.clf()

        dV = Vgrid - Vtens
        plt.pcolormesh(gridpts[i], gridpts[j], dV, cmap='RdBu')
        plt.xlabel('q%d' % i)
        plt.ylabel('q%d' % j)
        plt.colorbar()
        plt.savefig('diff_%d_%d.png' % (i, j))
        plt.clf()

        print("squared error: ", np.sum(dV**2))

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

            if "cores" in f:
                del f["cores"]
            for i, core in enumerate(self.cores):
                f.create_dataset("cores/%d" % i, data = self.cores[i])
    
    def ReadCoreTensors(self, IntsFile = None):
        if IntsFile is None:
            IntsFile = self.IntsFile

        with h5py.File(IntsFile, "r") as f:
            self.core_tensors = [f["core_tensors/%d" % i][()] for i in range(self.Nm)]
            self.cores = [f["cores/%d" % i][()] for i in range(self.Nm)]

    def kernel(self, x0 = None):
        self.Timer.start(0)
        self.CalcNM(x0 = x0)
        self.Timer.stop(0)

        if self.loc_method is not None:
            mlo = NMBoys(self)
            mol_lo = mlo.kernel()
            self.nm.freqs = mol_lo.Frequencies / constants.AU_TO_INVCM
            self.nm.Frequencies = mol_lo.Frequencies
            self.Frequencies = mol_lo.Frequencies
            self.nm.nm_coeff = mlo.Q_loc

        if self.order_method is not None:
            if self.order_method.upper() == 'DIST':
                self.distance_matrix = mlo.get_dist_matrix()
                path, path_length = traveling_salesman_brute_force(self.distance_matrix)
                self.nm.freqs = self.nm.freqs[path]
                self.nm.nm_coeff = self.nm.nm_coeff[:, :, path]
                self.nm.Frequencies = self.nm.Frequencies[path]
                self.Frequencies = self.Frequencies[path]
                print(path)
                print(self.Frequencies)

            elif self.order_method.upper() == 'COORD_ENTROPY':
                self.distance_matrix = np.zeros((self.nm.freqs.shape[0], self.nm.freqs.shape[0]))
                gridpts, coeff = self.get_heg(self.ngridpts)
                for i in range(self.nm.freqs.shape[0]):
                    for j in range(i+1, self.nm.freqs.shape[0]):
                        vgrid = np.array([[self.nm.potential_2mode(i, j, qi, qj) for qj in gridpts[j]] for qi in gridpts[i]])
                        u, s, vh = np.linalg.svd(vgrid)
                        self.distance_matrix[i, j] = -np.sum(s * np.log(s))
                        self.distance_matrix[j, i] = self.distance_matrix[i, j]
                x = fiedler_vector(self.distance_matrix)
                path = np.argsort(x)
                self.nm.freqs = self.nm.freqs[path]
                self.nm.nm_coeff = self.nm.nm_coeff[:, :, path]
                self.nm.Frequencies = self.nm.Frequencies[path]
                self.Frequencies = self.Frequencies[path]
                print(path)
                print(self.Frequencies)

        if not self.ReadTensors:
            self.CalcTT(tt_method = self.tt_method, rank = self.rank, tci_tol = self.tci_tol)
        else:
            self.ReadCoreTensors()
        print(self)

        self.Timer.report(self.TimerNames)

if __name__ == "__main__":
    from vstr.examples.potentials.h2o_n.h2o_n_pot import calc_h2o_n_pot
    from pyscf import gto, scf

    mass_h = 1836.152697 # in au
    mass_o = 32810.46286
    
    '''
    natoms = 3
    mass = [mass_h, mass_h, mass_o]
    def potential_cart_h2o(coords):
        #if np.array(coords).ndim == 3:
        #    return calc_h2o_pot(coords, len(coords))
        #else:
        #    return calc_h2o_pot([coords], 1)[0]
        return calc_h2o_n_pot(coords)[0]
    x0 = np.array([
        [0, -0.757, 0.587],
        [0,  0.757, 0.587],
        [0,  0    , 0    ]]) 
    '''


    natoms = 6
    mass = [mass_h, mass_h, mass_h, mass_h, mass_o, mass_o]
    def potential_cart_h2o(coords):
        return calc_h2o_n_pot(coords)[0]
    x0 = np.array([
      [-0.50113604,  -0.00000008,  -0.09491969],
      [-1.79430900,   0.00000593,  -0.89798097],
      [1.78424679,   0.75999552,   0.49046049],
      [1.78424142,  -0.76001421,   0.49044679],
      [-1.46030412,  -0.00000130,  -0.00111733],
      [1.46515384,  -0.00000381,   0.001]])

    mol = TCIMolecule(potential_cart_h2o, natoms, mass, ngridpts = 5, loc_method = 'boys', order_method = 'coord_entropy')
    x0 *= A2B

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

    print("8")
    mCI = TCIVHCI(mol, MaxTotalQuanta = 8)
    mCI.kernel(doVHCI = False)
