import sys
import numpy as np
import scipy
import numdifftools as nd
import h5py
from vstr.utils import init_funcs, constants
from vstr.ff.force_field import ScaleFC_me
from vstr.cpp_wrappers.vhci_jf.vhci_jf_functions import VCISparseHamNMode
from vstr.spectra.dipole import GetDipole
from pyscf import gto, scf, cc

import tntorch as tn
import torch

A2B = 1.88973
AU2CM = 219474.63 
AMU2AU = 1822.888486209

def ho_3d(x):
    w = [0.00751569, 0.01746802, 0.01797638]
    D = 0.1
    a = 1
    #w = [1642.62705755, 3825.7245683,  3929.27450625]
    return (0.5 * (w[1]**2 * x[0, 1]**2 + w[2]**2 * x[0, 2]**2) * 1836.15267389
        + D * (1 - np.exp(-a * x[0, 0]))**2)
#            + 0.25 * 100 * x[0, 1] * x[0, 1] * x[0, 2] * x[0, 2]) # + 0.001 * x[0, 0] * x[0, 1] * x[0, 2] + 0.01 * x[0, 0]**2 * x[0, 1]**2 * x[0, 2]**2)

def PyPotential(x, pymol, Method = 'rhf'):
    new_mol = pymol.copy()
    new_mol.unit = 'B'
    atom = []
    for i in range(pymol.natm):
        atom.append([pymol.atom_symbol(i), (x[i, 0], x[i, 1], x[i, 2])])
    new_mol.atom = atom
    new_mol.build()
    mf = scf.RHF(new_mol)
    mf.kernel(verbose = 0, max_cycle = 10000)
    E = mf.e_tot
    if Method == 'ccsd(t)':
        mcc = cc.CCSD(mf).run(verbose = 0)
        de = mcc.ccsd_t()
        E = mcc.e_tot + de
    elif Method == 'ccsd':
        mcc = cc.CCSD(mf).run(verbose = 0)
        E = mcc.e_tot
    return E

def PyDipole(x, pymol, Method = 'rhf', ReturnE = False):
    new_mol = pymol.copy()
    new_mol.unit = 'B'
    atom = []
    for i in range(pymol.natm):
        atom.append([pymol.atom_symbol(i), (x[i, 0], x[i, 1], x[i, 2])])
    new_mol.atom = atom
    new_mol.build()
    mf = scf.RHF(new_mol)
    mf.kernel(verbose = 0, max_cycle = 10000)
    return GetDipole(new_mol, mf, Method = Method, ReturnE = ReturnE)

class Molecule():
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
        self.Order = 2
        self.OrderPlus = None
        self.calc_integrals = True
        self.calc_dipole = False
        self.usePyPotDip = False
        self.IntsFile = "./ints.h5"
        self.ReadInt = False
        self.ReadDip = False
        self.ReadGeom = False
        self.doGeomOpt = True
        self.doShiftPotential = True

        self.__dict__.update(kwargs)

    def __str__(self):
        '''
        print("_______________________")
        print("|                     |")
        print("|      Molecule       |")
        print("|_____________________|")
        print("")
        print("Number of Atoms       :", self.natoms)
        print("Normal Mode Frequences:")
        for w in self.Frequencies:
            print("  ", w)
        print("Order                 :", self.Order)
        print("Number of Grid Points :", self.ngridpts)
        print("Dipole                :", self.calc_dipole)
        print("Intergral Path        :", self.IntsFile)
        print("Minimum Energy        :", self.V0, "cm-1")
        print("Minimum Geometry (A)  :\n", self.x0)
        '''

        mol_str = ""
        mol_str += "_______________________\n"
        mol_str += "|                     |\n"
        mol_str += "|      Molecule       |\n"
        mol_str += "|_____________________|\n\n"
        mol_str += "Number of Atoms       : %d\n" % self.natoms
        mol_str += "Normal Mode Frequences:\n"
        for w in self.Frequencies:
            mol_str += "    %.6f\n" % w
        mol_str += "Order                 : %d\n" % self.Order
        mol_str += "Number of Grid Points : %d\n" % self.ngridpts
        mol_str += "Dipole                : %s\n" % self.calc_dipole
        mol_str += "Intergral Path        : %s\n" % self.IntsFile
        mol_str += "Minimum Energy        : %.6f cm-1\n" % self.V0
        mol_str += "Minimum Geometry (A)  :\n"
        for i in range(self.x0.shape[0]):
            mol_str += "    %.6f    %.6f    %.6f\n" % (self.x0[i, 0], self.x0[i, 1], self.x0[i, 2])
        return mol_str

    def init_pypotential(self, pymol, Method = 'rhf'):
        def pypot(x):
            return PyPotential(x, pymol, Method = Method)
        self.potential_cart = pypot

    def init_pydipole(self, pymol, Method = 'rhf', doPotDip = False):
        self.calc_dipole = True
        if doPotDip:
            self.calc_integrals = False
            self.usePyPotDip = True
        def pydip(x):
            return PyDipole(x, pymol, Method = Method, ReturnE = doPotDip)
        self.dipole_cart = pydip

    def init_dipole(self, dipole_cart, doPotDip = False):
        self.calc_dipole = True
        self.dipole_cart = dipole_cart
        if doPotDip:
            self.calc_integrals = False
            self.usePyPotDip = True

    '''
    def init_pypotdip(self, pymol, Method = 'rhf'):
        self.usePyPotDip = True
        def pydip(x):
            return PyDipole(x, pymol, Method = Method, ReturnE = True)
        self.dipole_cart = pydip
    '''

    def calc_com(self, x):
        com = np.zeros(3)
        for i in range(self.natoms):
            com += self.mass[i] * x[i]
        com /= np.sum(self.mass)
        return com

    def _potential(self, x):
        """
        Calculate the potential with 3N-dimensional argument.
        """
        return self.potential_cart(x.reshape((self.natoms,3)))
    
    def _dipole(self, x):
        return self.dipole_cart(x.reshape((self.natoms, 3)))

    def _nmode_potential(self, x):
        V = 0.0
        V += self.nm.V0
        q = self.nm._cart2normal(x)

        if self.Order >= 1:
            for i in range(self.nm.nmodes):
                V += self.nm.potential_1mode(i, q[i])
        if self.Order >= 2:
            for i in range(self.nm.nmodes):
                for j in range(i + 1, self.nm.nmodes):
                    V += self.nm.potential_2mode(i, j, q[i], q[j])
        if self.Order >= 3:
            for i in range(self.nm.nmodes):
                for j in range(i + 1, self.nm.nmodes):
                    for k in range(j + 1, self.nm.nmodes):
                        V += self.nm.potential_3mode(i, j, k, q[i], q[j], q[k])
        return V
    
    def inverse_moment_of_inertia_cart(self, x):
        I = np.zeros((3, 3))
        for i in range(self.natoms):
            r = x[i] - self.CoM
            m = self.mass[i]
            I += m * (np.dot(r, r) * np.eye(3) - np.outer(r, r))
        return np.linalg.inv(I).reshape(-1)

    def _inverse_moment_of_inertia(self, x):
        return self.inverse_moment_of_inertia_cart(x.reshape((self.natoms, 3)))

    def ShiftPotential(self, V0):
        orig_pot = self.potential_cart
        def shifted_pot(coord):
            return orig_pot(coord) - V0
        self.potential_cart = shifted_pot

    def CalcNM(self, x0 = None):
        self.nm = NormalModes(self)
        if self.ReadGeom:
            self.ReadGeometry()
        else:
            self.nm.kernel(x0 = x0, doGeomOpt = self.doGeomOpt)
        self.Nm = self.nm.freqs.shape[0]
        #debug!!
        #c = np.zeros((self.natoms * 3, self.nm.nmodes))
        #c[:3,:3] = np.eye(3)
        #self.nm.nm_coeff = c.reshape((self.natoms, 3, self.nm.nmodes))
        #print(self.nm.freqs)
        #self.potential_cart = ho_3d
        #self.nm.x0 = np.zeros_like(self.nm.x0)
        # debug!!
        self.Frequencies = self.nm.freqs * constants.AU_TO_INVCM
        self.x0 = self.nm.x0 / constants.ANGSTROM_TO_AU
        self.V0 = self.nm.V0 * constants.AU_TO_INVCM
        if self.doShiftPotential:
            self.ShiftPotential(self.nm.V0)
            self.nm.V0 = 0
            self.V0 = 0
        self.mu0 = self.nm.mu0

        self.CoM = self.calc_com(self.nm.x0)

    def CalcNModePotential(self, Order = None, OrderPlus = None):
        if Order is None:
            Order = self.Order
        if OrderPlus is None:
            OrderPlus = self.OrderPlus
        self.nmode = NModePotential(self.nm)
        self.ints = [np.asarray([]), np.asarray([[[[[[]]]]]]), np.asarray([[[[[[[[[]]]]]]]]])]
        if self.ReadInt:
            self.ReadIntegrals()
        else:
            for i in range(Order):
                self.ints[i] = self.nmode.get_ints(i+1, ngridpts = self.ngridpts, onemode_coeff = self.onemode_coeff)
                if i == 0:
                    for j in range(self.Nm):
                        OMBasis = init_funcs.InitGridBasis([self.Frequencies[j]], [self.ngridpts])[0]
                        OMH = VCISparseHamNMode(OMBasis, OMBasis, [self.Frequencies[j]], 0.0, [self.ints[0][j].tolist()], self.ints[1].tolist(), self.ints[2].tolist(), True)
                        e, v = np.linalg.eigh(OMH.todense())
                        self.onemode_coeff.append(v)
                        self.onemode_eig.append(e)
            if OrderPlus is not None:
                for i in range(Order, OrderPlus):
                    if i == 2:
                        self.Divergent3Modes = self.FindDivergent3Modes()
                    self.ints[i] = self.nmode.get_ints(OrderPlus, ngridpts = self.ngridpts, onemode_coeff = self.onemode_coeff, modes = self.Divergent3Modes)
    
    def CalcNModeDipole(self, Order = None):
        if Order is None:
            Order = self.Order
        
        # Did not calculate one-mode integrals yet. Need to do that now
        if len(self.onemode_coeff) == 0:
            self.nmode = NModePotential(self.nm)
            ints1 = self.nmode.get_ints(1, ngridpts = self.ngridpts, onemode_coeff = self.onemode_coeff)
            for j in range(self.Nm):
                OMBasis = init_funcs.InitGridBasis([self.Frequencies[j]], [self.ngridpts])[0]
                OMH = VCISparseHamNMode(OMBasis, OMBasis, [self.Frequencies[j]], self.V0, [self.ints[0][j].tolist()], [[[[[[]]]]]], [[[[[[[[[]]]]]]]]], True)
                e, v = np.linalg.eigh(OMH.todense())
                self.onemode_coeff.append(v)
                self.onemode_eig.append(e)

        self.dip_ints = [np.asarray([[]] * 3), np.asarray([[[[[[[]]]]]]] * 3), np.asarray([[[[[[[[[[]]]]]]]]]] * 3)]
        if self.ReadDip:
            self.ReadDipoles()
        else:
            for i in range(Order):
                self.dip_ints[i] = self.nmode.get_dipole_ints(i + 1, ngridpts = self.ngridpts, onemode_coeff = self.onemode_coeff, usePyPotDip = self.usePyPotDip)

        # separate dipole and energy as necessary
        if self.usePyPotDip and not self.ReadDip:
            self.ints = [np.asarray([]), np.asarray([[[[[[]]]]]]), np.asarray([[[[[[[[[]]]]]]]]])]
            for n in range(Order):
                if n + 1 == 1:
                    self.ints[n] = np.empty(self.Nm, dtype=object)
                    for i in range(self.Nm):
                        self.ints[n][i] = self.dip_ints[n][0][i] * constants.AU_TO_INVCM
                    self.dip_ints[n] = np.delete(self.dip_ints[n], (0), axis = 0)
                if n + 1 == 2:
                    self.ints[n] = np.empty((self.Nm, self.Nm), dtype=object)
                    for i in range(self.Nm):
                        for j in range(self.Nm):
                            self.ints[n][i, j] = self.dip_ints[n][0][i, j] * constants.AU_TO_INVCM
                    self.dip_ints[n] = np.delete(self.dip_ints[n], (0), axis = 0)
                if n + 1 == 3:
                    self.ints[n] = np.empty((self.Nm, self.Nm, self.Nm), dtype=object)
                    for i in range(self.Nm):
                        for j in range(self.Nm):
                            for k in range(self.Nm):
                                self.ints[n][i, j, k] = self.dip_ints[n][0][i, j, k] * constants.AU_TO_INVCM
                    self.dip_ints[n] = np.delete(self.dip_ints[n], (0), axis = 0)

    def CalcNModeInvInertia(self, Order = None):
        if Order is None:
            Order = self.Order
        
        # Did not calculate one-mode integrals yet. Need to do that now
        if len(self.onemode_coeff) == 0:
            self.nmode = NModePotential(self.nm)
            ints1 = self.nmode.get_ints(1, ngridpts = self.ngridpts, onemode_coeff = self.onemode_coeff)
            for j in range(self.Nm):
                OMBasis = init_funcs.InitGridBasis([self.Frequencies[j]], [self.ngridpts])[0]
                OMH = VCISparseHamNMode(OMBasis, OMBasis, [self.Frequencies[j]], self.V0, [self.ints[0][j].tolist()], [[[[[[]]]]]], [[[[[[[[[]]]]]]]]], True)
                e, v = np.linalg.eigh(OMH.todense())
                self.onemode_coeff.append(v)
                self.onemode_eig.append(e)

        self.inv_inertia_ints = [np.asarray([[]] * 3), np.asarray([[[[[[[]]]]]]] * 3), np.asarray([[[[[[[[[[]]]]]]]]]] * 3)]
        if self.ReadInvInertia:
            self.ReadInvInertia()
        else:
            for i in range(Order):
                self.inv_inertia_ints[i] = self.nmode.get_inv_inertia_ints(i + 1, ngridpts = self.ngridpts, onemode_coeff = self.onemode_coeff)

    def InitializeBasis(self):
        self.FullBasis = init_funcs.InitTruncatedBasis(self.Nm, self.Frequencies, self.FullMaxQuanta, self.FullTotalQuanta)
        self.Basis = init_funcs.InitTruncatedBasis(self.Nm, self.Frequencies, self.InitMaxQuanta, self.InitTotalQuanta)

    def ReorganizeBasis(self):
        IndexBasis = []
        for B in self.Basis:
            IndexBasis.append(self.FullBasis.index(B))
        IndexOther = [i for i in range(len(self.FullBasis)) if i not in IndexBasis]
        FullBasis = [self.FullBasis[i] for i in IndexBasis] + [self.FullBasis[i] for i in IndexOther]
        self.FullBasis = FullBasis
        return FullBasis, IndexBasis, IndexOther

    def SaveIntegrals(self, IntsFile = None):
        if IntsFile is None:
            IntsFile = self.IntsFile

        if self.OrderPlus is None:
            MaxOrder = self.Order
        else:
            MaxOrder = self.OrderPlus
        
        with h5py.File(IntsFile, "a") as f:
            if "ints" in f:
                del f["ints"]
            g = f.create_group("ints")
            g1 = g.create_group("1")
            for i in range(self.Nm):
                g1.create_dataset("%d" % (i + 1), data = self.ints[0][i])
            if MaxOrder >= 2:
                g2 = g.create_group("2")
                for i in range(self.Nm):
                    for j in range(self.Nm):
                        g2.create_dataset("%d_%d" % (i + 1, j + 1), data = self.ints[1][i, j])
                if MaxOrder >= 3:
                    g3 = g.create_group("3")
                    for i in range(self.Nm):
                        for j in range(self.Nm):
                            for k in range(self.Nm):
                                g3.create_dataset("%d_%d_%d" %(i + 1, j + 1, k + 1), data = self.ints[2][i, j, k])

            if "onemode_coeff" in f:
                del f["onemode_coeff"]
            f1c = f.create_group("onemode_coeff")
            for i in range(self.Nm):
                f1c.create_dataset("%d" % (i + 1), data = self.onemode_coeff[i])
            if "onemode_eig" in f:
                del f["onemode_eig"]
            f1e = f.create_group("onemode_eig")
            for i in range(self.Nm):
                f1e.create_dataset("%d" % (i + 1), data = self.onemode_eig[i])
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

    def SaveDipoles(self, IntsFile = None):
        if IntsFile is None:
            IntsFile = self.IntsFile
        cart_coord = ['x', 'y', 'z']
        
        with h5py.File(IntsFile, "a") as f:
            if "dip_ints" in f:
                del f["dip_ints"]
            g = f.create_group("dip_ints")
            g1 = g.create_group("1")
            for x in range(3):
                g1x = g1.create_group(cart_coord[x])
                for i in range(self.Nm):
                    g1x.create_dataset("%d" % (i + 1), data = self.dip_ints[0][x, i])
            if self.Order >= 2:
                g2 = g.create_group("2")
                for x in range(3):
                    g2x = g2.create_group(cart_coord[x])
                    for i in range(self.Nm):
                        for j in range(self.Nm):
                            g2x.create_dataset("%d_%d" % (i + 1, j + 1), data = self.dip_ints[1][x, i, j])
                if self.Order >= 3:
                    g3 = g.create_group("3")
                    for x in range(3):
                        g3x = g3.create_group(cart_coord[x])
                        for i in range(self.Nm):
                            for j in range(self.Nm):
                                for k in range(self.Nm):
                                    g3x.create_dataset("%d_%d_%d" %(i + 1, j + 1, k + 1), data = self.dip_ints[2][x, i, j, k])

            if "onemode_coeff" in f:
                del f["onemode_coeff"]
            f1c = f.create_group("onemode_coeff")
            for i in range(self.Nm):
                f1c.create_dataset("%d" % (i + 1), data = self.onemode_coeff[i])
            if "onemode_eig" in f:
                del f["onemode_eig"]
            f1e = f.create_group("onemode_eig")
            for i in range(self.Nm):
                f1e.create_dataset("%d" % (i + 1), data = self.onemode_eig[i])
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

    def SaveInvInertia(self, IntsFile = None):
        if IntsFile is None:
            IntsFile = self.IntsFile
        cart_coord = ['xx', 'xy', 'xz', 'yx', 'yy', 'yz', 'zx', 'zy', 'zz']

        with h5py.File(IntsFile, "a") as f:
            if "inv_inertia_ints" in f:
                del f["inv_inertia_ints"]
            g = f.create_group("inv_inertia_ints")
            if self.Order >= 2:
                g2 = g.create_group("2")
                for x in range(9):
                    g2x = g2.create_group(cart_coord[x])
                    for i in range(self.Nm):
                        for j in range(self.Nm):
                             g2x.create_dataset("%d_%d" % (i + 1, j + 1), data = self.inv_inertia_ints[1][x, i, j])
                if self.Order >= 3:
                    g3 = g.create_group("3")
                    for x in range(9):
                        g3x = g3.create_group(cart_coord[x])
                        for i in range(self.Nm):
                            for j in range(self.Nm):
                                for k in range(self.Nm):
                                    g3x.create_dataset("%d_%d_%d" %(i + 1, j + 1, k + 1), data = self.inv_inertia_ints[2][x, i, j, k])

    def ReadGeometry(self, IntsFile = None):
        if IntsFile is None:
            IntsFile = self.IntsFile

        with h5py.File(IntsFile, "r") as f:
            self.nm.x0 = f["x0"][()]
            self.nm.nm_coeff = f["nm_coeff"][()]
            self.nm.freqs = f["freq"][()] / constants.AU_TO_INVCM
            self.nm.V0 = f["V0"][()]
            self.nm.mu0 = f["mu0"][()]

    def ReadIntegrals(self, IntsFile = None):
        if IntsFile is None:
            IntsFile = self.IntsFile
        MaxOrder = self.Order
        if self.OrderPlus is not None:
            MaxOrder = self.OrderPlus

        with h5py.File(IntsFile, "r") as f:
            for n in range(MaxOrder):
                if n == 0:
                    self.ints[n] = np.empty(self.Nm, dtype = object)
                    for i in range(self.Nm):
                        self.ints[n][i] = f["ints/%d/%d" % (n + 1, i + 1)][()]
                if n == 1:
                    self.ints[n] = np.empty((self.Nm, self.Nm), dtype = object)
                    for i in range(self.Nm):
                        for j in range(self.Nm):
                            self.ints[n][i, j] = f["ints/%d/%d_%d" % (n + 1, i + 1, j + 1)][()]
                if n == 2:
                    self.ints[n] = np.empty((self.Nm, self.Nm, self.Nm), dtype = object)
                    for i in range(self.Nm):
                        for j in range(self.Nm):
                            for k in range(self.Nm):
                                self.ints[n][i, j, k] = f["ints/%d/%d_%d_%d" % (n + 1, i + 1, j + 1, k + 1)][()]

            self.onemode_eig = []
            for i in range(self.Nm):
                self.onemode_eig.append(f["onemode_eig/%d" % (i + 1)][()])
            self.onemode_coeff = []
            for i in range(self.Nm):
                self.onemode_coeff.append(f["onemode_coeff/%d" % (i + 1)][()])

    def ReadDipoles(self, IntsFile = None):
        if IntsFile is None:
            IntsFile = self.IntsFile
        cart_coord = ['x', 'y', 'z']

        with h5py.File(IntsFile, "r") as f:
            for n in range(self.Order):
                if n == 0:
                    self.dip_ints[n] = np.empty((3, self.Nm), dtype = object)
                    for x in range(3):
                        for i in range(self.Nm):
                            self.dip_ints[n][x, i] = f["dip_ints/%d/%s/%d" % (n + 1, cart_coord[x], i + 1)][()]
                if n == 1:
                    self.dip_ints[n] = np.empty((3, self.Nm, self.Nm), dtype = object)
                    for x in range(3):
                        for i in range(self.Nm):
                            for j in range(self.Nm):
                                self.dip_ints[n][x, i, j] = f["dip_ints/%d/%s/%d_%d" % (n + 1, cart_coord[x], i + 1, j + 1)][()]
                if n == 2:
                    self.dip_ints[n] = np.empty((3, self.Nm, self.Nm, self.Nm), dtype = object)
                    for x in range(3):
                        for i in range(self.Nm):
                            for j in range(self.Nm):
                                for k in range(self.Nm):
                                    self.dip_ints[n][x, i, j, k] = f["dip_ints/%d/%s/%d_%d_%d" % (n + 1, cart_coord[x], i + 1, j + 1, k + 1)][()]

    def ReadInvInertia(self, IntsFile = None):
        if IntsFile is None:
            IntsFile = self.IntsFile
        cart_coord = ['xx', 'xy', 'xz', 'yx', 'yy', 'yz', 'zx', 'zy', 'zz']

        with h5py.File(IntsFile, "r") as f:
            for n in range(self.Order):
                if n == 1:
                    self.inv_intertia_ints[n] = np.empty((9, self.Nm, self.Nm), dtype = object)
                    for x in range(9):
                        for i in range(self.Nm):
                            for j in range(self.Nm):
                                self.inv_inertia_ints[n][x, i, j] = f["inv_inertia_ints/%d/%s/%d_%d" % (n + 1, cart_coord[x], i + 1, j + 1)][()]
                if n == 2:
                    self.dip_ints[n] = np.empty((9, self.Nm, self.Nm, self.Nm), dtype = object)
                    for x in range(9):
                        for i in range(self.Nm):
                            for j in range(self.Nm):
                                for k in range(self.Nm):
                                    self.inv_inertia_ints[n][x, i, j, k] = f["inv_inertia_ints/%d/%s/%d_%d_%d" % (n + 1, cart_coord[x], i + 1, j + 1, k + 1)][()]

    def ZeroCoupling(self, FreqTol = 1000, NCutoff = 0):
        LowFreqModes = np.where(self.Frequencies < FreqTol)[0]
        HighFreqModes = np.where(self.Frequencies > FreqTol)[0]
        for i in LowFreqModes:
            for j in HighFreqModes:
                '''
                D = np.diagonal(np.diagonal(np.diagonal(self.ints[1][i, j])))
                self.ints[1][i, j] = np.zeros_like(self.ints[1][i, j])
                self.ints[1][j, i] = np.zeros_like(self.ints[1][j, i])
                np.fill_diagonal(self.ints[1][i, j], D)
                '''
                for n in range(self.ints[1][i, j].shape[0]):
                    for m in range(self.ints[1][i, j].shape[3]):
                        '''
                        self.ints[1][i, j][n, :, :, m] = np.zeros_like(self.ints[1][i, j][n, :, :, m])
                        self.ints[1][j, i][n, :, :, m] = np.zeros_like(self.ints[1][j, i][n, :, :, m])
                        self.ints[1][i, j][:, n, m, :] = np.zeros_like(self.ints[1][i, j][:, n, m, :])
                        self.ints[1][j, i][:, n, m, :] = np.zeros_like(self.ints[1][j, i][:, n, m, :])
                        '''
                        self.ints[1][i, j][n, :, :, m] *= np.tri(self.ints[1][i, j].shape[1], k = NCutoff)
                        self.ints[1][i, j][n, :, :, m] *= np.tri(self.ints[1][i, j].shape[1], k = NCutoff).T
                        self.ints[1][j, i][n, :, :, m] *= np.tri(self.ints[1][j, i].shape[1], k = NCutoff)
                        self.ints[1][j, i][n, :, :, m] *= np.tri(self.ints[1][j, i].shape[1], k = NCutoff).T
                        self.ints[1][i, j][:, n, m, :] *= np.tri(self.ints[1][i, j].shape[1], k = NCutoff)
                        self.ints[1][i, j][:, n, m, :] *= np.tri(self.ints[1][i, j].shape[1], k = NCutoff).T
                        self.ints[1][j, i][:, n, m, :] *= np.tri(self.ints[1][j, i].shape[1], k = NCutoff)
                        self.ints[1][j, i][:, n, m, :] *= np.tri(self.ints[1][j, i].shape[1], k = NCutoff).T

    def ZeroCoupling2(self, FreqTol1 = 1200, FreqTol2 = 2000):
        LowFreqModes = np.where(self.Frequencies < FreqTol1)[0]
        MidFreqModes = np.where((self.Frequencies > FreqTol1) & (self.Frequencies < FreqTol2))[0]
        HighFreqModes = np.where(self.Frequencies > FreqTol2)[0]
        for i in LowFreqModes:
            for j in HighFreqModes:
                if i != j:
                    self.ints[1][i, j] = np.zeros_like(self.ints[1][i, j])
                    self.ints[1][j, i] = np.zeros_like(self.ints[1][j, i])
        for i in MidFreqModes:
            for j in HighFreqModes:
                if i != j:
                    self.ints[1][i, j] = np.zeros_like(self.ints[1][i, j])
                    self.ints[1][j, i] = np.zeros_like(self.ints[1][j, i])
        for i in LowFreqModes:
            for j in MidFreqModes:
                if i != j:
                    self.ints[1][i, j] = np.zeros_like(self.ints[1][i, j])
                    self.ints[1][j, i] = np.zeros_like(self.ints[1][j, i])
    
    def ZeroDivergentCoupling(self, FreqTol = 1000):
        Divergent3Modes = self.FindDivergent3Modes()
        LowFreqModes = np.where(self.Frequencies < FreqTol)[0]
        HighFreqModes = np.where(self.Frequencies > FreqTol)[0]
        for i in LowFreqModes:
            for j in HighFreqModes:
                for Triplet in Divergent3Modes:
                    if i in Triplet and j in Triplet:
                        self.ints[1][i, j] = np.zeros_like(self.ints[1][i, j])
                        self.ints[1][j, i] = np.zeros_like(self.ints[1][j, i])
                        break

    def ScanPES(self, Modes, Range = None, PlottedModes = None, NPoints = 10, SavePES = None):
        if Range is None:
            Range = [[-50, 50]] * len(Modes)
        ModePts = []
        for r in Range:
            ModePts.append(np.linspace(r[0], r[1], NPoints))

        xs = np.empty([NPoints] * len(Modes), dtype = object)
        PESNMode = np.empty([NPoints] * len(Modes), dtype = object)
        PESTrue = np.empty([NPoints] * len(Modes), dtype = object)
        qind = [0] * len(Modes)

        while(True):
            q = np.zeros(self.Nm)
            print(qind, flush = True)
            for m in range(len(Modes)):
                q[Modes[m]] = ModePts[m][qind[m]]
            x = self.nm._normal2cart(q)

            xs[tuple(qind)] = x
            PESTrue[tuple(qind)] = self.potential_cart(x)
            PESNMode[tuple(qind)] = self._nmode_potential(x)

            qind[0] = qind[0] + 1
            for m in range(len(Modes) - 1):
                if qind[m] == NPoints:
                    qind[m] = 0
                    qind[m + 1] = qind[m + 1] + 1
            if qind[-1] == NPoints:
                break

        if SavePES is not None:
            with h5py.File(SavePES, "a") as f:
                f.create_dataset("xs", data = xs)
                f.create_dataset("pes_nmode", data = PESNMode)
                f.create_dataset("pes_true", data = PESTrue)
        if PlottedModes is not None:
            midpt = NPoints // 2
            if len(PlottedModes) == 2:
                xax = ModePts[PlottedModes[0]]
                yax = ModePts[PlottedModes[1]]
            plot_pes_nmode = np.zeros((NPoints, NPoints))
            plot_pes_true = np.zeros((NPoints, NPoints))
            for i in range(NPoints):
                for j in range(NPoints):
                    ind = [midpt] * len(Modes)
                    ind[PlottedModes[0]] = i
                    ind[PlottedModes[1]] = j
                    plot_pes_nmode[i, j] = PESNMode[tuple(ind)]
                    plot_pes_true[i, j] = PESTrue[tuple(ind)]
            import matplotlib.pyplot as plt
            Y, X = np.meshgrid(xax, yax)
            plt.pcolormesh(X, Y, plot_pes_true, cmap = 'RdBu')
            plt.colorbar()
            plt.title('True PES')
            plt.savefig("pes_true.png")
            plt.pcolormesh(X, Y, plot_pes_nmode, cmap = 'RdBu')
            plt.title('NMode PES')
            plt.savefig("pes_nmode.png")

    def FindDivergent3Modes(self):
        DivergentTriplets = []
        V0 = self._nmode_potential(self.nm.x0)
        for i in range(self.Nm):
            for j in range(i + 1, self.Nm):
                for k in range(j + 1, self.Nm):
                    q = np.zeros(self.Nm)
                    q[i] = 1
                    q[j] = 1
                    q[k] = 1
                    q = q * 50
                    x = self.nm._normal2cart(q)
                    if (self._nmode_potential(x) - V0) < 0:
                        print("Divergent triplet found: ", i, j, k, flush=True)
                        DivergentTriplets.append((i, j, k))
        return DivergentTriplets

    def kernel(self, x0 = None):
        self.CalcNM(x0 = x0)
        if self.calc_integrals:
            print("Calculating integrals...", flush = True)
            self.CalcNModePotential(OrderPlus = self.OrderPlus)
        if self.calc_dipole:
            if not self.calc_integrals:
                self.CalcNModePotential(Order = 1)
            print("Calculating dipoles...", flush = True)
            self.CalcNModeDipole()
        print(self)

class NormalModes():

    def __init__(self, mol):
        self.mol = mol
        self.V0 = 0
        self.nm_coeff = None
        self.freqs = None

        self.nmodes = 3*self.mol.natoms - 6

    def kernel(self, x0=None, doGeomOpt = True):
        """
        Calculate the normal modes of the potential defined in mol.
    
        Parameters:
        x0 ((natoms,3) ndarray): Guess of the minimum of the PES
        """

        natoms = self.mol.natoms

        if x0 is None:
            x0 = np.random.random((natoms,3))

        pes = self.mol._potential
        if doGeomOpt:
            print("Beginning geometry optimization...", flush = True)
            result = scipy.optimize.minimize(pes, x0, tol = 1e-12)
            print("...geometry optimization complete.", flush = True)
            self.x0 = result.x.reshape((natoms,3))
        else:
            self.x0 = x0.reshape((natoms,3))
        self.V0 = pes(self.x0.reshape(-1))
        try:
            self.mu0 = self.mol._dipole(self.x0)
        except:
            self.mu0 = np.asarray([0.0, 0.0, 0.0])
        try:
            self.B0 = self.mol._inv_moment_of_inertia(self.x0)
        except:
            self.B0 = np.zeros((9))
        
        hess0 = self.hessian(self.x0).reshape((natoms*3, natoms*3))

        e, v = scipy.linalg.eigh(hess0)
        e, v = e[6:], v[:,6:] # remove translations and rotations
        v = v.reshape((natoms,3,self.nmodes))
        self.freqs, self.nm_coeff = np.sqrt(e), v
        if self.mol.LowFrequencyCutoff is not None:
            Cutoff = self.mol.LowFrequencyCutoff / constants.AU_TO_INVCM
            Keep = self.freqs > Cutoff
            self.freqs = self.freqs[Keep]
            self.nm_coeff = self.nm_coeff[:, :, Keep]
            self.nmodes = self.nm_coeff.shape[2]
        
        return self.freqs, self.nm_coeff

    def gradient(self, x):
        pes = self.mol._potential
        grad = nd.Gradient(pes)(x.reshape(-1))
        return grad.reshape((self.mol.natoms,3))

    def hessian(self, x):
        """
        Calculate the mass-weighted Hessian.
        """
        natoms = self.mol.natoms
        mass = self.mol.mass
        pes = self.mol._potential
        hess = nd.Hessian(pes)(x.reshape(-1))
        hess = hess.reshape((natoms,3,natoms,3))
        for i in range(natoms):
            for j in range(natoms):
                hess[i,:,j,:] /= np.sqrt(mass[i]*mass[j])
        return hess

    def get_ff(self, dx = 1e-4, Order = 4):
        """
        Calculate the cubic and quartic derivatives of the potential with respect to normal modes
        """
        natoms = self.mol.natoms
        mass = self.mol.mass
        pes = self.mol._potential
        freq_cm = self.freqs * constants.AU_TO_INVCM
        def hess(x):
            C = np.einsum('nij,n->nij', self.nm_coeff, 1/np.sqrt(mass))
            C = C.reshape(3 * natoms, self.nmodes)
            return np.einsum('ki,kl,lj->ij', C, nd.Hessian(pes)(x.reshape(-1)), C)
        H0 = hess(self.x0)
        V3 = []
        V4 = []
        for i in range(self.nmodes):
            q = np.zeros(self.nmodes)
            q[i] = dx
            x = self._normal2cart(q)
            Hpi = hess(x)

            q = np.zeros(self.nmodes)
            q[i] = -dx
            x = self._normal2cart(q)
            Hmi = hess(x)

            d3V = (Hpi - Hmi) / (2 * dx)
            d4Vd = (Hpi + Hmi - 2 * H0) / (dx**2)
            for j in range(i, self.nmodes):
                for k in range(j, self.nmodes):
                    Vijk = ScaleFC_me(d3V[j, k], freq_cm, [i, j, k])
                    Viijk = ScaleFC_me(d4Vd[j, k], freq_cm, [i, i, j, k])
                    if abs(Vijk) > 1.0:
                        V3.append((Vijk, [i, j, k]))
                    if abs(Viijk) > 1.0:
                        V4.append((Viijk, [i, i, j, k]))
                
            for j in range(i + 1, self.nmodes):
                q = np.zeros(self.nmodes)
                q[i] = dx
                q[j] = dx
                x = self._normal2cart(q)
                Hpipj = hess(x)

                q = np.zeros(self.nmodes)
                q[i] = dx
                q[j] = -dx
                x = self._normal2cart(q)
                Hpimj = hess(x)

                q = np.zeros(self.nmodes)
                q[i] = -dx
                q[j] = dx
                x = self._normal2cart(q)
                Hmipj = hess(x)

                q = np.zeros(self.nmodes)
                q[i] = -dx
                q[j] = -dx
                x = self._normal2cart(q)
                Hmimj = hess(x)

                d4V = (Hpipj + Hmimj - Hmipj - Hpimj) / (4 * dx**2)
                for k in range(j, self.nmodes):
                    for l in range(k, self.nmodes):
                        Vijkl = ScaleFC_me(d4V[k, l], freq_cm, [i, j, k, l])
                        if abs(Vijkl) > 1.0:
                            V4.append((Vijkl, [i, j, k, l]))
        return V3, V4

    def _normal2cart(self, q):
        return self.x0 + np.einsum('n,ndi,i->nd',1/np.sqrt(self.mol.mass),self.nm_coeff,q)

    def _cart2normal(self, x):
        return np.einsum('ndi,n,nd->i', self.nm_coeff, np.sqrt(self.mol.mass), x-self.x0)

    def potential_1mode(self, i, qi):
        """
        Calculate the 1-mode potential.
        
        Parameters:
        i (int): the mode index
        qi (float): the normal mode coordinate
        """
        q = np.zeros(self.nmodes)
        q[i] = qi
        x = self._normal2cart(q)
        return self.mol.potential_cart(x) - self.V0

    def potential_2mode(self, i, j, qi, qj):
        q = np.zeros(self.nmodes)
        q[i] = qi
        q[j] = qj
        x = self._normal2cart(q)
        return (self.mol.potential_cart(x)
                - self.potential_1mode(i,qi) - self.potential_1mode(j,qj) - self.V0)

    def potential_3mode(self, i, j, k, qi, qj, qk):
        q = np.zeros(self.nmodes)
        q[i] = qi
        q[j] = qj
        q[k] = qk
        x = self._normal2cart(q)
        return (self.mol.potential_cart(x) - self.potential_2mode(i, j, qi, qj) - self.potential_2mode(i, k, qi, qk) - self.potential_2mode(j, k, qj, qk) - self.potential_1mode(i, qi) - self.potential_1mode(j, qj) - self.potential_1mode(k, qk) - self.V0)

    def dipole_1mode(self, i, qi):
        """
        Calculate the 1-mode potential.
        
        Parameters:
        i (int): the mode index
        qi (float): the normal mode coordinate
        """
        q = np.zeros(self.nmodes)
        q[i] = qi
        x = self._normal2cart(q)
        return self.mol.dipole_cart(x) - self.mu0

    def dipole_2mode(self, i, j, qi, qj):
        q = np.zeros(self.nmodes)
        q[i] = qi
        q[j] = qj
        x = self._normal2cart(q)
        return (self.mol.dipole_cart(x)
                - self.dipole_1mode(i,qi) - self.dipole_1mode(j,qj) - self.mu0)

    def dipole_3mode(self, i, j, k, qi, qj, qk):
        q = np.zeros(self.nmodes)
        q[i] = qi
        q[j] = qj
        q[k] = qk
        x = self._normal2cart(q)
        return (self.mol.dipole_cart(x) - self.dipole_2mode(i, j, qi, qj) - self.dipole_2mode(i, k, qi, qk) - self.dipole_2mode(j, k, qj, qk) - self.dipole_1mode(i, qi) - self.dipole_1mode(j, qj) - self.dipole_1mode(k, qk) - self.mu0)

    def inv_moment_of_inertia_1mode(self, i, qi):
        """
        Calculate the 1-mode potential.
        
        Parameters:
        i (int): the mode index
        qi (float): the normal mode coordinate
        """
        q = np.zeros(self.nmodes)
        q[i] = qi
        x = self._normal2cart(q)
        return self.mol.inv_moment_of_inertia_cart(x) - self.B0

    def inv_moment_of_inertia_2mode(self, i, j, qi, qj):
        q = np.zeros(self.nmodes)
        q[i] = qi
        q[j] = qj
        x = self._normal2cart(q)
        return (self.mol.inv_moment_of_inertia_cart(x)
                - self.inv_moment_of_inertia_1mode(i,qi) - self.inv_moment_of_inertia_1mode(j,qj) - self.B0)

    def inv_moment_of_inertia_3mode(self, i, j, k, qi, qj, qk):
        q = np.zeros(self.nmodes)
        q[i] = qi
        q[j] = qj
        q[k] = qk
        x = self._normal2cart(q)
        return (self.mol.inv_moment_of_inertia_cart(x) - self.inv_moment_of_inertia_2mode(i, j, qi, qj) - self.inv_moment_of_inertia_2mode(i, k, qi, qk) - self.inv_moment_of_inertia_2mode(j, k, qj, qk) - self.inv_moment_of_inertia_1mode(i, qi) - self.inv_moment_of_inertia_1mode(j, qj) - self.inv_moment_of_inertia_1mode(k, qk) - self.B0)

    def potential_nm(self, q):
        x = self._normal2cart(q)
        return self.mol.potential_cart(x)

    def potential_nm_i(self, *argv):
        q = np.zeros(self.nmodes)
        for i in range(self.nmodes):
            q[i] = argv[i]
        return self.potential_nm(q)

    def potential_nm_torch(self, *argv):
        argv_np = [arg.detach().numpy() for arg in argv]
        potential_nm_vec = np.vectorize(self.potential_nm_i)
        return torch.tensor(potential_nm_vec(*argv_np))

    '''
    def get_heg(self, ngridpts = None):
        nmode = NModePotential(self)
        if ngridpts is None:
            ngridpts = [12]*self.nmodes
        elif isinstance(ngridpts, (int, np.integer)):
            ngridpts = [ngridpts]*self.nmodes
        self.gridpts, self.dvr_coeff = nmode.get_heg(ngridpts)
    '''

    def do_tci(self, gridpts):
        gridpts_torch = [torch.tensor(grid) for grid in gridpts]
        t = tn.cross(function = self.potential_nm_torch, domain = gridpts_torch)
        return t.cores

def get_qmat_ho(omega, nmax):
    qmat = np.zeros((nmax,nmax))
    for n in range(nmax-1):
        qmat[n,n+1] = qmat[n+1,n] = np.sqrt((n+1)/(2*omega))
    return qmat


def get_tmat_ho(omega, nmax):
    tmat = np.zeros((nmax,nmax))
    for n in range(nmax):
        tmat[n,n] = omega/2*(n+0.5)
        if n < nmax-2:
            tmat[n,n+2] = tmat[n+2,n] = -omega/4*np.sqrt((n+1)*(n+2))
    return tmat


class HEG1Mode():

    def __init__(self, nm):
        self.nm = nm

    def kernel(self, ngridpts=None):
        """
        Calculate the HEG grid points and functions by diagonalizing the Q
        matrix in a basis of harmonic oscillator eigenfunctions.

        Parameters:
        ngridpts (int or list of int): The number of harmonic oscillator
            eigenfunctions (per mode) to use when diagonalizing Q.
        """
        if ngridpts is None:
            ngridpts = [12]*self.nm.nmodes
        elif isinstance(ngridpts, (int, np.integer)):
            ngridpts = [ngridpts]*self.nm.nmodes 

        self.ngridpts = ngridpts
        gridpts = list()
        coeff = list()
        for i in range(self.nm.nmodes):
            nho = ngridpts[i]
            omega = self.nm.freqs[i]
            qmat = get_qmat_ho(omega, nho)
            q, u = scipy.linalg.eigh(qmat)
            gridpts.append(q)
            coeff.append(u)

        self.gridpts, self.coeff = gridpts, coeff
        return self.gridpts, self.coeff


class HEGOpt1Mode():
    def __init__(self, heg0):
        self.heg0 = heg0
        self.nm = heg0.nm

    def kernel(self, ngridpts):
        if ngridpts is None:
            ngridpts = [8]*self.nm.nmodes
        elif isinstance(ngridpts, (int, np.integer)):
            ngridpts = [ngridpts]*self.nm.nmodes 

        self.ngridpts = ngridpts

        gridpts = list()
        coeff = list()
        for i in range(self.nm.nmodes):
            nho = self.heg0.ngridpts[i]
            omega = self.nm.freqs[i]
            q0, v0 = self.heg0.gridpts[i], self.heg0.coeff[i]
            tmat = get_tmat_ho(omega, nho)
            vgrid = np.array([self.nm.potential_1mode(i,qi) for qi in q0]) # this should be vectorized
            vmat = np.dot(v0, np.dot(vgrid, v0.T))
            e, v = scipy.linalg.eigh(tmat+vmat)
            v = v[:,:ngridpts[i]]
            qmat0 = get_qmat_ho(omega, nho)
            qmat = np.dot(v.T, np.dot(qmat0, v))
            q, u = scipy.linalg.eigh(qmat)
            gridpts.append(q)
            coeff.append(u)

        self.gridpts, self.coeff = gridpts, coeff
        return self.gridpts, self.coeff
        

class NModePotential():
    def __init__(self, nm):
        self.nm = nm

    def get_ints(self, nmode, ngridpts=None, optimized=False, ngridpts0=None, onemode_coeff = None, modes = None):
        print("Calculating n-Mode integrals for n =", nmode, flush = True) 
        if optimized is False:
            if ngridpts is None:
                ngridpts = [12]*self.nm.nmodes
            elif isinstance(ngridpts, (int, np.integer)):
                ngridpts = [ngridpts]*self.nm.nmodes

            gridpts, coeff = self.get_heg(ngridpts)

        else:
            if ngridpts is None:
                ngridpts = [8]*self.nm.nmodes
            elif isinstance(ngridpts, (int, np.integer)):
                ngridpts = [ngridpts]*self.nm.nmodes

            if ngridpts0 is None:
                ngridpts0 = [12]*self.nm.nmodes
            elif isinstance(ngridpts0, (int, np.integer)):
                ngridpts0 = [ngridpts0]*self.nm.nmodes

            gridpts, coeff = self.get_heg(ngridpts, optimized, ngridpts0)

        nmodes = self.nm.nmodes
        if modes is None:
            if nmode == 1:
                ints = np.empty(nmodes, dtype=object)
                for i in range(nmodes):
                    vgrid = np.array([self.nm.potential_1mode(i,qi) for qi in gridpts[i]]) # this should be vectorized
                    vi = np.dot(coeff[i], np.dot(np.diag(vgrid), coeff[i].T))
                    ints[i] = vi * constants.AU_TO_INVCM
            elif nmode == 2:
                ints = np.empty((nmodes,nmodes), dtype=object)
                for i in range(nmodes):
                    Ci = coeff[i].T @ onemode_coeff[i]
                    for j in range(nmodes):
                        Cj = coeff[j].T @ onemode_coeff[j]
                        vgrid = np.array([[self.nm.potential_2mode(i,j,qi,qj) for qj in gridpts[j]] for qi in gridpts[i]]) # this should be vectorized
                        vij = np.einsum('gp,hq,gh,gr,hs->pqrs', Ci, Cj, vgrid, Ci, Cj, optimize=True)
                        ints[i, j] = vij * constants.AU_TO_INVCM

            elif nmode == 3:
                ints = np.empty((nmodes,nmodes,nmodes), dtype=object)
                for i in range(nmodes):
                    Ci = coeff[i].T @ onemode_coeff[i]
                    for j in range(nmodes):
                        Cj = coeff[j].T @ onemode_coeff[j]
                        for k in range(nmodes):
                            Ck = coeff[k].T @ onemode_coeff[k]
                            vgrid = np.array([[[self.nm.potential_3mode(i,j,k,qi,qj,qk) for qk in gridpts[k]] for qj in gridpts[j]] for qi in gridpts[i]])
                            vijk = np.einsum('gp,hq,fr,ghf,gs,ht,fu->pqrstu', Ci, Cj, Ck, vgrid, Ci, Cj, Ck, optimize=True)
                            ints[i, j, k] = vijk * constants.AU_TO_INVCM
        else:
            if nmode == 3:
                ints = np.empty((nmodes,nmodes,nmodes), dtype=object)
                for i in range(nmodes):
                    Ci = coeff[i].T @ onemode_coeff[i]
                    for j in range(i, nmodes):
                        Cj = coeff[j].T @ onemode_coeff[j]
                        for k in range(j, nmodes):
                            Ck = coeff[k].T @ onemode_coeff[k]
                            if (i, j, k) in modes:
                                vgrid = np.array([[[self.nm.potential_3mode(i,j,k,qi,qj,qk) for qk in gridpts[k]] for qj in gridpts[j]] for qi in gridpts[i]])
                                vijk = np.einsum('gp,hq,fr,ghf,gs,ht,fu->pqrstu', Ci, Cj, Ck, vgrid, Ci, Cj, Ck, optimize=True)
                                ints[i, j, k] = vijk * constants.AU_TO_INVCM
                                ints[j, i, k] = vijk * constants.AU_TO_INVCM
                                ints[i, k, j] = vijk * constants.AU_TO_INVCM
                                ints[k, i, j] = vijk * constants.AU_TO_INVCM
                                ints[j, k, i] = vijk * constants.AU_TO_INVCM
                                ints[k, j, i] = vijk * constants.AU_TO_INVCM
                            else:
                                ints[i, j, k] = np.zeros((ngridpts[i], ngridpts[j], ngridpts[k], ngridpts[i], ngridpts[j], ngridpts[k]))
                                ints[j, i, k] = np.zeros((ngridpts[i], ngridpts[j], ngridpts[k], ngridpts[i], ngridpts[j], ngridpts[k]))
                                ints[i, k, j] = np.zeros((ngridpts[i], ngridpts[j], ngridpts[k], ngridpts[i], ngridpts[j], ngridpts[k]))
                                ints[k, i, j] = np.zeros((ngridpts[i], ngridpts[j], ngridpts[k], ngridpts[i], ngridpts[j], ngridpts[k]))
                                ints[j, k, i] = np.zeros((ngridpts[i], ngridpts[j], ngridpts[k], ngridpts[i], ngridpts[j], ngridpts[k]))
                                ints[k, j, i] = np.zeros((ngridpts[i], ngridpts[j], ngridpts[k], ngridpts[i], ngridpts[j], ngridpts[k]))

        return ints

    def get_dipole_ints(self, nmode, ngridpts=None, optimized=False, ngridpts0=None, onemode_coeff = None, usePyPotDip = False):
        print("Calculating n-Mode dipole integrals for n =", nmode, flush = True) 
        if optimized is False:
            if ngridpts is None:
                ngridpts = [12]*self.nm.nmodes
            elif isinstance(ngridpts, (int, np.integer)):
                ngridpts = [ngridpts]*self.nm.nmodes

            gridpts, coeff = self.get_heg(ngridpts)

        else:
            if ngridpts is None:
                ngridpts = [8]*self.nm.nmodes
            elif isinstance(ngridpts, (int, np.integer)):
                ngridpts = [ngridpts]*self.nm.nmodes

            if ngridpts0 is None:
                ngridpts0 = [12]*self.nm.nmodes
            elif isinstance(ngridpts0, (int, np.integer)):
                ngridpts0 = [ngridpts0]*self.nm.nmodes

            gridpts, coeff = self.get_heg(ngridpts, optimized, ngridpts0)

        nmodes = self.nm.nmodes
        if nmode == 1:
            ints = np.empty((3, nmodes), dtype=object)
            if usePyPotDip:
                ints = np.empty((4, nmodes), dtype=object)
            for i in range(nmodes):
                vgrid = np.array([self.nm.dipole_1mode(i,qi) for qi in gridpts[i]]) # this should be vectorized
                Ci = coeff[i].T @ onemode_coeff[i]
                vi = np.einsum('gp,gx,gq->xpq', Ci, vgrid, Ci, optimize=True)
                ints[0, i] = vi[0]
                ints[1, i] = vi[1]
                ints[2, i] = vi[2]
                if usePyPotDip:
                    ints[3, i] = vi[3]
        elif nmode == 2:
            ints = np.empty((3, nmodes, nmodes), dtype=object)
            if usePyPotDip:
                ints = np.empty((4, nmodes, nmodes), dtype=object)
            for i in range(nmodes):
                Ci = coeff[i].T @ onemode_coeff[i]
                for j in range(nmodes):
                    Cj = coeff[j].T @ onemode_coeff[j]
                    vgrid = np.array([[self.nm.dipole_2mode(i,j,qi,qj) for qj in gridpts[j]] for qi in gridpts[i]]) # this should be vectorized
                    vij = np.einsum('gp,hq,ghx,gr,hs->xpqrs', Ci, Cj, vgrid, Ci, Cj, optimize=True)
                    ints[0, i, j] = vij[0]
                    ints[1, i, j] = vij[1]
                    ints[2, i, j] = vij[2]
                    if usePyPotDip:
                        ints[3, i, j] = vij[3]

        elif nmode == 3:
            ints = np.empty((3, nmodes, nmodes, nmodes), dtype=object)
            if usePyPotDip:
                ints = np.empty((4, nmodes, nmodes, nmodes), dtype=object)
            for i in range(nmodes):
                Ci = coeff[i].T @ onemode_coeff[i]
                for j in range(nmodes):
                    Cj = coeff[j].T @ onemode_coeff[j]
                    for k in range(nmodes):
                        Ck = coeff[k].T @ onemode_coeff[k]
                        vgrid = np.array([[[self.nm.dipole_3mode(i,j,k,qi,qj,qk) for qk in gridpts[k]] for qj in gridpts[j]] for qi in gridpts[i]])
                        vijk = np.zeros((ngridpts[i], ngridpts[j], ngridpts[k], ngridpts[i], ngridpts[j], ngridpts[k]))
                        vijk = np.einsum('gp,hq,fr,ghfx,gs,ht,fu->xpqrstu', Ci, Cj, Ck, vgrid, Ci, Cj, Ck, optimize=True)
                        ints[0, i, j, k] = vijk[0]
                        ints[1, i, j, k] = vijk[1]
                        ints[2, i, j, k] = vijk[2]
                        if usePyPotDip:
                            ints[3, i, j, k] = vijk[3]

        return ints

    def get_inv_inertia_ints(self, nmode, ngridpts=None, optimized=False, ngridpts0=None, onemode_coeff = None):
        print("Calculating n-Mode inverse inertia integrals for n =", nmode, flush = True) 
        if optimized is False:
            if ngridpts is None:
                ngridpts = [12]*self.nm.nmodes
            elif isinstance(ngridpts, (int, np.integer)):
                ngridpts = [ngridpts]*self.nm.nmodes

            gridpts, coeff = self.get_heg(ngridpts)

        else:
            if ngridpts is None:
                ngridpts = [8]*self.nm.nmodes
            elif isinstance(ngridpts, (int, np.integer)):
                ngridpts = [ngridpts]*self.nm.nmodes

            if ngridpts0 is None:
                ngridpts0 = [12]*self.nm.nmodes
            elif isinstance(ngridpts0, (int, np.integer)):
                ngridpts0 = [ngridpts0]*self.nm.nmodes

            gridpts, coeff = self.get_heg(ngridpts, optimized, ngridpts0)

        nmodes = self.nm.nmodes
        if nmode == 1:
            ints = np.empty((9, nmodes), dtype=object)
            for i in range(nmodes):
                vgrid = np.array([self.nm.inv_moment_of_inertia_1mode(i,qi) for qi in gridpts[i]]) # this should be vectorized
                Ci = coeff[i].T @ onemode_coeff[i]
                vi = np.einsum('gp,gx,gq->xpq', Ci, vgrid, Ci, optimize=True)
                ints[0, i] = vi[0]
                ints[1, i] = vi[1]
                ints[2, i] = vi[2]
                ints[3, i] = vi[3]
                ints[4, i] = vi[4]
                ints[5, i] = vi[5]
                ints[6, i] = vi[6]
                ints[7, i] = vi[7]
                ints[8, i] = vi[8]

        elif nmode == 2:
            ints = np.empty((3, nmodes, nmodes), dtype=object)
            for i in range(nmodes):
                Ci = coeff[i].T @ onemode_coeff[i]
                for j in range(nmodes):
                    Cj = coeff[j].T @ onemode_coeff[j]
                    vgrid = np.array([[self.nm.inv_moment_of_inertia_2mode(i,j,qi,qj) for qj in gridpts[j]] for qi in gridpts[i]]) # this should be vectorized
                    vij = np.einsum('gp,hq,ghx,gr,hs->xpqrs', Ci, Cj, vgrid, Ci, Cj, optimize=True)
                    ints[0, i, j] = vij[0]
                    ints[1, i, j] = vij[1]
                    ints[2, i, j] = vij[2]
                    ints[3, i, j] = vij[3]
                    ints[4, i, j] = vij[4]
                    ints[5, i, j] = vij[5]
                    ints[6, i, j] = vij[6]
                    ints[7, i, j] = vij[7]
                    ints[8, i, j] = vij[8]

        elif nmode == 3:
            ints = np.empty((3, nmodes, nmodes, nmodes), dtype=object)
            for i in range(nmodes):
                Ci = coeff[i].T @ onemode_coeff[i]
                for j in range(nmodes):
                    Cj = coeff[j].T @ onemode_coeff[j]
                    for k in range(nmodes):
                        Ck = coeff[k].T @ onemode_coeff[k]
                        vgrid = np.array([[[self.nm.inv_moment_of_inertia_3mode(i,j,k,qi,qj,qk) for qk in gridpts[k]] for qj in gridpts[j]] for qi in gridpts[i]])
                        vijk = np.zeros((ngridpts[i], ngridpts[j], ngridpts[k], ngridpts[i], ngridpts[j], ngridpts[k]))
                        vijk = np.einsum('gp,hq,fr,ghfx,gs,ht,fu->xpqrstu', Ci, Cj, Ck, vgrid, Ci, Cj, Ck, optimize=True)
                        ints[0, i, j, k] = vijk[0]
                        ints[1, i, j, k] = vijk[1]
                        ints[2, i, j, k] = vijk[2]
                        ints[3, i, j, k] = vijk[3]
                        ints[4, i, j, k] = vijk[4]
                        ints[5, i, j, k] = vijk[5]
                        ints[6, i, j, k] = vijk[6]
                        ints[7, i, j, k] = vijk[7]
                        ints[8, i, j, k] = vijk[8]

        return ints


    def get_heg(self, ngridpts, optimized=False, ngridpts0=None):
        """
        Calculate the HEG grid points and functions by diagonalizing the Q
        matrix in a basis of harmonic oscillator eigenfunctions.

        Parameters:
        ngridpts (int or list of int): The number of harmonic oscillator
            eigenfunctions (per mode) to use when diagonalizing Q.
        """
        self.ngridpts = ngridpts
        gridpts = list()
        coeff = list()
        for i in range(self.nm.nmodes):
            if optimized: nho = ngridpts0[i]
            else: nho = ngridpts[i]
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

        self.gridpts, self.coeff = gridpts, coeff
        return self.gridpts, self.coeff 

if __name__ == '__main__':
    from vstr.examples.potentials.h2o.h2o_pot import calc_h2o_pot
    def pot_cart(coords):
        if np.array(coords).ndim == 3:
            return calc_h2o_pot(coords, len(coords))
        else:
            return calc_h2o_pot([coords], 1)[0]

    x0 = np.array([
        [0, -0.757, 0.587],
        [0,  0.757, 0.587],
        [0,  0    , 0    ]]) 
    x0 *= constants.ANGSTROM_TO_AU
    mass_h = 1836.152697 # in au
    mass_o = 29148.94642

    mass = [mass_h, mass_h, mass_o]

    pymol = gto.M()
    pymol.atom = 'H 0 -.757, .587; H 0 .757 .587; O 0 0 0'
    pymol.basis = 'sto-3g'
    pymol.build()

    vmol = Molecule(pot_cart, x0.shape[0], mass, ngridpts = 2, Order = 2)
    vmol.IntsFile = './ints.h5'
    vmol.calc_dipole = True
    vmol.init_pypotential(pymol)
    vmol.init_pydipole(pymol, doPotDip = True)
    vmol.kernel(x0 = x0)
    vmol.SaveIntegrals()
    vmol.SaveDipoles()

    print(vmol.ints[1][0, 1])

    '''
    vmol_read = Molecule(pot_cart, x0.shape[0], mass, ngridpts=2, Order = 3)
    vmol_read.ReadInt = True
    vmol_read.ReadDip = True
    vmol_read.calc_dipole = True
    vmol_read.IntsFile = './ints.h5'
    vmol_read.kernel(x0 = x0)
    
    print(vmol_read.ints[0])
    print(vmol_read.dip_ints[0])
    '''
