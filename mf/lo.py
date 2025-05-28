import numpy as np
import math
from pyscf.lo import Boys
from pyscf.soscf import ciah
from pyscf import lib
from pyscf.lib import logger
from pyscf import __config__
from pyscf.lo import orth, cholesky_mos

class NMOptimizer():
    def __init__(self, mol, maxiter = 1000, **kwargs):
        self.mol = mol
        self.nm_coeff = mol.nm.nm_coeff
        self.Q_loc = self.nm_coeff.copy()
        self.U = np.eye(mol.nm.nmodes)
        self.coords = mol.x0
        self.nm = mol.nm
        self.maxiter = maxiter
        self.tol = 1e-6

        self.nmodes = self.nm_coeff.shape[2]
        self.triu_idx = np.triu_indices(self.nmodes, k = 1)
        self.nidx = len(self.triu_idx[0])

        self.__dict__.update(kwargs)

    def kernel(self):
        K = self.get_K()
        cost_old = self.cost_function_i(0, 0, K)
        for i in range(self.maxiter):
            for j in range(self.nidx):
                p = self.triu_idx[0][j]
                q = self.triu_idx[1][j]
                uj1, uj2 = self.calc_angle(K, p, q)
                cost_new1 = self.cost_function_i(uj1, j, K)
                cost_new2 = self.cost_function_i(uj2, j, K)
                if cost_new1 > cost_new2:
                    uj = uj1
                    cost_new = cost_new1
                else:
                    uj = uj2
                    cost_new = cost_new2
                self.update_Q(uj, j)
                self.update_U(uj, j)
                K = self.get_K(self.Q_loc)
                print("Iteration: %d, Mode: %d, Cost: %f" % (i, j, cost_new), flush = True)
                #if abs(cost_new - cost_old) < self.tol:
                #    return self.Q_loc
                #cost_old = cost_new
            if abs(cost_new - cost_old) < self.tol:
                break
            cost_old = cost_new
            if i == self.maxiter - 1:
                print("Warning: Maximum number of iterations reached")

        # Reorder local modes
        S = np.zeros((self.nmodes, self.nmodes))
        np.fill_diagonal(S, self.mol.Frequencies)
        SLoc = self.U.T @ S @ self.U
        idx = np.argsort(SLoc.diagonal())
        self.U = self.U[:, idx]
        self.Q_loc = self.Q_loc[:, :, idx]
        self.Frequencies = SLoc.diagonal()[idx]
        self.mol.Frequencies = SLoc.diagonal()[idx]
        self.mol.nm.nm_coeff = self.Q_loc
        return self.mol

    def cost_function(self, U):
        pass

    def cost_function_i(self, u, i, K):
        U = np.eye(self.nmodes) 
        U[self.triu_idx[0][i], self.triu_idx[0][i]] = np.cos(u)
        U[self.triu_idx[1][i], self.triu_idx[1][i]] = np.cos(u)
        U[self.triu_idx[0][i], self.triu_idx[1][i]] = -np.sin(u)
        U[self.triu_idx[1][i], self.triu_idx[0][i]] = np.sin(u)
        return np.einsum('sr,tr,ur,vr,stuv->', U, U, U, U, K, optimize = True)

    def get_K(self, QLoc = None):
        pass

    def update_Q(self, u, i):
        Ui = np.eye(self.nmodes)
        Ui[self.triu_idx[0][i], self.triu_idx[0][i]] = np.cos(u)
        Ui[self.triu_idx[1][i], self.triu_idx[1][i]] = np.cos(u)
        Ui[self.triu_idx[0][i], self.triu_idx[1][i]] = -np.sin(u)
        Ui[self.triu_idx[1][i], self.triu_idx[0][i]] = np.sin(u)
        self.Q_loc = np.einsum('kp,nxk->nxp', Ui, self.Q_loc, optimize = True)

    def update_U(self, u, i):
        Ui = np.eye(self.nmodes)
        Ui[self.triu_idx[0][i], self.triu_idx[0][i]] = np.cos(u)
        Ui[self.triu_idx[1][i], self.triu_idx[1][i]] = np.cos(u)
        Ui[self.triu_idx[0][i], self.triu_idx[1][i]] = -np.sin(u)
        Ui[self.triu_idx[1][i], self.triu_idx[0][i]] = np.sin(u)
        self.U = self.U @ Ui

    def calc_angle(self, K, p, q):
        A = K[p, q, p, q] - (K[p, p, p, p] + K[q, q, q, q] - 2 * K[p, p, q, q]) / 4
        B = K[p, p, p, q] - K[q, q, q, p]
        
        a1 = math.asin( B / np.sqrt(A**2 + B**2)) / 4
        a2 = math.acos(-A / np.sqrt(A**2 + B**2)) / 4

        assert (abs(a1) < np.pi / 4) or (abs(a2) < np.pi / 4)
        return a1, a2


def RCenter(mLO, U):
    QNew = np.einsum('kp,nxk->nxp', U, mLO.nm_coeff, optimize = True)
    C = (QNew * QNew).sum(axis = 1)
    return QNew, C #, np.einsum('np,nx->xp', C, mLO.coords, optimize = True)

def ContractU(mLO, Us):
    U = np.eye(mLO.nm_coeff.shape[2])
    for i in range(mLO.nidx):
        u = np.eye(mLO.nm_coeff.shape[2])
        u[mLO.triu_idx[0][i], mLO.triu_idx[0][i]] = np.cos(Us[i])
        u[mLO.triu_idx[1][i], mLO.triu_idx[1][i]] = np.cos(Us[i])
        u[mLO.triu_idx[0][i], mLO.triu_idx[1][i]] = np.sin(Us[i])
        u[mLO.triu_idx[1][i], mLO.triu_idx[0][i]] = -np.sin(Us[i])
        U = U @ u
    return U

def cost_function_boys(mLO, U = None):
    if U is None:
        U = np.eye(mLO.nm_coeff.shape[2])
    QNew, C = mLO.RCenter(U)
    RC = np.einsum('np,nx->px', C, mLO.coords, optimize = True)
    return (RC * RC).sum()

def get_K_boys(mLO, QLoc = None):
    if QLoc is None:
        QLoc = mLO.Q_loc
    K1 = np.einsum('kx,kys,kyt->xst', mLO.coords, QLoc, QLoc, optimize = True)
    return np.einsum('xst,xuv->stuv', K1, K1, optimize = True)

def get_dist_matrix_boys(mLO):
    QNew, C = mLO.RCenter(mLO.U)
    RC = np.einsum('np,nx->px', C, mLO.coords, optimize = True)
    M = np.zeros((mLO.nmodes, mLO.nmodes))
    for i in range(mLO.nmodes):
        for j in range(mLO.nmodes):
            M[i, j] = np.linalg.norm(RC[i] - RC[j])
    return M

class NMBoys(NMOptimizer):
    RCenter = RCenter
    cost_function = cost_function_boys
    get_K = get_K_boys
    get_dist_matrix = get_dist_matrix_boys


if __name__ == '__main__':
    from vstr.examples.potentials.h2o_n.h2o_n_pot import calc_h2o_n_pot
    from vstr.utils import constants 
    from vstr.nmode.mol import Molecule

    def pot_cart(coords):
        return calc_h2o_n_pot(coords)[0]

    x0 = np.array([
        [0, -0.757, 0.587],
        [0,  0.757, 0.587],
        [0,  0    , 0    ]]) 
    x0 *= constants.ANGSTROM_TO_AU
    mass_h = 1836.152697 # in au
    mass_o = 29148.94642

    mass = [mass_h, mass_h, mass_o]

    vmol = Molecule(pot_cart, x0.shape[0], mass, ngridpts = 4, Order = 1)
    vmol.ReadGeom = True
    vmol.ReadInt = True
    vmol.kernel(x0 = x0)
    #vmol.SaveIntegrals()
    mlo = NMBoys(vmol)

    print("Normal Modes")
    print(mlo.nm_coeff)
    mol_lo = mlo.kernel()
    print("Local Modes")
    print(mlo.Q_loc.reshape(9, 3))
    print(np.einsum('nxk,kp->nxp', mlo.nm_coeff, mlo.U).reshape(9,3))

    S = np.zeros((3, 3))
    S[0, 0] = vmol.Frequencies[0]
    S[1, 1] = vmol.Frequencies[1]
    S[2, 2] = vmol.Frequencies[2]
    print(mlo.U.T @ S @ mlo.U)

