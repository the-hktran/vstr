import numpy as np
from pyscf.lo import Boys

def RCenter(mLO, U):
    QNew = np.einsum('pk,nxk->nxp', U, mLO.nm_coeff, optimize = True)
    C = QNew.sum(axis = 1)
    return QNew, C #, np.einsum('np,nx->xp', C, mLO.coords, optimize = True)

def gen_g_hop(mLO, U):
    g = mLO.get_grad(U)

    K, P = g.shape
    h = np.zeros((K, P, K, P))

    QNew, C = mLO.RCenter(U)
    RR = mLO.coords @ mLO.coords.T

    h = 4 * np.einsum('pq,nm,nxk,nxp,myp,myl->kplq', np.eye(P), RR, mLO.nm_coeff, QNew, 2 * QNew, mLO.nm_coeff, optimize = True) + 4 * np.einsum('pq,nm,nxk,mp,nxl->kplq', np.eye(P), RR, mLO.nm_coeff, C, mLO.nm_coeff, optimize = True)

    h_diag = np.diag(h.reshape((K * P, K * P))).reshape(K, P)

    def h_op(x):
        return np.einsum('kplq,lq->kp', h, x, optimize = True)
    
    return g, h_op, h_diag
    
def get_grad(mLO, U = None):
    if U is None:
        U = np.eye(mLO.nm_coeff.shape[2])
    QNew, C = mLO.RCenter(U)
    RR = mLO.coords @ mLO.coords.T
    return 4 * np.einsum('nm,mp,nxp,nxk->kp', RR, C, QNew, mLO.nm_coeff, optimize = True)

def cost_function(mLO, U = None):
    if U is None:
        U = np.eye(mLO.nm_coeff.shape[2])
    QNew, C = mLO.RCenter(U)
    RC = np.einsum('np,nx->px', C, mLO.coords, optimize = True)
    return (RC * RC).sum()

class NMBoys(Boys):
    RCenter = RCenter
    gen_g_hop = gen_g_hop
    get_grad = get_grad
    cost_function = cost_function

    def __init__(self, mol, **kwargs):
        Boys.__init__(self, mol, mol.nm.nm_coeff, **kwargs)
        self.nm_coeff = mol.nm.nm_coeff
        self.coords = mol.x0

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

    vmol = Molecule(pot_cart, x0.shape[0], mass, ngridpts = 8, Order = 2)
    vmol.kernel()

    mlo = NMBoys(vmol)
    mlo.kernel()

