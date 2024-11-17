import numpy as np
from pyscf.lo import Boys
from pyscf.soscf import ciah
from pyscf import lib
from pyscf.lib import logger
from pyscf import __config__
from pyscf.lo import orth, cholesky_mos

def kernel(localizer, mo_coeff=None, callback=None, verbose=None):
    from pyscf.tools import mo_mapping
    if mo_coeff is not None:
        localizer.mo_coeff = np.asarray(mo_coeff, order='C')
    if localizer.mo_coeff.shape[1] <= 1:
        return localizer.mo_coeff

    if localizer.verbose >= logger.WARN:
        localizer.check_sanity()
    localizer.dump_flags()

    cput0 = (logger.process_clock(), logger.perf_counter())
    log = logger.new_logger(localizer, verbose=verbose)

    if localizer.conv_tol_grad is None:
        conv_tol_grad = np.sqrt(localizer.conv_tol*.1)
        log.info('Set conv_tol_grad to %g', conv_tol_grad)
    else:
        conv_tol_grad = localizer.conv_tol_grad

    if mo_coeff is None:
        if getattr(localizer, 'mol', None) and localizer.mol.natm == 0:
            # For customized Hamiltonian
            u0 = localizer.get_init_guess('random')
        else:
            #u0 = localizer.get_init_guess(localizer.init_guess)
            u0 = np.eye(localizer.mo_coeff.shape[0])
    else:
        u0 = localizer.get_init_guess(None)

    rotaiter = ciah.rotate_orb_cc(localizer, u0, conv_tol_grad, verbose=log)
    u, g_orb, stat = next(rotaiter)
    cput1 = log.timer('initializing CIAH', *cput0)

    tot_kf = stat.tot_kf
    tot_hop = stat.tot_hop
    conv = False
    e_last = 0
    for imacro in range(localizer.max_cycle):
        norm_gorb = np.linalg.norm(g_orb)
        u0 = lib.dot(u0, u)
        e = localizer.cost_function(u0)
        e_last, de = e, e-e_last

        log.info('macro= %d  f(x)= %.14g  delta_f= %g  |g|= %g  %d KF %d Hx' %(
                 imacro+1, e, de, norm_gorb, stat.tot_kf+1, stat.tot_hop))
        cput1 = log.timer('cycle= %d'%(imacro+1), *cput1)

        if (norm_gorb < conv_tol_grad and abs(de) < localizer.conv_tol
                and stat.tot_hop < localizer.ah_max_cycle):
            conv = True

        if callable(callback):
            callback(locals())

        if conv:
            break

        u, g_orb, stat = rotaiter.send(u0)
        tot_kf += stat.tot_kf
        tot_hop += stat.tot_hop

    rotaiter.close()
    #log.info
    print('macro X = %d  f(x)= %.14g  |g|= %g  %d intor %d KF %d Hx' %(
             imacro+1, e, norm_gorb,
             (imacro+1)*2, tot_kf+imacro+1, tot_hop))
# Sort the localized orbitals, to make each localized orbitals as close as
# possible to the corresponding input orbitals
    sorted_idx = mo_mapping.mo_1to1map(u0)
    localizer.mo_coeff = lib.dot(localizer.mo_coeff, u0[:,sorted_idx])
    return localizer.mo_coeff

def RCenter(mLO, U):
    QNew = np.einsum('kp,nxk->nxp', U, mLO.nm_coeff, optimize = True)
    C = (QNew * QNew).sum(axis = 1)
    return QNew, C #, np.einsum('np,nx->xp', C, mLO.coords, optimize = True)

def gen_g_hop(mLO, U):
    print(U)
    g = mLO.get_grad(U)
    g = mLO.pack_uniq_var(g)

    K = mLO.nm_coeff.shape[0]
    P = mLO.nm_coeff.shape[1]
    h = np.zeros((K, P, K, P))

    QNew, C = mLO.RCenter(U)
    RR = mLO.coords @ mLO.coords.T

    h = 4 * np.einsum('pq,nm,nxk,nxp,myp,myl->kplq', np.eye(P), RR, mLO.nm_coeff, QNew, 2 * QNew, mLO.nm_coeff, optimize = True) + 4 * np.einsum('pq,nm,nxk,mp,nxl->kplq', np.eye(P), RR, mLO.nm_coeff, C, mLO.nm_coeff, optimize = True)
    h_diag = mLO.pack_uniq_var(np.diag(h.reshape(K * P, K * P)).reshape(K, P))

    def h_op(x):
        x = mLO.unpack_uniq_var(x)
        return mLO.pack_uniq_var(np.einsum('kplq,lq->kp', h, x, optimize = True))
    
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

    kernel = kernel

    def __init__(self, mol, **kwargs):
        mol.stdout = None
        mol.verbose = 0
        mol.natm = mol.natoms
        mol.mo_coeff = mol.nm.nm_coeff

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

    vmol = Molecule(pot_cart, x0.shape[0], mass, ngridpts = 4, Order = 1)
    vmol.ReadGeom = True
    vmol.ReadInt = True
    vmol.kernel(x0 = x0)
    #vmol.SaveIntegrals()
    mlo = NMBoys(vmol)

    '''
    # numerical derivative
    U = np.eye(3)
    G = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            dU = np.zeros((3, 3))
            dU[i, j] = 1e-2
            G[i, j] = (mlo.cost_function(U + dU) - mlo.cost_function(U - dU)) / 2e-2
    print(G)
    print(mlo.get_grad(U))
    # numerical hessian
    H = np.zeros((3, 3, 3, 3))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    dU1 = np.zeros((3, 3))
                    dU2 = np.zeros((3, 3))
                    dU1[i, j] = 1e-2
                    dU2[k, l] = 1e-2
                    H[i, j, k, l] = (mlo.cost_function(U + dU1 + dU2) - mlo.cost_function(U + dU1 - dU2) - mlo.cost_function(U - dU1 + dU2) + mlo.cost_function(U - dU1 - dU2)) / 4e-4
    print(H)
    '''
    print("Normal Modes")
    print(mlo.nm_coeff)
    lo_coeff = mlo.kernel()
    print("Local Modes")
    print(lo_coeff)
