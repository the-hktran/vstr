import numpy as np
from pyscf import gto, scf, hessian, cc
from pyscf.cc import ccsd_t_lambda_slow as ccsd_t_lambda
from pyscf.cc import ccsd_t_rdm_slow as ccsd_t_rdm
from vstr.ff.force_field import PerturbCoord
from vstr.ff.normal_modes import CoordToAtom, AtomToCoord
from vstr.cpp_wrappers.vhci_jf.vhci_jf_functions import FConst
from vstr.utils import constants

def MakeDipoleList(mu_raw):
    DipoleSurface = []
    for m in mu_raw:
        d = FConst(m[0] * np.sqrt(2**len(m[1])), m[1], True)
        DipoleSurface.append(d)
    return DipoleSurface

def DipoleNorm(mf):
    return np.sqrt((mf.dip_moment()**2).sum())

def EthansDipole(mf):
    ao_dip = mf.mol.intor_symmetric('int1e_irp')
    e1_dip = np.einsum('xij,ji->x', ao_dip, mf.make_rdm1())

    charges = mol.atom_charges()
    coords = mol.atom_coords()
    nucl_dip = np.einsum('i,ix->x', charges, coords)
    mol_dip = e1_dip #nucl_dip - e1_dip

    from pyscf.data import nist
    return mol_dip *2.541746 

def ScaleDipole(mu, Freq, Modes):
    n = len(Modes)
    ScaledMu = mu / ((2 * np.pi)**(n / 2) * constants.AMU_TO_ME**(n / 2) * constants.C_AU**(n / 2) * constants.BOHR_TO_CM**(n / 2))
    if Freq is not None:
        for i in Modes:
            ScaledMu = ScaledMu / np.sqrt(Freq[i])
    # Units are now scaled by unitless factor. Dipole derivative has the same units as dipole.
    return ScaledMu

def GetDipole(mol, new_mf, Method = 'rhf'):
    if Method == 'rhf':
        return new_mf.dip_moment()
    elif Method == 'ccsd':
        mcc = cc.CCSD(mf)
        mcc.kernel()
        P = mcc.make_rdm1(ao_repr=True)
        return new_mf.dip_moment(mol, P)
    elif Method == 'ccsd_t':
        mcc = cc.CCSD(mf)
        mcc.kernel()
        eris = mcc.ao2mo()
        conv, l1, l2 = ccsd_t_lambda.kernel(mcc, eris, mcc.t1, mcc.t2)
        P = ccsd_t_rdm.make_rdm1(mcc, mcc.t1, mcc.t2, l1, l2, eris, ao_repr = True)
        return new_new.dip_moment(mol, P)
    else:
        raise RuntimeError("Unrecognized method for dipole calculation")

def PerturbDipole(X0, Modes, Coords, dx, atom0, mol, Method = 'rhf'):
    X = PerturbCoord(X0, Modes, Coords, dx)
    atom = CoordToAtom(atom0, X)
    mol.atom = atom
    mol.build()
    new_mf = scf.RHF(mol)
    new_mf.kernel(verbose = 0)
    #return DipoleNorm(new_mf)
    return GetDipole(mol, new_mf, Method = Method)

def GetDipoleSurface(mf, Coords, Freq = None, Order = 1, dx = 1e-1, Method = 'rhf'):
    DipoleX = []
    DipoleY = []
    DipoleZ = []
    NCoord = Coords.shape[1]
    X0 = AtomToCoord(mf) # in Bohr
    atom0 = mf.mol._atom.copy()
    mu0 = GetDipole(mol, mf, Method = Method)
    new_mol = mf.mol.copy()
    new_mol.unit = 'B'

    DipoleX.append([(mu0[0], [])])
    DipoleY.append([(mu0[1], [])])
    DipoleZ.append([(mu0[2], [])])

    mu1X = []
    mu1Y = []
    mu1Z = []
    mu2X = []
    mu2Y = []
    mu2Z = []
    mu3X = []
    mu3Y = []
    mu3Z = []
    mu4X = []
    mu4Y = []
    mu4Z = []

    if Order >= 1:
        for i in range(NCoord):
            mu_p = PerturbDipole(X0, [[i, 1]], Coords, dx, atom0, new_mol, Method = Method)
            mu_m = PerturbDipole(X0, [[i, -1]], Coords, dx, atom0, new_mol, Method = Method)

            mu_i = (mu_p - mu_m) / (2 * dx)
            mu_i = ScaleDipole(mu_i, Freq, [i])
            mu1X.append((mu_i[0], [i]))
            mu1Y.append((mu_i[1], [i]))
            mu1Z.append((mu_i[2], [i]))

            if Order >= 2:
                mu_ii = (mu_p + mu_m - 2 * mu0) / (dx * dx)
                mu_ii = ScaleDipole(mu_ii, Freq, [i, i])
                mu2X.append((mu_ii[0], [i, i]))
                mu2Y.append((mu_ii[1], [i, i]))
                mu2Z.append((mu_ii[2], [i, i]))

            if Order >= 3:
                mu_pp = PerturbDipole(X0, [[i, 2]], Coords, dx, atom0, new_mol, Method = Method)
                mu_mm = PerturbDipole(X0, [[i, -2]], Coords, dx, atom0, new_mol, Method = Method)
                mu_iii = (mu_pp - 2 * mu_p + 2 * mu_m - mu_mm) / (2 * dx * dx * dx)
                mu_iii = ScaleDipole(mu_iii, Freq, [i, i, i])
                mu3X.append((mu_iii[0], [i, i, i]))
                mu3Y.append((mu_iii[1], [i, i, i]))
                mu3Z.append((mu_iii[2], [i, i, i]))
                if Order >= 4:
                    mu_iiii = (mu_pp - 4 * mu_p + 6 * mu0 - 4 * mu_m + mu_mm) / (dx**4)
                    mu_iiii = ScaleDipole(mu_iiii, Freq, [i, i, i, i])
                    mu4X.append((mu_iiii[0], [i, i, i, i]))
                    mu4Y.append((mu_iiii[1], [i, i, i, i]))
                    mu4Z.append((mu_iiii[2], [i, i, i, i]))

        DipoleX.append(mu1X)
        DipoleY.append(mu1Y)
        DipoleZ.append(mu1Z)
    
    if Order >= 2:
        for i in range(NCoord):
            for j in range(i + 1, NCoord):
                mu_pp = PerturbDipole(X0, [[i, 1], [j, 1]], Coords, dx, atom0, new_mol, Method = Method)
                mu_pm = PerturbDipole(X0, [[i, 1], [j, -1]], Coords, dx, atom0, new_mol, Method = Method)
                mu_mp = PerturbDipole(X0, [[i, -1], [j, 1]], Coords, dx, atom0, new_mol, Method = Method)
                mu_mm = PerturbDipole(X0, [[i, -1], [j, -1]], Coords, dx, atom0, new_mol, Method = Method)
                mu_ij = (mu_pp + mu_mm - mu_pm - mu_mp) / (4 * dx * dx)
                mu_ij = ScaleDipole(mu_ij, Freq, [i, j])
                mu2X.append((mu_ij[0], [i, j]))
                mu2Y.append((mu_ij[1], [i, j]))
                mu2Z.append((mu_ij[2], [i, j]))
        DipoleX.append(mu2X)
        DipoleY.append(mu2Y)
        DipoleZ.append(mu2Z)

    if Order >= 3:
        for i in range(NCoord):
            # iii done before
            for j in range(i + 1, NCoord):
                mu_ppp = PerturbDipole(X0, [[i, 2], [j, 1]], Coords, dx, atom0, new_mol, Method = Method)
                mu_mmp = PerturbDipole(X0, [[i, -2], [j, 1]], Coords, dx, atom0, new_mol, Method = Method)
                mu_00p = PerturbDipole(X0, [[j, 1]], Coords, dx, atom0, new_mol, Method = Method)
                mu_ppm = PerturbDipole(X0, [[i, 2], [j, -1]], Coords, dx, atom0, new_mol, Method = Method)
                mu_mmm = PerturbDipole(X0, [[i, -2], [j, -1]], Coords, dx, atom0, new_mol, Method = Method)
                mu_00m = PerturbDipole(X0, [[j, -1]], Coords, dx, atom0, new_mol, Method = Method)
                mu_iij = (mu_ppp + mu_mmp - 2 * mu_00p - mu_ppm - mu_mmm + 2 * mu_00m) / (8 * dx**3)
                mu_iij = ScaleDipole(mu_iij, Freq, [i, i, j])
                mu3X.append((mu_iij[0], [i, i, j]))
                mu3Y.append((mu_iij[1], [i, i, j]))
                mu3Z.append((mu_iij[2], [i, i, j]))
                for k in range(j + 1, NCoord):
                    mu_ppp = PerturbDipole(X0, [[i, 1], [j, 1], [k, 1]], Coords, dx, atom0, new_mol, Method = Method)
                    mu_ppm = PerturbDipole(X0, [[i, 1], [j, 1], [k, -1]], Coords, dx, atom0, new_mol, Method = Method)
                    mu_pmp = PerturbDipole(X0, [[i, 1], [j, -1], [k, 1]], Coords, dx, atom0, new_mol, Method = Method)
                    mu_pmm = PerturbDipole(X0, [[i, 1], [j, -1], [k, -1]], Coords, dx, atom0, new_mol, Method = Method)
                    mu_mpp = PerturbDipole(X0, [[i, -1], [j, 1], [k, 1]], Coords, dx, atom0, new_mol, Method = Method)
                    mu_mpm = PerturbDipole(X0, [[i, -1], [j, 1], [k, -1]], Coords, dx, atom0, new_mol, Method = Method)
                    mu_mmp = PerturbDipole(X0, [[i, -1], [j, -1], [k, 1]], Coords, dx, atom0, new_mol, Method = Method)
                    mu_mmm = PerturbDipole(X0, [[i, -1], [j, -1], [k, -1]], Coords, dx, atom0, new_mol, Method = Method)
                    mu_ijk = (mu_ppp - mu_ppm - mu_pmp - mu_mpp + mu_mmp + mu_mpm + mu_pmm - mu_mmm) / (8 * dx**3)
                    mu_ijk = ScaleDipole(mu_ijk, Freq, [i, j, k])
                    mu3X.append((mu_ijk[0], [i, j, k]))
                    mu3Y.append((mu_ijk[1], [i, j, k]))
                    mu3Z.append((mu_ijk[2], [i, j, k]))
        DipoleX.append(mu3X)
        DipoleY.append(mu3Y)
        DipoleZ.append(mu3Z)

    if Order >= 4:
        for i in range(NCoord):
        # iiii terms done in linear part
            for j in range(i + 1, NCoord):
                #iiij terms
                mu_pppp = PerturbDipole(X0, [[i, 3], [j, 1]], Coords, dx, atom0, new_mol, Method = Method)
                mu_p00p = PerturbDipole(X0, [[i, 1], [j, 1]], Coords, dx, atom0, new_mol, Method = Method)
                mu_mmmp = PerturbDipole(X0, [[i, -3], [j, 1]], Coords, dx, atom0, new_mol, Method = Method)
                mu_m00p = PerturbDipole(X0, [[i, -1], [j, 1]], Coords, dx, atom0, new_mol, Method = Method)
                mu_pppm = PerturbDipole(X0, [[i, 3], [j, -1]], Coords, dx, atom0, new_mol, Method = Method)
                mu_p00m = PerturbDipole(X0, [[i, 1], [j, -1]], Coords, dx, atom0, new_mol, Method = Method)
                mu_mmmm = PerturbDipole(X0, [[i, -3], [j, -1]], Coords, dx, atom0, new_mol, Method = Method)
                mu_m00m = PerturbDipole(X0, [[i, -1], [j, -1]], Coords, dx, atom0, new_mol, Method = Method)
                mu_iiij = (mu_pppp - 3 * mu_p00p + 3 * mu_m00p - mu_mmmp - mu_pppm + 3 * mu_p00m - 3 * mu_m00m + mu_mmmm) / (16 * dx**4)
                mu_iiij = ScaleDipole(mu_ij, Freq, [i, i, i, j])
                mu4X.append((mu_iiij[0], [i, i, i, j]))
                mu4Y.append((mu_iiij[1], [i, i, i, j]))
                mu4Z.append((mu_iiij[2], [i, i, i, j]))
                for k in range(j + 1, NCoord):
                    #iijk terms
                    mu_pppp = PerturbDipole(X0, [[i, 2], [j, 1], [k, 1]], Coords, dx, atom0, new_mol, Method = Method)
                    mu_mmpp = PerturbDipole(X0, [[i, -2], [j, 1], [k, 1]], Coords, dx, atom0, new_mol, Method = Method)
                    mu_00pp = PerturbDipole(X0, [[j, 1], [k, 1]], Coords, dx, atom0, new_mol, Method = Method)
                    mu_ppmp = PerturbDipole(X0, [[i, 2], [j, -1], [k, 1]], Coords, dx, atom0, new_mol, Method = Method)
                    mu_mmmp = PerturbDipole(X0, [[i, -2], [j, -1], [k, 1]], Coords, dx, atom0, new_mol, Method = Method)
                    mu_00mp = PerturbDipole(X0, [[j, -1], [k, 1]], Coords, dx, atom0, new_mol, Method = Method)
                    mu_pppm = PerturbDipole(X0, [[i, 2], [j, 1], [k, -1]], Coords, dx, atom0, new_mol, Method = Method)
                    mu_mmpm = PerturbDipole(X0, [[i, -2], [j, 1], [k, -1]], Coords, dx, atom0, new_mol, Method = Method)
                    mu_00pm = PerturbDipole(X0, [[j, 1], [k, -1]], Coords, dx, atom0, new_mol, Method = Method)
                    mu_ppmm = PerturbDipole(X0, [[i, 2], [j, -1], [k, -1]], Coords, dx, atom0, new_mol, Method = Method)
                    mu_mmmm = PerturbDipole(X0, [[i, -2], [j, -1], [k, -1]], Coords, dx, atom0, new_mol, Method = Method)
                    mu_00mm = PerturbDipole(X0, [[j, -1], [k, -1]], Coords, dx, atom0, new_mol, Method = Method)
                    mu_iijk = (mu_pppp + mu_mmpp - 2 * mu_00pp - mu_ppmp - mu_mmmp + 2 * mu_00mp - mu_pppm - mu_mmpm + 2 * mu_00pm + mu_ppmm + mu_mmmm - 2 * mu_00mm) / (16 * dx**4)
                    mu_iijk = ScaleDipole(mu_iijk, Freq, [i, i, j, k])
                    mu4X.append((mu_iijk[0], [i, i, j, k]))
                    mu4Y.append((mu_iijk[1], [i, i, j, k]))
                    mu4Z.append((mu_iijk[2], [i, i, j, k]))
                    for l in range(k + 1, NCoord):
                        #ijkl terms
                        mu_pppp = PerturbDipole(X0, [[i, 1], [j, 1], [k, 1], [l, 1]], Coords, dx, atom0, new_mol, Method = Method)
                        mu_mppp = PerturbDipole(X0, [[i, -1], [j, 1], [k, 1], [l, 1]], Coords, dx, atom0, new_mol, Method = Method)
                        mu_pmpp = PerturbDipole(X0, [[i, 1], [j, -1], [k, 1], [l, 1]], Coords, dx, atom0, new_mol, Method = Method)
                        mu_mmpp = PerturbDipole(X0, [[i, -1], [j, -1], [k, 1], [l, 1]], Coords, dx, atom0, new_mol, Method = Method)
                        mu_ppmp = PerturbDipole(X0, [[i, 1], [j, 1], [k, -1], [l, 1]], Coords, dx, atom0, new_mol, Method = Method)
                        mu_mpmp = PerturbDipole(X0, [[i, -1], [j, 1], [k, -1], [l, 1]], Coords, dx, atom0, new_mol, Method = Method)
                        mu_pmmp = PerturbDipole(X0, [[i, 1], [j, -1], [k, -1], [l, 1]], Coords, dx, atom0, new_mol, Method = Method)
                        mu_mmmp = PerturbDipole(X0, [[i, -1], [j, -1], [k, -1], [l, 1]], Coords, dx, atom0, new_mol, Method = Method)
                        mu_pppm = PerturbDipole(X0, [[i, 1], [j, 1], [k, 1], [l, -1]], Coords, dx, atom0, new_mol, Method = Method)
                        mu_mppm = PerturbDipole(X0, [[i, -1], [j, 1], [k, 1], [l, -1]], Coords, dx, atom0, new_mol, Method = Method)
                        mu_pmpm = PerturbDipole(X0, [[i, 1], [j, -1], [k, 1], [l, -1]], Coords, dx, atom0, new_mol, Method = Method)
                        mu_mmpm = PerturbDipole(X0, [[i, -1], [j, -1], [k, 1], [l, -1]], Coords, dx, atom0, new_mol, Method = Method)
                        mu_ppmm = PerturbDipole(X0, [[i, 1], [j, 1], [k, -1], [l, -1]], Coords, dx, atom0, new_mol, Method = Method)
                        mu_mpmm = PerturbDipole(X0, [[i, -1], [j, 1], [k, -1], [l, -1]], Coords, dx, atom0, new_mol, Method = Method)
                        mu_pmmm = PerturbDipole(X0, [[i, 1], [j, -1], [k, -1], [l, -1]], Coords, dx, atom0, new_mol, Method = Method)
                        mu_mmmm = PerturbDipole(X0, [[i, -1], [j, -1], [k, -1], [l, -1]], Coords, dx, atom0, new_mol, Method = Method)
                        mu_ijkl = (mu_pppp - mu_mppp - mu_pmpp + mu_mmpp - mu_ppmp + mu_mpmp + mu_pmmp - mu_mmmp - mu_pppm + mu_mppm + mu_pmpm - mu_mmpm + mu_ppmm - mu_mpmm - mu_pmmm + mu_mmmm) / (16 * dx**4)
                        mu_ijkl = ScaleDipole(mu_ijkl, Freq, [i, j, k, l])
                        mu4X.append((mu_ijkl[0], [i, j, k, l]))
                        mu4Y.append((mu_ijkl[1], [i, j, k, l]))
                        mu4Z.append((mu_ijkl[2], [i, j, k, l]))
        DipoleX.append(mu4X)
        DipoleY.append(mu4Y)
        DipoleZ.append(mu4Z)

    return [DipoleX, DipoleY, DipoleZ]

if __name__ == "__main__":
    from vstr.ff.normal_modes import GetNormalModes

    mol = gto.M()
    mol.atom = 'O 0 0 0; H 0 0 1; H 0 1 0'
    
    '''
    C     
    H  1 1.121896     
    O  1 1.212288 2 120.679502     
    C  1 1.514668 2 114.798218 3 -179.988745 0     
    H  4 1.101458 1 110.154420 2 -180.000000 0     
    H  4 1.105672 1 109.602564 2 -58.727388 0     
    H  4 1.105673 1 109.602602 2 58.729508 0   
    '''
    mol.basis='sto-3g'
    mol.build()
    mf = scf.RHF(mol)
    mf.kernel()

    w, C = GetNormalModes(mf, Method = 'rhf')

    mu = GetDipoleSurface(mf, C, dx = 1e-3, Order = 2, Method = 'rhf')
    print(mu)
    mu = GetDipoleSurface(mf, C, dx = 1e-3, Order = 2, Method = 'ccsd')
    print(mu)
    
