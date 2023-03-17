import numpy as np
from pyscf import gto, scf, hessian
from vstr.ff.force_field import PerturbCoord
from vstr.ff.normal_modes import CoordToAtom, AtomToCoord
from vstr.cpp_wrappers.vhci_jf.vhci_jf_functions import FConst

def MakeDipoleList(mu_raw):
    DipoleSurface = []
    for m in mu_raw:
        d = FConst(m[0], m[1], True)
        DipoleSurface.append(d)
    return DipoleSurface

def DipoleNorm(mf):
    return np.sqrt((mf.dip_moment()**2).sum())

def PerturbDipole(X0, Modes, Coords, dx, atom0, mol, Method = 'rhf'):
    X = PerturbCoord(X0, Modes, Coords, dx)
    atom = CoordToAtom(atom0, X)
    mol.atom = atom
    mol.build()
    new_mf = scf.RHF(mol)
    new_mf.kernel(verbose = 0)
    return DipoleNorm(new_mf)

def GetDipoleSurface(mf, Coords, Order = 1, dx = 1e-1):
    Dipole = []
    NCoord = Coords.shape[1]
    X0 = AtomToCoord(mf) # in Bohr
    atom0 = mf.mol._atom.copy()
    mu0 = DipoleNorm(mf)
    new_mol = mf.mol.copy()
    new_mol.unit = 'B'

    Dipole.append([(mu0, [])])

    mu1 = []
    mu2 = []
    mu3 = []
    mu4 = []

    if Order >= 1:
        for i in range(NCoord):
            mu_p = PerturbDipole(X0, [[i, 1]], Coords, dx, atom0, new_mol)
            mu_m = PerturbDipole(X0, [[i, -1]], Coords, dx, atom0, new_mol)

            mu_i = (mu_p - mu_m) / (2 * dx)
            mu1.append((mu_i, [i]))

            if Order >= 2:
                mu_ii = (mu_p + mu_m - 2 * mu0) / (dx * dx)
                mu2.append((mu_ii, [i, i]))

            if Order >= 3:
                mu_pp = PerturbDipole(X0, [[i, 2]], Coords, dx, atom0, new_mol)
                mu_mm = PerturbDipole(X0, [[i, -2]], Coords, dx, atom0, new_mol)
                mu_iii = (mu_pp - 2 * mu_p + 2 * mu_m - mu_mm) / (2 * dx * dx * dx)
                mu3.append((mu_iii, [i, i, i]))
                if Order >= 4:
                    mu_iiii = (mu_pp - 4 * mu_p + 6 * mu0 - 4 * mu_m + mu_mm) / (dx**4)
                    mu4.append((mu_iiii, [i, i, i, i]))

        Dipole.append(mu1)
    
    if Order >= 2:
        for i in range(NCoord):
            for j in range(i + 1, NCoord):
                mu_pp = PerturbDipole(X0, [[i, 1], [j, 1]], Coords, dx, atom0, new_mol)
                mu_pm = PerturbDipole(X0, [[i, 1], [j, -1]], Coords, dx, atom0, new_mol)
                mu_mp = PerturbDipole(X0, [[i, -1], [j, 1]], Coords, dx, atom0, new_mol)
                mu_mm = PerturbDipole(X0, [[i, -1], [j, -1]], Coords, dx, atom0, new_mol)
                mu_ij = (mu_pp + mu_mm - mu_pm - mu_mp) / (4 * dx * dx)
                mu2.append((mu_ij, [i, j]))
        Dipole.append(mu2)

    if Order >= 3:
        for i in range(NCoord):
            # iii done before
            for j in range(i + 1, NCoord):
                mu_ppp = PerturbDipole(X0, [[i, 2], [j, 1]], Coords, dx, atom0, new_mol)
                mu_mmp = PerturbDipole(X0, [[i, -2], [j, 1]], Coords, dx, atom0, new_mol)
                mu_00p = PerturbDipole(X0, [[j, 1]], Coords, dx, atom0, new_mol)
                mu_ppm = PerturbDipole(X0, [[i, 2], [j, -1]], Coords, dx, atom0, new_mol)
                mu_mmm = PerturbDipole(X0, [[i, -2], [j, -1]], Coords, dx, atom0, new_mol)
                mu_00m = PerturbDipole(X0, [[j, -1]], Coords, dx, atom0, new_mol)
                mu_iij = (mu_ppp + mu_mmp - 2 * mu_00p - mu_ppm - mu_mmm + 2 * mu_00m) / (8 * dx**3)
                mu3.append((mu_iij, [i, i, j]))
                for k in range(j + 1, NCoord):
                    mu_ppp = PerturbDipole(X0, [[i, 1], [j, 1], [k, 1]], Coords, dx, atom0, new_mol)
                    mu_ppm = PerturbDipole(X0, [[i, 1], [j, 1], [k, -1]], Coords, dx, atom0, new_mol)
                    mu_pmp = PerturbDipole(X0, [[i, 1], [j, -1], [k, 1]], Coords, dx, atom0, new_mol)
                    mu_pmm = PerturbDipole(X0, [[i, 1], [j, -1], [k, -1]], Coords, dx, atom0, new_mol)
                    mu_mpp = PerturbDipole(X0, [[i, -1], [j, 1], [k, 1]], Coords, dx, atom0, new_mol)
                    mu_mpm = PerturbDipole(X0, [[i, -1], [j, 1], [k, -1]], Coords, dx, atom0, new_mol)
                    mu_mmp = PerturbDipole(X0, [[i, -1], [j, -1], [k, 1]], Coords, dx, atom0, new_mol)
                    mu_mmm = PerturbDipole(X0, [[i, -1], [j, -1], [k, -1]], Coords, dx, atom0, new_mol)
                    mu_ijk = (mu_ppp - mu_ppm - mu_pmp - mu_mpp + mu_mmp + mu_mpm + mu_pmm - mu_mmm) / (8 * dx**3)
                    mu3.append((mu_ijk, [i, j, k]))
        Dipole.append(mu3)

    if Order >= 4:
        for i in range(NCoord):
        # iiii terms done in linear part
            for j in range(i + 1, NCoord):
                #iiij terms
                mu_pppp = PerturbDipole(X0, [[i, 3], [j, 1]], Coords, dx, atom0, new_mol)
                mu_p00p = PerturbDipole(X0, [[i, 1], [j, 1]], Coords, dx, atom0, new_mol)
                mu_mmmp = PerturbDipole(X0, [[i, -3], [j, 1]], Coords, dx, atom0, new_mol)
                mu_m00p = PerturbDipole(X0, [[i, -1], [j, 1]], Coords, dx, atom0, new_mol)
                mu_pppm = PerturbDipole(X0, [[i, 3], [j, -1]], Coords, dx, atom0, new_mol)
                mu_p00m = PerturbDipole(X0, [[i, 1], [j, -1]], Coords, dx, atom0, new_mol)
                mu_mmmm = PerturbDipole(X0, [[i, -3], [j, -1]], Coords, dx, atom0, new_mol)
                mu_m00m = PerturbDipole(X0, [[i, -1], [j, -1]], Coords, dx, atom0, new_mol)
                mu_iiij = (mu_pppp - 3 * mu_p00p + 3 * mu_m00p - mu_mmmp - mu_pppm + 3 * mu_p00m - 3 * mu_m00m + mu_mmmm) / (16 * dx**4)
                mu4.append((mu_iiij, [i, i, i, j]))
                for k in range(j + 1, NCoord):
                    #iijk terms
                    mu_pppp = PerturbDipole(X0, [[i, 2], [j, 1], [k, 1]], Coords, dx, atom0, new_mol)
                    mu_mmpp = PerturbDipole(X0, [[i, -2], [j, 1], [k, 1]], Coords, dx, atom0, new_mol)
                    mu_00pp = PerturbDipole(X0, [[j, 1], [k, 1]], Coords, dx, atom0, new_mol)
                    mu_ppmp = PerturbDipole(X0, [[i, 2], [j, -1], [k, 1]], Coords, dx, atom0, new_mol)
                    mu_mmmp = PerturbDipole(X0, [[i, -2], [j, -1], [k, 1]], Coords, dx, atom0, new_mol)
                    mu_00mp = PerturbDipole(X0, [[j, -1], [k, 1]], Coords, dx, atom0, new_mol)
                    mu_pppm = PerturbDipole(X0, [[i, 2], [j, 1], [k, -1]], Coords, dx, atom0, new_mol)
                    mu_mmpm = PerturbDipole(X0, [[i, -2], [j, 1], [k, -1]], Coords, dx, atom0, new_mol)
                    mu_00pm = PerturbDipole(X0, [[j, 1], [k, -1]], Coords, dx, atom0, new_mol)
                    mu_ppmm = PerturbDipole(X0, [[i, 2], [j, -1], [k, -1]], Coords, dx, atom0, new_mol)
                    mu_mmmm = PerturbDipole(X0, [[i, -2], [j, -1], [k, -1]], Coords, dx, atom0, new_mol)
                    mu_00mm = PerturbDipole(X0, [[j, -1], [k, -1]], Coords, dx, atom0, new_mol)
                    mu_iijk = (mu_pppp + mu_mmpp - 2 * mu_00pp - mu_ppmp - mu_mmmp + 2 * mu_00mp - mu_pppm - mu_mmpm + 2 * mu_00pm + mu_ppmm + mu_mmmm - 2 * mu_00mm) / (16 * dx**4)
                    mu4.append((mu_iijk, [i, i, j, k]))
                    for l in range(k + 1, NCoord):
                        #ijkl terms
                        mu_pppp = PerturbDipole(X0, [[i, 1], [j, 1], [k, 1], [l, 1]], Coords, dx, atom0, new_mol)
                        mu_mppp = PerturbDipole(X0, [[i, -1], [j, 1], [k, 1], [l, 1]], Coords, dx, atom0, new_mol)
                        mu_pmpp = PerturbDipole(X0, [[i, 1], [j, -1], [k, 1], [l, 1]], Coords, dx, atom0, new_mol)
                        mu_mmpp = PerturbDipole(X0, [[i, -1], [j, -1], [k, 1], [l, 1]], Coords, dx, atom0, new_mol)
                        mu_ppmp = PerturbDipole(X0, [[i, 1], [j, 1], [k, -1], [l, 1]], Coords, dx, atom0, new_mol)
                        mu_mpmp = PerturbDipole(X0, [[i, -1], [j, 1], [k, -1], [l, 1]], Coords, dx, atom0, new_mol)
                        mu_pmmp = PerturbDipole(X0, [[i, 1], [j, -1], [k, -1], [l, 1]], Coords, dx, atom0, new_mol)
                        mu_mmmp = PerturbDipole(X0, [[i, -1], [j, -1], [k, -1], [l, 1]], Coords, dx, atom0, new_mol)
                        mu_pppm = PerturbDipole(X0, [[i, 1], [j, 1], [k, 1], [l, -1]], Coords, dx, atom0, new_mol)
                        mu_mppm = PerturbDipole(X0, [[i, -1], [j, 1], [k, 1], [l, -1]], Coords, dx, atom0, new_mol)
                        mu_pmpm = PerturbDipole(X0, [[i, 1], [j, -1], [k, 1], [l, -1]], Coords, dx, atom0, new_mol)
                        mu_mmpm = PerturbDipole(X0, [[i, -1], [j, -1], [k, 1], [l, -1]], Coords, dx, atom0, new_mol)
                        mu_ppmm = PerturbDipole(X0, [[i, 1], [j, 1], [k, -1], [l, -1]], Coords, dx, atom0, new_mol)
                        mu_mpmm = PerturbDipole(X0, [[i, -1], [j, 1], [k, -1], [l, -1]], Coords, dx, atom0, new_mol)
                        mu_pmmm = PerturbDipole(X0, [[i, 1], [j, -1], [k, -1], [l, -1]], Coords, dx, atom0, new_mol)
                        mu_mmmm = PerturbDipole(X0, [[i, -1], [j, -1], [k, -1], [l, -1]], Coords, dx, atom0, new_mol)
                        mu_ijkl = (mu_pppp - mu_mppp - mu_pmpp + mu_mmpp - mu_ppmp + mu_mpmp + mu_pmmp - mu_mmmp - mu_pppm + mu_mppm + mu_pmpm - mu_mmpm + mu_ppmm - mu_mpmm - mu_pmmm + mu_mmmm) / (16 * dx**4)
                        mu4.append((mu_ijkl, [i, j, k, l]))
        Dipole.append(mu4)

    return Dipole

if __name__ == "__main__":
    from vstr.ff.normal_modes import GetNormalModes

    mol = gto.M()
    mol.atom = '''
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

    mu = GetDipoleSurface(mf, C, Order = 2)
    print(mu)
    
