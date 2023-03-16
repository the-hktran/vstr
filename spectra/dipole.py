import numpy as np
from pyscf import gto, scf, hessian
from vstr.ff.force_field import PerturbCoord
from vstr.ff.normal_modes import CoordToAtom, AtomToCoord
from vstr.cpp_wrappers.vhci_jf.vhci_jf_functions import FConst

def MakeDipoleList(mu_raw):
    DipoleSurface = []
    for m in mu_raw:
        d = FConst(m[0], m[1], False)
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

def GetDipoleSurface(mf, Coords, Order = 1, dx = 1e-4):
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

    if Order >= 1:
        for i in range(NCoord):
            mu_p = PerturbDipole(X0, [[i, 1]], Coords, dx, atom0, new_mol)
            mu_m = PerturbDipole(X0, [[i, -1]], Coords, dx, atom0, new_mol)

            mu_i = (mu_p - mu_m) / (2 * dx)
            mu1.append((mu_i, [i]))

            if Order >= 2:
                mu_ii = (mu_p + mu_m - 2 * mu0) / (dx * dx)
                mu2.append((mu_ii, [i, i]))
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
    
