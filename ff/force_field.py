import numpy as np
from pyscf import gto, scf, hessian
from vstr.ff.normal_modes import AtomToCoord, CoordToAtom, GetHessian
from vstr.utils import constants

'''
Gets the force field constants from numerical derivatives of the hessian in the provided cooordinates.
The FF constants are stored in 
Cubic = [(Vijk, [i, j, k]), ...]
...
V = [Cubic, Quartic, ..]
'''

def PerturbCoord(X0, Modes, Coords, dx):
    X = X0.copy()
    for M in Modes:
        X += M[1] * Coords[:, M[0]] * dx
    return X

def CoordToHessian(X, atom0, mol, Method = 'rhf'):
    atom = CoordToAtom(atom0, X)
    mol.atom = atom
    mol.build()
    new_mf = scf.RHF(mol)
    new_mf.kernel(verbose = 0)
    return GetHessian(new_mf, Method = Method, MassWeighted = False)

def PerturbHessian(X0, Modes, Coords, dx, atom0, mol, Method = 'rhf'):
    X = PerturbCoord(X0, Modes, Coords, dx)
    H = CoordToHessian(X, atom0, mol, Method = Method)
    return Coords.T @ H @ Coords

'''
Scales force field constant to wavenumbers: Assumes energy in hartrees, mass weighted normal modes in bohr * sqrt(amu)
and frequencies in cm-1
'''
def ScaleFC(FC, Freqs, Modes):
    n = len(Modes)
    ScaledFC = FC / (constants.AMU_TO_ME**(n / 2) * constants.C_AU**(n / 2 + 1) * constants.BOHR_TO_CM**(n / 2 + 1))
    for i in Modes:
        ScaledFC = ScaledFC / np.sqrt(Freqs[i])
    return ScaledFC

def GetFF(mf, Coords, Freqs, Order = 4, Method = 'rhf', dx = 1e-2, tol = 1e-4):
    V = []
    NCoord = Coords.shape[1]
    X0 = AtomToCoord(mf) # in Bohr
    atom0 = mf.mol._atom.copy()
    H0 = GetHessian(mf, Method = Method, MassWeighted = True)
    H0 = Coords.T @ H0 @ Coords
    new_mol = mf.mol.copy()
    new_mol.unit = 'B'
    
    V3 = []
    V4 = []
    V5 = []
    V6 = []

    if Order >= 3:
        for i in range(NCoord):
            Hpi = PerturbHessian(X0, [[i, 1]], Coords, dx, atom0, new_mol, Method = Method)
            Hmi = PerturbHessian(X0, [[i, -1]], Coords, dx, atom0, new_mol, Method = Method)

            dHdxi = (Hpi - Hmi) / (2 * dx)

            for j in range(i, NCoord):
                for k in range(j, NCoord):
                    Hijk = ScaleFC(dHdxi[j, k], Freqs, [i, j, k])
                    if abs(Hijk) > tol:
                        V3.append((Hijk, [i, j, k]))

            if Order >= 4:
                # Can do iijk here as well
                d2Hdxidxi = (Hpi + Hmi - 2 * H0) / (dx * dx)
                for j in range(i, NCoord):
                    for k in range(j, NCoord):
                        Hiijk = ScaleFC(d2Hdxidxi[j, k], Freqs, [i, i, j, k])
                        if abs(Hiijk) > tol:
                            V4.append((Hiijk, [i, i, j, k]))
        V.append(V3)

    if Order >= 4:
        for i in range(NCoord):
            # iijk terms are handled in the cubic part
            for j in range(i + 1, NCoord):
                Hpipj = PerturbHessian(X0, [[i, 1], [j, 1]], Coords, dx, atom0, new_mol, Method = Method)
                Hmimj = PerturbHessian(X0, [[i, -1], [j, -1]], Coords, dx, atom0, new_mol, Method = Method)
                Hpimj = PerturbHessian(X0, [[i, 1], [j, -1]], Coords, dx, atom0, new_mol, Method = Method)
                Hmipj = PerturbHessian(X0, [[i, -1], [j, 1]], Coords, dx, atom0, new_mol, Method = Method)
                
                Hij = (Hpipj + Hmimj - Hpimj - Hmipj) / (4 * dx * dx)
                
                for k in range(j, NCoord):
                    for l in range(k, NCoord):
                        Hijkl = ScaleFC(Hij[k, l], Freqs, [i, j, k, l])
                        if abs(Hijkl) > tol:
                            V4.append((Hijkl, [i, j, k, l]))
        V.append(V4)

    if Order >= 5:
        for i in range(NCoord):
            #iiijk terms
            Hpipipi = PerturbHessian(X0, [[i, 3]], Coords, dx, atom0, new_mol, Method = Method)
            Hpi = PerturbHessian(X0, [[i, 1]], Coords, dx, atom0, new_mol, Method = Method)
            Hmi = PerturbHessian(X0, [[i, -1]], Coords, dx, atom0, new_mol, Method = Method)
            Hmimimi = PerturbHessian(X0, [[i, -3]], Coords, dx, atom0, new_mol, Method = Method)
            
            Hiii = (Hpipipi - 3 * Hpi + 3 * Hmi - Hmimimi) / (8 * dx * dx * dx)
            for j in range(i, NCoord):
                for k in range(j, NCoord):
                    Hiiijk = ScaleFC(Hiii[j, k], Freqs, [i, i, i, j, k])
                    if abs(Hiiikl) > tol:
                        V5.append((Hiiikl, [i, i, i, j, k]))

            # Handle sextic iiiijk terms here
            if Order >= 6:
                Hpipi = PerturbHessian(X0, [[i, 2]], Coords, dx, atom0, new_mol, Method = Method)
                Hmimi = PerturbHessian(X0, [[i, -2]], Coords, dx, atom0, new_mol, Method = Method)

                Hiiii = (Hpipi - 4 * Hpi + 6 * H0 - 4 * Hmi + Hmimi) / (dx * dx * dx * dx)
                for j in range(i, NCoord):
                    for k in range(j, NCoord):
                        Hiiiijk = ScaleFC(Hiiii[i, j], Freqs, [i, i, i, i, j, k])
                        if abs(Hiiiijk) > tol:
                            V6.append((Hiiiijk, [i, i, i, i, j, k]))

            for j in range(i + 1, NCoord):
                # iijkl elements.
                Hpipipj = PerturbHessian(X0, [[i, 2], [j, 1]], Coords, dx, atom0, new_mol, Method = Method)
                Hmimipj = PerturbHessian(X0, [[i, -2], [j, 1]], Coords, dx, atom0, new_mol, Method = Method)
                Hpj = PerturbHessian(X0, [[j, 1]], Coords, dx, atom0, new_mol, Method = Method)
                Hpipimj = PerturbHessian(X0, [[i, 2], [j, -1]], Coords, dx, atom0, new_mol, Method = Method)
                Hmimimj = PerturbHessian(X0, [[i, -2], [j, -1]], Coords, dx, atom0, new_mol, Method = Method)
                Hmj = PerturbHessian(X0, [[j, -1]], Coords, dx, atom0, new_mol, Method = Method)
                
                Hiij = (Hpipipj + Hmimipj - 2 * Hpj - Hpipimj - Hmimimj + 2 * Hmj) / (8 * dx * dx * dx)
                for k in range(j, NCoord):
                    for l in range(k, NCoord):
                        Hiijkl = ScaleFC(Hiij[k, l], Freqs, [i, i, j, k, l])
                        if abs(Hiijkl) > tol:
                            V5.append((Hiijkl, [i, i, j, k, l]))

                for k in range(j + 1, NCoord):
                    Hpipjpk = PerturbHessian(X0, [[i, 1], [j, 1], [k, 1]], Coords, dx, atom0, new_mol, Method = Method)
                    Hpipjmk = PerturbHessian(X0, [[i, 1], [j, 1], [k, -1]], Coords, dx, atom0, new_mol, Method = Method)
                    Hpimjpk = PerturbHessian(X0, [[i, 1], [j, -1], [k, 1]], Coords, dx, atom0, new_mol, Method = Method)
                    Hmipjpk = PerturbHessian(X0, [[i, -1], [j, 1], [k, 1]], Coords, dx, atom0, new_mol, Method = Method)
                    Hpimjmk = PerturbHessian(X0, [[i, 1], [j, -1], [k, -1]], Coords, dx, atom0, new_mol, Method = Method)
                    Hmipjmk = PerturbHessian(X0, [[i, -1], [j, 1], [k, -1]], Coords, dx, atom0, new_mol, Method = Method)
                    Hmimjpk = PerturbHessian(X0, [[i, -1], [j, -1], [k, 1]], Coords, dx, atom0, new_mol, Method = Method)
                    Hmimjmk = PerturbHessian(X0, [[i, -1], [j, -1], [k, -1]], Coords, dx, atom0, new_mol, Method = Method)

                    Hijk = (Hpipjpk - Hpipjmk - Hpimjpk - Hmipjpk + Hmimjpk + Hmipjmk + Hpimjmk - Hmimjmk) / (8 * dx * dx * dx)
                    for l in range(k, NCoord):
                        for p in range(l, NCoord):
                            Hijklp = ScaleFC(Hijklp, Freqs, [i, j, k, l, p])
                            if abs(Hijklp) > tol:
                                V5.append((Hijklp, [i, j, k, l, p]))
        V.append(V5)
    
    if Order >= 6:
        for i in range(NCoord):
            # iiiijk terms handled in quintic derivatives
            for j in range(i + 1, NCoord):
                #iiijkl terms
                Hpipipipj = PerturbHessian(X0, [[i, 3], [j, 1]], Coords, dx, atom0, new_mol, Method = Method)
                Hpipj = PerturbHessian(X0, [[i, 1], [j, 1]], Coords, dx, atom0, new_mol, Method = Method)
                Hmipj = PerturbHessian(X0, [[i, -1], [j, 1]], Coords, dx, atom0, new_mol, Method = Method)
                Hmimimipj = PerturbHessian(X0, [[i, -3], [j, 1]], Coords, dx, atom0, new_mol, Method = Method)
                Hpipipimj = PerturbHessian(X0, [[i, 3], [j, -1]], Coords, dx, atom0, new_mol, Method = Method)
                Hpimj = PerturbHessian(X0, [[i, 1], [j, -1]], Coords, dx, atom0, new_mol, Method = Method)
                Hmimj = PerturbHessian(X0, [[i, -1], [j, -1]], Coords, dx, atom0, new_mol, Method = Method)
                Hmimimimj = PerturbHessian(X0, [[i, -3], [j, -1]], Coords, dx, atom0, new_mol, Method = Method)
                
                Hiiij = (Hpipipipj - 3 * Hpipj + 3 * Hmipj - Hmimimipj - Hpipipimj + 3 * Hpimj - 3 * Hmipj + Hmimimimj) / (16 * dx**4)
                for k in range(j, NCoord):
                    for l in range(k, NCoord):
                        Hiiijkl = ScaleFC(Hiiijkl, Freqs, [i, i, i, j, k, l])
                        if abs(Hiiijkl) > tol:
                            V6.append((Hiiijkl, [i, i, i, j, k, l]))

                for k in range(j + 1, NCoord):
                    Hpipipjpk = PerturbHessian(X0, [[i, 2], [j, 1], [k, 1]], Coords, dx, atom0, new_mol, Method = Method)
                    Hmimipjpk = PerturbHessian(X0, [[i, -2], [j, 1], [k, 1]], Coords, dx, atom0, new_mol, Method = Method)
                    Hpjpk = PerturbHessian(X0, [[j, 1], [k, 1]], Coords, dx, atom0, new_mol, Method = Method)
                    Hpipimjpk = PerturbHessian(X0, [[i, 2], [j, -1], [k, 1]], Coords, dx, atom0, new_mol, Method = Method)
                    Hmimimjpk = PerturbHessian(X0, [[i, -2], [j, -1], [k, 1]], Coords, dx, atom0, new_mol, Method = Method)
                    Hmjpk = PerturbHessian(X0, [[j, -1], [k, 1]], Coords, dx, atom0, new_mol, Method = Method)
                    Hpipipjmk = PerturbHessian(X0, [[i, 2], [j, 1], [k, -1]], Coords, dx, atom0, new_mol, Method = Method)
                    Hmimipjmk = PerturbHessian(X0, [[i, -2], [j, 1], [k, -1]], Coords, dx, atom0, new_mol, Method = Method)
                    Hpjmk = PerturbHessian(X0, [[j, 1], [k, -1]], Coords, dx, atom0, new_mol, Method = Method)
                    Hpipimjmk = PerturbHessian(X0, [[i, 2], [j, -1], [k, -1]], Coords, dx, atom0, new_mol, Method = Method)
                    Hmimimjmk = PerturbHessian(X0, [[i, -2], [j, -1], [k, -1]], Coords, dx, atom0, new_mol, Method = Method)
                    Hmjmk = PerturbHessian(X0, [[j, -1], [k, -1]], Coords, dx, atom0, new_mol, Method = Method)

                    Hiijk = (Hpipipjpk + Hmimipjpk - 2 * Hpjpk - Hpipimjpk - Hmimimjpk + 2 * Hmjpk - Hpipipjmk - Hmimipjmk + 2 * Hpjmk + Hpipimjmk + Hmimimjmk - 2 * Hmjmk) / (16 * dx**4)
                    for l in range(k, NCoord):
                        for p in range(l, NCoord):
                            Hiijklp = ScaleFC(Hiijk[l, p], Freqs, [i, i, j, k, l, p])
                            if abs(Hiijklp) > tol:
                                V6.append((Hiijklp, [i, i, j, k, l, p]))

                    for l in range(k + 1, NCoord):
                        Hpppp = PerturbHessian(X0, [[i, 1], [j, 1], [k, 1], [l, 1]], Coords, dx, atom0, new_mol, Method = Method)
                        Hmppp = PerturbHessian(X0, [[i, -1], [j, 1], [k, 1], [l, 1]], Coords, dx, atom0, new_mol, Method = Method)
                        Hpmpp = PerturbHessian(X0, [[i, 1], [j, -1], [k, 1], [l, 1]], Coords, dx, atom0, new_mol, Method = Method)
                        Hmmpp = PerturbHessian(X0, [[i, -1], [j, -1], [k, 1], [l, 1]], Coords, dx, atom0, new_mol, Method = Method)
                        Hppmp = PerturbHessian(X0, [[i, 1], [j, 1], [k, -1], [l, 1]], Coords, dx, atom0, new_mol, Method = Method)
                        Hmpmp = PerturbHessian(X0, [[i, -1], [j, 1], [k, -1], [l, 1]], Coords, dx, atom0, new_mol, Method = Method)
                        Hpmmp = PerturbHessian(X0, [[i, 1], [j, -1], [k, -1], [l, 1]], Coords, dx, atom0, new_mol, Method = Method)
                        Hmmmp = PerturbHessian(X0, [[i, -1], [j, -1], [k, -1], [l, 1]], Coords, dx, atom0, new_mol, Method = Method)
                        Hpppm = PerturbHessian(X0, [[i, 1], [j, 1], [k, 1], [l, -1]], Coords, dx, atom0, new_mol, Method = Method)
                        Hmppm = PerturbHessian(X0, [[i, -1], [j, 1], [k, 1], [l, -1]], Coords, dx, atom0, new_mol, Method = Method)
                        Hpmpm = PerturbHessian(X0, [[i, 1], [j, -1], [k, 1], [l, -1]], Coords, dx, atom0, new_mol, Method = Method)
                        Hmmpm = PerturbHessian(X0, [[i, -1], [j, -1], [k, 1], [l, -1]], Coords, dx, atom0, new_mol, Method = Method)
                        Hppmm = PerturbHessian(X0, [[i, 1], [j, 1], [k, -1], [l, -1]], Coords, dx, atom0, new_mol, Method = Method)
                        Hmpmm = PerturbHessian(X0, [[i, -1], [j, 1], [k, -1], [l, -1]], Coords, dx, atom0, new_mol, Method = Method)
                        Hpmmm = PerturbHessian(X0, [[i, 1], [j, -1], [k, -1], [l, -1]], Coords, dx, atom0, new_mol, Method = Method)
                        Hmmmm = PerturbHessian(X0, [[i, -1], [j, -1], [k, -1], [l, -1]], Coords, dx, atom0, new_mol, Method = Method)

                        Hijkl = (Hpppp - Hmppp - Hpmpp + Hmmpp - Hppmp + Hmpmp + Hpmmp - Hmmmp - Hpppm + Hmppm + Hpmpm - Hmmpm + Hppmm - Hmpmm - Hpmmm + Hmmmm) / (16 * dx**4)
                        for p in range(l, NCoord):
                            for q in range(p, NCoord):
                                Hijklpq = ScaleFC(Hijkl[p, q], Freqs, [i, j, k, l, p, q])
                                if abs(Hijklpq) > tol:
                                    V6.append((Hijklpq, [i, j, k, l, p, q]))
        V.append(V6)
    
    return V

if __name__ == "__main__":
    from vstr.ff.normal_modes import GetNormalModes

    mol = gto.M()
    mol.fromfile("h2o.xyz")
    mol.basis='sto-3g'
    mol.build()
    mf = scf.RHF(mol)
    mf.kernel()

    w, C = GetNormalModes(mf)
    print(w, C)
    V = GetFF(mf, C, w)
    print(V)
