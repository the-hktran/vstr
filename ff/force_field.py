import numpy as np
from pyscf import gto, scf, hessian
from vstr.ff.normal_modes import AtomToCoord, CoordToAtom, GetHessian

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

def GetFF(mf, Coords, Order = 4, Method = 'rhf', dx = 1e-2):
    V = []
    NCoord = Coords.shape[1]
    X0 = AtomToCoord(mf) # in Bohr
    atom0 = mf.mol._atom.copy()
    H0 = GetHessian(mf, Method = Method, MassWeighted = False)
    new_mol = mf.mol.copy()
    new_mol.unit = 'B'
    
    V3 = []
    V4 = []
    if Order >= 3:
        for i in range(NCoord):
            X0pi = X0 + Coords[:, i] * dx
            X0mp = X0 - Coords[:, i] * dx

            atompi = CoordToAtom(atom0, X0pi)
            atommi = CoordToAtom(atom0, X0mi)

            new_mol.atom = atompi
            new_mol.build()
            new_mf = scf.RHF(new_mol)
            new_mf.kernel()
            Hpi = GetHessian(new_mf, Method = Method, MassWeighted = False)

            new_mol.atom = atommi
            new_mol.build()
            new_mf = scf.RHF(new_mol)
            new_mf.kernel()
            Hmi = GetHessian(new_mf, Method = Method, MassWeighted = False)

            dHdxi = (Hpi - Hmi) / (2 * dx)

            for j in range(i, NCoord):
                for k in range(j, NCoord):
                    V3.append((dHdxi[j, k], [i, j, k]))

            if Order >= 4:
                # Can do iijk here as well
                d2Hdxidxi = (Hpi + Hmi - 2 * H0) / (dx * dx)
                for j in range(i, NCoord):
                    for k in range(j, NCoord):
                        V4.append((d2Hdxidxi[j, k], [i, i, j, k]))
        V.append(V3)

    if Order >= 4:
        for i in range(NCoord):
            # iijk terms are handled in the cubic part
            for j in range(i + 1, NCoord):
                X0pipj = X0 + Coords[:, i] * dx + Coords[:, j] * dx
                X0mimj = X0 - Coords[:, i] * dx - Coords[:, j] * dx
                X0pimj = X0 + Coords[:, i] * dx - Coords[:, j] * dx
                X0mipj = X0 - Coords[:, i] * dx + Coords[:, j] * dx
                
                atompipj = CoordToAtom(atom0, X0pipj)
                atommimj = CoordToAtom(atom0, X0mimj)
                atompimj = CoordToAtom(atom0, X0pimj)
                atommipj = CoordToAtom(atom0, X0mipj)
                
                new_mol.atom = atompipj
                new_mol.build()
                new_mf = scf.RHF(new_mol)
                new_mf.kernel()
                Hpipj = GetHessian(new_mf, Method = Method, MassWeighted = False)
                
                new_mol.atom = atommimj
                new_mol.build()
                new_mf = scf.RHF(new_mol)
                new_mf.kernel()
                Hmimj = GetHessian(new_mf, Method = Method, MassWeighted = False)
                
                new_mol.atom = atompimj
                new_mol.build()
                new_mf = scf.RHF(new_mol)
                new_mf.kernel()
                Hpimj = GetHessian(new_mf, Method = Method, MassWeighted = False)

                new_mol.atom = atommipj
                new_mol.build()
                new_mf = scf.RHF(new_mol)
                new_mf.kernel()
                Hmipj = GetHessian(new_mf, Method = Method, MassWeighted = False)

                Hij = (Hpipj + Hmimj - Hpimj - Hmipj) / (4 * dx * dx)
                
                for k in range(j, NCoord):
                    for l in range(k, NCoord):
                        V4.append((Hij[k, l], [i, j, k, l]))
        V.append(V4)

    if Order >= 5:
        for i in range(NCoord):
            #iiijk terms
            X0pipipi = X0 + 3 * Coords[:, i] * dx
            X0pi = X0 + Coords[:, i] * dx
            X0mi = X0 - Coords[:, i] * dx
            X0mimimi = X0 - 3 * Coords[:, i] * dx

            atompipipi = CoordToAtom(atom0, X0pipipi)
            atommimimi = CoordToAtom(atom0, X0mimimi)
            atompi = CoordToAtom(atom0, X0pi)
            atommi = CoordToAtom(atom0, X0mi)

            new_mol.atom = atompipipi
            new_mol.build()
            new_mf = scf.RHF(new_mol)
            new_mf.kernel()
            Hpipipi = GetHessian(new_mf, Method = Method, MassWeighted = False)

            new_mol.atom = atommimimi
            new_mol.build()
            new_mf = scf.RHF(new_mol)
            new_mf.kernel()
            Hmimimi = GetHessian(new_mf, Method = Method, MassWeighted = False)

            new_mol.atom = atompi
            new_mol.build()
            new_mf = scf.RHF(new_mol)
            new_mf.kernel()
            Hpi = GetHessian(new_mf, Method = Method, MassWeighted = False)

            new_mol.atom = atommi
            new_mol.build()
            new_mf = scf.RHF(new_mol)
            new_mf.kernel()
            Hmi = GetHessian(new_mf, Method = Method, MassWeighted = False)

            Hiii = (Hpipipi - 3 * Hpi + 3 * Hmi - Hmimimi) / (8 * dx * dx * dx)
            for j in range(i, NCoord):
                for k in range(j, NCoord):
                    V5.append((Hiii[j, k], [i, i, i, j, k]))

            # Handle sextic iiiijk terms here
            if Order >= 6:
                X0pipi = X0 + 2 * Coords[:, i] * dx
                X0mimi = X0 - 2 * Coords[:, i] * dx
                
                atompipi = CoordToAtom(atom0, X0pipi)
                atommimi = CoordToAtom(atom0, X0mimi)

                new_mol.atom = atompipi
                new_mol.build()
                new_mf = scf.RHF(new_mol)
                new_mf.kernel()
                Hpipi = GetHessian(new_mf, Method = Method, MassWeighted = False)

                new_mol.atom = atommimi
                new_mol.build()
                new_mf = scf.RHF(new_mol)
                new_mf.kernel()
                Hmimi = GetHessian(new_mf, Method = Method, MassWeighted = False)

                Hiiii = (Hpipi - 4 * Hpi + 6 * H0 - 4 * Hmi + Hmimi) / (dx * dx * dx * dx)
                for j in range(i, NCoord):
                    for k in range(j, NCoord):
                        V6.append((Hiiii[j, k], [i, i, i, i, j, k]))

            for j in range(i + 1, NCoord):
                # iijkl elements.
                X0pipipj = X0 + 2 * Coords[:, i] * dx + Coords[:, j] * dx
                X0mimipj = X0 - 2 * Coords[:, i] * dx + Coords[:, j] * dx
                X0pj = X0 + Coords[:, j] * dx
                X0pipimj = X0 + 2 * Coords[:, i] * dx - Coords[:, j] * dx
                X0mimimj = X0 - 2 * Coords[:, i] * dx - Coords[:, j] * dx
                X0mj = X0 - Coords[:, j] * dx

                atompipipj = CoordToAtom(atom0, X0pipipj)
                atommimipj = CoordToAtom(atom0, X0mimipj)
                atompj = CoordToAtom(atom0, X0pj)
                atompipimj = CoordToAtom(atom0, X0pipimj)
                atommimimj = CoordToAtom(atom0, X0mimimj)
                atommj = CoordToAtom(atom0, X0mj)

                new_mol.atom = atompipipj
                new_mol.build()
                new_mf = scf.RHF(new_mol)
                new_mf.kernel()
                Hpipipj = GetHessian(new_mf, Method = Method, MassWeighted = False)

                new_mol.atom = atommimipj
                new_mol.build()
                new_mf = scf.RHF(new_mol)
                new_mf.kernel()
                Hmimipj = GetHessian(new_mf, Method = Method, MassWeighted = False)

                new_mol.atom = atompj
                new_mol.build()
                new_mf = scf.RHF(new_mol)
                new_mf.kernel()
                Hpj = GetHessian(new_mf, Method = Method, MassWeighted = False)
                
                new_mol.atom = atompipimj
                new_mol.build()
                new_mf = scf.RHF(new_mol)
                new_mf.kernel()
                Hpipimj = GetHessian(new_mf, Method = Method, MassWeighted = False)
                
                new_mol.atom = atommimimj
                new_mol.build()
                new_mf = scf.RHF(new_mol)
                new_mf.kernel()
                Hmimimj = GetHessian(new_mf, Method = Method, MassWeighted = False)
                
                new_mol.atom = atommj
                new_mol.build()
                new_mf = scf.RHF(new_mol)
                new_mf.kernel()
                Hmj = GetHessian(new_mf, Method = Method, MassWeighted = False)

                Hiij = (Hpipipj + Hmimipj - 2 * Hpj - Hpipimj - Hmimimj + 2 * Hmj) / (8 * dx * dx * dx)
                for k in range(j, NCoord):
                    for l in range(k, NCoord):
                        V5.append((Hiij[k, l], [i, i, j, k, l]))

                for k in range(j + 1, NCoord):
                    X0pipjpk = X0 + Coords[:, i] * dx + Coords[:, j] * dx + Coords[:, k] * dx
                    X0pipjmk = X0 + Coords[:, i] * dx + Coords[:, j] * dx - Coords[:, k] * dx
                    X0pimjpk = X0 + Coords[:, i] * dx - Coords[:, j] * dx + Coords[:, k] * dx
                    X0mipjpk = X0 - Coords[:, i] * dx + Coords[:, j] * dx + Coords[:, k] * dx
                    X0pimjmk = X0 + Coords[:, i] * dx - Coords[:, j] * dx - Coords[:, k] * dx
                    X0mimjpk = X0 - Coords[:, i] * dx - Coords[:, j] * dx + Coords[:, k] * dx
                    X0mipjmk = X0 - Coords[:, i] * dx + Coords[:, j] * dx - Coords[:, k] * dx
                    X0mimjmk = X0 - Coords[:, i] * dx - Coords[:, j] * dx - Coords[:, k] * dx

                    atompipjpk = CoordToAtom(atom0, X0pipjpk)
                    atompipjmk = CoordToAtom(atom0, X0pipjmk)
                    atompimjpk = CoordToAtom(atom0, X0pimjpk)
                    atommipjpk = CoordToAtom(atom0, X0mipjpk)
                    atommimjpk = CoordToAtom(atom0, X0mimjpk)
                    atommipjmk = CoordToAtom(atom0, X0mipjmk)
                    atompimjmk = CoordToAtom(atom0, X0pimjmk)
                    atommimjmk = CoordToAtom(atom0, X0mimjmk)
 
                    new_mol.atom = atompipjpk
                    new_mol.build()
                    new_mf = scf.RHF(new_mol)
                    new_mf.kernel()
                    Hpipjpk = GetHessian(new_mf, Method = Method, MassWeighted = False)
                           
                    new_mol.atom = atompipjmk
                    new_mol.build()
                    new_mf = scf.RHF(new_mol)
                    new_mf.kernel()
                    Hpipjmk = GetHessian(new_mf, Method = Method, MassWeighted = False)

                    new_mol.atom = atompimjpk
                    new_mol.build()
                    new_mf = scf.RHF(new_mol)
                    new_mf.kernel()
                    Hpimjpk = GetHessian(new_mf, Method = Method, MassWeighted = False)

                    new_mol.atom = atommipjpk
                    new_mol.build()
                    new_mf = scf.RHF(new_mol)
                    new_mf.kernel()
                    Hmipjpk = GetHessian(new_mf, Method = Method, MassWeighted = False)

                    new_mol.atom = atommimjpk
                    new_mol.build()
                    new_mf = scf.RHF(new_mol)
                    new_mf.kernel()
                    Hmimjpk = GetHessian(new_mf, Method = Method, MassWeighted = False)
                   
                    new_mol.atom = atommipjmk
                    new_mol.build()
                    new_mf = scf.RHF(new_mol)
                    new_mf.kernel()
                    Hmipjmk = GetHessian(new_mf, Method = Method, MassWeighted = False)

                    new_mol.atom = atompimjmk
                    new_mol.build()
                    new_mf = scf.RHF(new_mol)
                    new_mf.kernel()
                    Hpimjmk = GetHessian(new_mf, Method = Method, MassWeighted = False)
                   
                    new_mol.atom = atommimjmk
                    new_mol.build()
                    new_mf = scf.RHF(new_mol)
                    new_mf.kernel()
                    Hmimjmk = GetHessian(new_mf, Method = Method, MassWeighted = False)

                    Hijk = (Hpipjpk - Hpipjmk - Hpimjpk - Hmipjpk + Hmimjpk + Hmipjmk + Hpimjmk - Hmimjmk) / (8 * dx * dx * dx)
                    for l in range(k, NCoord):
                        for p in range(l, NCoord):
                            V5.append((Hijk[l, p], [i, j, k, l, p]))
        V.append(V5)
    
    if Order >= 6:
        for i in range(NCoord):
            # iiiijk terms handled in quintic derivatives
            for j in range(i + 1, NCoord):
                #iiijkl terms
                X0pipipipj = X0 + 3 * Coords[:, i] * dx + Coords[:, j] * dx
                X0pipj = X0 + Coords[:, i] * dx + Coords[:, j] * dx
                X0mipj = X0 - Coords[:, i] * dx + Coords[:, j] * dx
                X0mimimipj = X0 - 3 * Coords[:, i] * dx + Coords[:, j] * dx
                X0pipipimj = X0 + 3 * Coords[:, i] * dx - Coords[:, j] * dx
                X0pimj = X0 + Coords[:, i] * dx - Coords[:, j] * dx
                X0mimj = X0 - Coords[:, i] * dx - Coords[:, j] * dx
                X0mimimimj = X0 - 3 * Coords[:, i] * dx - Coords[:, j] * dx

                Hpipipipj = CoordToHessian(X0pipipipj, atom0, new_mol, Method = Method)
                Hpipj = CoordToHessian(X0pipj, atom0, new_mol, Method = Method)
                Hmipj = CoordToHessian(X0mipj, atom0, new_mol, Method = Method)
                Hmimimipj = CoordToHessian(X0mimimipj, atom0, new_mol, Method = Method)
                Hpipipimj = CoordToHessian(X0pipipimj, atom0, new_mol, Method = Method)
                Hpimj = CoordToHessian(X0pimj, atom0, new_mol, Method = Method)
                Hmimj = CoordToHessian(X0mimj, atom0, new_mol, Method = Method)
                Hmimimimj = CoordToHessian(X0mimimimj, atom0, new_mol, Method = Method)

                Hiiij = (Hpipipipj - 3 * Hpipj + 3 Hmipj - Hmimimipj - Hpipipimj + 3 * Hpimj - 3 * Hmipj + Hmimimimj) / (16 * dx**4)
                for k in range(j, NCoord):
                    for l in range(k, NCoord):
                        V6.append((Hiiij[k, l], [i, i, i, j, k, l]))

                for k in range(j + 1, NCoord):
                    Xpipipjpk = PerturbCoord(X0, [[i, 2], [j, 1], [k, 1]], Coords, dx)
                    Xmimipjpk = PerturbCoord(X0, [[i, -2], [j, 1], [k, 1]], Coords, dx)
                    Xpjpk = PerturbCoord(X0, [[j, 1], [k, 1]], Coords, dx)
                    Xpipimjpk = PerturbCoord(X0, [[i, 2], [j, -1], [k, 1]], Coords, dx)
                    Xmimimjpk = PerturbCoord(X0, [[i, -2], [j, -1], [k, 1]], Coords, dx)
                    Xmjpk = PerturbCoord(X0, [[j, -1], [k, 1]], Coords, dx)
                    Xpipipjmk = PerturbCoord(X0, [[i, 2], [j, 1], [k, -1]], Coords, dx)
                    Xmimipjmk = PerturbCoord(X0, [[i, -2], [j, 1], [k, -1]], Coords, dx)
                    Xpjmk = PerturbCoord(X0, [[j, 1], [k, -1]], Coords, dx)
                    Xpipimjmk = PerturbCoord(X0, [[i, 2], [j, -1], [k, -1]], Coords, dx)
                    Xmimimjmk = PerturbCoord(X0, [[i, -2], [j, -1], [k, -1]], Coords, dx)
                    Xmjmk = PerturbCoord(X0, [[j, -1], [k, -1]], Coords, dx)

                    Hpipipjpk = CoordToHessian(X0pipipjpk, atom0, new_mol, Method = Method)
                    Hmimipjpk = CoordToHessian(X0mimipjpk, atom0, new_mol, Method = Method)
                    Hpjpk = CoordToHessian(X0pjpk, atom0, new_mol, Method = Method)
                    Hpipimjpk = CoordToHessian(X0pipimjpk, atom0, new_mol, Method = Method)
                    Hmimimjpk = CoordToHessian(X0mimimjpk, atom0, new_mol, Method = Method)
                    Hmjpk = CoordToHessian(X0mjpk, atom0, new_mol, Method = Method)
                    Hpipipjmk = CoordToHessian(X0pipipjmk, atom0, new_mol, Method = Method)
                    Hmimipjmk = CoordToHessian(X0mimipjmk, atom0, new_mol, Method = Method)
                    Hpjmk = CoordToHessian(X0pjmk, atom0, new_mol, Method = Method)
                    Hpipimjmk = CoordToHessian(X0pipimjmk, atom0, new_mol, Method = Method)
                    Hmimimjmk = CoordToHessian(X0mimimjmk, atom0, new_mol, Method = Method)
                    Hmjmk = CoordToHessian(X0mjmk, atom0, new_mol, Method = Method)

                    Hiijk = (Hpipipjpk + Hmimipjpk - 2 * Hpjpk - Hpipimjpk - Hmimimjpk + 2 * Hmjpk - Hpipipjmk - Hmimipjmk + 2 * Hpjmk + Hpipimjmk + Hmimimjmk - 2 * Hmjmk) / (16 * dx**4)
                    for l in range(k, NCoord):
                        for p in range(l, NCoord):
                            V6.append((Hiijk[l, p], [i, i, j, k, l, p]))

                    for l in range(k + 1, NCoord):
                        Xpppp = PerturbCoord(X0, [[i, 1], [j, 1], [k, 1], [l, 1]], Coords, dx)
                        Xmppp = PerturbCoord(X0, [[i, -1], [j, 1], [k, 1], [l, 1]], Coords, dx)
                        Xpmpp = PerturbCoord(X0, [[i, 1], [j, -1], [k, 1], [l, 1]], Coords, dx)
                        Xmmpp = PerturbCoord(X0, [[i, -1], [j, -1], [k, 1], [l, 1]], Coords, dx)
                        Xppmp = PerturbCoord(X0, [[i, 1], [j, 1], [k, -1], [l, 1]], Coords, dx)
                        Xmpmp = PerturbCoord(X0, [[i, -1], [j, 1], [k, -1], [l, 1]], Coords, dx)
                        Xpmmp = PerturbCoord(X0, [[i, 1], [j, -1], [k, -1], [l, 1]], Coords, dx)
                        Xmmmp = PerturbCoord(X0, [[i, -1], [j, -1], [k, -1], [l, 1]], Coords, dx)
                        Xpppm = PerturbCoord(X0, [[i, 1], [j, 1], [k, 1], [l, -1]], Coords, dx)
                        Xmppm = PerturbCoord(X0, [[i, -1], [j, 1], [k, 1], [l, -1]], Coords, dx)
                        Xpmpm = PerturbCoord(X0, [[i, 1], [j, -1], [k, 1], [l, -1]], Coords, dx)
                        Xmmpm = PerturbCoord(X0, [[i, -1], [j, -1], [k, 1], [l, -1]], Coords, dx)
                        Xppmm = PerturbCoord(X0, [[i, 1], [j, 1], [k, -1], [l, -1]], Coords, dx)
                        Xmpmm = PerturbCoord(X0, [[i, -1], [j, 1], [k, -1], [l, -1]], Coords, dx)
                        Xpmmm = PerturbCoord(X0, [[i, 1], [j, -1], [k, -1], [l, -1]], Coords, dx)
                        Xmmmm = PerturbCoord(X0, [[i, -1], [j, -1], [k, -1], [l, -1]], Coords, dx)

                        Hpppp = CoordToHessian(Xpppp, atom0, new_mol, Method = Method)
                        Hmppp = CoordToHessian(Xmppp, atom0, new_mol, Method = Method)
                        Hpmpp = CoordToHessian(Xpmpp, atom0, new_mol, Method = Method)
                        Hmmpp = CoordToHessian(Xmmpp, atom0, new_mol, Method = Method)
                        Hppmp = CoordToHessian(Xppmp, atom0, new_mol, Method = Method)
                        Hmpmp = CoordToHessian(Xmpmp, atom0, new_mol, Method = Method)
                        Hpmmp = CoordToHessian(Xpmmp, atom0, new_mol, Method = Method)
                        Hmmmp = CoordToHessian(Xmmmp, atom0, new_mol, Method = Method)
                        Hpppm = CoordToHessian(Xpppm, atom0, new_mol, Method = Method)
                        Hmppm = CoordToHessian(Xmppm, atom0, new_mol, Method = Method)
                        Hpmpm = CoordToHessian(Xpmpm, atom0, new_mol, Method = Method)
                        Hmmpm = CoordToHessian(Xmmpm, atom0, new_mol, Method = Method)
                        Hppmm = CoordToHessian(Xppmm, atom0, new_mol, Method = Method)
                        Hmpmm = CoordToHessian(Xmpmm, atom0, new_mol, Method = Method)
                        Hpmmm = CoordToHessian(Xpmmm, atom0, new_mol, Method = Method)
                        Hmmmm = CoordToHessian(Xmmmm, atom0, new_mol, Method = Method)

                        Hijkl = (Hpppp - Hmppp - Hpmpp + Hmmpp - Hppmp + Hmpmp + Hpmmp - Hmmmp - Hpppm + Hmppm + Hpmpm - Hmmpm + Hppmm - Hmpmm - Hpmmm + Hmmmm) / (16 * dx**4)
                        for p in range(l, NCoord):
                            for q in range(p, NCoord):
                                V6.append((Hijkl[p, q], [i, j, k, l, p, q]))
        V.append(V6)
    
    return V

if __name__ == "__main__":
    from vstr.ff.normal_modes import GetNormalModes
