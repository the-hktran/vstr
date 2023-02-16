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
    if Order <= 3:
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

            if Order <= 4:
                # Can do iijk here as well
                d2Hdxidxi = (Hpi + Hmi - 2 * H0) / (dx * dx)
                for j in range(i, NCoord):
                    for k in range(j, NCoord):
                        V4.append((d2Hdxidxi[j, k], [i, i, j, k]))
        V.append(V3)

    if Order <= 4:
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





if __name__ == "__main__":
    from vstr.ff.normal_modes import GetNormalModes
