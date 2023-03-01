from pyscf import gto, scf, hessian
import numpy as np
#from pyscf.data import nist
from vstr.utils import constants

def AtomToCoord(mf):
    '''
    X = np.zeros((3 * mf.mol.natm))
    for i in range(mf.mol.natm):
        X[i * 3:i * 3 + 3] = mf.mol._atom[i][1]
    '''
    X = mf.mol.atom_coords(unit='Bohr').reshape(3 * mf.mol.natm)
    return X

def CoordToAtom(atom, X):
    atom_new = atom.copy()
    for i in range(len(atom)):
        atom_new[i] = (atom_new[i][0], list(X[i*3:i*3 + 3]))
    return atom_new

def GetNumHessian(mf, Coords = None, Method = 'rhf', dx = 1e-2, MassWeighted = True, isotope_avg = True):
    if Coords is None:
        Coords = np.eye(mf.mol.natm * 3)
    Mass = mf.mol.atom_mass_list(isotope_avg=isotope_avg)
    NCoords = Coords.shape[1]
    X0 = AtomToCoord(mf)
    atom0 = mf.mol._atom.copy()
    H = np.zeros((NCoords, NCoords))
    E0 = mf.e_tot
    for i in range(NCoords):
        X0pipi = X0 + 2 * Coords[:, i] * dx
        X0mimi = X0 - 2 * Coords[:, i] * dx

        atompipi = CoordToAtom(atom0, X0pipi)
        atommimi = CoordToAtom(atom0, X0mimi)

        new_mol = mf.mol.copy()
        new_mol.atom = atompipi
        new_mol.unit = 'B'
        new_mol.build()
        new_mf = scf.RHF(new_mol)
        new_mf.kernel(verbose=0)
        Epipi = new_mf.e_tot

        new_mol = mf.mol.copy()
        new_mol.atom = atommimi
        new_mol.unit = 'B'
        new_mol.build()
        new_mf = scf.RHF(new_mol)
        new_mf.kernel(verbose=0)
        Emimi = new_mf.e_tot

        H[i, i] = (Epipi + Emimi - 2 * E0) / (4 * dx * dx)
        if MassWeighted:
            H[i, i] = H[i, i] / Mass[i // 3]

        for j in range(i + 1, NCoords):
            print(i, j)
            X0pipj = X0 + Coords[:, i] * dx + Coords[:, j] * dx
            X0mimj = X0 - Coords[:, i] * dx - Coords[:, j] * dx
            X0pimj = X0 + Coords[:, i] * dx - Coords[:, j] * dx
            X0mipj = X0 - Coords[:, i] * dx + Coords[:, j] * dx
            x = X0pipi - X0
            atompipj = CoordToAtom(atom0, X0pipj)
            atommimj = CoordToAtom(atom0, X0mimj)
            atompimj = CoordToAtom(atom0, X0pimj)
            atommipj = CoordToAtom(atom0, X0mipj)
            
            new_mol = mf.mol.copy()
            new_mol.atom = atompipj
            new_mol.unit = 'B'
            new_mol.build()
            new_mf = scf.RHF(new_mol)
            new_mf.kernel(verbose=0)
            Epipj = new_mf.e_tot
            
            new_mol = mf.mol.copy()
            new_mol.atom = atommimj
            new_mol.unit = 'B'
            new_mol.build()
            new_mf = scf.RHF(new_mol)
            new_mf.kernel(verbose=0)
            Emimj = new_mf.e_tot
            
            new_mol = mf.mol.copy()
            new_mol.atom = atompimj
            new_mol.unit = 'B'
            new_mol.build()
            new_mf = scf.RHF(new_mol)
            new_mf.kernel(verbose=0)
            Epimj = new_mf.e_tot

            new_mol = mf.mol.copy()
            new_mol.atom = atommipj
            new_mol.unit = 'B'
            new_mol.build()
            new_mf = scf.RHF(new_mol)
            new_mf.kernel(verbose=0)
            Emipj = new_mf.e_tot

            H[i, j] = (Epipj + Emimj - Epimj - Emipj) / (4 * dx * dx)
            if MassWeighted:
                H[i, j] = H[i, j] / (np.sqrt(Mass[i // 3] * Mass[j // 3]))
            H[j, i] = H[i, j]

    return H

def GetHessian(mf, Method = 'rhf', MassWeighted = False, isotope_avg=True):
    mass = mf.mol.atom_mass_list(isotope_avg=isotope_avg)
    if Method == 'rhf':
        HRaw = hessian.RHF(mf).kernel()
        if MassWeighted:
            H = np.einsum('pqxy,p,q->pqxy', HRaw, mass**-0.5, mass**-0.5)
        else:
            H = HRaw
        H = H.transpose(0, 2, 1, 3).reshape(mf.mol.natm * 3, mf.mol.natm * 3)
        '''
        H = np.zeros((mf.mol.natm * 3, mf.mol.natm * 3))
        for i in range(mf.mol.natm):
            I = list(range(3 * i, 3 * i + 3))
            for j in range(mf.mol.natm):
                J = list(range(3 * j, 3 * j + 3))
                H[np.ix_(I, J)] = (HRaw[i, j] / (np.sqrt(mass[i] * mass[j])))
        '''
        return H
    else:
        raise ValueError("No hessian method available for that method")

def GetNormalModes(mf, H = None, Method = 'rhf', tol = 1e-1):
    if H is None:
        H = GetHessian(mf, Method = Method, MassWeighted = True)
    w, C = np.linalg.eigh(H)
    print(w)
    C = C[:, np.where(w > tol)[0]]
    w = w[np.where(w > tol)[0]]
    w = np.sqrt(w / constants.AMU_TO_ME) / (2 * np.pi * constants.C_AU * constants.BOHR_TO_CM)

    # Need to mass weight C
    Mass = mf.mol.atom_mass_list(isotope_avg=True)
    for i in range(mf.mol.natm):
        C[(3 * i):(3 * i + 3), :] /= np.sqrt(Mass[i])
    return w, C

if __name__ == '__main__':
    mol = gto.M()
    mol.fromfile("h2o.xyz")
    mol.basis='sto-3g'
    mol.build()
    mf = scf.RHF(mol)
    mf.kernel()

    #pyH = GetHessian(mf)# hessian.RHF(mf).kernel()
    #myH = GetNumHessian(mf)
    #print(pyH)
    #print(myH)

    w, C = GetNormalModes(mf)
    
    HA = GetHessian(mf, MassWeighted = True)
    #HN = GetNumHessian(mf, MassWeight=True)
    #print(HA-HN)
    #print(((HA-HN)**2).sum())
    #w, C = GetNormalModes(mf, H = HN)
    #CMW = C.copy()
    #for j in range(mf.mol.natm):
    #    CMW[(3 * j):(3 * j + 3),:] /= np.sqrt(mf.mol.atom_mass_list(isotope_avg=True)[j])
    HN = GetNumHessian(mf, Coords = C, MassWeighted = False)
    #print(HA)
    #Htest = GetNumHessian(mf, Coords = C, MassWeight=False)
    #print(C.T @ HA @ C)
    print(HN)
    #print(C.T @ HN @ C)
    #print(Htest)
