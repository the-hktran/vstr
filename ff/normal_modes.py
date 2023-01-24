from pyscf import gto, scf, hessian
import numpy as np

def GetHessian(mf, Method = 'rhf'):
    if Method == 'rhf':
        HRaw = hessian.RHF(mf).kernel()
        H = np.zeros((mf.mol.natm * 3, mf.mol.natm * 3))
        for i in range(mf.mol.natm):
            I = list(range(3 * i, 3 * i + 3))
            for j in range(mf.mol.natm):
                J = list(range(3 * j, 3 * j + 3))
                H[np.ix_(I, J)] = (HRaw[i, j] / (np.sqrt(mf.mol.atom_mass_list()[i] * mf.mol.atom_mass_list()[j])))
        return H
    else:
        raise ValueError("No hessian method available for that method")

def GetNormalModes(mf, Method = 'rhf'):
    H = GetHessian(mf, Method = Method)
    print(np.allclose(H, H.T))
    w, C = np.linalg.eigh(H)
    print(w)
    return w, C

if __name__ == '__main__':
    mol = gto.M()
    mol.fromfile("h2o.xyz")
    mf = scf.RHF(mol)
    mf.kernel()
    GetNormalModes(mf)
