import numpy as np
from pyscf import gto, scf, hessian
from vstr.ff.normal_modes import AtomToCoord, CoordToAtom, GetHessian, GetNumHessian
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
    mol.verbose = 0
    mol.build()
    new_mf = scf.RHF(mol)
    new_mf.kernel()
    return GetHessian(new_mf, Method = Method, MassWeighted = False)

def PerturbHessian(X0, Modes, Coords, dx, atom0, mol, Method = 'rhf'):
    X = PerturbCoord(X0, Modes, Coords, dx)
    H = CoordToHessian(X, atom0, mol, Method = Method)
    return Coords.T @ H @ Coords

def PerturbEnergy(X0, Modes, Coords, dx, atom0, mol, Method = 'rhf'):
    X = PerturbCoord(X0, Modes, Coords, dx)
    atom = CoordToAtom(atom0, X)
    mol.atom = atom
    mol.build()
    new_mf = scf.RHF(mol)
    new_mf.kernel(verbose = 0)
    return new_mf.e_tot


'''
Scales force field constant to wavenumbers: Assumes energy in hartrees, mass weighted normal modes in bohr * sqrt(amu)
and frequencies in cm-1
'''
def ScaleFC(FC, Freqs, Modes):
    n = len(Modes)
    ScaledFC = FC / ((2 * np.pi)**(n / 2) * constants.AMU_TO_ME**(n / 2) * constants.C_AU**(n / 2) * constants.BOHR_TO_CM**(n / 2)) * constants.HARTREE_TO_INVCM #/ constants.C_AU**((n-1)/4)
    for i in Modes:
        ScaledFC = ScaledFC / np.sqrt(Freqs[i])
    return ScaledFC

def GetFF(mf, Coords, Freqs, Order = 4, Method = 'rhf', dx = 1e-4, tol = 1.0, QuarticMin = None):
    V = []
    NCoord = Coords.shape[1]
    X0 = AtomToCoord(mf) # in Bohr
    atom0 = mf.mol._atom.copy()
    H0 = GetHessian(mf, Method = Method, MassWeighted = False)
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
                            # We want to screen terms like iiii and iijk but we are okay with terms like iijj
                            if QuarticMin is not None:
                                if Hiijk < QuarticMin and (j != k or (i == j and j == k)):
                                    continue
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
                            if QuarticMin is not None:
                                if Hijkl < QuarticMin:
                                    continue
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
                    if abs(Hiiijk) > tol:
                        V5.append((Hiiijk, [i, i, i, j, k]))

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
                            Hijklp = ScaleFC(Hijk[l, p], Freqs, [i, j, k, l, p])
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
                        Hiiijkl = ScaleFC(Hiiij[k, l], Freqs, [i, i, i, j, k, l])
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

def GetCHARMMFF(HessDir, H0, Freqs, dx, BaseName = 'mode', Order = 4, Semidiagonal = True, tol = 1.0):
    '''
    args
        HessDir (string): Directory where the hessians are stored
        BaseName (string): Name of hessian files following format BaseName+n.hess
        Order (int): Order of expansion
        Semidiagonal (bool): Keep only semidiagonal terms in even order expansions
    '''
    from vstr.utils.charmm_tools import ReadHessian
    V = [] 
    V3 = []
    V4 = []
    V5 = []
    V6 = []

    # read first line of file
    with open(HessDir + '/' + BaseName + '.hess', 'r') as f:
        NAtom = int(f.readline())
    NCoord = 3 * NAtom - 6

    if Order >= 3:
        for i in range(NCoord):
            Hpi = ReadHessian(HessDir + '/' + BaseName + '+' + str(i) + '.hess', NAtom)
            Hmi = ReadHessian(HessDir + '/' + BaseName + '-' + str(i) + '.hess', NAtom)

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
        if Semidiagonal:
            V.append(V4)

    if Order >= 4 and not Semidiagonal:
        for i in range(NCoord):
            # iijk terms are handled in the cubic part
            for j in range(i + 1, NCoord):
                Hpipj = ReadHessian(HessDir + '/' + BaseName + '+' + str(i) + '+' + str(j) + '.hess', NAtom)
                Hmimj = ReadHessian(HessDir + '/' + BaseName + '-' + str(i) + '-' + str(j) + '.hess', NAtom)
                Hpimj = ReadHessian(HessDir + '/' + BaseName + '+' + str(i) + '-' + str(j) + '.hess', NAtom)
                Hmipj = ReadHessian(HessDir + '/' + BaseName + '-' + str(i) + '+' + str(j) + '.hess', NAtom)
                
                Hij = (Hpipj + Hmimj - Hpimj - Hmipj) / (4 * dx * dx)
                
                for k in range(j, NCoord):
                    for l in range(k, NCoord):
                        Hijkl = ScaleFC(Hij[k, l], Freqs, [i, j, k, l])
                        if abs(Hijkl) > tol:
                            V4.append((Hijkl, [i, j, k, l]))
        V.append(V4)

    return V

def MakeMatrix(VList, N):
    V3 = np.zeros((N, N, N))
    for v in VList[0]:
        I = v[1]
        from itertools import permutations
        Is = list(permutations(I))
        for i in Is:
            V3[i] = v[0]
    V4 = np.zeros((N, N, N, N))
    for v in VList[1]:
        I = v[1]
        Is = list(permutations(I))
        for i in Is:
            V4[i] = v[0]
    return V3, V4

def MakeInputFile(Vs, Freqs, InpFile):
    f = open(InpFile, "w")
    f.write("Comment:\n")
    f.write("HCI_Eps: 0.1\n")
    f.write("NStates: 100\n")
    f.write("PT2: 1\n")
    f.write("PT2_Eps: 0.001\n")
    f.write("SPT2_Eps: -1 100 50\n")
    f.write("Max_quanta: 8\n")
    f.write("Modes: " + str(Freqs.shape[0]) + "\n")
    for i, w in enumerate(Freqs):
        f.write(str(i) + " " + str(w) + " " + str(20) + "\n")
    FCNum = 0
    for V in Vs:
        FCNum += len(V)
    f.write("Force_constants: " + str(FCNum) + "\n")
    for V in Vs:
        for v in V:
            Qs = " "
            for j in v[1]:
                Qs += str(j) + " "
            f.write(str(len(v[1])) + Qs + str(v[0]) + "\n")

def PruneVs(Vs, Max = None, type = []):
    PrunedV = Vs.copy()
    if Max is not None:
        PrunedV = []
        for V in Vs:
            Vi = []
            for v in V:
                if abs(v[0]) < Max:
                    Vi.append(v)
            PrunedV.append(Vi)

    for t in type:
        if t.upper() == 'SEMIDIAGONAL':
            SQFF = []
            for v in Vs[1]:
                if v[1][0] == v[1][1] and v[1][2] == v[1][3]:
                    SQFF.append(v)
            PrunedV[1] = SQFF
        if t.upper() == 'POSITIVE':
            V4 = []
            for v in Vs[1]:
                if v[1][0] == v[1][1]:
                    #iiii, iijj, iijk
                    if (v[1][1] == v[1][2] and v[1][1] == v[1][3]) or v[1][2] == v[1][3] or (v[1][1] != v[1][2] and v[1][2] != v[1][3]):
                        if v[0] > -20:
                            V4.append(v)
                    else:
                        V4.append(v)
                else:
                    V4.append(v)
            PrunedV[1] = V4
    return PrunedV

if __name__ == "__main__":
    from vstr.ff.normal_modes import GetNormalModes

    mol = gto.M()
    mol.atom ='''
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
    print(w, C)
    #np.save("normal_modes", C)
    '''
    I = np.eye(mol.natm * 3)
    V = GetFF(mf, I, w, Order = 5, Method = 'rhf')
    #print(V)
    V3, V4 = MakeMatrix(V, mol.natm * 3)
    V3NM = np.einsum('ijk,ia,jb,kc->abc', V3, C, C, C)
    V4NM = np.einsum('ijkl,ia,jb,kc,ld->abcd', V4, C, C, C, C)
    print(V3NM)
    print(V4NM)
    '''

    V = GetFF(mf, C, w, Order = 4, dx = 1e-1, Method = 'rhf')
    print(V) 
    MakeInputFile(V, w, "ch3cho.inp")
    #V = GetFF(mf, C, w, Order = 6, dx = 1e-2, Method = 'rhf')
    #print(V)
    #V = GetFF(mf, C, w, Order = 6, dx = 1e-3, Method = 'rhf')
    #print(V)
    #V = GetFF(mf, C, w, Order = 6, dx = 1e-4, Method = 'rhf')
    #print(V)
    #from vstr.vhci.vhci import VHCI
    #mVCI = VHCI(w, V, MaxQuanta = 10, MaxTotalQuanta = 4, eps1=1.0, Nstates = 10)
    #mVCI.kernel()
