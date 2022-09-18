import numpy as np
from functools import reduce

def FormW(mVHCI, V):
    Order = len(V.shape)
    Modes = mVHCI.w.shape[0]
    ws = []
    for i in range(Order):
        ws.append(mVHCI.w)
    wProd = reduce(np.multiply, np.ix_(*ws))
    wProd *= 2**Order
    wProd = wProd**-0.5

    return V * wProd

def ConnectedD(D, Qs):
    Conn0 = [D.copy()]
    Conn1 = []
    for q in Qs:
        for Con0 in Conn0:
            tmpC = Con0.copy()
            tmpC[q] += 1
            Conn1.append(tmpC)
            tmpC = Con0.copy()
            tmpC[q] -= 1
            if tmpC[q] < 0:
                continue
            Conn1.append(tmpC)
        Conn0 = Conn1.copy()
        Conn1 = []

    # Now remove duplicates
    Conn1 = []
    for C0 in Conn0:
        if C0 not in Conn1:
            Conn1.append(C0)
    
    return Conn1

def ScreenD(mVHCI, Ws = None, C = None, eps = 1):
    if Ws is None:
        Ws = mVHCI.Ws
    if C is None:
        C = mVHCI.Cs[0]
    NewDs = []
    for W in Ws:
        WShape = W.shape
        WN = np.prod(W.shape)
        WFlat = W.reshape(WN)
        WSortedIndex = np.flip(np.argsort(abs(WFlat)))
        CSortedIndex = np.flip(np.argsort(abs(C)))
        for i in WSortedIndex:
            for n in CSortedIndex:
                if abs(WFlat[i] * C[n]) > eps:
                    Qs = np.unravel_index(i, WShape)
                    AddedD = ConnectedD(mVHCI.Ds[n], Qs)
                    NewDs += AddedD
                else:
                    break
    NewDsUnique = []
    for D in NewDs:
        if D not in NewDsUnique and D not in mVHCI.Ds:
            NewDsUnique.append(D)
    return NewDsUnique, len(NewDsUnique)

def HCIStep(mVHCI, eps = 1):
    NewDs, NAdded = mVHCI.ScreenD(Ws = mVHCI.Ws, C = mVHCI.Cs[0], eps = eps)
    mVHCI.Ds += NewDs
    return NAdded

def HCI(mVHCI):
    NAdded = len(mVHCI.Ds)
    it = 1
    while float(NAdded) / float(len(mVHCI.Ds)) > mVHCI.tol:
        NAdded = mVHCI.HCIStep(eps = mVHCI.eps1)
        mVHCI.Diagonalize()
        print("VHCI Iteration", it, "complete with", NAdded, "new configurations.")
        it += 1
        if it > mVHCI.MaxIter:
            raise RuntimeError("VHCI did not converge.")

def PT2(mVHCI, NStates = 1):
    dEs_PT2 = [None] * NStates
    for i in range(NStates):
        PertDs, NPert = mVHCI.ScreenD(Ws = mVHCI.Ws, C = mVHCI.Cs[i], eps = mVHCI.eps2)
        H_MN = mVHCI.HamV(Ds = mVHCI.Ds, DLefts = PertDs) # OD Hamiltonian between original and perturbative space
        E_M = mVHCI.HamV(Ds = PertDs, DiagonalOnly = True)
        Hc2 = H_MN @ mVHCI.Cs[i]
        Hc2 = Hc2 * Hc2
        Denom = (mVHCI.Es[i] - E_M)**(-1)
        dEs_PT2[i] = np.einsum('m,m->', Hc2, Denom)
    mVHCI.dEs_PT2 = dEs_PT2

def Diagonalize(mVHCI):
    mVHCI.H = mVHCI.HamV()
    mVHCI.Es, mVHCI.Cs = np.linalg.eigh(mVHCI.H)

'''
Forms the explicit vibrational Hamiltonian in the basis given by Ds
'''
def HamV(mVHCI, Ds = None, DLefts = None, DiagonalOnly = False):
    if Ds is None:
        Ds = mVHCI.Ds
    OffDiagonal = True 
    if DLefts is None:
        # DLefts is the Bras, if off diagonal elements are desired.
        OffDiagonal = False
        DLefts = Ds
    N = len(Ds)
    NL = len(DLefts)
    H = np.zeros((NL, N))
    if not OffDiagonal:
        for i in range(N):
            for j, n in enumerate(Ds[i]):
                H[i, i] += (n + 0.5) * mVHCI.w[j] # HO diagonal elements
    if DiagonalOnly:
        return H.diagonal()

    # Now we loop through each V
    for i, D in enumerate(Ds):
        for W in mVHCI.Ws:
            # Go through each derivative
            WShape = W.shape
            WN = np.prod(W.shape)
            WFlat = W.reshape(WN)
            for iw in range(WN):
                WElement = WFlat[iw]
                if abs(WElement) < 1e-12:
                    continue
                QIndices = np.unravel_index(iw, WShape)
                # Determine the connected determinants for this derivative
                Conn0 = ConnectedD(D, QIndices)
                # Loop through the connected determinants and add in matrix elements
                for DConn in Conn0:
                    try:
                        j = DLefts.index(DConn)
                    except:
                        continue
                    if not OffDiagonal and j < i:
                        continue
                    Hji = 1
                    for Di, Dv in enumerate(D):
                        if D[Di] != DConn[Di]:
                            n = max(D[Di], DConn[Di])
                            Hji *= np.sqrt(n)
                    Hji *= WElement
                    H[j, i] += Hji
                    if not OffDiagonal:
                        H[i, j] += Hji
    if not OffDiagonal:
        assert(np.allclose(H, H.T))
    return H


            
'''
Class that handles VHCI
'''
class VHCI:
    FormW = FormW
    HamV = HamV
    HCI = HCI
    Diagonalize = Diagonalize
    HCIStep = HCIStep
    ConnectedD = ConnectedD
    ScreenD = ScreenD
    PT2 = PT2

    def __init__(self, w, Vs, Ds):
        self.w = w # Harmonic Frequencies
        self.Vs = Vs
        self.Ws = []
        self.Ds = Ds
        self.eps1 = 0.1
        self.eps2 = 0.01
        self.tol = 0.01
        for V in self.Vs:
            self.Ws.append(self.FormW(V))
        self.MaxIter = 1000

        # Initialize the Energies and Coefficients
        self.Diagonalize()

    def kernel(self):
        pass

if __name__ == "__main__":
    V2 = np.asarray([[0, 1], [1, 0]])
    Ds = [[0, 0], [0, 1], [0, 2], [1, 0], [1, 1], [1, 2], [2, 0], [2, 1], [2, 2]]
    w = np.asarray([1, 2])

    mVHCI = VHCI(w, [V2], Ds)
    print(mVHCI.Es)
    mVHCI.HCI()
    print(mVHCI.Es)
    mVHCI.PT2(NStates = 1)
