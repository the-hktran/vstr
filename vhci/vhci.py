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
        if C0 not in Conn:
            Conn1.append(C0)
    
    return Conn1

def ScreenD(mVHCI, Ws = None, C = None, eps = 1e-3):
    if Ws is None:
        Ws = mVHCI.Ws
    if C is None:
        C = mVHCI.C
    NewDs = []
    NAdded = 0
    for W in Ws:
        WShape = W.shape
        WN = np.prod(W.shape)
        WFlat = W.reshape(WN)
        WSortedIndex = np.argsort(abs(WFlat))
        CSortedIndex = np.argsort(abs(C))
        for i in WSortedIndex:
            for n in CSortedIndex:
                if abs(WFlat[i] * mVHCI.C[n]) > eps:
                    Qs = np.unravel_index(i, WShape)
                    AddedD = ConnectedD([mVHCI.Ds[n].copy()], Qs)
                    NewDs += AddedD
                    NAdded += len(AddedD)
                else:
                    break
    return NewDs, NAdded

def HCIStep(mVHCI, eps = 1e-2):
    NewDs, NAdded = mVHCI.ScreenD(Ws = mVHCI.Ws, C = mVHCI.Cs[0], eps = 1e-3)
    mVHCI.Ds += NewDs
    return NAdded

def HCI(mVHCI, tol = 0.01):
    NAdded = len(mVHCI.Ds)
    it = 1
    while float(NAdded) / float(len(mVHCI.Ds)) > tol:
        NAdded = mVHCI.HCIStep()
        mVHCI.Diagonalize()
        print("VHCI Iteration", it, "complete with", NAdded, "new configurations.")
        it += 1
        if it > mVHCI.MaxIter:
            raise RuntimeError("VHCI did not converge.")

def PT2(mVHCI, NStates, eps = 1e-4):
    for i in range(NStates):
        PertDs = mVHCI.ScreenD(Ws = mVHCI.Ws, C = mVHCI.Cs[i], eps = eps)
        PertH = mVHCI.HamV(Ds = PertDs) # Hamiltonian in the perturbative space
        Hc2 = 
    

def Diagonalize(mVHCI):
    mVHCI.H = mVHCI.HamV()
    mVHCI.Es, mVHCI.Cs = np.linalg.eigh(mVHCI.H)
    mVHCI.E = Es[0]
    mVHCI.C = mVHCI.Cs[0]

'''
Forms the explicit vibrational Hamiltonian in the basis given by Ds
'''
def HamV(mVHCI, Ds = None, DLefts = None):
    if Ds = None:
        Ds = mVHCI.Ds
    N = len(Ds)
    H = np.zeros((N, N))
    for i in range(N):
        for j, n in enumerate(Ds[i]):
            H[i, i] += (n + 0.5) * mVHCI.w[j] # HO diagonal elements

    # Now we loop through each V
    for i, D in enumerate(Ds):
        for W in mVHCI.Ws:
            # Go through each derivative
            WShape = W.shape
            WN = np.prod(W.shape)
            WFlat = W.reshape(WN)
            ConnectedD = []
            for iw in range(WN):
                WElement = WFlat[iw]
                if abs(WElement) < 1e-12:
                    continue
                QIndices = np.unravel_index(iw, WShape)
                # Determine the connected determinants for this derivative
                Conn0 = ConnectedD([D.copy()], QIndices)
                # Loop through the connected determinants and add in matrix elements
                for DConn in Conn0:
                    try:
                        j = Ds.index(DConn)
                    except:
                        continue
                    Hij = 1
                    for Di, Dv in enumerate(D):
                        if D[Di] != DConn[Di]:
                            n = max(D[Di], DConn[Di])
                            Hij *= np.sqrt(n)
                    Hij *= WElement
                    H[i, j] += Hij
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

    def __init__(self, Vs, w, Ds):
        self.Vs = Vs # PES Derivatives
        self.w = w # Harmonic Frequencies
        self.Ws = []
        self.Ds = Ds
        for V in self.Vs:
            self.Ws.append(self.FormW(V))
        self.MaxIter = 1000

    def kernel(self):
        pass

if __name__ == "__main__":
    V2 = np.asarray([[0, 1], [1, 0]])
    Ds = [[0, 0], [0, 1], [0, 2], [1, 0], [1, 1], [1, 2], [2, 0], [2, 1], [2, 2]]
    w = np.asarray([1, 1])

    mVHCI = VHCI([V2], w, Ds)
    mVHCI.HamV()
    print(mVHCI.H)
