import numpy as np
from functools import reduce
import itertools
import math

def FormW(mVHCI, V):
    print('form W...')
    Order = len(V.shape)
    Modes = mVHCI.w.shape[0]
    ws = []
    for i in range(Order):
        ws.append(mVHCI.w)
    wProd = reduce(np.multiply, np.ix_(*ws))
    wProd *= 2**Order
    wProd = wProd**-0.5
    print('...done')

    return V * wProd

def ConnectedBasis(mVHCI, B, Qs):
#    print('form connections...')
    Conn0 = [B.copy()]
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
#    print('...done')
    return Conn1

def ScreenBasis(mVHCI, Ws = None, C = None, eps = 1):
    print('screen basis...')
    if Ws is None:
        Ws = mVHCI.Ws
    if C is None:
        C = mVHCI.Cs[0]
    NewBasis = []
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
                    AddedB = mVHCI.ConnectedBasis(mVHCI.Basis[n], Qs)
                    NewBasis += AddedB
                else:
                    break
    UniqueBasis = []
    for B in NewBasis:
        if B not in UniqueBasis and B not in mVHCI.Basis:
            UniqueBasis.append(B)
    print('..done')
    return UniqueBasis, len(UniqueBasis)

def HCIStep(mVHCI, eps = 1):
    NewBasis, NAdded = mVHCI.ScreenBasis(Ws = mVHCI.Ws, C = mVHCI.Cs[0], eps = eps)
    mVHCI.Basis += NewBasis
    return NAdded

def HCI(mVHCI):
    NAdded = len(mVHCI.Basis)
    it = 1
    while float(NAdded) / float(len(mVHCI.Basis)) > mVHCI.tol:
        NAdded = mVHCI.HCIStep(eps = mVHCI.eps1)
        mVHCI.Diagonalize()
        print("VHCI Iteration", it, "complete with", NAdded, "new configurations.")
        it += 1
        if it > mVHCI.MaxIter:
            raise RuntimeError("VHCI did not converge.")

def PT2(mVHCI, NStates = 1):
    dEs_PT2 = [None] * NStates
    for i in range(NStates):
        PertBasis, NPert = mVHCI.ScreenBasis(Ws = mVHCI.Ws, C = mVHCI.Cs[i], eps = mVHCI.eps2)
        H_MN = mVHCI.HamV(Basis = mVHCI.Basis, BasisBras = PertBasis) # OD Hamiltonian between original and perturbative space
        E_M = mVHCI.HamV(Basis = PertBasis, DiagonalOnly = True)
        Hc2 = H_MN @ mVHCI.Cs[i]
        Hc2 = Hc2 * Hc2
        Denom = (mVHCI.Es[i] - E_M)**(-1)
        dEs_PT2[i] = np.einsum('m,m->', Hc2, Denom)
    mVHCI.dEs_PT2 = dEs_PT2

def Diagonalize(mVHCI):
    mVHCI.H = mVHCI.HamV()
    print('diagonalize...')
    mVHCI.Es, mVHCI.Cs = np.linalg.eigh(mVHCI.H)
    print('...done')

'''
Forms the explicit vibrational Hamiltonian in the basis given by self.Basis
'''
def HamV(mVHCI, Basis = None, BasisBras = None, DiagonalOnly = False):
    if Basis is None:
        Basis = mVHCI.Basis
    OffDiagonal = True 
    if BasisBras is None:
        # BasisBras is the Bras, if off diagonal elements are desired.
        OffDiagonal = False
        BasisBras = Basis
    N = len(Basis)
    NL = len(BasisBras)
    H = np.zeros((NL, N))
    BasisBrasDict = {str(B): i for B, i in enumerate(BasisBras)}
    if not OffDiagonal:
        for i in range(N):
            for j, n in enumerate(Basis[i]):
                H[i, i] += (n + 0.5) * mVHCI.w[j] # HO diagonal elements
    if DiagonalOnly:
        return H.diagonal()

    print('forming H...')
    # Now we loop through each V
    for i, B in enumerate(Basis):
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
                Conn0 = mVHCI.ConnectedBasis(B, QIndices)
                # Loop through the connected determinants and add in matrix elements
                for BConn in Conn0:
                    try:
                        j = BasisBrasDict[str(BConn)] #BasisBras.index(BConn)
                        print(j)
                    except:
                        continue
                    if not OffDiagonal and j < i:
                        continue
                    Hji = 1
                    for Bi, Bv in enumerate(B):
                        if B[Bi] != BConn[Bi]:
                            n = max(B[Bi], BConn[Bi])
                            Hji *= np.sqrt(n)
                    Hji *= WElement
                    H[j, i] += Hji
                    if not OffDiagonal:
                        H[i, j] += Hji
    if not OffDiagonal:
        assert(np.allclose(H, H.T))
    print('...done')
    return H

'''
Defines basis set with a max quanta designation for each mode
    MaxQuanta: (list) Maximum quanta in each mode
    return: List of lists, each list is a HO basis function
'''
def InitTruncatedBasis(mVHCI, MaxQuanta):
    QuantaList = []
    for n in MaxQuanta:
        QuantaList.append(list(range(n)))
    Basis = list(itertools.product(*QuantaList))
    for i in range(len(Basis)):
        Basis[i] = list(Basis[i])
    return Basis

            
'''
Class that handles VHCI
'''
class VHCI:
    FormW = FormW
    HamV = HamV
    HCI = HCI
    Diagonalize = Diagonalize
    HCIStep = HCIStep
    ConnectedBasis = ConnectedBasis
    ScreenBasis = ScreenBasis
    PT2 = PT2
    InitTruncatedBasis = InitTruncatedBasis

    def __init__(self, w, Vs, MaxQuanta = 2, **kwargs):
        self.w = w # Harmonic Frequencies
        self.Vs = Vs # Derivatives of the PES
        self.Ws = [] # Scaled derivatives to be filled later
        print(len(Vs))
        for V in self.Vs:
            self.Ws.append(self.FormW(V))
        if isinstance(MaxQuanta, int):
            MaxQuanta = [MaxQuanta] * len(w)
        print('init basis...')
        self.Basis = self.InitTruncatedBasis(MaxQuanta)
        print('...done')
        self.eps1 = 0.1
        self.eps2 = 0.01
        self.tol = 0.01
        self.MaxIter = 1000

        self.__dict__.update(kwargs)

        # Initialize the Energies and Coefficients
        self.Diagonalize()

    def kernel(self):
        pass

if __name__ == "__main__":
    '''
    V2 = np.asarray([[0, 1], [1, 0]])
    Ds = [[0, 0], [0, 1], [0, 2], [1, 0], [1, 1], [1, 2], [2, 0], [2, 1], [2, 2]]
    w = np.asarray([1, 2])

    mVHCI = VHCI(w, [V2], MaxQuanta = 3)
    print(mVHCI.Es)
    mVHCI.HCI()
    print(mVHCI.Es)
    mVHCI.PT2(NStates = 1)
    '''

    from read_jf_input import Read
    w, MaxQuanta, Vs, eps1, eps2, NStates = Read('C2H4_quartic.inp')
    print(w)
    mVHCI = VHCI(np.asarray(w), Vs, MaxQuanta = 3, eps1 = eps1, eps2 = eps2)
    mVHCI.HCI()
    print(mVHCI.Es)
