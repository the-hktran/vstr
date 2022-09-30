import numpy as np
from vstr import utils
from functools import reduce
import itertools
import math

def FormW(mVHCI, V):
    # This part is for when the frequency appears in the denominator, but I guess that is usually included
    '''
    Modes = mVHCI.w.shape[0]
    ws = []
    for i in range(Order):
        ws.append(mVHCI.w)
    wProd = reduce(np.multiply, np.ix_(*ws))
    wProd *= 2**Order
    wProd = wProd**-0.5
    '''
    W = []
    Order = len(V[0][1])
    for v in V:
        WElement = v[0] / np.sqrt(2**Order)
        QIndUniq = list(set(list(v[1])))
        QIndCount = []
        for k in QIndUniq:
            QIndCount.append(v[1].count(k))
        DerivCoeff = 1.0
        for k in QIndCount:
            DerivCoeff /= math.factorial(k)
        WElement *= DerivCoeff
        W.append((WElement, v[1]))
    return W

def FormWSD(mVHCI):
    WSD = [[], []]
    for W in mVHCI.Ws:
        if len(W[0][1]) == 3:
            for i in range(len(mVHCI.w)):
                Wi = 0.0
                for w in W:
                    # Two cases, Wiii and Wijj
                    if w[1].count(i) == 1:
                        Wi += 2.0 * w[0]
                    elif w[1].count(i) == 3:
                        Wi += 3.0 * w[0]
                if abs(Wi) > 1e-12:
                    WSD[0].append((Wi, [i]))
        elif len(W[0][1]) == 4:
            for i in range(len(mVHCI.w)):
                for j in range(i, len(mVHCI.w)):
                    Wij = 0.0
                    for w in W:
                    # Four cases, Wiiii, Wiikk, Wijjj, Wijkk
                        if w[1].count(i) == 1 and w[1].count(j) == 1 and i != j:
                            Wij += 2.0 * w[0]
                        elif w[1].count(i) == 2 and len(list(set(w[1]))) == 2 and i == j:
                            Wij += 2.0 * w[0]
                        elif (w[1].count(i) == 1 and w[1].count(j) == 3) or (w[1].count(i) == 3 and w[1].count(j) == 1):
                            Wij += 3.0 * w[0]
                        elif w[1].count(i) == 4 and i == j:
                            Wij += 4.0 * w[0]
                    if abs(Wij) > 1e-12:
                        WSD[1].append((Wij, [i, j]))
    mVHCI.WSD = WSD

'''
Determines all possible basis funtions connected to the given coupling element
and returns the connected basis function with the coefficient as a list of tuples 
of a list and float. Each outer list element is a basis function connected to
the given basis function. The inner list is the basis function and the float
is the associated raising and lowering coefficient.
'''
def ConnectedBasis(mVHCI, B, Qs):
    Conn0 = [(B.copy(), 1.0)]
    Conn1 = []
    # Loop through each derivative and generate the raised and lowered basis function for each.
    # After each derivative, recursively replace the list of basis functions
    for q in Qs:
        for Con0, Coeff in Conn0:
            tmpC = Con0.copy()
            tmpC[q] += 1
            tmpCoeff = Coeff
            tmpCoeff *= np.sqrt(tmpC[q])
            Conn1.append((tmpC, tmpCoeff))
            tmpC = Con0.copy()
            tmpCoeff = Coeff
            tmpCoeff *= np.sqrt(tmpC[q])
            tmpC[q] -= 1
            if tmpC[q] < 0:
                continue
            Conn1.append((tmpC, tmpCoeff))
        Conn0 = Conn1.copy()
        Conn1 = []
    
    return Conn0

def ScreenBasis(mVHCI, Ws = None, C = None, eps = 0.01):
    if Ws is None:
        Ws = mVHCI.Ws + mVHCI.WSD
    if C is None:
        C = mVHCI.Cs[0]
    NewBasis = []
    for W in Ws:
        WElement = []
        for w in W:
            WElement.append(w[0])
        WElement = np.asarray(WElement)
        WSortedIndex = np.flip(np.argsort(abs(WElement)))
        CSortedIndex = np.flip(np.argsort(abs(C)))
        for i in WSortedIndex:
            #if abs(WElement[i] * C[CSortedIndex[0]]) > eps:
                for n in CSortedIndex:
                    if abs(WElement[i] * C[n]) > eps:
                        AddedB = mVHCI.ConnectedBasis(mVHCI.Basis[n], W[i][1])
                        for AddB in AddedB:
                            NewBasis += [AddB[0]]
                    else:
                        break
            #else:
            #    break
    UniqueBasis = []
    for B in NewBasis:
        if B not in UniqueBasis and B not in mVHCI.Basis:
            UniqueBasis.append(B)
    return UniqueBasis, len(UniqueBasis)

def HCIStep(mVHCI, eps = 0.01):
    # We use the maximum Cn from the first NStates states as our C vector
    NewBasis, NAdded = mVHCI.ScreenBasis(Ws = mVHCI.Ws + mVHCI.WSD, C = abs(mVHCI.Cs[:, :mVHCI.NStates]).max(axis = 1), eps = eps)
    mVHCI.Basis += NewBasis
    return NAdded

def HCI(mVHCI):
    NAdded = len(mVHCI.Basis)
    it = 1
    while float(NAdded) / float(len(mVHCI.Basis)) > mVHCI.tol:
        NAdded = mVHCI.HCIStep(eps = mVHCI.eps1)
        mVHCI.Diagonalize()
        print("VHCI Iteration", it, "complete with", NAdded, "new configurations and a total of", len(mVHCI.Basis))
        it += 1
        if it > mVHCI.MaxIter:
            raise RuntimeError("VHCI did not converge.")

def PT2(mVHCI):
    dEs_PT2 = [None] * mVHCI.NStates
    PertBasis, NPert = mVHCI.ScreenBasis(Ws = mVHCI.Ws + mVHCI.WSD, C = abs(mVHCI.Cs[:, :mVHCI.NStates]).max(axis = 1), eps = mVHCI.eps2)
    print("Perturbative Space contains", NPert, "basis states.")
    H_MN = mVHCI.HamV(Basis = mVHCI.Basis, BasisBras = PertBasis) # OD Hamiltonian between original and perturbative space
    E_M = np.diag(mVHCI.HamV(Basis = PertBasis))
    Hc = H_MN @ mVHCI.Cs[:, :mVHCI.NStates]
    for n in range(mVHCI.NStates):
        Denom = (mVHCI.Es[n] - E_M)**(-1)
        dEs_PT2[n] = np.einsum('m,m->', Hc[:, n]**2., Denom)
    mVHCI.dEs_PT2 = dEs_PT2

def Diagonalize(mVHCI):
    mVHCI.H = mVHCI.HamV()
    mVHCI.Es, mVHCI.Cs = np.linalg.eigh(mVHCI.H)

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
    BasisBrasDict = {str(B): i for i, B in enumerate(BasisBras)}
    if not OffDiagonal:
        for i in range(N):
            for j, n in enumerate(Basis[i]):
                H[i, i] += (n + 0.5) * mVHCI.w[j] # HO diagonal elements

    # Now we loop through each V
    for i, B in enumerate(Basis):
        for W in mVHCI.Ws:
            # Go through each derivative
            for w in W:
                WElement = w[0]
                if abs(WElement) < 1e-12:
                    continue # If there are zeros, skip them
                QIndices = w[1]

                # Determine the connected determinants for this derivative
                Conn0 = mVHCI.ConnectedBasis(B, QIndices)

                # Loop through the connected determinants and add in matrix elements
                for BConn in Conn0:
                    try:
                        j = BasisBrasDict[str(BConn[0])] #BasisBras.index(BConn)
                    except:
                        continue
                    if DiagonalOnly:
                        if j != i:
                            continue
                    if not OffDiagonal and j > i:
                        continue
                    Hji = (BConn[1] * WElement)
                    H[j, i] += Hji
                    if not OffDiagonal and i != j:
                        H[i, j] += Hji
    if not OffDiagonal:
        assert(np.allclose(H, H.T))
    return H

'''
Defines basis set with a max quanta designation for each mode
    MaxQuanta: (list) Maximum quanta in each mode
    return: List of lists, each list is a HO basis function
'''
def InitTruncatedBasis(mVHCI, MaxQuanta, MaxTotalQuanta = None):
    QuantaList = []
    for n in MaxQuanta:
        QuantaList.append(list(range(n)))
    Basis = list(itertools.product(*QuantaList))
    for i in range(len(Basis)):
        Basis[i] = list(Basis[i])
    if MaxTotalQuanta is not None:
        NewBasis = []
        for B in Basis:
            if sum(B) < MaxTotalQuanta + 1:
                NewBasis.append(B)
        Basis = NewBasis
    print("Initial basis functions are:\n", Basis)
    return Basis

            
'''
Class that handles VHCI
'''
class VHCI:
    FormW = FormW
    FormWSD = FormWSD
    HamV = HamV
    HCI = HCI
    Diagonalize = Diagonalize
    HCIStep = HCIStep
    ConnectedBasis = ConnectedBasis
    ScreenBasis = ScreenBasis
    PT2 = PT2
    InitTruncatedBasis = InitTruncatedBasis

    def __init__(self, w, Vs, MaxQuanta = 2, MaxTotalQuanta = 2, NStates = 10, **kwargs):
        self.w = w # Harmonic Frequencies
        self.Vs = Vs # Derivatives of the PES: Format: List of lists, the first index is the order and each inner list is a list of indices
        self.Ws = [] # Scaled derivatives to be filled later
        for V in self.Vs:
            self.Ws.append(self.FormW(V))
        self.FormWSD()
        if isinstance(MaxQuanta, int):
            MaxQuanta = [MaxQuanta] * len(w)
        self.MaxTotalQuanta = MaxTotalQuanta
        self.Basis = self.InitTruncatedBasis(MaxQuanta, MaxTotalQuanta = MaxTotalQuanta)
        self.eps1 = 0.1
        self.eps2 = 0.01
        self.tol = 0.01
        self.MaxIter = 1000
        self.NStates = NStates

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

    from vstr.utils.read_jf_input import Read
    w, MaxQuanta, MaxTotalQuanta, Vs, eps1, eps2, NStates = Read('CLO2.inp')
    mVHCI = VHCI(np.asarray(w), Vs, MaxQuanta = MaxQuanta, MaxTotalQuanta = MaxTotalQuanta, eps1 = eps1, eps2 = eps2, NStates = NStates)
    #mVHCI.Diagonalize()
    mVHCI.HCI()
    print(mVHCI.Es[:NStates])
    mVHCI.PT2()
    print(mVHCI.Es[:NStates])
    #print(mVHCI.dEs_PT2)
