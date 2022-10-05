import numpy as np
from vstr import utils
#from vstr.vhci.vci_headers import VDeriv
from vstr.cpp_wrappers.vhci_functions import VDerivCPP as VDeriv
from vstr.cpp_wrappers.vhci_functions import FormBasisConnectionsCPP, HamVCPP
from functools import reduce
import itertools
import math
from scipy import sparse

def FormW(mVHCI, V):
    # This part is for when the frequency appears in the denominator, but I guess that is usually included
    Ws = []
    #Order = len(V[0][1])
    for v in V:
        #WElement = v[0] / np.sqrt(2**Order)
        #QIndUniq = list(set(list(v[1])))
        #QIndCount = []
        #for k in QIndUniq:
        #    QIndCount.append(v[1].count(k))
        #DerivCoeff = 1.0
        #for k in QIndCount:
        #    DerivCoeff /= math.factorial(k)
        #WElement *= DerivCoeff
        W = VDeriv(v[0], v[1], True)
        Ws.append(W)
    return Ws

def FormWSD(mVHCI):
    WSD = [[], []]
    for Wp in mVHCI.Ws:
        if Wp[0].Order == 3:
            for i in range(len(mVHCI.w)):
                Wi = 0.0
                for W in Wp:
                    # Two cases, Wiii and Wijj
                    if W.QIndices.count(i) == 1:
                        Wi += 2.0 * W.W
                    elif W.QIndices.count(i) == 3:
                        Wi += 3.0 * W.W
                if abs(Wi) > 1e-12:
                    mW = VDeriv(Wi, [i], False)
                    WSD[0].append(mW)
        elif Wp[0].Order == 4:
            for i in range(len(mVHCI.w)):
                for j in range(i, len(mVHCI.w)):
                    Wij = 0.0
                    for W in Wp:
                    # Four cases, Wiiii, Wiikk, Wijjj, Wijkk
                        if W.QIndices.count(i) == 1 and W.QIndices.count(j) == 1 and i != j:
                            Wij += 2.0 * W.W
                        elif W.QIndices.count(i) == 2 and len(W.QUnique) == 2 and i == j:
                            Wij += 2.0 * W.W
                        elif (W.QIndices.count(i) == 1 and W.QIndices.count(j) == 3) or (W.QIndices.count(i) == 3 and W.QIndices.count(j) == 1):
                            Wij += 3.0 * W.W
                        elif W.QIndices.count(i) == 4 and i == j:
                            Wij += 4.0 * W.W
                    if abs(Wij) > 1e-12:
                        mW = VDeriv(Wij, [i, j], False)
                        WSD[1].append(mW)
    mVHCI.WSD = WSD

'''
Determines all possible basis funtions connected to the given coupling element
and returns the connected basis function with the coefficient as a list of tuples 
of a list and float. Each outer list element is a basis function connected to
the given basis function. The inner list is the basis function and the float
is the associated raising and lowering coefficient.
'''
def ConnectedBasis(mVHCI, B, Qs):
    Conn0 = [[B.copy(), 1.0]]
    Conn1 = []
    # Loop through each derivative and generate the raised and lowered basis function for each.
    # After each derivative, recursively replace the list of basis functions
    for q in Qs:
        for Con0, Coeff in Conn0:
            tmpC = Con0.copy()
            tmpC[q] += 1
            tmpCoeff = Coeff
            tmpCoeff *= np.sqrt(tmpC[q])
            Conn1.append([tmpC, tmpCoeff])
            tmpC = Con0.copy()
            tmpCoeff = Coeff
            tmpCoeff *= np.sqrt(tmpC[q])
            tmpC[q] -= 1
            if tmpC[q] < 0:
                continue
            Conn1.append([tmpC, tmpCoeff])
        Conn0 = Conn1.copy()
        Conn1 = []
    
    return Conn0

def FormBasisConnections(mVHCI, Basis):
    BasisConnections = []
    for B in Basis:
        ConnB = []
        for p in range(len(mVHCI.Ws)):
            for W in mVHCI.Ws[p]:
                ConnBByW = mVHCI.ConnectedBasis(B, W.QIndices)
                for i, conn in enumerate(ConnBByW):
                    ConnBByW[i][1] = conn[1] * W.W
                ConnB += ConnBByW
        BasisConnections.append(ConnB)
    return BasisConnections

def ScreenBasis(mVHCI, Ws = None, C = None, eps = 0.01):
    if Ws is None:
        Ws = mVHCI.Ws + mVHCI.WSD
    if C is None:
        C = mVHCI.Cs[0]
    NewBasis = []
    for Wp in Ws:
        WElement = []
        for W in Wp:
            WElement.append(W.W)
        WElement = np.asarray(WElement)
        WSortedIndex = np.flip(np.argsort(abs(WElement)))
        CSortedIndex = np.flip(np.argsort(abs(C)))
        for i in WSortedIndex:
            if abs(WElement[i] * C[CSortedIndex[0]]) > eps:
                for n in CSortedIndex:
                    if abs(WElement[i] * C[n]) > eps:
                        AddedB = mVHCI.ConnectedBasis(mVHCI.Basis[n], Wp[i].QIndices)
                        for AddB in AddedB:
                            NewBasis += [AddB[0]]
                    else:
                        break
            else:
                break
    UniqueBasis = []
    for B in NewBasis:
        if B not in UniqueBasis and B not in mVHCI.Basis:
            UniqueBasis.append(B)
    return UniqueBasis, len(UniqueBasis)

def HCIStep(mVHCI, eps = 0.01):
    # We use the maximum Cn from the first NStates states as our C vector
    mVHCI.OldBasis = mVHCI.Basis.copy()
    mVHCI.OldBasisConn = mVHCI.BasisConn.copy()
    
    NewBasis, NAdded = mVHCI.ScreenBasis(Ws = mVHCI.Ws + mVHCI.WSD, C = abs(mVHCI.Cs[:, :mVHCI.NStates]).max(axis = 1), eps = eps)
    mVHCI.NewBasis = NewBasis.copy()
    NewBasisConn = FormBasisConnectionsCPP(mVHCI.Ws, NewBasis)
    #NewBasisConn = mVHCI.FormBasisConnections(NewBasis)
    mVHCI.NewBasisConn = NewBasisConn.copy()
    
    mVHCI.Basis += NewBasis
    mVHCI.BasisConn += NewBasisConn
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
    PertBasisConn = FormBasisConnectionsCPP(mVHCI.Ws, PertBasis) #mVHCI.FormBasisConnections(PertBasis)
    print("Perturbative Space contains", NPert, "basis states.")
    #H_MN = mVHCI.HamV(Basis = mVHCI.Basis, BasisConn = mVHCI.BasisConn, BasisBras = PertBasis) # OD Hamiltonian between original and perturbative space
    H_MN = HamVCPP(mVHCI.Basis, mVHCI.BasisConn, PertBasis, mVHCI.w, mVHCI.Ws, False, True) # OD Hamiltonian between original and perturbative space
    #E_M = np.diag(mVHCI.HamV(Basis = PertBasis, BasisConn = PertBasisConn))
    E_M = np.diag(HamVCPP(PertBasis, PertBasisConn, PertBasis, mVHCI.w, mVHCI.Ws, True, False))
    Hc = H_MN @ mVHCI.Cs[:, :mVHCI.NStates]
    for n in range(mVHCI.NStates):
        Denom = (mVHCI.Es[n] - E_M)**(-1)
        dEs_PT2[n] = np.einsum('m,m->', Hc[:, n]**2., Denom)
    mVHCI.dEs_PT2 = dEs_PT2

def Diagonalize(mVHCI):
    #if mVHCI.NewBasis is None:
    mVHCI.H = HamVCPP(mVHCI.Basis, mVHCI.BasisConn, mVHCI.Basis, mVHCI.w, mVHCI.Ws, False, False) #mVHCI.HamV()
    '''
    else:
        HOld = mVHCI.H
        DimOld = HOld.shape[0]
        mVHCI.H = np.zeros((len(mVHCI.Basis), len(mVHCI.Basis)))
        HON = HamVCPP(mVHCI.OldBasis, mVHCI.OldBasisConn, mVHCI.NewBasis, mVHCI.w, mVHCI.Ws, False, True)
        HNN = HamVCPP(mVHCI.NewBasis, mVHCI.NewBasisConn, mVHCI.NewBasis, mVHCI.w, mVHCI.Ws, False, False)
        mVHCI.H[:DimOld, :DimOld] = HOld
        mVHCI.H[DimOld:, :DimOld] = HON
        mVHCI.H[:DimOld, DimOld:] = HON.T
        mVHCI.H[DimOld:, DimOld:] = HNN
    '''
    mVHCI.Es, mVHCI.Cs = np.linalg.eigh(mVHCI.H)

def SparseDiagonalize(mVHCI):
    mVHCI.H = mVHCI.SparseHamV()
    mVHCI.Es, mVHCI.Cs = sparse.linalg.eigsh(mVHCI.H, k = mVHCI.NStates, which = 'SM')

'''
Forms the explicit vibrational Hamiltonian in the basis given by self.Basis
'''
def HamV(mVHCI, Basis = None, BasisConn = None, BasisBras = None, DiagonalOnly = False):
    if Basis is None:
        Basis = mVHCI.Basis
    if BasisConn is None:
        BasisConn = mVHCI.BasisConn
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
        for BConn in BasisConn[i]:
            try:
                j = BasisBrasDict[str(BConn[0])]
            except:
                continue
            if DiagonalOnly:
                if j != i:
                    continue
            if not OffDiagonal and j > i:
                continue
            H[j, i] += BConn[1]
            if not OffDiagonal and i != j:
                H[i, j] += BConn[1]
    if not OffDiagonal:
        assert(np.allclose(H, H.T))
    return H

def SparseHamV(mVHCI, Basis = None, BasisConn = None, BasisBras = None, DiagonalOnly = False):
    if Basis is None:
        Basis = mVHCI.Basis
    if BasisConn is None:
        BasisConn = mVHCI.BasisConn
    OffDiagonal = True 
    if BasisBras is None:
        # BasisBras is the Bras, if off diagonal elements are desired.
        OffDiagonal = False
        BasisBras = Basis
    N = len(Basis)
    NL = len(BasisBras)
    HijDict = {}

    def AddElement(i, j, Hij, Dict):
        Key = str(i) + ',' + str(j)
        try:
            Dict[Key] += Hij
        except:
            Dict.update({Key : Hij})

    BasisBrasDict = {str(B): i for i, B in enumerate(BasisBras)}
    if not OffDiagonal:
        for i in range(N):
            for j, n in enumerate(Basis[i]):
                AddElement(i, i, (n + 0.5) * mVHCI.w[j], HijDict) # HO diagonal elements

    # Now we loop through each V
    for i, B in enumerate(Basis):
        for BConn in BasisConn[i]:
            try:
                j = BasisBrasDict[str(BConn[0])]
            except:
                continue
            if DiagonalOnly:
                if j != i:
                    continue
            if not OffDiagonal and j > i:
                continue
            AddElement(j, i, BConn[1], HijDict)
            if not OffDiagonal and i != j:
                AddElement(i, j, BConn[1], HijDict)

    # Turn dictionary into lists of indices and values
    HElem = []
    I = []
    J = []
    for i in range(len(Basis)):
        for j in range(len(BasisBras)):
            try:
                Hij = HijDict[str(i) + ',' + str(j)]
                I.append(i)
                J.append(j)
                HElem.append(Hij)
            except:
                continue
    H = sparse.csr_matrix((HElem, (I, J)), shape = (NL, N))
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
    SparseHamV = SparseHamV
    HCI = HCI
    Diagonalize = Diagonalize
    SparseDiagonalize = SparseDiagonalize
    HCIStep = HCIStep
    ConnectedBasis = ConnectedBasis
    FormBasisConnections = FormBasisConnections
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
        self.BasisConn = FormBasisConnectionsCPP(self.Ws, self.Basis) #self.FormBasisConnections(self.Basis)
        self.NewBasis = None
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
    mVHCI.Diagonalize()
    mVHCI.HCI()
    print(mVHCI.Es[:NStates])
    mVHCI.PT2()
    print(mVHCI.Es[:NStates])
    #print(mVHCI.dEs_PT2)
