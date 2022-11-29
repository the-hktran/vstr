from vstr.cpp_wrappers.vhci_jf.vhci_jf_functions import WaveFunction, FConst, HOFunc # classes from JF's code
import itertools
import numpy as np

def FormW(V):
    Ws = []
    for v in V:
        W = FConst(v[0], v[1], True);
        Ws.append(W)
    return Ws

def InitTruncatedBasis(NModes, Frequencies, MaxQuanta, MaxTotalQuanta = None):
    Basis = []
    Bs = []
    B0 = [0] * NModes
    Bs.append(B0)
    Basis.append(B0)
    if MaxTotalQuanta is None:
        MaxTotalQuanta = max(MaxQuanta)
    for m in range(MaxTotalQuanta):
        BNext = []
        for B in Bs:
            for i in range(len(B)):
                NewB = B.copy()
                NewB[i] += 1
                if (NewB[i] < MaxQuanta[i]):
                    if NewB not in BNext and NewB not in Basis:
                        BNext.append(NewB)
        Basis = Basis + BNext
        Bs = BNext.copy()
    
    print("Initial basis functions are:\n", Basis)

    BasisWF = []
    for B in Basis:
        WF = WaveFunction(B, Frequencies)
        BasisWF.append(WF)
    return BasisWF

def InitGridBasis(Frequencies, MaxQuanta, ListOnly = False):
    def IncrementBasis(B, Max, Mode):
        if B[Mode] >= Max[Mode] - 1:
            B[Mode] = 0
            Mode += 1
            return IncrementBasis(B, Max, Mode)
        else:
            B[Mode] += 1
            return B
    
    NModes = len(MaxQuanta)
    NBas = np.prod(MaxQuanta)
    B0 = np.zeros((NModes), dtype = np.uint64)
    BasisList = []
    BasisList.append(B0)
    for i in range(NBas - 1):
        Bi = BasisList[-1].copy()
        Bi = IncrementBasis(Bi, MaxQuanta, 0)
        BasisList.append(Bi)
    if ListOnly:
        return BasisList
    Basis = []
    for B in BasisList:
        Basis.append(WaveFunction(B, Frequencies))
    return Basis, BasisList
