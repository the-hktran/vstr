import numpy as np
from vstr import utils
from vstr.cpp_wrappers.vhci_jf.vhci_jf_functions import WaveFunction, FConst, HOFunc # classes from JF's code
from vstr.cpp_wrappers.vhci_jf.vhci_jf_functions import GenerateHamV, GenerateSparseHamV, AddStatesHB, HeatBath_Sort_FC, DoPT2, DoSPT2, VCIHamFromVSCF, VCISparseHamFromVSCF, AddStatesHBWithMax
from vstr.utils.perf_utils import TIMER
from functools import reduce
import itertools
import math
from scipy import sparse

def FormW(mVHCI, V):
    Ws = []
    for v in V:
        W = FConst(v[0], v[1], True);
        Ws.append(W)
    return Ws

def FormWSD(mVHCI):
    WSD = [[], []]
    W3 = mVHCI.Potential[0]
    W4 = mVHCI.Potential[1]
    for i in range(mVHCI.NModes):
        Wi = 0.0
        for W in W3:
            # Two cases, Wiii and Wijj
            if W.QIndices.count(i) == 1:
                Wi += 2.0 * W.fc
            elif W.QIndices.count(i) == 3:
                Wi += 3.0 * W.fc
        if abs(Wi) > 1e-12:
            mW = FConst(Wi, [i], False)
            WSD[0].append(mW)
    for i in range(mVHCI.NModes):
        for j in range(i, mVHCI.NModes):
            Wij = 0.0
            for W in W4:
            # Four cases, Wiiii, Wiikk, Wijjj, Wijkk
                if W.QIndices.count(i) == 1 and W.QIndices.count(j) == 1 and i != j:
                    Wij += 2.0 * W.fc
                elif W.QIndices.count(i) == 2 and len(W.QUnique) == 2 and i == j:
                    Wij += 2.0 * W.fc
                elif (W.QIndices.count(i) == 1 and W.QIndices.count(j) == 3) or (W.QIndices.count(i) == 3 and W.QIndices.count(j) == 1):
                    Wij += 3.0 * W.fc
                elif W.QIndices.count(i) == 4 and i == j:
                    Wij += 4.0 * W.fc
            if abs(Wij) > 1e-12:
                mW = FConst(Wij, [i, j], False)
                WSD[1].append(mW)
    mVHCI.PotentialSD = WSD

def ScreenBasis(mVHCI, Ws = None, C = None, eps = 0.01):
    if Ws is None:
        Ws = mVHCI.PotentialListFull
    if C is None:
        C = mVHCI.C[0]
    UniqueBasis = AddStatesHBWithMax(mVHCI.Basis, Ws, C, eps, mVHCI.MaxQuanta)
    return UniqueBasis, len(UniqueBasis)

def HCIStep(mVHCI, eps = 0.01):
    mVHCI.Timer.start(2)
    NewBasis, NAdded = mVHCI.ScreenBasis(Ws = mVHCI.PotentialListFull, C = abs(mVHCI.C[:, :mVHCI.NStates]).max(axis = 1), eps = eps)
    mVHCI.Basis += NewBasis
    mVHCI.Timer.stop(2)
    return NAdded

def HCI(mVHCI):
    NAdded = len(mVHCI.Basis)
    it = 1
    while (float(NAdded) / float(len(mVHCI.Basis))) > mVHCI.tol:
        NAdded = mVHCI.HCIStep(eps = mVHCI.eps1)
        print("VHCI Iteration", it, "complete with", NAdded, "new configurations and a total of", len(mVHCI.Basis))
        mVHCI.SparseDiagonalize()
        it += 1
        if it > mVHCI.MaxIter:
            raise RuntimeError("VHCI did not converge.")

def PT2(mVHCI, doStochastic = False):
    assert(mVHCI.eps2 < mVHCI.eps1)
    if doStochastic:
        if mVHCI.eps3 < 0:
            mVHCI.dE_PT2, mVHCI.sE_PT2 = DoSPT2(mVHCI.C, mVHCI.E, mVHCI.Basis, mVHCI.PotentialListFull, mVHCI.PotentialList, mVHCI.Potential[0], mVHCI.Potential[1], mVHCI.Potential[2], mVHCI.Potential[3], mVHCI.eps2, mVHCI.NStates, mVHCI.NWalkers, mVHCI.NSamples, False, mVHCI.eps3)
        else:
            assert (mVHCI.eps3 < mVHCI.eps2)
            mVHCI.dE_PT2, mVHCI.sE_PT2 = DoSPT2(mVHCI.C, mVHCI.E, mVHCI.Basis, mVHCI.PotentialListFull, mVHCI.PotentialList, mVHCI.Potential[0], mVHCI.Potential[1], mVHCI.Potential[2], mVHCI.Potential[3], mVHCI.eps2, mVHCI.NStates, mVHCI.NWalkers, mVHCI.NSamples, True, mVHCI.eps3)
    else:
        mVHCI.dE_PT2 = DoPT2(mVHCI.C, mVHCI.E, mVHCI.Basis, mVHCI.PotentialListFull, mVHCI.PotentialList, mVHCI.Potential[0], mVHCI.Potential[1], mVHCI.Potential[2], mVHCI.Potential[3], mVHCI.eps2, mVHCI.NStates)
    mVHCI.E_HCI_PT2 = mVHCI.E_HCI + mVHCI.dE_PT2

def Diagonalize(mVHCI):
    mVHCI.Timer.start(1)
    H = VCIHamFromVSCF(mVHCI.Basis, mVHCI.Frequencies, mVHCI.PotentialList, mVHCI.ModalCs, mVHCI.GenericV)
    mVHCI.Timer.stop(1)
    mVHCI.Timer.start(0)
    mVHCI.E, mVHCI.C = np.linalg.eigh(H)
    mVHCI.Timer.stop(0)
    mVHCI.E_HCI = mVHCI.E[:mVHCI.NStates].copy()

def SparseDiagonalize(mVHCI):
    H = VCISparseHamFromVSCF(mVHCI.Basis, mVHCI.Frequencies, mVHCI.PotentialList, mVHCI.ModalCs, mVHCI.GenericV)
    mVHCI.Timer.stop(1)
    mVHCI.Timer.start(0)
    mVHCI.E, mVHCI.C = sparse.linalg.eigsh(H, k = mVHCI.NStates, which = 'SM')
    mVHCI.Timer.stop(0)
    mVHCI.E_HCI = mVHCI.E[:mVHCI.NStates].copy()

def InitTruncatedBasis(mVHCI, MaxQuanta, MaxTotalQuanta = None):
    Basis = []
    Bs = []
    B0 = [0] * mVHCI.NModes
    Bs.append(B0)
    Basis.append(B0)
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
    
    # We need to translate this into the wavefunction object
    BasisWF = []
    for B in Basis:
        WF = WaveFunction(B, mVHCI.Frequencies)
        BasisWF.append(WF)
    return BasisWF

def TranslateBasisToString(B):
    BString = ""
    for j, HO in enumerate(B.Modes):
        if HO.Quanta == 1:
            BString += 'w%d + ' % (j + 1)
        elif HO.Quanta > 1:
            BString += '%dw%d + ' % (HO.Quanta, j + 1)
    return BString[:-3]

'''
This function prints the linear combination of the most important product states
'''
def LCLine(mVHCI, n, thr = 1e-2):
    def BasisToString(B):
        BString = '|'
        for HO in B.Modes:
            BString += str(HO.Quanta)
        BString += '>'
        return BString
    
    LC = ''
    for i in range(mVHCI.C.shape[0]):
        if abs(mVHCI.C[i, n]) > thr:
            LC += str(mVHCI.C[i, n])
            LC += BasisToString(mVHCI.Basis[i])
            LC += ' + '
    return LC[:-3]
    

def PrintResults(mVHCI):
    if mVHCI.dE_PT2 is None:
        FinalE = mVHCI.E_HCI
    else:
        FinalE = mVHCI.E_HCI_PT2
    for n in range(FinalE.shape[0]):
        MaxBasis = np.argmax(abs(mVHCI.C[:, n]))
        BString = TranslateBasisToString(mVHCI.Basis[MaxBasis])
        Outline = '{:.8f}\t'.format(FinalE[n])
        if mVHCI.sE_PT2 is not None:
            Outline += '+/- {:.8E}\t'.format(mVHCI.sE_PT2[n])
        Outline += '\t%s' % (BString)
        LCString = mVHCI.LCLine(n)
        Outline += '\t%s' % (LCString)
        print(Outline)
            
'''
Class that handles VHCI
'''
class VCI:
    FormW = FormW
    FormWSD = FormWSD
    HCI = HCI
    Diagonalize = Diagonalize
    SparseDiagonalize = SparseDiagonalize
    HCIStep = HCIStep
    ScreenBasis = ScreenBasis
    PT2 = PT2
    InitTruncatedBasis = InitTruncatedBasis
    PrintResults = PrintResults
    LCLine = LCLine

    def __init__(self, mVSCF, MaxTotalQuanta = 2, NStates = 10, **kwargs):
        self.mVSCF = mVSCF
        self.ModalCs = mVSCF.Cs
        self.Frequencies = mVSCF.Frequencies # 1D array of all harmonic frequencies.
        self.NModes = mVSCF.Frequencies.shape[0]
        self.MaxQuanta = mVSCF.MaxQuanta
        self.Potential = mVSCF.Potential
        self.FormWSD()

        self.PotentialListFull = []
        for Wp in self.PotentialSD:
            self.PotentialListFull += Wp
        self.PotentialList = []
        for Wp in self.Potential:
            self.PotentialList += Wp
            self.PotentialListFull += Wp
        self.PotentialListFull = HeatBath_Sort_FC(self.PotentialListFull) # Only need to sort these once
        self.GenericV = mVSCF.AnharmTensor

        self.MaxTotalQuanta = MaxTotalQuanta
        for Nm in self.MaxQuanta:
            assert(Nm >= self.MaxTotalQuanta)
        self.Basis = utils.init_funcs.InitTruncatedBasis(self.NModes, self.Frequencies, self.MaxQuanta, MaxTotalQuanta = self.MaxTotalQuanta)
        self.eps1 = 0.1 # HB epsilon
        self.eps2 = 0.01 # PT2/SPT2 epsilon
        self.eps3 = -1.0 # SSPT2 epsilon, < 0 means do not do semi-stochastic
        self.tol = 0.01
        self.MaxIter = 1000
        self.NStates = NStates
        self.NWalkers = 200
        self.NSamples = 50
        self.dE_PT2 = None
        self.sE_PT2 = None

        self.__dict__.update(kwargs)

        self.Timer = TIMER(4)
        self.TimerNames = ['Diagonalize', 'Form Hamiltonian', 'Screen Basis', 'PT2 correction']

        # Initialize the Energies and Coefficients
        self.Timer.start(0)
        self.Diagonalize()
        self.Timer.stop(0)

    def kernel(self, doVHCI = True):
        self.Timer.start(0)
        self.Diagonalize()
        self.Timer.stop(0)

        if doVHCI:
            self.HCI()
        self.Timer.report(self.TimerNames)

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
    from vstr.mf.vscf import VSCF
    from vstr.vhci.vhci import VHCI
    w, MaxQuanta, MaxTotalQuanta, Vs, eps1, eps2, eps3, NWalkers, NSamples, NStates = Read('CLO2.inp')
    mf = VSCF(w, Vs, MaxQuanta = MaxQuanta, NStates = NStates, SlowV = False, ModeOcc = [0, 0, 0])
    mf.SCF(DoDIIS = False)
    mVCI = VCI(mf, MaxTotalQuanta, eps1 = eps1, eps2 = eps2, eps3 = eps3, NWalkers = NWalkers, NSamples = NSamples, NStates = NStates)
    #mVCI.Diagonalize()
    mVCI.kernel()
    mVHCI = VHCI(np.asarray(w), Vs, MaxQuanta = MaxQuanta, MaxTotalQuanta = MaxTotalQuanta, eps1 = eps1, eps2 = eps2, eps3 = eps3, NWalkers = NWalkers, NSamples = NSamples, NStates = NStates)
    mVHCI.HCI()
    #mVHCI.Diagonalize()
    #print(mVHCI.E_HCI)
    #print(mVHCI.E[:NStates])
    #mVHCI.PT2(doStochastic = True)
    #print(mVHCI.E[:NStates])
    #print(mVHCI.dE_PT2)
    #print(mVHCI.E_HCI_PT2)
    #mVHCI.PrintResults()
