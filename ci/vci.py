import numpy as np
from vstr import utils
from vstr.cpp_wrappers.vhci_jf.vhci_jf_functions import WaveFunction, FConst, HOFunc # classes from JF's code
from vstr.cpp_wrappers.vhci_jf.vhci_jf_functions import HeatBath_Sort_FC, DoPT2FromVSCF, DoSPT2FromVSCF, VCIHamFromVSCF, VCISparseHamFromVSCF, AddStatesHB, AddStatesHBWithMax, AddStatesHBFromVSCF, ContractedAnharmonicPotential, ContractedHOTerms
from vstr.utils.perf_utils import TIMER
from functools import reduce
import itertools
import math
from scipy import sparse

def ReadBasisFromFile(mVHCI, FileName):
    mVHCI.Basis = []
    with open(FileName + "_basis.chk") as CHKFile:
        for Line in CHKFile:
            B = Line.split()
            B = [int(b) for b in B]
            mVHCI.Basis.append(WaveFunction(B, mVHCI.Frequencies))
    mVHCI.E = np.load(FileName + "_E.npy")
    mVHCI.C = np.load(FileName + "_C.npy")

def SaveBasisToFile(mVHCI, FileName):
    CHKFile = open(FileName + "_basis.chk", 'w')
    for B in mVHCI.Basis:
        Line = ""
        for i in range(mVHCI.NModes):
            Line += str(B.Modes[i].Quanta) + " "
        CHKFile.write(Line[:-1])
        CHKFile.write('\n')
    CHKFile.close()

    np.save(mVHCI.CHKFile + "_E", mVHCI.E)
    np.save(mVHCI.CHKFile + "_C", mVHCI.C)

def SumQuanta(mVHCI, Basis = None):
    if Basis is None:
        Basis = mVHCI.Basis
    SummedQuanta = [0] * mVHCI.NModes
    for B in Basis:
        for i in range(mVHCI.NModes):
            SummedQuanta[i] += B.Modes[i].Quanta
    return SummedQuanta

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
    
    if mVHCI.HBMethod == 'exact':
        UniqueBasis = AddStatesHBFromVSCF(mVHCI.Basis, Ws, C, eps, mVHCI.Ys)
    elif mVHCI.HBMethod == 'ho_max':
        UniqueBasis = AddStatesHBWithMax(mVHCI.Basis, Ws, C, eps, mVHCI.MaxQuanta, mVHCI.HighestQuanta)
    elif mVHCI.HBMethod == 'ho_orig':
        UniqueBasis = AddStatesHB(mVHCI.Basis, Ws, C, eps)
    return UniqueBasis, len(UniqueBasis)

def HCIStep(mVHCI, eps = 0.01):
    mVHCI.Timer.start(2)
    NewBasis, NAdded = mVHCI.ScreenBasis(Ws = mVHCI.PotentialListFull, C = abs(mVHCI.C[:, :mVHCI.NStates]).max(axis = 1), eps = eps)
    mVHCI.Basis += NewBasis
    mVHCI.Timer.stop(2)
    return NAdded, NewBasis

def HCI(mVHCI):
    NAdded = len(mVHCI.Basis)
    it = 1
    while (float(NAdded) / float(len(mVHCI.Basis))) > mVHCI.tol:
        NAdded, mVHCI.NewBasis = mVHCI.HCIStep(eps = mVHCI.eps1)
        print("VHCI Iteration", it, "complete with", NAdded, "new configurations and a total of", len(mVHCI.Basis), flush = True)
        mVHCI.SparseDiagonalize()
        if mVHCI.PrintHCISteps:
            mVHCI.PrintResults()
        if mVHCI.SaveToFile:
            mVHCI.SaveBasisToFile(mVHCI.CHKFile)
        it += 1
        if it > mVHCI.MaxIter:
            raise RuntimeError("VHCI did not converge.")
    mVHCI.NewBasis = None
    mVHCI.H = None

def PT2(mVHCI, doStochastic = False):
    assert(mVHCI.eps2 < mVHCI.eps1)
    if doStochastic:
        if mVHCI.eps3 < 0:
            mVHCI.Timer.start(4)
            mVHCI.dE_PT2, mVHCI.sE_PT2 = DoSPT2FromVSCF(mVHCI.C, mVHCI.E, mVHCI.Basis, mVHCI.PotentialListFull, mVHCI.PotentialList, mVHCI.eps2, mVHCI.NStates, mVHCI.NWalkers, mVHCI.NSamples, mVHCI.Ys, mVHCI.Xs, False, mVHCI.eps3)
            mVHCI.Timer.stop(4)
        else:
            mVHCI.Timer.start(5)
            assert (mVHCI.eps3 < mVHCI.eps2)
            mVHCI.dE_PT2, mVHCI.sE_PT2 = DoSPT2FromVSCF(mVHCI.C, mVHCI.E, mVHCI.Basis, mVHCI.PotentialListFull, mVHCI.PotentialList, mVHCI.eps2, mVHCI.NStates, mVHCI.NWalkers, mVHCI.NSamples, mVHCI.Ys, mVHCI.Xs, True, mVHCI.eps3)
            mVHCI.Timer.stop(5)
    else:
        mVHCI.Timer.start(3)
        mVHCI.dE_PT2 = DoPT2FromVSCF(mVHCI.C, mVHCI.E, mVHCI.Basis, mVHCI.PotentialListFull, mVHCI.PotentialList, mVHCI.eps2, mVHCI.NStates, mVHCI.Ys, mVHCI.Xs)
        mVHCI.Timer.stop(3)
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
    mVHCI.Timer.start(1)
    if mVHCI.H is None:
        mVHCI.H = VCISparseHamFromVSCF(mVHCI.Basis, mVHCI.Basis, mVHCI.Frequencies, mVHCI.PotentialList, mVHCI.Ys, mVHCI.Xs, True)
    else:
        if len(mVHCI.NewBasis) != 0:
            HIJ = VCISparseHamFromVSCF(mVHCI.Basis[:-len(mVHCI.NewBasis)], mVHCI.NewBasis, mVHCI.Frequencies, mVHCI.PotentialList, mVHCI.Ys, mVHCI.Xs, False)
            HJJ = VCISparseHamFromVSCF(mVHCI.NewBasis, mVHCI.NewBasis, mVHCI.Frequencies, mVHCI.PotentialList, mVHCI.Ys, mVHCI.Xs, True)
            mVHCI.H = sparse.hstack([mVHCI.H, HIJ])
            mVHCI.H = sparse.vstack([mVHCI.H, sparse.hstack([HIJ.transpose(), HJJ])])
    mVHCI.Timer.stop(1)
    mVHCI.Timer.start(0)
    mVHCI.E, mVHCI.C = sparse.linalg.eigsh(mVHCI.H, k = mVHCI.NStates, which = 'SM')
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
        print(Outline, flush = True)

def PrintParameters(mVHCI):
    print("______________________________")
    print("|                            |")
    print("|       PARAMETER LIST       |")
    print("|____________________________|")
    print("")
    print("eps_1          :", mVHCI.eps1)
    print("eps_2          :", mVHCI.eps2)
    print("eps_3          :", mVHCI.eps3)
    print("NStates        :", mVHCI.NStates)
    print("NWalkers       :", mVHCI.NWalkers)
    print("NSamples       :", mVHCI.NSamples)
    print("HB Criterion   :", mVHCI.HBMethod)
    print("", flush = True)

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
    SumQuanta = SumQuanta
    PrintResults = PrintResults
    PrintParameters = PrintParameters
    LCLine = LCLine
    SaveBasisToFile = SaveBasisToFile
    ReadBasisFromFile = ReadBasisFromFile

    def __init__(self, mVSCF, MaxTotalQuanta = 2, NStates = 10, **kwargs):
        self.mVSCF = mVSCF
        self.ModalCs = mVSCF.Cs
        self.Frequencies = mVSCF.Frequencies # 1D array of all harmonic frequencies.
        self.NModes = mVSCF.Frequencies.shape[0]
        self.MaxQuanta = mVSCF.MaxQuanta
        self.Potential= mVSCF.Potential
        self.FormWSD()

        self.PotentialListFull = []
        for Wp in self.PotentialSD:
            self.PotentialListFull += Wp
        self.PotentialList = []
        for Wp in self.Potential:
            self.PotentialListFull += Wp
            self.PotentialList += Wp
        self.PotentialListFull = HeatBath_Sort_FC(self.PotentialListFull) # Only need to sort these once
        self.GenericV = mVSCF.AnharmTensor

        self.MaxTotalQuanta = MaxTotalQuanta
        self.IncludeSqrt = True
        for Nm in self.MaxQuanta:
            assert(Nm >= self.MaxTotalQuanta)
        self.H = None
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
        self.HBMethod = 'exact' # ['ho_orig', 'ho_max', 'exact']

        self.CHKFile = None
        self.ReadFromFile = False
        self.SaveToFile = False
        self.PrintHCISteps = False

        self.__dict__.update(kwargs)

        self.Timer = TIMER(6)
        self.TimerNames = ['Diagonalize', 'Form Hamiltonian', 'Screen Basis', 'PT2 correction', 'SPT2 correction', 'SSPT2 correction']

    def kernel(self, doVCI = True, doVHCI = True, doPT2 = False, doSPT2 = False, ComparePT2 = False):
        if self.SaveToFile or self.ReadFromFile:
            assert(self.CHKFile is not None)

        if self.ReadFromFile:
            self.ReadBasisFromFile(self.CHKFile)
        else:
            self.Basis = utils.init_funcs.InitTruncatedBasis(self.NModes, self.Frequencies, self.MaxQuanta, MaxTotalQuanta = self.MaxTotalQuanta)

        self.Ys = ContractedAnharmonicPotential(self.ModalCs, self.GenericV)
        self.Xs = ContractedHOTerms(self.ModalCs, self.Frequencies)
        
        self.PrintParameters()

        if doVCI:
            self.Timer.start(0)
            self.SparseDiagonalize()
            print("===== VSCF-VCI RESULTS =====", flush = True)
            self.PrintResults()
            print("")
            self.Timer.stop(0)
        
        if doVHCI:
            print("===== VSCF-VHCI RESULTS =====", flush = True)
            self.HCI()
            self.PrintResults()
            print("")

        if doPT2:
            print("===== VSCF-VHCI+PT2 RESULTS =====", flush = True)
            self.PT2(doStochastic = False)
            self.PrintResults()
            print("")
        if doSPT2:
            print("===== VSCF-VHCI+SPT2 RESULTS =====", flush = True)
            self.PT2(doStochastic = True)
            self.PrintResults()
            print("")
        if ComparePT2:
            print("===== VSCF-VHCI+SSPT2 RESULTS =====", flush = True)
            eps = self.eps2
            self.eps2 *= 10
            self.eps3 = eps
            self.PrintParameters()
            self.PT2(doStochastic = True)
            self.PrintResults()
            print("")
        self.Timer.report(self.TimerNames)

    @property
    def AverageQuanta(self):
        if self.IncludeSqrt:
            return [Q / len(self.Basis) for Q in self.SummedQuanta]
        else:
            return [0] * self.NModes
    @property
    def HighestQuanta(self):
        HQ = [0] * self.NModes
        for B in self.Basis:
            for n in range(self.NModes):
                if B.Modes[n].Quanta > HQ[n]:
                    HQ[n] = B.Modes[n].Quanta
        return HQ

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
    #w, MaxQuanta, MaxTotalQuanta, Vs, eps1, eps2, eps3, NWalkers, NSamples, NStates = Read('../examples/jf_input/CLO2.inp')
    w, MaxQuanta, MaxTotalQuanta, Vs, eps1, eps2, eps3, NWalkers, NSamples, NStates = Read('CLO2.inp')
    mf = VSCF(w, Vs, MaxQuanta = MaxQuanta, NStates = NStates)
    mf.kernel()
    mVCI = VCI(mf, MaxTotalQuanta, eps1 = eps1, eps2 = eps2, eps3 = eps3, NWalkers = NWalkers, NSamples = NSamples, NStates = NStates)
    #mVCI.CHKFile = "chk"
    #mVCI.ReadFromFile=True
    #mVCI.SaveToFile=True
    mVCI.kernel(doVHCI = True, doPT2 = True, doSPT2 = True, ComparePT2 = True)
    #mVCI.PrintResults()

    '''
    mVHCI = VHCI(np.asarray(w), Vs, MaxQuanta = MaxQuanta, MaxTotalQuanta = MaxTotalQuanta, eps1 = eps1, eps2 = eps2, eps3 = eps3, NWalkers = NWalkers, NSamples = NSamples, NStates = NStates)
    mVHCI.kernel(doVHCI = True, doPT2 = False, doSPT2 = False, ComparePT2 = False)
    #mVHCI.PT2(doStochastic = True)
    '''
