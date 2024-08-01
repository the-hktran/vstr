import numpy as np
from vstr.utils import init_funcs
from vstr.utils.perf_utils import TIMER
from vstr.utils import constants
from vstr.cpp_wrappers.vhci_jf.vhci_jf_functions import WaveFunction, FConst, HOFunc # classes from JF's code
from vstr.cpp_wrappers.vhci_jf.vhci_jf_functions import GenerateHamV, GenerateSparseHamV, GenerateSparseHamVOD, GenerateHamAnharmV, AddStatesHB, AddStatesHBWithMax, AddStatesHBFromVSCF, HeatBath_Sort_FC, DoPT2, DoSPT2, AddStatesHBStoreCoupling, VCISparseHamNMode, VCISparseHamNModeFromOM, ConnectedStatesCIPSI, AddStatesCIPSI, AddStatesHB2Mode, VCISparseT
from functools import reduce
import itertools
import math
from scipy import sparse
import numdifftools as nd

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

def MakeAnharmTensor(mVHCI, PotentialList = None):
    if PotentialList is None:
        PotentialList = mVHCI.PotentialList
    AnharmTensor = []
    
    MaxQ = [mVHCI.MaxQuanta[0]]
    Freq = [1.0]
    GenericBasis = init_funcs.InitGridBasis(Freq, MaxQ)[0]

    AnharmTensor.append(np.zeros((0,0)))
    for p in range(1, 7):
        W = FConst(1.0, [0] * p, False)
        W = [W]
        W.append(FConst(0.0, [0] * 6, False))
        CubicFC = []
        QuarticFC = []
        QuinticFC = []
        SexticFC = []
        if p == 1 or p == 3:
            CubicFC = W.copy()
        elif p == 2 or p == 4:
            QuarticFC = W.copy()
        elif p == 5:
            QuinticFC = W.copy()
        elif p == 6:
            SexticFC = W.copy()
        H = GenerateHamAnharmV(GenericBasis, Freq, W, CubicFC, QuarticFC, QuinticFC, SexticFC)
        AnharmTensor.append(H)
    AnharmTensor[0] = AnharmTensor[1]
    return AnharmTensor

def FormW(mVHCI, V, Scaled = True):
    Ws = []
    for v in V:
        W = FConst(v[0], v[1], Scaled);
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

    if mVHCI.HBMethod == 'orig':
        UniqueBasis = AddStatesHB(mVHCI.Basis, Ws, C, eps)
    elif mVHCI.HBMethod == 'max':
        UniqueBasis = AddStatesHBWithMax(mVHCI.Basis, Ws, C, eps, mVHCI.MaxQuanta, mVHCI.HighestQuanta)
    elif mVHCI.HBMethod == 'exact' or mVHCI.HBMethod.upper() == 'QFF':
        UniqueBasis = AddStatesHBFromVSCF(mVHCI.Basis, Ws, C, eps, mVHCI.Ys)
    elif mVHCI.HBMethod == 'coupling':
        UniqueBasis = AddStatesHBStoreCoupling(mVHCI.Basis, Ws, C, eps, mVHCI.Ys)
        return UniqueBasis, len(UniqueBasis[0])
    elif mVHCI.HBMethod.upper() == '2MODE':
        UniqueBasis = AddStatesHB2Mode(mVHCI.Basis, mVHCI.mol.ints[1], mVHCI.Sorted2Mode, C, eps, True) 
    elif mVHCI.HBMethod.upper() == 'CIPSI':
        ConnectedBasis = ConnectedStatesCIPSI(mVHCI.Basis, mVHCI.MaxQuanta, mVHCI.mol.Order)
        UniqueBasis = AddStatesCIPSI(mVHCI.Basis, ConnectedBasis, C, mVHCI.E, mVHCI.Frequencies, mVHCI.mol.V0, mVHCI.mol.ints[0], mVHCI.mol.ints[1], mVHCI.mol.ints[2], eps)

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
    #mVHCI.H = None
    mVHCI.NewBasis = None

def PT2(mVHCI, doStochastic = False):
    assert(mVHCI.eps2 < mVHCI.eps1)
    if doStochastic:
        if mVHCI.eps3 < 0:
            mVHCI.Timer.start(4)
            mVHCI.dE_PT2, mVHCI.sE_PT2 = DoSPT2(mVHCI.C, mVHCI.E, mVHCI.Basis, mVHCI.PotentialListFull, mVHCI.PotentialList, mVHCI.Potential[0], mVHCI.Potential[1], mVHCI.Potential[2], mVHCI.Potential[3], mVHCI.Ys, mVHCI.eps2, mVHCI.NStatesPT2, mVHCI.NWalkers, mVHCI.NSamples, False, mVHCI.eps3)
            mVHCI.Timer.stop(4)
        else:
            mVHCI.Timer.start(5)
            assert (mVHCI.eps3 < mVHCI.eps2)
            mVHCI.dE_PT2, mVHCI.sE_PT2 = DoSPT2(mVHCI.C, mVHCI.E, mVHCI.Basis, mVHCI.PotentialListFull, mVHCI.PotentialList, mVHCI.Potential[0], mVHCI.Potential[1], mVHCI.Potential[2], mVHCI.Potential[3], mVHCI.Ys, mVHCI.eps2, mVHCI.NStatesPT2, mVHCI.NWalkers, mVHCI.NSamples, True, mVHCI.eps3)
            mVHCI.Timer.stop(5)
    else:
        mVHCI.Timer.start(3)
        mVHCI.dE_PT2 = DoPT2(mVHCI.C, mVHCI.E, mVHCI.Basis, mVHCI.PotentialListFull, mVHCI.PotentialList, mVHCI.Potential[0], mVHCI.Potential[1], mVHCI.Potential[2], mVHCI.Potential[3], mVHCI.Ys, mVHCI.eps2, mVHCI.NStatesPT2)
        mVHCI.Timer.stop(3)
    mVHCI.E_HCI_PT2 = mVHCI.E_HCI[:mVHCI.NStatesPT2] + mVHCI.dE_PT2

def Diagonalize(mVHCI):
    mVHCI.Timer.start(1)
    H = GenerateHamV(mVHCI.Basis, mVHCI.Frequencies, mVHCI.PotentialList, mVHCI.Potential[0], mVHCI.Potential[1], mVHCI.Potential[2], mVHCI.Potential[3])
    mVHCI.Timer.stop(1)
    mVHCI.Timer.start(0)
    mVHCI.E, mVHCI.C = np.linalg.eigh(H)
    mVHCI.Timer.stop(0)
    mVHCI.E_HCI = mVHCI.E[:mVHCI.NStates].copy()

def SparseDiagonalize(mVHCI):
    mVHCI.Timer.start(1)
    if mVHCI.H is None:
        mVHCI.H = GenerateSparseHamV(mVHCI.Basis, mVHCI.Frequencies, mVHCI.PotentialList, mVHCI.Potential[0], mVHCI.Potential[1], mVHCI.Potential[2], mVHCI.Potential[3])
    else:
        if len(mVHCI.NewBasis) != 0:
            HIJ = GenerateSparseHamVOD(mVHCI.Basis[:-len(mVHCI.NewBasis)], mVHCI.NewBasis, mVHCI.Frequencies, mVHCI.PotentialList, mVHCI.Potential[0], mVHCI.Potential[1], mVHCI.Potential[2], mVHCI.Potential[3])
            HJJ = GenerateSparseHamV(mVHCI.NewBasis, mVHCI.Frequencies, mVHCI.PotentialList, mVHCI.Potential[0], mVHCI.Potential[1], mVHCI.Potential[2], mVHCI.Potential[3])
            mVHCI.H = sparse.hstack([mVHCI.H, HIJ])
            mVHCI.H = sparse.vstack([mVHCI.H, sparse.hstack([HIJ.transpose(), HJJ])])
    mVHCI.Timer.stop(1)
    mVHCI.Timer.start(0)
    mVHCI.E, mVHCI.C = sparse.linalg.eigsh(mVHCI.H, k = mVHCI.NStates, which = 'SA')
    mVHCI.Timer.stop(0)
    mVHCI.E_HCI = mVHCI.E[:mVHCI.NStates].copy()

def SparseDiagonalizeNMode(mVHCI):
    mVHCI.Timer.start(1)
    if mVHCI.H is None:
        if mVHCI.mol.use_onemode_states:
            mVHCI.H = VCISparseHamNModeFromOM(mVHCI.Basis, mVHCI.Basis, mVHCI.Frequencies, mVHCI.mol.V0, mVHCI.mol.onemode_eig, mVHCI.mol.ints[1].tolist(), mVHCI.mol.ints[2].tolist(), True)
        else:
            mVHCI.H = VCISparseHamNMode(mVHCI.Basis, mVHCI.Basis, mVHCI.Frequencies, mVHCI.mol.V0, mVHCI.mol.ints[0].tolist(), mVHCI.mol.ints[1].tolist(), mVHCI.mol.ints[2].tolist(), True)
    else:
        if len(mVHCI.NewBasis) != 0:
            if mVHCI.mol.use_onemode_states:
                HIJ = VCISparseHamNModeFromOM(mVHCI.Basis[:-len(mVHCI.NewBasis)], mVHCI.NewBasis, mVHCI.Frequencies, mVHCI.mol.V0, mVHCI.mol.onemode_eig, mVHCI.mol.ints[1].tolist(), mVHCI.mol.ints[2].tolist(), False)
                HJJ = VCISparseHamNModeFromOM(mVHCI.NewBasis, mVHCI.NewBasis, mVHCI.Frequencies, mVHCI.mol.V0, mVHCI.mol.onemode_eig, mVHCI.mol.ints[1].tolist(), mVHCI.mol.ints[2].tolist(), True)
            else:
                HIJ = VCISparseHamNMode(mVHCI.Basis[:-len(mVHCI.NewBasis)], mVHCI.NewBasis, mVHCI.Frequencies, mVHCI.mol.V0, mVHCI.mol.ints[0].tolist(), mVHCI.mol.ints[1].tolist(), mVHCI.mol.ints[2].tolist(), False)
                HJJ = VCISparseHamNMode(mVHCI.NewBasis, mVHCI.NewBasis, mVHCI.Frequencies, mVHCI.mol.V0, mVHCI.mol.ints[0].tolist(), mVHCI.mol.ints[1].tolist(), mVHCI.mol.ints[2].tolist(), True)
            mVHCI.H = sparse.hstack([mVHCI.H, HIJ])
            mVHCI.H = sparse.vstack([mVHCI.H, sparse.hstack([HIJ.transpose(), HJJ])])
    mVHCI.Timer.stop(1)
    mVHCI.Timer.start(0)
    mVHCI.E, mVHCI.C = sparse.linalg.eigsh(mVHCI.H, k = mVHCI.NStates, which = 'SA')
    mVHCI.Timer.stop(0)
    mVHCI.E_HCI = mVHCI.E[:mVHCI.NStates].copy()

def VCISparseHamTCI(Basis1, Basis2, Frequencies, V0, CoreTensors, OffDiagonal):
    T = VCISparseT(Basis1, Basis2, Frequencies, OffDiagonal)
    N1 = len(Basis1)
    N2 = len(Basis2)
    V = sparse.lil_matrix((N1, N2))
    thr = 1e-12
    if OffDiagonal:
        for i, B1 in enumerate(Basis1):
            for j, B2 in enumerate(Basis2):
                G = CoreTensors[0][:, :, B1.Modes[0].Quanta, B2.Modes[0].Quanta]
                for m in range(Frequencies.shape[0] - 1):
                    G = G @ CoreTensors[m + 1][:, :, B1.Modes[m + 1].Quanta, B2.Modes[m + 1].Quanta]
                Vij = G[0, 0]
                if abs(Vij) > thr:
                    V[i, j] = Vij
    else:
        for i in range(N1):
            for j in range(i, N2):
                G = CoreTensors[0][:, :, Basis1[i].Modes[0].Quanta, Basis2[j].Modes[0].Quanta]
                for m in range(Frequencies.shape[0] - 1):
                    G = G @ CoreTensors[m + 1][:, :, Basis1[i].Modes[m + 1].Quanta, Basis2[j].Modes[m + 1].Quanta]
                Vij = G[0, 0]
                if abs(Vij) > thr:
                    V[i, j] = Vij
                    V[j, i] = Vij
    V = V * constants.AU_TO_INVCM
    return T + V

def SparseDiagonalizeTCI(mVHCI):
    mVHCI.Timer.start(1)
    if mVHCI.H is None:
        mVHCI.H = VCISparseHamTCI(mVHCI.Basis, mVHCI.Basis, mVHCI.Frequencies, mVHCI.mol.V0, mVHCI.mol.core_tensors, False)
    else:
        if len(mVHCI.NewBasis) != 0:
            HIJ = VCISparseHamTCI(mVHCI.Basis[:-len(mVHCI.NewBasis)], mVHCI.NewBasis, mVHCI.Frequencies, mVHCI.mol.V0, mVHCI.mol.core_tensors, True)
            HJJ = VCISparseHamTCI(mVHCI.NewBasis, mVHCI.NewBasis, mVHCI.Frequencies, mVHCI.mol.V0, mVHCI.mol.core_tensors, False)
            mVHCI.H = sparse.hstack([mVHCI.H, HIJ])
            mVHCI.H = sparse.vstack([mVHCI.H, sparse.hstack([HIJ.transpose(), HJJ])])
    mVHCI.Timer.stop(1)
    mVHCI.Timer.start(0)
    mVHCI.E, mVHCI.C = sparse.linalg.eigsh(mVHCI.H, k = mVHCI.NStates, which = 'SA')
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
    
    # We need to translate this into the wavefunction object
    BasisWF = []
    for B in Basis:
        WF = WaveFunction(B, mVHCI.Frequencies)
        BasisWF.append(WF)
    return BasisWF

def InitC(mVHCI):
    from vstr.harmonic.harm_analysis import HarmonicEnergy
    mVHCI.C = np.eye(len(mVHCI.Basis))
    mVHCI.E = []
    ZPE = 0.0
    for w in mVHCI.Frequencies:
        ZPE += 0.5 * w
    for B in mVHCI.Basis:
        b = []
        for i in range(len(B.Modes)):
            b.append(B.Modes[i].Quanta)
        mVHCI.E.append(HarmonicEnergy(mVHCI.Frequencies, b, ZPE = ZPE))

def ExpectedQ(mVHCI, State = 0):
    MaxQ = [mVHCI.MaxQuanta[0]]
    Freq = [1.0]
    GenericBasis = init_funcs.InitGridBasis(Freq, MaxQ)[0]

    ExpQ = []
    NModes = mVHCI.Frequencies.shape[0]
    for i in range(NModes):
        W = FConst(1.0, [i], False)
        W = [W]
        W.append(FConst(0.0, [0] * 6, False))
        CubicFC = []
        QuarticFC = []
        QuinticFC = []
        SexticFC = []
        CubicFC = W.copy()
        Q = GenerateSparseHamV(mVHCI.Basis, mVHCI.Frequencies, W, CubicFC, QuarticFC, QuinticFC, SexticFC)
        QBar = mVHCI.C[:, State].T @ Q @ mVHCI.C[:, State]
        ExpQ.append(QBar)
    return np.array(ExpQ)

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
def LCLine(mVHCI, n, thr = 2.5e-1, C = None):
    if C is None:
        C = mVHCI.C
    def BasisToString(B):
        BString = '|'
        for HO in B.Modes:
            BString += str(HO.Quanta)
        BString += '>'
        return BString
    
    LC = ''
    for i in range(C.shape[0]):
        if abs(C[i, n]) > thr:
            LC += str(C[i, n])
            LC += BasisToString(mVHCI.Basis[i])
            LC += ' + '
    return LC[:-3]
    

def PrintResults(mVHCI, thr = 2.5e-1):
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
        LCString = mVHCI.LCLine(n, thr = thr)
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
class VHCI:
    FormW = FormW
    FormWSD = FormWSD
    MakeAnharmTensor = MakeAnharmTensor
    HCI = HCI
    Diagonalize = Diagonalize
    SparseDiagonalize = SparseDiagonalize
    HCIStep = HCIStep
    ScreenBasis = ScreenBasis
    PT2 = PT2
    InitTruncatedBasis = InitTruncatedBasis
    InitC = InitC
    ExpectedQ = ExpectedQ
    PrintResults = PrintResults
    LCLine = LCLine
    PrintParameters = PrintParameters
    ReadBasisFromFile = ReadBasisFromFile
    SaveBasisToFile = SaveBasisToFile

    def __init__(self, Frequencies, UnscaledPotential, MaxQuanta = 2, MaxTotalQuanta = 2, NStates = 10, **kwargs):
        self.Frequencies = Frequencies # 1D array of all harmonic frequencies.
        self.NModes = Frequencies.shape[0]
        self.UnscaledPotential = UnscaledPotential
        self.Potential = [[]] * 4 # Cubic, quartic, quintic, sextic
        for V in self.UnscaledPotential:
            Wp = self.FormW(V)
            self.Potential[Wp[0].Order - 3] = Wp
        self.FormWSD()

        self.PotentialListFull = []
        for Wp in self.PotentialSD:
            self.PotentialListFull += Wp
        self.PotentialList = []
        for Wp in self.Potential:
            self.PotentialList += Wp
            self.PotentialListFull += Wp
        self.PotentialListFull = HeatBath_Sort_FC(self.PotentialListFull) # Only need to sort these once
        self.MaxQuanta = MaxQuanta
        if isinstance(MaxQuanta, int):
            self.MaxQuanta = [MaxQuanta] * self.NModes
        
        self.MaxTotalQuanta = MaxTotalQuanta
        self._HighestQuanta = [MaxTotalQuanta] * self.NModes
        self.H = None
        self.eps1 = 0.1 # HB epsilon
        self.eps2 = 0.01 # PT2/SPT2 epsilon
        self.eps3 = -1.0 # SSPT2 epsilon, < 0 means do not do semi-stochastic
        self.tol = 0.01
        self.MaxIter = 1000
        self.NStates = NStates
        self.NStatesPT2 = NStates
        self.NWalkers = 200
        self.NSamples = 50
        self.dE_PT2 = None
        self.sE_PT2 = None
        self.HBMethod = 'exact' #['orig', 'max', 'exact']

        self.CHKFile = None
        self.ReadFromFile = False
        self.SaveToFile = False
        self.PrintHCISteps = False

        self.__dict__.update(kwargs)

        self.Timer = TIMER(6)
        self.TimerNames = ['Diagonalize', 'Form Hamiltonian', 'Screen Basis', 'PT2 correction', 'SPT2 correction', 'SSPT2 correction']

    def kernel(self, doVCI = True, doVHCI = True, doPT2 = False, doSPT2 = False, ComparePT2 = False):
        assert(self.HBMethod in ['orig', 'max', 'exact'])
        if self.HBMethod == 'exact':
            self.Ys = [self.MakeAnharmTensor()] * self.NModes
            #Is = []
            #for n in self.MaxQuanta:
            #    Is.append(np.eye(n))
            #self.Xs = ContractedHOTerms(Is, self.Frequencies)

        if self.SaveToFile or self.ReadFromFile:
            assert(self.CHKFile is not None)

        if self.ReadFromFile:
            self.ReadBasisFromFile(self.CHKFile)
        else:
            self.Basis = init_funcs.InitTruncatedBasis(self.NModes, self.Frequencies, self.MaxQuanta, MaxTotalQuanta = self.MaxTotalQuanta)

        self.PrintParameters()

        if doVCI:
            self.Timer.start(0)
            self.SparseDiagonalize()
            print("===== VCI RESULTS =====", flush = True)
            print("Initial Basis Size:", len(self.Basis))
            self.PrintResults()
            print("")
            self.Timer.stop(0)
        
        if doVHCI:
            print("===== VHCI RESULTS =====", flush = True)
            self.HCI()
            self.PrintResults()
            print("")

        if doPT2:
            print("===== VHCI+PT2 RESULTS =====", flush = True)
            self.PT2(doStochastic = False)
            self.PrintResults()
            print("")
        if doSPT2:
            print("===== VHCI+SPT2 RESULTS =====", flush = True)
            self.PT2(doStochastic = True)
            self.PrintResults()
            print("")
        if ComparePT2:
            print("===== VHCI+SSPT2 RESULTS =====", flush = True)
            eps = self.eps2
            self.eps2 *= 10
            self.eps3 = eps
            self.PrintParameters()
            self.PT2(doStochastic = True)
            self.PrintResults()
            print("")
        self.Timer.report(self.TimerNames)


    @property
    def SummedQuanta(self):
        SummedQuanta = [0] * self.NModes
        for B in self.Basis:
            for n in range(self.NModes):
                SummedQuanta[n] += B.Modes[n].Quanta
        return SummedQuanta

    @property
    def AverageQuanta(self):
        SummedQuanta = self.SummedQuanta
        return [Q / len(self.Basis) for Q in SummedQuanta]

    @property
    def HighestQuanta(self):
        try:
            for B in self.NewBasis:
                for n in range(self.NModes):
                    if B.Modes[n].Quanta > self._HighestQuanta[n]:
                        self._HighestQuanta[n] = B.Modes[n].Quanta
        except:
            return self._HighestQuanta
        return self._HighestQuanta
            
class NModeVHCI(VHCI):
    SparseDiagonalize = SparseDiagonalizeNMode
    
    def __init__(self, mol, NStates = 10, **kwargs):
        self.mol = mol
        self.Frequencies = mol.Frequencies # 1D array of all harmonic frequencies.
        self.NModes = self.Frequencies.shape[0]
        self.MaxQuanta = [mol.ngridpts] * self.NModes
        
        self.MaxTotalQuanta = 5
        for m in self.MaxQuanta:
            assert(m <= mol.ngridpts)
        self.H = None
        self.eps1 = 0.1 # HB epsilon
        self.eps2 = 0.01 # PT2/SPT2 epsilon
        self.eps3 = -1.0 # SSPT2 epsilon, < 0 means do not do semi-stochastic
        self.tol = 0.01
        self.MaxIter = 1000
        self.NStates = NStates
        self.NStatesPT2 = NStates
        self.NWalkers = 200
        self.NSamples = 50
        self.dE_PT2 = None
        self.sE_PT2 = None
        self.HBMethod = 'qff' #['qff', '2mode']

        self.CHKFile = None
        self.ReadFromFile = False
        self.SaveToFile = False
        self.PrintHCISteps = False

        self.__dict__.update(kwargs)

        self.Timer = TIMER(6)
        self.TimerNames = ['Diagonalize', 'Form Hamiltonian', 'Screen Basis', 'PT2 correction', 'SPT2 correction', 'SSPT2 correction']

    def kernel(self, doVCI = True, doVHCI = True, doPT2 = False, doSPT2 = False, ComparePT2 = False):
        if self.HBMethod.upper() == 'QFF':
            V3, V4 = self.mol.nm.get_ff()
            self.Potential = [[]] * 4 # Cubic, quartic, quintic, sextic
            for V in [V3, V4]:
                Wp = self.FormW(V)
                self.Potential[Wp[0].Order - 3] = Wp
            self.FormWSD()

            self.PotentialListFull = []
            for Wp in self.PotentialSD:
                self.PotentialListFull += Wp
            self.PotentialList = []
            for Wp in self.Potential:
                self.PotentialList += Wp
                self.PotentialListFull += Wp
            self.PotentialListFull = HeatBath_Sort_FC(self.PotentialListFull) # Only need to sort these once
            self.Ys = [self.MakeAnharmTensor()] * self.NModes
        else:
            self.PotentialListFull = []

        if self.HBMethod.upper() == '2MODE':
            self.Sorted2Mode = np.empty((self.mol.Nm, self.mol.Nm, self.mol.ngridpts, self.mol.ngridpts, self.mol.ngridpts**2, 2), dtype = np.int8)
            for i in range(self.mol.Nm):
                for j in range(self.mol.Nm):
                    for ni in range(self.mol.ngridpts):
                        for nj in range(self.mol.ngridpts):
                            Sorted = np.argsort(-abs(self.mol.ints[1][i, j][ni, nj].reshape(-1)))
                            Sorted = np.unravel_index(Sorted, (self.mol.ngridpts, self.mol.ngridpts))
                            Sorted = np.vstack((Sorted[0], Sorted[1]))
                            self.Sorted2Mode[i, j, ni, nj] = Sorted.T

        if self.SaveToFile or self.ReadFromFile:
            assert(self.CHKFile is not None)

        if self.ReadFromFile:
            self.ReadBasisFromFile(self.CHKFile)
        else:
            self.Basis = init_funcs.InitTruncatedBasis(self.NModes, self.Frequencies, self.MaxQuanta, MaxTotalQuanta = self.MaxTotalQuanta)

        self.PrintParameters()

        if doVCI:
            self.Timer.start(0)
            self.SparseDiagonalize()
            print("===== VCI RESULTS =====", flush = True)
            print("Initial Basis Size:", len(self.Basis))
            self.PrintResults()
            print("")
            self.Timer.stop(0)
        
        if doVHCI:
            print("===== VHCI RESULTS =====", flush = True)
            self.HCI()
            self.PrintResults()
            print("")

        if doPT2:
            print("===== VHCI+PT2 RESULTS =====", flush = True)
            self.PT2(doStochastic = False)
            self.PrintResults()
            print("")
        if doSPT2:
            print("===== VHCI+SPT2 RESULTS =====", flush = True)
            self.PT2(doStochastic = True)
            self.PrintResults()
            print("")
        if ComparePT2:
            print("===== VHCI+SSPT2 RESULTS =====", flush = True)
            eps = self.eps2
            self.eps2 *= 10
            self.eps3 = eps
            self.PrintParameters()
            self.PT2(doStochastic = True)
            self.PrintResults()
            print("")
        self.Timer.report(self.TimerNames)

class TCIVHCI(VHCI):
    SparseDiagonalize = SparseDiagonalizeTCI
    def __init__(self, mol, NStates = 10, **kwargs):
        self.mol = mol
        self.Frequencies = mol.Frequencies # 1D array of all harmonic frequencies.
        self.NModes = self.Frequencies.shape[0]
        self.MaxQuanta = [mol.ngridpts] * self.NModes
        
        self.MaxTotalQuanta = 5
        for m in self.MaxQuanta:
            assert(m <= mol.ngridpts)
        self.H = None
        self.eps1 = 0.1 # HB epsilon
        self.eps2 = 0.01 # PT2/SPT2 epsilon
        self.eps3 = -1.0 # SSPT2 epsilon, < 0 means do not do semi-stochastic
        self.tol = 0.01
        self.MaxIter = 1000
        self.NStates = NStates
        self.NStatesPT2 = NStates
        self.NWalkers = 200
        self.NSamples = 50
        self.dE_PT2 = None
        self.sE_PT2 = None
        self.HBMethod = 'pass' #['qff', '2mode']

        self.CHKFile = None
        self.ReadFromFile = False
        self.SaveToFile = False
        self.PrintHCISteps = False

        self.__dict__.update(kwargs)

        self.Timer = TIMER(6)
        self.TimerNames = ['Diagonalize', 'Form Hamiltonian', 'Screen Basis', 'PT2 correction', 'SPT2 correction', 'SSPT2 correction']
    
    def kernel(self, doVCI = True, doVHCI = True, doPT2 = False, doSPT2 = False, ComparePT2 = False):
        if self.HBMethod.upper() == 'MAXTENSOR':
            self.MaxTensors = [abs(G.reshape(G.shape[0] * G.shape[1], G.shape[2], G.shape[3])).max(axis = 0) for G in self.mol.core_tensors]
            self.MSortedIndices = [np.argsort(-abs(M), axis = 0) for M in self.MaxTensors]

        if self.SaveToFile or self.ReadFromFile:
            assert(self.CHKFile is not None)

        if self.ReadFromFile:
            self.ReadBasisFromFile(self.CHKFile)
        else:
            self.Basis = init_funcs.InitTruncatedBasis(self.NModes, self.Frequencies, self.MaxQuanta, MaxTotalQuanta = self.MaxTotalQuanta)

        self.PrintParameters()

        if doVCI:
            self.Timer.start(0)
            self.SparseDiagonalize()
            print("===== VCI RESULTS =====", flush = True)
            print("Initial Basis Size:", len(self.Basis))
            self.PrintResults()
            print("")
            self.Timer.stop(0)
        
        if doVHCI:
            print("===== VHCI RESULTS =====", flush = True)
            self.HCI()
            self.PrintResults()
            print("")

        if doPT2:
            print("===== VHCI+PT2 RESULTS =====", flush = True)
            self.PT2(doStochastic = False)
            self.PrintResults()
            print("")
        if doSPT2:
            print("===== VHCI+SPT2 RESULTS =====", flush = True)
            self.PT2(doStochastic = True)
            self.PrintResults()
            print("")
        if ComparePT2:
            print("===== VHCI+SSPT2 RESULTS =====", flush = True)
            eps = self.eps2
            self.eps2 *= 10
            self.eps3 = eps
            self.PrintParameters()
            self.PT2(doStochastic = True)
            self.PrintResults()
            print("")
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
    w, MaxQuanta, MaxTotalQuanta, Vs, eps1, eps2, eps3, NWalkers, NSamples, NStates = Read('CLO2.inp')
    mVHCI = VHCI(np.asarray(w), Vs, MaxQuanta = MaxQuanta, MaxTotalQuanta = MaxTotalQuanta, eps1 = eps1, eps2 = eps2, eps3 = eps3, NWalkers = NWalkers, NSamples = NSamples, NStates = NStates, NStatesPT2 = 2)
    mVHCI.kernel(doPT2 = True, ComparePT2 = True)
    #print(mVHCI.E[:NStates])
    #mVHCI.PT2(doStochastic = True)
    #print(mVHCI.E[:NStates])
    #print(mVHCI.dE_PT2)
    #print(mVHCI.E_HCI_PT2)
