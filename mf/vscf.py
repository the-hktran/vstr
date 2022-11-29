import numpy as np
from vstr.utils.init_funcs import FormW, InitTruncatedBasis, InitGridBasis
from vstr.cpp_wrappers.vhci_jf.vhci_jf_functions import WaveFunction, GenerateHam0V, GenerateSparseHamAnharmV, GetVEffCPP, MakeCTensorCPP, GetVEffFASTCPP
from vstr.utils.perf_utils import TIMER

def InitCs(mVSCF):
    Cs = []
    for Mode in range(mVSCF.NModes):
        C = np.eye(mVSCF.MaxQuanta[Mode])
        Cs.append(C)
    return Cs

'''
Gets the index of one modal basis function tensor product with the basis set of the other modals
i \otimes J \otimes K \otimes L ...
'''
def GetModalSlices(mVSCF):
    # In this simple implementation, this simply means all the basis functions with i in the
    # place of the mode we are interested in.
    ModalSlices = []
    for N in range(mVSCF.NModes):
        BasisByMode = []
        for i in range(mVSCF.MaxQuanta[N]):
            BasisByModeByIndex = []
            for j, B in enumerate(mVSCF.Basis):
                if B.Modes[N].Quanta == i:
                    BasisByModeByIndex.append(j)
            BasisByMode.append(BasisByModeByIndex)
        ModalSlices.append(BasisByMode)
    return ModalSlices

'''
This function collects all the nonzero anharmonic terms, for each anharmonic force constant
'''
def MakeAnharmTensorSLOW(mVSCF):
    AnharmTensor = []
    for W in mVSCF.PotentialList:
        CubicFC = []
        QuarticFC = []
        QuinticFC = []
        SexticFC = []
        if W.Order == 3:
            CubicFC.append(W)
        elif W.Order == 4:
            QuarticFC.append(W)
        elif W.Order == 5:
            QuinticFC.append(W)
        elif W.Order == 6:
            SexticFC.append(W)
        H = GenerateSparseHamAnharmV(mVSCF.Basis, mVSCF.Frequencies, mVSCF.PotentialList, CubicFC, QuarticFC, QuinticFC, SexticFC)
        AnharmTensor.append(H)
    return AnharmTensor

def MakeAnharmTensor(mVSCF):
    AnharmTensor = []
    RestrictedBases = []
    QUniques = []
    for W in mVSCF.PotentialList:
        # First we need to make a basis restricted only to relevant modes. These need to be saved to generate the new C tensors
        mVSCF.Timer.start(0)
        RestrictedMaxQ = [1] * mVSCF.NModes
        QUniques.append(W.QUnique)
        for i in W.QUnique:
            RestrictedMaxQ[i] = mVSCF.MaxQuanta[i]
        RestrictedBasis = InitGridBasis(mVSCF.Frequencies, RestrictedMaxQ)[0]
        RestrictedBases.append(RestrictedBasis)
        mVSCF.Timer.stop(0)
        mVSCF.Timer.start(1)
        CubicFC = []
        QuarticFC = []
        QuinticFC = []
        SexticFC = []
        if W.Order == 3:
            CubicFC.append(W)
        elif W.Order == 4:
            QuarticFC.append(W)
        elif W.Order == 5:
            QuinticFC.append(W)
        elif W.Order == 6:
            SexticFC.append(W)
        H = GenerateSparseHamAnharmV(RestrictedBasis, mVSCF.Frequencies, mVSCF.PotentialList, CubicFC, QuarticFC, QuinticFC, SexticFC)
        AnharmTensor.append(H)
        mVSCF.Timer.stop(1)
    return AnharmTensor, RestrictedBases, QUniques

'''
def GetVEffByMode(mVSCF, Mode):
    VEff = np.zeros((mVSCF.MaxQuanta[Mode], mVSCF.MaxQuanta[Mode]))
    for n in range(VEff.shape[0]):
        for m in range(n, VEff.shape[1]):
            Vnm = 0.0
            for H in mVSCF.AnharmTensor:
                I, J = H.nonzero()
                for i in I:
                    if i not in mVSCF.ModalSlices[Mode][n]:
                        continue
                    BraModeBasis = []
                    for WF in mVSCF.Basis[i].Modes:
                        BraModeBasis.append(WF.Quanta)
                    for j in J:
                        if j not in VSCF.ModalSlices[Mode][m]:
                            continue
                        # Now we know we have a nonzero element. We need to associate this with a set of modal basis functions.
                        KetModeBasis = []
                        for WF in mVSCF.Basis[j].Modes:
                            KetModeBasis.append(WF.Quanta)
                        VContracted = H[i, j]
                        for a, A in enumerate(BraModeBasis):
                            if a == Mode:
                                continue
                            VContracted *= mVSCF.CByMode[a][A]
                        for b, B in enumerate(KetModeBasis):
                            if b == Mode:
                                continue
                            VContracted *= mVSCF.CByMode[b][B]
                        Vnm += VContracted
            VEff[n, m] = Vnm
            VEff[m, n] = Vnm
    return VEff

def GetVEff(mVSCF):
    VEff = []
    for N in range(mVSCF.NModes):
        VEffN = mVSCF.GetVEffByMode(N)
        VEff.append(VEffN)
    return VEff
'''

'''
Forms the Harmonic Oscillator Hamiltonian for each mode. This plays the role of
the one electron integrals in electronic structure. This is stored and is later
used to form the contracted "h" part of the Fock matrix.
'''
def MakeHOHam(mVSCF):
    HamHO = []
    for Mode in range(mVSCF.NModes):
        hM = np.zeros((mVSCF.MaxQuanta[Mode], mVSCF.MaxQuanta[Mode]))
        for i in range(mVSCF.MaxQuanta[Mode]):
            hM[i, i] = (i + 0.5) * mVSCF.Frequencies[Mode]
        HamHO.append(hM)
    return HamHO

'''
Takes a stored HO Hamiltonian and contracts it into the modal basis
'''
def GetHCore(mVSCF, Cs = None, HamHO = None):
    if Cs is None:
        Cs = mVSCF.Cs
    if HamHO is None:
        HamHO = mVSCF.HamHO
    hs = []
    for Mode in range(mVSCF.NModes):
        hs.append(Cs[Mode].T @ HamHO[Mode] @ Cs[Mode])
    return hs

def CalcESCF(mVSCF, ModeOcc = None, V0 = None):
    if ModeOcc is None:
        ModeOcc = mVSCF.ModeOcc
    if V0 is None:
        if mVSCF.SlowV:
            V0 = mVSCF.GetVEffSLOW(ModeOcc = ModeOcc, FirstV = True)[0]
        else:
            V0 = mVSCF.GetVEff(ModeOcc = ModeOcc, FirstV = True)[0]
    
    DC = (mVSCF.Cs[0][:, ModeOcc[0]].T @ V0 @ mVSCF.Cs[0][:, ModeOcc[0]])
    E_SCF = 0.0
    for Mode, E in enumerate(mVSCF.Es):
        E_SCF += E[ModeOcc[Mode]]
    return E_SCF - (mVSCF.NModes - 1) * DC

'''
Takes stored anharmonic Hamiltonians and contracts it into the modal basis
'''
def GetVEffSLOW(mVSCF, ModeOcc = None, FirstV = False):
    if ModeOcc is None:
        ModeOcc = mVSCF.ModeOcc
    return GetVEffCPP(mVSCF.AnharmTensor, mVSCF.Basis, mVSCF.Cs, mVSCF.ModalSlices, mVSCF.MaxQuanta, ModeOcc, FirstV)

def GetRestrictedSlice(RestrictedBasis, Ind, Mode):
    RSlice = []
    for i in range(len(RestrictedBasis)):
        if RestrictedBasis[i].Modes[Mode].Quanta == Ind:
            RSlice.append(i)
    return RSlice

'''
def GetVEff(mVSCF, ModeOcc = None, FirstV = False):
    if ModeOcc is None:
        ModeOcc = mVSCF.ModeOcc
    VEff = []
    NModes = mVSCF.NModes
    if FirstV:
        NModes = 1
    for m in range(mVSCF.NModes):
        Vm = np.zeros((mVSCF.MaxQuanta[m], mVSCF.MaxQuanta[m]))
        for q in range(len(mVSCF.AnharmTensor)):
            Vmq = np.zeros_like(Vm)
            if m not in mVSCF.QUniques[q]:
                CTensor = MakeCTensorCPP(mVSCF.Cs, mVSCF.QUniques[q], [ModeOcc], mVSCF.RestrictedBases[q])
                Vmq = np.eye(mVSCF.MaxQuanta[m]) * (CTensor.T @ mVSCF.AnharmTensor[q].dot(CTensor))[0, 0]
            else:
                Qs = mVSCF.QUniques[q].copy()
                Qs.remove(m)
                for i in range(Vm.shape[0]):
                    RSliceI = GetRestrictedSlice(mVSCF.RestrictedBases[q], i, m)
                    RBasisI = [B for i, B in enumerate(mVSCF.RestrictedBases[q]) if i in RSliceI]
                    CTensorI = MakeCTensorCPP(mVSCF.Cs, Qs, [ModeOcc], RBasisI)
                    for j in range(i, Vm.shape[1]):
                        if j == i:
                            RSliceJ = RSliceI.copy()
                        else:
                            RSliceJ = GetRestrictedSlice(mVSCF.RestrictedBases[q], j, m)
                        VReduced = mVSCF.AnharmTensor[q][np.ix_(RSliceI, RSliceJ)]
                        RBasisJ = [B for i, B in enumerate(mVSCF.RestrictedBases[q]) if i in RSliceJ]
                        CTensorJ = MakeCTensorCPP(mVSCF.Cs, Qs, [ModeOcc], RBasisJ)
                        Vmq[i, j] = (CTensorI.T @ VReduced.dot(CTensorJ))[0, 0]
                        Vmq[j, i] = Vmq[i, j]
            Vm += Vmq
        VEff.append(Vm)
    return VEff
'''

def GetVEff(mVSCF, ModeOcc = None, FirstV = False):
    if ModeOcc is None:
        ModeOcc = mVSCF.ModeOcc
    return GetVEffFASTCPP(mVSCF.AnharmTensor, mVSCF.RestrictedBases, mVSCF.QUniques, mVSCF.Cs, mVSCF.MaxQuanta, ModeOcc, FirstV)

def GetFock(mVSCF, hs = None, Cs = None, CalcE = False):
    if Cs is None:
        Cs = mVSCF.Cs
    if hs is None:
        hs = mVSCF.HamHO

    if mVSCF.SlowV:
        Vs = mVSCF.GetVEffSLOW()
    else:
        Vs = mVSCF.GetVEff()

    # Save the double counting term while we have the effective potentials
    if CalcE:
        mVSCF.ESCF = mVSCF.CalcESCF(ModeOcc = mVSCF.ModeOcc, V0 = Vs[0])
    Fs = []
    for Mode in range(mVSCF.NModes):
        Fs.append(hs[Mode] + Vs[Mode])
    return Fs

def FockError(F, C):
    return C @ C.T @ F - F @ C @ C.T

def StoreFock(mVSCF):
    # First, update the list of Fock and Error matrices
    mVSCF.AllFs.append(mVSCF.Fs)
    Errs = []
    for Mode in range(mVSCF.NModes):
        Err = FockError(mVSCF.Fs[Mode], mVSCF.Cs[Mode])
        Errs.append(Err)
    mVSCF.AllErrs.append(Errs)

    # Restrict the size of the space
    if len(mVSCF.AllFs) > mVSCF.DIISSpace:
        mVSCF.AllFs = mVSCF.AllFs[-mVSCF.DIISSpace:]
        mVSCF.AllErrs = mVSCF.AllErrs[-mVSCF.DIISSpace:]
                
def DIISUpdate(mVSCF):
    for Mode in range(mVSCF.NModes):
        B = np.ones((mVSCF.DIISSpace + 1, mVSCF.DIISSpace + 1))
        for i in range(mVSCF.DIISSpace):
            for j in range(mVSCF.DIISSpace):
                B[i, j] = (mVSCF.AllErrs[i][Mode] * mVSCF.AllErrs[j][Mode]).sum()
                B[j, i] = B[i, j]
        B[-1, -1] = 0.0
        x = np.zeros((mVSCF.DIISSpace + 1, 1))
        x[-1, 0] = 1.0
        Coeff = np.linalg.solve(B, x)
        NewF = np.zeros(mVSCF.Fs[Mode].shape)
        for i in range(mVSCF.DIISSpace):
            NewF += Coeff[i] * mVSCF.AllFs[i][Mode]
        mVSCF.Fs[Mode] = NewF 
        

def SCFIteration(mVSCF, It, DoDIIS = True):
    mVSCF.Timer.start(2)
    mVSCF.Fs = mVSCF.GetFock(CalcE = True)
    mVSCF.Timer.stop(2)
    mVSCF.Timer.start(3)
    if DoDIIS:
        mVSCF.StoreFock()
        if It > mVSCF.DIISStart:
            mVSCF.DIISUpdate() # Replaces Fock matrices with DIIS updated fock matrices
    mVSCF.Timer.stop(3)
    mVSCF.Timer.start(4)
    COld = mVSCF.Cs.copy()
    mVSCF.Cs = []
    mVSCF.Es = []
    for F in mVSCF.Fs:
        E, C = np.linalg.eigh(F)
        mVSCF.Es.append(E)
        mVSCF.Cs.append(C)
    SCFErr = 0.0
    for Mode in range(mVSCF.NModes):
        SCFErr = ((abs(COld[Mode]) - abs(mVSCF.Cs[Mode]))**2).sum()
    mVSCF.Timer.stop(4)
    return SCFErr

def SCF(mVSCF, DoDIIS = True, tol = 1e-8, etol = 1e-6):
    if DoDIIS:
        mVSCF.AllFs = []
        mVSCF.AllErrs = []

    mVSCF.Es = []
    for HHO in mVSCF.HamHO:
        mVSCF.Es.append(np.diag(HHO))
    ConvErr = 1
    It = 1
    EnergyErr = 1
    while(ConvErr > tol or EnergyErr > etol):
        EnergyErr = mVSCF.ESCF
        ConvErr = mVSCF.SCFIteration(It, DoDIIS = DoDIIS)
        EnergyErr = abs(EnergyErr - mVSCF.ESCF)
        print("VSCF Iteration %d complete with an SCF error of %.12f/%.12f and SCF Energy of %.6f" % (It, ConvErr, EnergyErr, mVSCF.ESCF))
        It += 1
        if It > mVSCF.MaxIterations:
            raise RuntimeError("Maximum number of SCF iterations reached without convergence.")
    mVSCF.Timer.report(mVSCF.TimerNames)

def LCLine(mVSCF, ModeOcc, thr = 1e-2):
    def BasisToString(B):
        BString = '|'
        for Q in B:
            BString += str(Q)
        BString += '>'
        return BString

    LC = ''
    for B in mVSCF.BasisList:
        Coeff = 1.0
        for i, n in enumerate(B):
            Coeff *= mVSCF.Cs[i][n, ModeOcc[i]]
            if abs(Coeff) < thr:
                break
        if abs(Coeff) > thr:
            LC += str(Coeff) + BasisToString(B) + ' + '
    return LC[:-3]

def LowestStates(mVSCF, NStates, MaxQuanta = None):
    LStates = []
    LStates.append([0] * mVSCF.NModes)
    EMax = 1e10 #np.sum([mVSCF.MaxQuanta[i] * e for i, e in enumerate(mVSCF.Frequencies)])
    for n in range(NStates):
        NextMode = 0
        EOld = EMax.copy()
        for B in LStates:
            for m in range(mVSCF.NModes):
                if MaxQuanta is not None and B[m] >= MaxQuanta[m] - 1
                BTestIncr = B.copy()
                BTestIncr[m] += 1
                if BTestIncr in LStates:
                    continue
                E = np.sum([BTestIncr[i] * e for i, e in enumerate(mVSCF.Frequencies)])
                if E < EOld:
                    NextMode = m
                    BSave = B.copy()
                    EOld = E.copy()
        BSave[NextMode] += 1
        LStates.append(BSave)
    return LStates


def PrintResults(mVSCF, NStates = None):
    FinalE = []
    if NStates is None:
        NStates = np.prod(mVSCF.MaxQuanta)
    EnergyBList = mVSCF.LowestStates(NStates)
    for B in  EnergyBList:
        FinalE.append(mVSCF.CalcESCF(ModeOcc = B))
    SortedInd = np.argsort(np.asarray(FinalE))

    mVSCF.BasisList = InitGridBasis(mVSCF.Frequencies, mVSCF.MaxQuanta, ListOnly = True)
    
    for i in SortedInd:
        OutLine = '{:.8f}\t'.format(FinalE[i])
        ModeLabel = ""
        for j, n in enumerate(EnergyBList[i]):
            if n > 1:
                ModeLabel += str(n)
            if n > 0:
                ModeLabel += 'w{} + '.format(j + 1)
        ModeLabel = ModeLabel[:-3]
        OutLine += ModeLabel
        LCString = mVSCF.LCLine(EnergyBList[i])
        OutLine += '\t%s' % (LCString)
        print(OutLine)

class VSCF:
    InitCs = InitCs
    GetModalSlices = GetModalSlices
    MakeAnharmTensor = MakeAnharmTensor
    MakeAnharmTensorSLOW = MakeAnharmTensorSLOW
    MakeHOHam = MakeHOHam
   
    GetHCore = GetHCore
    GetVEff = GetVEff
    GetVEffSLOW = GetVEffSLOW
    GetFock = GetFock
    StoreFock = StoreFock
    DIISUpdate = DIISUpdate
    SCF = SCF
    SCFIteration = SCFIteration
    LowestStates = LowestStates
    CalcESCF = CalcESCF

    PrintResults = PrintResults
    LCLine = LCLine
    def __init__(self, Frequencies, UnscaledPotential, MaxQuanta = 2, NStates = 10, SlowV = False, **kwargs):
        self.Frequencies = Frequencies
        self.NModes = self.Frequencies.shape[0]

        self.Timer = TIMER(5)
        self.Timer.set_overhead(self.Timer.estimate_overhead())
        self.TimerNames = ["Basis Init", "Anharm Pot Init", "Fock Generation", "DIIS", "Solve Fock"]

        self.Potential = [[]] * 4
        for V in UnscaledPotential:
            Wp = FormW(V)
            self.Potential[Wp[0].Order - 3] = Wp
        self.PotentialList = []
        for Wp in self.Potential:
            self.PotentialList += Wp

        if isinstance(MaxQuanta, int):
            self.MaxQuanta = [MaxQuanta] * self.NModes
        else:
            self.MaxQuanta = MaxQuanta

        self.SlowV = SlowV
        if self.SlowV:
            self.Timer.start(0)
            self.Basis, self.BasisList = InitGridBasis(self.Frequencies, MaxQuanta)
            self.ModalSlices = self.GetModalSlices()
            self.Timer.stop(0)
        self.HamHO = self.MakeHOHam()
        if self.SlowV:
            self.Timer.start(1)
            self.AnharmTensor = self.MakeAnharmTensorSLOW()
            self.Timer.stop(1)
        else:
            self.AnharmTensor, self.RestrictedBases, self.QUniques = self.MakeAnharmTensor()
        self.Cs = self.InitCs()
        self.ModeOcc = [0] * self.NModes
        self.ESCF = 0.0

        self.DoDIIS = False 
        self.DIISSpace = 5
        self.DIISStart = 10
        self.MaxIterations = 1000
        self.__dict__.update(kwargs)

if __name__ == "__main__":
    from vstr.utils.read_jf_input import Read
    w, MaxQuanta, MaxTotalQuanta, Vs, eps1, eps2, eps3, NWalkers, NSamples, NStates = Read('CLO2.inp')
    mf = VSCF(w, Vs, MaxQuanta = MaxQuanta, NStates = NStates, SlowV = False, ModeOcc = [0, 0, 0])
    mf.SCF(DoDIIS = False)
    print(mf.CalcESCF())
    mf.PrintResults(NStates = 10)
