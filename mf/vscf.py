import numpy as np
from vstr.utils.init_funcs import FormW, InitTruncatedBasis, InitGridBasis
from vstr.cpp_wrappers.vhci_jf.vhci_jf_functions import WaveFunction, FConst, GenerateHam0V, GenerateSparseHamAnharmV, GenerateHamAnharmV, GetVEffCPP
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

def MakeAnharmTensor(mVSCF, PotentialList = None):
    if PotentialList is None:
        PotentialList = mVSCF.PotentialList
    AnharmTensor = []
    FCs = []
    QUniques = []
    QPowers = []
        
    mVSCF.Timer.start(0)
    MaxQ = [mVSCF.MaxQuanta[0]]
    Freq = [1.0]
    GenericBasis = InitGridBasis(Freq, MaxQ)[0]
    mVSCF.Timer.stop(0)

    mVSCF.Timer.start(1)
    for W in PotentialList:
        FCs.append(W.fc)
        QUniques.append(W.QUnique)
        QPowers.append(W.QPowers)
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
        H = GenerateSparseHamAnharmV(GenericBasis, Freq, W, CubicFC, QuarticFC, QuinticFC, SexticFC)
        AnharmTensor.append(H)
    AnharmTensor[0] = AnharmTensor[1]
    mVSCF.Timer.stop(1)
    return AnharmTensor, FCs, QUniques, QPowers

def UpdateFC(mVSCF, Frequencies, PotentialList):
    FCs = []
    QUniques = []
    QPowers = []
    
    mVSCF.Potential = [[], [], [], []]
    for W in PotentialList:
        FCs.append(W.fc)
        QUniques.append(W.QUnique)
        QPowers.append(W.QPowers)
        if W.Order == 3:
            mVSCF.Potential[0].append(W)
        elif W.Order == 2 or W.Order == 4:
            mVSCF.Potential[1].append(W)
        elif W.Order == 5:
            mVSCF.Potential[2].append(W)
        elif W.Order == 6:
            mVSCF.Potential[3].append(W)

    mVSCF.FCs = FCs
    mVSCF.QUniques = QUniques
    mVSCF.QPowers = QPowers
    mVSCF.Frequencies = Frequencies
    mVSCF.HamHO = mVSCF.MakeHOHam()

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
    mVSCF.Es = []
    for Mode in range(mVSCF.NModes):
        E = []
        hM = np.zeros((mVSCF.MaxQuanta[Mode], mVSCF.MaxQuanta[Mode]))
        for i in range(mVSCF.MaxQuanta[Mode]):
            hM[i, i] = (i + 0.5) * mVSCF.Frequencies[Mode]
            E.append(hM[i, i])
        HamHO.append(hM)
        mVSCF.Es.append(E)
    return HamHO

'''
Takes a stored HO Hamiltonian and contracts it into the modal basis
'''
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
'''
def CalcESCF(mVSCF, ModeOcc = None, V0 = None):
    if ModeOcc is None:
        ModeOcc = mVSCF.ModeOcc
    if V0 is None:
        V0 = mVSCF.GetVEff(ModeOcc1 = ModeOcc, FirstV = True)[0]
    DC = (mVSCF.Cs[0][:, ModeOcc[0]].T @ V0 @ mVSCF.Cs[0][:, ModeOcc[0]])
    E_SCF = 0.0
    for Mode, E in enumerate(mVSCF.Es):
        E_SCF += E[ModeOcc[Mode]]
    return E_SCF - (mVSCF.NModes - 1) * DC

'''
Takes stored anharmonic Hamiltonians and contracts it into the modal basis
'''
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
'''
def GetVEff(mVSCF, ModeOcc = None, FirstV = False):
    if ModeOcc is None:
        ModeOcc = mVSCF.ModeOcc
    return GetVEffFASTCPP(mVSCF.AnharmTensor, mVSCF.RestrictedBases, mVSCF.QUniques, mVSCF.Cs, mVSCF.MaxQuanta, ModeOcc, FirstV)
'''
def GetVEff(mVSCF, ModeOcc1 = None, ModeOcc2 = None, FirstV = False):
    if ModeOcc1 is None:
        ModeOcc1 = mVSCF.ModeOcc
    if ModeOcc2 is None:
        ModeOcc2 = ModeOcc1
    return GetVEffCPP(mVSCF.AnharmTensor, mVSCF.FCs, mVSCF.QUniques, mVSCF.QPowers, mVSCF.Cs, mVSCF.MaxQuanta, ModeOcc1, ModeOcc2, FirstV)


def GetFock(mVSCF, hs = None, Cs = None, CalcE = False):
    if Cs is None:
        Cs = mVSCF.Cs
    if hs is None:
        hs = mVSCF.HamHO

    Vs = mVSCF.GetVEff()

    # Save the double counting term while we have the effective potentials
    if CalcE:
        mVSCF.ESCF = mVSCF.CalcESCF(ModeOcc = mVSCF.ModeOcc, V0 = Vs[0])
    Fs = []
    for Mode in range(mVSCF.NModes):
        Fs.append(hs[Mode] + Vs[Mode])
    return Fs

def GetFockNMode(mVSCF, ModeOcc1 = None, ModeOcc2 = None, Cs = None, CalcE = False):
    if ModeOcc1 is None:
        ModeOcc1 = mVSCF.ModeOcc
    if ModeOcc2 is None:
        ModeOcc2 = ModeOcc1

    if Cs is None:
        Cs = mVSCF.Cs

    ints = mVSCF.ContractIntegrals(Cs)[0]
    Vmf = np.zeros((mVSCF.MaxQuanta[0], mVSCF.MaxQuanta[0]))
    for n in range(mVSCF.mol.Order):
        if n == 1:
            for j in range(1, mVSCF.NModes):
                Vmf += ints[n][0, j][:, ModeOcc1[j], :, ModeOcc2[j]]
        if n == 2:
            for j in range(1, mVSCF.NModes):
                for k in range(j + 1, mVSCF.NModes):
                    Vmf += ints[n][0, j, k][:, ModeOcc1[j], ModeOcc1[k], :, ModeOcc2[j], ModeOcc2[k]]

    Fs = []
    for i in range(mVSCF.NModes):
        F = np.zeros((mVSCF.MaxQuanta[i], mVSCF.MaxQuanta[i]))
        np.fill_diagonal(F, mVSCF.mol.onemode_eig[i])
        #F = ints[0][i]
        #Vmf += (F - np.diag(np.diag(F)))
        '''
        if CalcE:
            V0 = mVSCF.mol.V0
            if mVSCF.mol.Order > 0:
                for j in range(mVSCF.NModes):
                    if j != i:
                        if ModeOcc1[j] == ModeOcc2[j]:
                            V0 += ints[0][j][ModeOcc1[j], ModeOcc2[j]]
            if mVSCF.mol.Order > 1:
                for j in range(mVSCF.NModes):
                    for k in range(mVSCF.NModes):
                        if j != i and k != i:
                            V0 += ints[1][j, k][ModeOcc1[j], ModeOcc1[k], ModeOcc2[j], ModeOcc2[k]]
            if mVSCF.mol.Order > 2:
                for j in range(mVSCF.NModes):
                    for k in range(mVSCF.NModes):
                        for l in range(mVSCF.NModes):
                            if j != i and k != i and l != i:
                                V0 += ints[2][j, k, l][ModeOcc1[j], ModeOcc1[k], ModeOcc1[l], ModeOcc2[j], ModeOcc2[k], ModeOcc2[l]]
            #F += np.eye(mVSCF.MaxQuanta[i]) * V0
        '''
        if mVSCF.mol.Order > 1:
            for j in range(mVSCF.NModes):
                if j != i:
                    F += ints[1][i, j][:, ModeOcc1[j], :, ModeOcc2[j]]
        if mVSCF.mol.Order > 2:
            for j in range(mVSCF.NModes):
                for k in range(j + 1, mVSCF.NModes):
                    if j != i and k != i:
                        F += ints[2][i, j, k][:, ModeOcc1[j], ModeOcc1[k], :, ModeOcc2[j], ModeOcc2[k]]
        #if i == 0:
        #    Vmf += np.eye(mVSCF.MaxQuanta[0]) * mVSCF.mol.V0
        Fs.append(F)
        if CalcE:
            mVSCF.ESCF = mVSCF.CalcESCF(ModeOcc = ModeOcc1, V0 = Vmf)

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
        if mVSCF.verbose > 0:
            print("VSCF Iteration %d complete with an SCF error of %.12f/%.12f and SCF Energy of %.6f" % (It, ConvErr, EnergyErr, mVSCF.ESCF), flush = True)
        It += 1
        if It > mVSCF.MaxIterations:
            raise RuntimeError("Maximum number of SCF iterations reached without convergence.")
    mVSCF.Converged = True

def SCFNMode(mVSCF, DoDIIS = True, tol = 1e-8, etol = 1e-6):
    if DoDIIS:
        mVSCF.AllFs = []
        mVSCF.AllErrs = []

    ConvErr = 1
    It = 1
    EnergyErr = 1
    while(ConvErr > tol or EnergyErr > etol):
        EnergyErr = mVSCF.ESCF
        ConvErr = mVSCF.SCFIteration(It, DoDIIS = DoDIIS)
        EnergyErr = abs(EnergyErr - mVSCF.ESCF)
        if mVSCF.verbose > 0:
            print("VSCF Iteration %d complete with an SCF error of %.12f/%.12f and SCF Energy of %.6f" % (It, ConvErr, EnergyErr, mVSCF.ESCF), flush = True)
        It += 1
        if It > mVSCF.MaxIterations:
            raise RuntimeError("Maximum number of SCF iterations reached without convergence.")
    mVSCF.Converged = True


def LCLine(mVSCF, ModeOcc, thr = 2.5e-1):
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
        EOld = EMax
        for B in LStates:
            for m in range(mVSCF.NModes):
                if MaxQuanta is not None and B[m] >= MaxQuanta[m] - 1:
                    continue
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

def MakeCoeffMatrix(mVSCF, NStates = None):
    if NStates is None:
        NStates = np.prod(mVSCF.MaxQuanta)

    mVSCF.MakeBasisList()

    VSCFBasis0 = mVSCF.LowestStates(NStates - 1, MaxQuanta = mVSCF.MaxQuanta)
    
    E = []
    for B in VSCFBasis0:
        E.append(mVSCF.CalcESCF(ModeOcc = B))
        if len(E) == NStates:
            break
    E = np.asarray(E)
    '''
    SortedInd = np.argsort(E)
    mVSCF.E = E[SortedInd]
    VSCFBasis = []
    for i in SortedInd:
        VSCFBasis.append(VSCFBasis0[i])
    '''
    VSCFBasis = VSCFBasis0
    mVSCF.E = E

    mVSCF.C = np.ones((len(mVSCF.BasisList), NStates))
    for n in range(NStates):
        for i, B in enumerate(mVSCF.BasisList):
            for m in range(mVSCF.NModes):
                mVSCF.C[i, n] *= mVSCF.Cs[m][B[m], VSCFBasis[n][m]]
        '''
        Cn = mVSCF.Cs[0][:, mVSCF.BasisList[n][0]]
        for i in range(1, mVSCF.NModes):
            Cn = np.kron(mVSCF.Cs[i][:, mVSCF.BasisList[n][i]], Cn)
        mVSCF.C[:, n] = Cn
        '''
def MakeBasisList(mVSCF):
    mVSCF.Basis, mVSCF.BasisList = InitGridBasis(mVSCF.Frequencies, mVSCF.MaxQuanta)

def PrintResults(mVSCF, NStates = None, PrintLC = False):
    FinalE = []
    if NStates is None:
        NStates = np.prod(mVSCF.MaxQuanta)
    #assert(NStates <= np.prod(mVSCF.MaxQuanta))
    EnergyBList = mVSCF.LowestStates(NStates, MaxQuanta = mVSCF.MaxQuanta)
    for B in  EnergyBList:
        FinalE.append(mVSCF.CalcESCF(ModeOcc = B))
    SortedInd = np.argsort(np.asarray(FinalE))

    if PrintLC:
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
        if PrintLC:
            LCString = mVSCF.LCLine(EnergyBList[i])
            OutLine += '\t%s' % (LCString)
        print(OutLine, flush = True)

def AnalyzeModals(mVSCF):
    print("|*****************************************|")
    print("|********** VSCF Modal Analysis **********|", flush = True)
    print("|*****************************************|")
    print("")

    for m in range(mVSCF.NModes):
        print(" == Modal", m + 1, "== ")
        Lines = [""] * (mVSCF.MaxQuanta[m] + 1)
        for i in range(mVSCF.MaxQuanta[m]):
            Lines[0] += '\tMO_%d' % (i + 1)
        for n in range(mVSCF.MaxQuanta[m]):
            Lines[n + 1] += 'HO_%d' % (n + 1)
            for nn in range(mVSCF.MaxQuanta[m]):
                Lines[n + 1] += '\t{:.6f}'.format(mVSCF.Cs[m][n, nn])
        for Line in Lines:
            print(Line, flush = True)
        print("")
    
def PrintPotential(mVSCF, Normalized = True):
    print("|*****************************************|")
    print("|********** Printing  Potential **********|", flush = True)
    print("|*****************************************|")
    print("")

    print("Normal Mode Frequencies:")
    for w in mVSCF.Frequencies:
        print(w)
    print("")
    print("Potential:")
    for FCp in mVSCF.Potential:
        for FC in FCp:
            V = FC.fc
            if Normalized:
                V *= np.sqrt(pow(2.0, FC.Order))
                for q in FC.QPowers:
                    V *= np.math.factorial(q)
            FCLine = str(FC.Order) + " "
            for q in FC.QIndices:
                FCLine += (str(q) + " ")
            FCLine += str(V)
            print(FCLine, flush = True)

class VSCF:
    InitCs = InitCs
    GetModalSlices = GetModalSlices
    MakeAnharmTensor = MakeAnharmTensor
    MakeHOHam = MakeHOHam
    UpdateFC = UpdateFC

    GetVEff = GetVEff
    GetFock = GetFock
    StoreFock = StoreFock
    DIISUpdate = DIISUpdate
    SCF = SCF
    SCFIteration = SCFIteration
    LowestStates = LowestStates
    CalcESCF = CalcESCF

    PrintResults = PrintResults
    LCLine = LCLine
    AnalyzeModals = AnalyzeModals
    PrintPotential = PrintPotential

    MakeBasisList = MakeBasisList
    MakeCoeffMatrix = MakeCoeffMatrix
    def __init__(self, Frequencies, UnscaledPotential, MaxQuanta = 2, NStates = 10, **kwargs):
        self.Frequencies = Frequencies
        self.NModes = self.Frequencies.shape[0]

        self.Timer = TIMER(5)
        self.Timer.set_overhead(self.Timer.estimate_overhead())
        self.TimerNames = ["Basis Init", "Anharm Pot Init", "Fock Generation", "DIIS", "Solve Fock"]

        self.Potential = [[], [], [], []]
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

        '''
        self.SlowV = SlowV
        self.FastType = FastType
        if self.SlowV:
            self.Timer.start(0)
            self.Basis, self.BasisList = InitGridBasis(self.Frequencies, MaxQuanta)
            self.ModalSlices = self.GetModalSlices()
            self.Timer.stop(0)
        '''
        self.HamHO = self.MakeHOHam()
        '''
        if self.SlowV:
            self.Timer.start(1)
            self.AnharmTensor = self.MakeAnharmTensorSLOW()
            self.Timer.stop(1)
        else:
            if self.FastType == 1:
                self.AnharmTensor, self.RestrictedBases, self.QUniques = self.MakeAnharmTensor()
            else:
        '''
        self.AnharmTensor, self.FCs, self.QUniques, self.QPowers = self.MakeAnharmTensor()
        # These needs to be saved for VCI
        #self.Potential = []
        self.PotentialList = []
        
        self.ModeOcc = [0] * self.NModes
        self.ESCF = 0.0
        self.Cs = self.InitCs()

        self.DoDIIS = False 
        self.DIISSpace = 5
        self.DIISStart = 10
        self.MaxIterations = 1000

        self.verbose = 2
        self.Converged = False

        self.__dict__.update(kwargs)

    def kernel(self, C0s = None):
        if C0s is not None:
            self.Cs = C0s
        else:
            self.Cs = self.InitCs()
        
        self.SCF(DoDIIS = self.DoDIIS)
        if self.verbose > 1:
            self.Timer.report(self.TimerNames)
        return self.ESCF

class NModeVSCF(VSCF):
    GetFock = GetFockNMode
    SCF = SCFNMode

    def __init__(self, mol, **kwargs):
        self.mol = mol
        self.Frequencies = mol.Frequencies
        self.NModes = self.Frequencies.shape[0]
        self.Es = mol.onemode_eig

        self.Timer = TIMER(5)
        self.Timer.set_overhead(self.Timer.estimate_overhead())
        self.TimerNames = ["Basis Init", "Anharm Pot Init", "Fock Generation", "DIIS", "Solve Fock"]
        self.MaxQuanta = [mol.ngridpts] * self.NModes
        self.ModeOcc = [0] * self.NModes
        self.ESCF = 0.0
        self.Cs = self.InitCs()

        self.DoDIIS = False
        self.DIISSpace = 5
        self.DIISStart = 10
        self.MaxIterations = 1000

        self.verbose = 2
        self.Converged = False

        self.__dict__.update(kwargs)

    def ContractIntegrals(self, Cs = None, withDipole = False):
        if Cs is None:
            Cs = self.Cs
        ints = [np.asarray([]), np.asarray([[[[[[]]]]]]), np.asarray([[[[[[[[[]]]]]]]]])]
        for n in range(self.mol.Order):
            if n == 0:
                ints[0] = np.empty(self.NModes, dtype = object)
                for i in range(self.NModes):
                    ints[0][i] = self.Cs[i].T @ np.diag(self.mol.onemode_eig[i]) @ self.Cs[i]
            if n == 1:
                ints[1] = np.empty((self.NModes, self.NModes), dtype = object)
                for i in range(self.NModes):
                    for j in range(self.NModes):
                        ints[1][i, j] = np.einsum('ip,jq,ijkl,kr,ls->pqrs', self.Cs[i], self.Cs[j], self.mol.ints[1][i, j], self.Cs[i], self.Cs[j], optimize = True)
            if n == 2:
                ints[2] = np.empty((self.NModes, self.NModes, self.NModes), dtype = object)
                for i in range(self.NModes):
                    for j in range(self.NModes):
                        for k in range(self.NModes):
                            ints[2][i, j, k] = np.einsum('ip,jq,kr,ijklmn,ls,mt,nu->pqrstu', self.Cs[i], self.Cs[j], self.Cs[k], self.mol.ints[2][i, j, k], self.Cs[i], self.Cs[j], self.Cs[k], optimize = True)
        
        dip_ints = None
        if withDipole:
            dip_ints = [np.asarray([[]] * 3), np.asarray([[[[[[[]]]]]]] * 3), np.asarray([[[[[[[[[[]]]]]]]]]] * 3)]
            for n in range(self.mol.Order):
                if n == 0:
                    dip_ints[0] = np.empty((3, self.NModes), dtype = object)
                    for i in range(self.NModes):
                        dip_ints[0][i] = np.einsum('ip,jq,xij->xpq', self.Cs[i], self.Cs[i], self.mol.dip_ints[i], optimize = True)
                if n == 1:
                    dip_ints[1] = np.empty((3, self.NModes, self.NModes), dtype = object)
                    for i in range(self.NModes):
                        for j in range(self.NModes):
                            dip_ints[1][i, j] = np.einsum('ip,jq,xijkl,kr,ls->xpqrs', self.Cs[i], self.Cs[j], self.mol.dip_ints[1][i, j], self.Cs[i], self.Cs[j], optimize = True)
                if n == 2:
                    dip_ints[2] = np.empty((3, self.NModes, self.NModes, self.NModes), dtype = object)
                    for i in range(self.NModes):
                        for j in range(self.NModes):
                            for k in range(self.NModes):
                                dip_ints[2][i, j, k] = np.einsum('ip,jq,kr,xijklm,ls,mt,nu->xpqrstu', self.Cs[i], self.Cs[j], self.Cs[k], self.mol.dip_ints[2][i, j, k], self.Cs[i], self.Cs[j], self.Cs[k], optimize = True)
        return ints, dip_ints

if __name__ == "__main__":
    from vstr.utils.read_jf_input import Read
    w, MaxQuanta, MaxTotalQuanta, Vs, eps1, eps2, eps3, NWalkers, NSamples, NStates = Read('CLO2.inp')
    print(MaxQuanta)
    mf = VSCF(w, Vs, MaxQuanta = MaxQuanta, NStates = NStates, SlowV = False)
    mf.SCF(DoDIIS = False)
    mf.AnalyzeModals()
    mf.PrintResults(NStates = 10)
    mf.MakeCoeffMatrix()
    '''
    mo = [0]*48
    #mo[-1]=2
    #print(mf.CalcESCF(ModeOcc = mo))
    mo[-1] = 3
    print(mf.CalcESCF(ModeOcc = mo))
    mf.ModeOcc = mo
    mf.SCF(DoDIIS=False)
    mf.PrintResults(NStates = 10)#mf.AnalyzeModals()
    '''
