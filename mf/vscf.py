import numpy as np
from vstr.utils.init_funcs import FormW, InitTruncatedBasis, InitGridBasis
from vstr.cpp_wrappers.vhci_jf.vhci_jf_functions import WaveFunction, GenerateHam0V, GenerateSparseHamAnharmV, GetVEffCPP

def GetModalBasis(mVSCF, Mode, Quanta, MBasis = None):
    if MBasis is None:
        MB = []
        for n in range(self.NModes):
            if n != Mode:
                for i in range(mVSCF.MaxQuanta[n]):
                    B = [0] * self.NModes
                    B[n] = i;
                    B[Mode] = Quanta
                    WF = WaveFunction(B, mVSCF.Frequencies)
                    MB.append(WF)
    else:
        MB = MBasis.copy()
        for N, WF in enumerate(MB):
            WF.Modes[Mode].Quanta = Quanta
            MB[N] = WF
    return MB

'''
def InitBasis(mVSCF):
    B0 = [0] * mVSCF.NModes
    Basis = []
    Basis.append(B0)
    for N in range(mVSCF.NModes):
        for n in range(mVSCF.MaxQuanta[N]):
            for B in Basis:
                BC = B.copy()
                BC[N] = n
                Basis.append(BC)
    return Basis
'''

def InitModalBasis(mVSCF):
    ModalBasis = []
    for N in range(mVSCF.NModes):
        MBasis = []
        for n in range(mVSCF.MaxQuanta[N]):
            B = [0] * mVSCF.NModes
            B[N] = n;
            MBasis.append(WaveFunction(B, mVSCF.Frequencies))
        ModalBasis.append(MBasis)
    return ModalBasis

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
def MakeAnharmTensor(mVSCF):
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
        print(H)
        AnharmTensor.append(H)
    return AnharmTensor

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
def MakeHOHam(mVSCF, ModalBasis = None):
    if ModalBasis is None:
        ModalBasis = mVSCF.ModalBasis
    HamHO = []
    for Mode in range(mVSCF.NModes):
        hM = GenerateHam0V(ModalBasis[Mode], mVSCF.Frequencies)
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

'''
Takes stored anharmonic Hamiltonians and contracts it into the modal basis
'''
def GetVEff(mVSCF):
    return GetVEffCPP(mVSCF.AnharmTensor, mVSCF.Basis, mVSCF.Cs, mVSCF.ModalSlices, mVSCF.MaxQuanta, mVSCF.ModeOcc)

def GetFock(mVSCF, hs = None, Cs = None):
    if Cs is None:
        Cs = mVSCF.Cs
    if hs is None:
        hs = mVSCF.GetHCore(Cs = Cs)

    print(mVSCF.Cs)
    Vs = mVSCF.GetVEff()
    print(Vs)
    print(diieie)
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
    mVSCF.Fs = mVSCF.GetFock()
    if DoDIIS:
        mVSCF.StoreFock()
        if It > mVSCF.DIISStart:
            mVSCF.DIISUpdate() # Replaces Fock matrices with DIIS updated fock matrices
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
    return SCFErr

def SCF(mVSCF, DoDIIS = True, tol = 1e-8):
    if DoDIIS:
        mVSCF.AllFs = []
        mVSCF.AllErrs = []

    ConvErr = 1
    It = 1
    while(ConvErr > tol):
        ConvErr = mVSCF.SCFIteration(It, DoDIIS = DoDIIS)
        print("VSCF Iteration %d complete with an SCF error of %.12f" % (It, ConvErr))
        It += 1
        if It > mVSCF.MaxIterations:
            raise RuntimeError("Maximum number of SCF iterations reached without convergence.")

class VSCF:
    InitModalBasis = InitModalBasis
    InitCs = InitCs
    GetModalSlices = GetModalSlices
    MakeAnharmTensor = MakeAnharmTensor
    MakeHOHam = MakeHOHam
   
    GetHCore = GetHCore
    GetVEff = GetVEff
    GetFock = GetFock
    StoreFock = StoreFock
    DIISUpdate = DIISUpdate
    SCF = SCF
    SCFIteration = SCFIteration

    def __init__(self, Frequencies, UnscaledPotential, MaxQuanta = 2, NStates = 10, **kwargs):
        self.Frequencies = Frequencies
        self.NModes = self.Frequencies.shape[0]

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

        self.Basis = InitGridBasis(self.Frequencies, MaxQuanta)
        self.ModalBasis = self.InitModalBasis()
        self.ModalSlices = self.GetModalSlices()
        self.HamHO = self.MakeHOHam()
        self.AnharmTensor = self.MakeAnharmTensor()
        self.Cs = self.InitCs()
        self.ModeOcc = [0] * self.NModes

        self.DoDIIS = True
        self.DIISSpace = 5
        self.DIISStart = 10
        self.MaxIterations = 1000
        self.__dict__.update(kwargs)

if __name__ == "__main__":
    from vstr.utils.read_jf_input import Read
    w, MaxQuanta, MaxTotalQuanta, Vs, eps1, eps2, eps3, NWalkers, NSamples, NStates = Read('test.inp')
    mf = VSCF(w, Vs, MaxQuanta = MaxQuanta, NStates = NStates)
    mf.SCF()
