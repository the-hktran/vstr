import numpy as np
from vstr.utils.init_functions import FormW, InitTruncatedBasis

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

def InitModalBasis(mVSCF)
    B0 = [0] * mVSCF.NModes
    Basis = []
    Basis.append(B0)
    for N in range(mVSCF.NModes):
        for n in range(mVSCF.MaxQuanta[N]):
            for B in Basis:
                BC = B.copy()
                B[N] = n
                Basis.append(B)
    return Basis

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
            for B in mVSCF.Basis:
                if B[N] == i
                    BasisByModeByIndex.append(B)
            BasisByMode.append(BasisByModeByIndex)
        ModalSlices.append(BasisByMode)
    return ModalSlices

'''
This function collects all the nonzero anharmonic terms, for each anharmonic force constant
'''
def MakeAnharmTensor(mVSCF)
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
        H = GenerateSparseHamV(mVSCF.Basis, mVSCF.Frequencies, mVSCF.PotentialList, CubicFC, QuarticFC, QuinticFC, SexticFC)
        AnharmTensor.append(H)
    return AnharmTensor

class VSCF:
    InitModalBasis = InitModalBasis
    GetModalSlices = GetModalSlices
    MakeAnharmTensor = MakeAnharmTensor
    
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
            MaxQuanta = [MaxQuanta] * self.NModes

        self.Basis = self.InitModalBasis()
        self.ModalSlices = self.GetModalSlices()
        self.AnharmTensor = self.MakeAnharmTensor()
