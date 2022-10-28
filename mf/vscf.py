import numpy as np
from vstr.utils.init_functions import FormW, InitTruncatedBasis

class VSCF:
    def __init__(self, Frequencies, UnscaledPotential, MaxQuanta = 2, MaxTotalQuanta = 2, NStates = 10, **kwargs):
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

