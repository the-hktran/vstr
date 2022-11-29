import numpy as np

def LowestStates(Frequencies, NStates, MaxQuanta = None):
    NModes = Frequencies.shape[0]
    LStates = []
    LStates.append([0] * NModes)
    EMax = 1e10 #np.sum([mVSCF.MaxQuanta[i] * e for i, e in enumerate(mVSCF.Frequencies)])
    for n in range(NStates):
        NextMode = 0
        EOld = EMax
        for B in LStates:
            for m in range(NModes):
                if MaxQuanta is not None and B[m] >= MaxQuanta[m] - 1:
                    continue
                BTestIncr = B.copy()
                BTestIncr[m] += 1
                if BTestIncr in LStates:
                    continue
                E = np.sum([BTestIncr[i] * e for i, e in enumerate(Frequencies)])
                if E < EOld:
                    NextMode = m
                    BSave = B.copy()
                    EOld = E.copy()
        BSave[NextMode] += 1
        LStates.append(BSave)
    return LStates

def ModeName(Basis):
    ModeLabel = ""
    for j, n in enumerate(Basis):
        if n > 1:
            ModeLabel += str(n)
        if n > 0:
            ModeLabel += 'w{} + '.format(j + 1)
    return ModeLabel[:-3]

def HarmonicEnergy(Frequencies, Basis, ZPE = 0.0):
    return np.sum([Basis[i] * e for i, e in enumerate(Frequencies)]) + ZPE

def HarmonicAnalysis(Frequencies, NStates):
    States = LowestStates(Frequencies, NStates)
    ZPE = np.sum(0.5 * Frequencies)
    for State in States:
        E = HarmonicEnergy(Frequencies, State, ZPE = ZPE)
        Label = ModeName(State)
        print(str(E), Label, flush = True)

if __name__ == '__main__':
    from vstr.utils.read_jf_input import Read
    w, MaxQuanta, MaxTotalQuanta, Vs, eps1, eps2, eps3, NWalkers, NSamples, NState = Read('../examples/jf_input/acetonitrile.inp')
    HarmonicAnalysis(w, 25)
