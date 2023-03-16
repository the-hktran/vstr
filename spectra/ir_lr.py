import numpy as np
from vstr.cpp_wrappers.vhci_jf.vhci_jf_functions import WaveFunction, FConst, HOFunc, GenerateHamV, GenerateSparseHamV, GenerateHamAnharmV, VCISparseHamFromVSCF
from vstr.spectra.dipole import GetDipoleSurface, MakeDipoleList
from scipy import sparse
import matplotlib.pyplot as plt

def GetTransitionDipoleMatrix(mIR):
    mIR.D = GenerateHamAnharmV(mIR.mVCI.Basis, list(mIR.mVCI.Frequencies), mIR.DipoleSurfaceList, mIR.DipoleSurface[3], mIR.DipoleSurface[4], mIR.DipoleSurface[5], mIR.DipoleSurface[6])
    mIR.D += np.eye(mIR.D.shape[0]) * mIR.DipoleSurface[0][0]

def GetTransitionDipoleMatrixFromVSCF(mIR):
    X0 = []
    for X in mIR.Xs:
        X0.append(np.zeros(X.shape))
    mIR.D = VCISparseHamFromVSCF(mIR.Basis, mIR.Basis, mIR.Frequencies, mIR.DipoleSurfaceList, mIR.Ys, X0, True).todense()
    mIR.D += np.eye(mIR.D.shape[0]) * mIR.DipoleSurface[0][0]

def GetAb(mIR, w, Basis = None):
    if Basis is None:
        H = mIR.mVCI.H
        C = mIR.mVCI.C
        E = mIR.mVCI.E
    else:
        H = GenerateSparseHamV(Basis, mIR.mVCI.Frequencies, mIR.mVCI.PotentialList, mIR.mVCI.Potential[0], mIR.mVCI.Potential[1], mIR.mVCI.Potential[2], mIR.mVCI.Potential[3])
        E, C = sparse.linalg.eigsh(H, k = mIR.mVCI.NStates, which = 'SM')

    A = np.eye(H.shape[0]) * (w + mIR.mVCI.E[0]) - H + np.eye(H.shape[0]) * mIR.eta * 1.j
    b = mIR.D @ C[:, 0]
    return A, b

def SolveAxb(A, b):
    return np.linalg.solve(A, b)

def Intensity(mIR, w):
    # Should define new basis with HCI and then solve VHCI here, be sure to update mVCI object
    A, b = mIR.GetAb(w)
    x = SolveAxb(A, b)
    return (1.j * np.dot(b, x)).imag

def PlotSpectrum(mIR, PlotName, XLabel = "Frequency", YLabel = "Intensity", Title = "IR Spectrum"):
    plt.scatter(mIR.ws, mIR.Is)
    plt.xlabel(XLabel)
    plt.ylabel(YLabel)
    plt.title(Title)
    plt.savefig(PlotName)

class LinearResponseIR:
    GetAb = GetAb
    GetTransitionDipoleMatrix = GetTransitionDipoleMatrix
    Intensity = Intensity
    PlotSpectrum = PlotSpectrum

    def __init__(self, mf, mVCI, FreqRange = [0, 5000], NPoints = 100, eta = 10, NormalModes = None, DipoleSurface = None, **kwargs):
        self.mf = mf
        self.mVCI = mVCI
        self.C = mVHCI.C
        self.E = mVHCI.E
        self.Frequencies = mVHCI.Frequencies
        self.Basis = mVHCI.Basis
        self.NormalModes = NormalModes
        self.FreqRange = FreqRange
        self.NPoints = NPoints
        self.eta = eta
        self.Order = 1
        self.DipoleSurface = DipoleSurface

        self.__dict__.update(kwargs)

    def kernel(self):
        # Make dipole surface
        if self.DipoleSurface is None:
            mu_raw = GetDipoleSurface(self.mf, self.NormalModes, Order = self.Order)
            self.DipoleSurface = []
            self.DipoleSurface.append(mu_raw[0][0])
            for n in range(1, self.Order + 1):
                Dn = MakeDipoleList(mu_raw[n])
                self.DipoleSurface.append(Dn)
            for n in range(self.Order + 1, 7):
                self.DipoleSurface.append([])
        self.DipoleSurfaceList = []
        for Dn in self.DipoleSurface[1:]:
            self.DipoleSurfaceList += Dn

        # Cubic/linear and quadratic/quartic terms must stack, and there must be something at the highest linear order, so we will make sure that happens here
        self.DipoleSurface[3] += self.DipoleSurface[1]
        self.DipoleSurface[4] += self.DipoleSurface[2]
        self.DipoleSurfaceList.append(FConst(0.0, [0] * 6, False))

        self.GetTransitionDipoleMatrix()

        self.ws = np.linspace(self.FreqRange[0], self.FreqRange[1], num = self.NPoints)
        I = []
        for w in self.ws:
            I.append(self.Intensity(w))
        self.Is = np.asarray(I)

if __name__ == "__main__":
    from vstr.ff.normal_modes import GetNormalModes
    from vstr.ff.force_field import GetFF
    from vstr.vhci.vhci import VHCI
    from pyscf import gto, scf

    mol = gto.M()
    mol.atom ='''
    O
    H 1 0.95
    H 1 0.95 2 104
    '''
    mol.basis = 'sto-3g'
    mol.build()
    mf = scf.RHF(mol)
    mf.kernel()

    w, NormalModes = GetNormalModes(mf)

    V = GetFF(mf, NormalModes, w, Order = 4)
    
    mVHCI = VHCI(w, V, MaxQuanta = 10, MaxTotalQuanta = 3, eps1 = 10, eps2 = 0.01, eps3 = -1, NWalkers = 50, NSamples = 50, NStates = 1)
    mVHCI.kernel()

    mIR = LinearResponseIR(mf, mVHCI, FreqRange = [0, 5000], NPoints = 100, eta = 10, NormalModes = NormalModes)
    mIR.kernel()
    print(mIR.ws)
    print(mIR.Is)
    mIR.PlotSpectrum("water_spectrum_lr.png")
