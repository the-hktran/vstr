import numpy as np
from vstr.cpp_wrappers.vhci_jf.vhci_jf_functions import WaveFunction, FConst, HOFunc, GenerateHamV, GenerateSparseHamV, GenerateSparseHamAnharmV, VCISparseHamFromVSCF
from vstr.spectra.dipole import GetDipoleSurface, MakeDipoleList
from scipy import sparse
import matplotlib.pyplot as plt

def GetTransitionDipoleMatrix(mIR, IncludeZeroth = False):
    mIR.D = GenerateSparseHamAnharmV(mIR.Basis, list(mIR.Frequencies), mIR.DipoleSurfaceList, mIR.DipoleSurface[3], mIR.DipoleSurface[4], mIR.DipoleSurface[5], mIR.DipoleSurface[6])
    if IncludeZeroth:
        for i in range(mIR.D.shape[0]):
            mIR.D[i, i] += mIR.DipoleSurface[0][0]

def GetTransitionDipoleMatrixFromVSCF(mIR, IncludeZeroth = False):
    X0 = []
    for X in mIR.Xs:
        X0.append(np.zeros(X.shape))
    mIR.D = VCISparseHamFromVSCF(mIR.Basis, mIR.Basis, mIR.Frequencies, mIR.DipoleSurfaceList, mIR.Ys, X0, True).todense()
    if IncludeZeroth:
        for i in range(mIR.D.shape[0]):
            mIR.D[i, i] += mIR.DipoleSurface[0][0]

def GetSpectralIntensities(mIR):
    CDC = mIR.C[:, 0].T @ mIR.D @ mIR.C[:, 1:]
    if CDC.ndim == 2:
        CDC = np.asarray(CDC).ravel()
    mIR.Intensities = CDC**2
    mIR.Excitations = []
    for w in mIR.E[1:]:
        mIR.Excitations.append(w - mIR.E[0])

def Lorentzian(x, x0, L):
    return 0.5 * L / (np.pi * ((x - x0)**2 + 0.25 * L**2))

def PlotSpectrum(mIR, PlotName, NPoints = 1000, L = 100, XLabel = "Frequency", YLabel = "Intensity", Title = "IR Spectrum", XMin = None, XMax = None):
    if XMin is None:
        XMin = 0
    if XMax is None:
        XMax = mIR.Excitations[-1] + 100
    X = np.linspace(XMin, XMax, num = NPoints)
    Y = []
    for x in X:
        y = 0
        for n in range(len(mIR.Excitations)):
            y += mIR.Intensities[n] * Lorentzian(x, mIR.Excitations[n], L = L)
        Y.append(y)
    mIR.ws = X
    mIR.Is = np.asarray(Y)
    mIR.Is = mIR.Is / max(mIR.Is) 
    plt.plot(mIR.ws, mIR.Is, linestyle = '-', marker = None)
    plt.xlabel(XLabel)
    plt.ylabel(YLabel)
    plt.title(Title)
    plt.savefig(PlotName)

def SaveSpectrum(mIR, SaveName):
    A = np.zeros((mIR.ws.shape[0], 2))
    A[:, 0] = mIR.ws
    A[:, 1] = mIR.Is
    np.savetxt(SaveName + ".csv", A, delimiter=",")

class IRSpectra:
    GetTransitionDipoleMatrix = GetTransitionDipoleMatrix
    GetSpectralIntensities = GetSpectralIntensities
    PlotSpectrum = PlotSpectrum
    SaveSpectrum = SaveSpectrum

    def __init__(self, mf, mVHCI, NormalModes = None, DipoleSurface = None, **kwargs):
        self.mf = mf
        self.C = mVHCI.C
        self.E = mVHCI.E
        self.Frequencies = mVHCI.Frequencies
        self.Basis = mVHCI.Basis
        self.NormalModes = NormalModes
        self.Order = 1
        self.DipoleSurface = DipoleSurface
        
        self.__dict__.update(kwargs)

    def kernel(self):
        if self.DipoleSurface is None:
            # Make dipole surface
            mu_raw = GetDipoleSurface(self.mf, self.NormalModes, Freq = self.Frequencies, Order = self.Order)
            self.DipoleSurface = []
            self.DipoleSurface.append(mu_raw[0][0])
            for n in range(1, self.Order + 1):
                Dn = MakeDipoleList(mu_raw[n])
                self.DipoleSurface.append(Dn)
        else:
            self.Order = len(self.DipoleSurface) - 1

        for n in range(self.Order + 1, 7):
            self.DipoleSurface.append([])
        self.DipoleSurfaceList = []
        for Dn in self.DipoleSurface[1:]:
            self.DipoleSurfaceList += Dn

        # Cubic/linear and quadratic/quartic terms must stack, and there must be something at the highest linear order, so we will make sure that happens here
        self.DipoleSurface[3] = self.DipoleSurface[1] + self.DipoleSurface[3]
        self.DipoleSurface[4] = self.DipoleSurface[2] + self.DipoleSurface[4]
        self.DipoleSurfaceList.append(FConst(0.0, [0] * 6, False))

        self.GetTransitionDipoleMatrix()
        self.GetSpectralIntensities()

class VSCFIRSpectra(IRSpectra):
    GetTransitionDipoleMatrix = GetTransitionDipoleMatrixFromVSCF

    def __init__(self, mf, mVHCI, NormalModes = None, **kwargs):
        IRSpectra.__init__(self, mf, mVHCI, NormalModes = NormalModes)
        self.__dict__.update(kwargs)

        self.Ys = mVHCI.Ys
        self.Xs = mVHCI.Xs

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
    
    mVHCI = VHCI(w, V, MaxQuanta = 10, MaxTotalQuanta = 1, eps1 = 1, eps2 = 0.01, eps3 = -1, NWalkers = 50, NSamples = 50, NStates = 3)
    mVHCI.kernel()
    #from scipy import sparse
    #mVHCI.E, mVHCI.C = np.linalg.eigh(mVHCI.H.todense())
    #mVHCI.E_HCI = mVHCI.E
    #mVHCI.PrintResults(thr=0)

    mIR = IRSpectra(mf, mVHCI, NormalModes = NormalModes, Order = 4)
    mIR.kernel()
    print(mIR.Intensities)
    print(mIR.Excitations)
    mIR.PlotSpectrum("water_spectrum.png")

    '''
    from vstr.mf.vscf import VSCF
    from vstr.ci.vci import VCI
    vmf = VSCF(w, V, MaxQuanta = 10, NStates = 10)
    vmf.kernel()

    mVCI = VCI(vmf, 3, eps1 = 10, eps2 = 0.01, eps3 = -1, NWalkers = 50, NSamples = 50, NStates = 10)
    mVCI.kernel()
    mVSCFIR = VSCFIRSpectra(mf, mVCI, NormalModes = NormalModes)
    mVSCFIR.kernel()
    print(mVSCFIR.Intensities)
    print(mVSCFIR.Excitations)
    '''
