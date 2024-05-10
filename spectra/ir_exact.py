import numpy as np
from vstr.cpp_wrappers.vhci_jf.vhci_jf_functions import WaveFunction, FConst, HOFunc, GenerateHamV, GenerateSparseHamV, GenerateSparseHamAnharmV, VCISparseHamFromVSCF
from vstr.spectra.dipole import GetDipoleSurface, MakeDipoleList
from vstr.utils.perf_utils import TIMER
from scipy import sparse
import matplotlib.pyplot as plt

def GetTransitionDipoleMatrix(mIR, IncludeZeroth = False):
    mIR.D = []
    for xi in range(3):
        mIR.D.append(GenerateSparseHamAnharmV(mIR.Basis, list(mIR.Frequencies), mIR.DipoleSurfaceList[xi], mIR.DipoleSurface[xi][3], mIR.DipoleSurface[xi][4], mIR.DipoleSurface[xi][5], mIR.DipoleSurface[xi][6]))
        if IncludeZeroth:
            for i in range(mIR.D[xi].shape[0]):
                mIR.D[xi][i, i] += mIR.DipoleSurface[xi][0][0]

def GetTransitionDipoleMatrixFromVSCF(mIR, IncludeZeroth = False):
    X0 = []
    for X in mIR.Xs:
        X0.append(np.zeros(X.shape))
    mIR.D = []
    for xi in range(3):
        mIR.D.append(VCISparseHamFromVSCF(mIR.Basis, mIR.Basis, mIR.Frequencies, mIR.DipoleSurfaceList[xi], mIR.Ys, X0, True))
        if IncludeZeroth:
            for i in range(mIR.D[xi].shape[0]):
                mIR.D[xi][i, i] += mIR.DipoleSurface[xi][0][0]

def GetSpectralIntensities(mIR):
    mIR.Intensities = [None] * 3
    for xi in range(3):
        if mIR.C is None:
            CDC = mIR.D[xi][0, 1:].todense()
        else:
            CDC = mIR.C[:, 0].T @ mIR.D[xi] @ mIR.C[:, 1:]
        if CDC.ndim == 2:
            CDC = np.asarray(CDC).ravel()
        mIR.Intensities[xi] = CDC**2
    mIR.Excitations = []
    for w in mIR.E[1:]:
        mIR.Excitations.append(w - mIR.E[0])

def Lorentzian(x, x0, L):
    #return 0.5 * L / (np.pi * ((x - x0)**2 + 0.25 * L**2))
    return L / (np.pi * ((x - x0)**2 + L**2))

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
            y += (mIR.Intensities[0][n] + mIR.Intensities[1][n] + mIR.Intensities[2][n]) * Lorentzian(x, mIR.Excitations[n], L = L)
        Y.append(y)
    mIR.ws = X
    mIR.Is = np.asarray(Y)
    if mIR.Normalize:
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
        self.Normalize = False
        
        self.__dict__.update(kwargs)

        self.Timer = TIMER(2)
        self.TimerNames = ["Generate Dipole Matrix", "Plot Spectrum"]

    def kernel(self):
        if self.DipoleSurface is None:
            mu_raw = GetDipoleSurface(self.mf, self.NormalModes, Freq = self.Frequencies, Order = self.Order, dx = 1e-1)
            self.DipoleSurface = []
            for x in range(3):
                DS = []
                DS.append(mu_raw[x][0][0])
                for n in range(1, self.Order + 1):
                    Dn = MakeDipoleList(mu_raw[x][n])
                    DS.append(Dn.copy())
                for n in range(self.Order + 1, 7):
                    DS.append([])
                self.DipoleSurface.append(DS.copy())
        
        self.DipoleSurfaceList = []
        for x in range(3):
            DSListX = []
            for Dn in self.DipoleSurface[x][1:]:
                DSListX += Dn
            self.DipoleSurfaceList.append(DSListX.copy())

        # Cubic/linear and quadratic/quartic terms must stack, and there must be something at the highest even order, so we will make sure that happens here
        for x in range(3):
            self.DipoleSurface[x][3] += self.DipoleSurface[x][1].copy()
            self.DipoleSurface[x][4] += self.DipoleSurface[x][2].copy()
            self.DipoleSurface[x][1] = []
            self.DipoleSurface[x][2] = []
            self.DipoleSurfaceList[x].append(FConst(0.0, [0] * 6, False))

        self.Timer.start(0)
        self.GetTransitionDipoleMatrix()
        self.Timer.stop(0)
        self.Timer.start(1)
        self.GetSpectralIntensities()
        self.Timer.stop(1)

class VSCFIRSpectra(IRSpectra):
    GetTransitionDipoleMatrix = GetTransitionDipoleMatrixFromVSCF

    def __init__(self, mf, mVHCI, NormalModes = None, **kwargs):
        IRSpectra.__init__(self, mf, mVHCI, NormalModes = NormalModes)
        self.__dict__.update(kwargs)

        self.Ys = mVHCI.Ys
        self.Xs = mVHCI.Xs

class IRSpectraNMode(IRSpectra):
    from vstr.spectra.ir_lr import GetTransitionDipoleMatrixNMode
    GetTransitionDipoleMatrix = GetTransitionDipoleMatrixNMode

    def __init__(self, mf, mVHCI, NormalModes = None, SpectralHBMethod = 2, **kwargs):
        IRSpectra.__init__(self, mf, mVHCI, NormalModes = NormalModes)
        self.__dict__.update(kwargs)
        self.mVCI = mVHCI
        self.mol = mVHCI.mol

        self.DipoleSurface = [[], [], [], []]
    
    def kernel(self):
        self.Timer.start(0)
        self.GetTransitionDipoleMatrix()
        self.Timer.stop(0)
        self.Timer.start(1)
        self.GetSpectralIntensities()
        self.Timer.stop(1)

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
    mol.basis = 'cc-pvdz'
    mol.build()
    mf = scf.RHF(mol)
    mf.kernel()

    w, NormalModes = GetNormalModes(mf)

    V = GetFF(mf, NormalModes, w, Order = 4)
    
    mVHCI = VHCI(w, V, MaxQuanta = 10, MaxTotalQuanta = 3, eps1 = 0.1, eps2 = 0.01, eps3 = -1, NWalkers = 50, NSamples = 50, NStates = 10)
    mVHCI.kernel()
    #from scipy import sparse
    #mVHCI.E, mVHCI.C = np.linalg.eigh(mVHCI.H.todense())
    #mVHCI.E_HCI = mVHCI.E
    #mVHCI.PrintResults(thr=0)

    mIR = IRSpectra(mf, mVHCI, NormalModes = NormalModes, Order = 2)
    mIR.kernel()
    print(mIR.Intensities)
    print(mIR.Excitations)
    mIR.PlotSpectrum("water_spectrum.png", XMin = 0, XMax = 5000)

    from vstr.mf.vscf import VSCF
    from vstr.ci.vci import VCI
    vmf = VSCF(w, V, MaxQuanta = 10, NStates = 100)
    vmf.kernel()
    #vmf.InitCs()
    #vmf.MakeCoeffMatrix(NStates = 6)
    #vmf.PrintResults(NStates = 100)

    '''
    mIR = IRSpectra(mf, vmf, NormalModes = NormalModes, Order = 2)
    mIR.kernel()
    mIR.PlotSpectrum("water_spectrum.png", XMin = 0, XMax = 5000)
    print(mIR.C.T @ mIR.D @ mIR.C)
    '''

    mVCI = VCI(vmf, 10, NStates = 100)
    mVCI.InitBasisAndC()#(Basis = vmf.Basis)
    mIR = VSCFIRSpectra(mf, mVCI, NormalModes = NormalModes, Order = 2)
    mIR.kernel()
    mIR.PlotSpectrum("water_spectrum.png", XMin = 0, XMax = 5000)

    '''
    mVCI = VCI(vmf, 3, eps1 = 10, eps2 = 0.01, eps3 = -1, NWalkers = 50, NSamples = 50, NStates = 10)
    mVCI.kernel()
    mVSCFIR = VSCFIRSpectra(mf, mVCI, NormalModes = NormalModes)
    mVSCFIR.kernel()
    print(mVSCFIR.Intensities)
    print(mVSCFIR.Excitations)
    '''
