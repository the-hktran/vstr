import numpy as np
from vstr.cpp_wrappers.vhci_jf.vhci_jf_functions import WaveFunction, FConst, HOFunc, GenerateHamV, GenerateSparseHamV, GenerateHamAnharmV, VCISparseHamFromVSCF, HeatBath_Sort_FC
from vstr.spectra.dipole import GetDipoleSurface, MakeDipoleList
from scipy import sparse
import matplotlib.pyplot as plt

def GetTransitionDipoleMatrix(mIR, IncludeZeroth = True):
    mIR.D = GenerateHamAnharmV(mIR.mVCI.Basis, list(mIR.mVCI.Frequencies), mIR.DipoleSurfaceList, mIR.DipoleSurface[3], mIR.DipoleSurface[4], mIR.DipoleSurface[5], mIR.DipoleSurface[6])
    if IncludeZeroth:
        mIR.D += np.eye(mIR.D.shape[0]) * mIR.DipoleSurface[0][0]
    else:
        np.fill_diagonal(mIR.D, 0)

def GetTransitionDipoleMatrixFromVSCF(mIR, IncludeZeroth = True):
    X0 = []
    for X in mIR.Xs:
        X0.append(np.zeros(X.shape))
    mIR.D = VCISparseHamFromVSCF(mIR.Basis, mIR.Basis, mIR.Frequencies, mIR.DipoleSurfaceList, mIR.Ys, X0, True).todense()
    if IncludeZeroth:
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

def ApproximateAInv(mIR, w, Order = 1):
    HOD = mIR.mVCI.H.todense()
    D = np.asarray(HOD.diagonal().copy() + (w + mIR.mVCI.E[0] + 1.j * mIR.eta)).ravel()
    DInv = D**(-1)
    np.fill_diagonal(HOD, 0)
    
    AInv = np.diag(DInv) + np.einsum('ij,i,j->ij', HOD, DInv, DInv, optimize = True)
    if Order >= 2:
        HDH = np.einsum('ij,j,jk->ik', HOD, DInv, HOD, optimize = True)
        AInv += 0.5 * np.einsum('ij,i,j->ij', HDH, DInv, DInv, optimize = True)
    return AInv

def SpectralScreenBasis(mIR, Type = 'ignore_H', eps = 0.01, InitState = None):
    if InitState is None:
        InitState = mIR.InitState
    return mIR.mVCI.ScreenBasis(Ws = mIR.DipoleSurfaceList, C = abs(mIR.mVCI.C[:, InitState]), eps = eps)

def SpectralHCIStep(mIR, w, eps = 0.01, InitState = 0):
    NewBasis, NAdded = mIR.SpectralScreenBasis(eps = eps, InitState = InitState)
    # This new basis set is pruned for all basis sets giving small denominators
    
def SpectralHCI(mIR, w):
    # This ignores the numerator
    # First, get new basis based on the dipole operator


def Intensity(mIR, w):
    # Should define new basis with HCI and then solve VHCI here, be sure to update mVCI object
    A, b = mIR.GetAb(w)
    x = SolveAxb(A, b)
    #AInv = mIR.ApproximateAInv(w)
    #AInvX = np.linalg.inv(A)
    #print(w, ((AInv - AInvX).conj() * (AInv - AInvX)).sum().real)
    #x = AInv @ b
    x = np.asarray(x).ravel()
    b = np.asarray(b).ravel()
    return (1.j * np.dot(b, x)).real / np.pi

def PlotSpectrum(mIR, PlotName, XLabel = "Frequency", YLabel = "Intensity", Title = "IR Spectrum"):
    plt.plot(mIR.ws, mIR.Is, linestyle = '-', marker = None)
    plt.xlabel(XLabel)
    plt.ylabel(YLabel)
    plt.title(Title)
    plt.savefig(PlotName)

def TestPowerSeries(mIR):
    H = mIR.mVCI.H.todense()
    n = H.shape[0]
    for w in (mIR.mVCI.E - mIR.mVCI.E[0]): #range(1000,16000,1):
        D = np.diag(H.diagonal().copy()) + np.eye(H.shape[0]) * (w + mIR.mVCI.E[0] + mIR.eta * 1.j)
        np.fill_diagonal(H, 0)
        DInv = np.linalg.inv(D)
        AInvX = np.linalg.inv(D - H)
        AInvP = DInv + DInv @ H @ DInv
        #print(AInvP)
        #print(AInvX)
        np.save("AInvP", AInvP)
        np.save("AInvX", AInvX)
        dA = AInvX - AInvP
        dA2 = (dA.conj() * dA).real
        plt.clf()
        plt.imshow(dA2, interpolation='nearest')
        plt.colorbar()
        plt.savefig("dA")
        err = (dA.conj() * dA).sum()
        if w > 1000:
            break
        print("err", w, err.real)

class LinearResponseIR:
    GetAb = GetAb
    GetTransitionDipoleMatrix = GetTransitionDipoleMatrix
    Intensity = Intensity
    ApproximateAInv = ApproximateAInv
    PlotSpectrum = PlotSpectrum

    TestPowerSeries = TestPowerSeries

    def __init__(self, mf, mVCI, FreqRange = [0, 5000], NPoints = 100, eta = 10, NormalModes = None, DipoleSurface = None, **kwargs):
        self.mf = mf
        self.mVCI = mVCI
        self.C = mVHCI.C
        self.E = mVHCI.E
        self.Frequencies = mVHCI.Frequencies
        self.Basis = mVHCI.Basis
        self.NormalModes = NormalModes
        self.InitState = 0
        self.FreqRange = FreqRange
        self.NPoints = NPoints
        self.eta = eta
        self.Order = 1
        self.DipoleSurface = DipoleSurface

        self.__dict__.update(kwargs)

    def kernel(self):
        # Make dipole surface
        if self.DipoleSurface is None:
            mu_raw = GetDipoleSurface(self.mf, self.NormalModes, Freq = self.Frequencies, Order = self.Order, dx = 1e-1)
            print(mu_raw)
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
        self.DipoleSurfaceList = HeatBath_Sort_FC(self.DipoleSurfaceList)

        self.GetTransitionDipoleMatrix(IncludeZeroth = False)
        np.fill_diagonal(self.D, 0)

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
    
    mVHCI = VHCI(w, V, MaxQuanta = 10, MaxTotalQuanta = 1, eps1 = 1, eps2 = 0.001, eps3 = -1, NWalkers = 50, NSamples = 50, NStates = 3)
    mVHCI.kernel()
    #mVHCI.E, mVHCI.C = np.linalg.eigh(mVHCI.H.todense())
    #mVHCI.E_HCI = mVHCI.E

    mIR = LinearResponseIR(mf, mVHCI, FreqRange = [0, 16000], NPoints = 1000, eta = 100, NormalModes = NormalModes, Order = 4)
    mIR.kernel()
    mIR.PlotSpectrum("water_spectrum_lr.png")
    mIR.TestPowerSeries()
