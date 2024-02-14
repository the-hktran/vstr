import numpy as np
from vstr.cpp_wrappers.vhci_jf.vhci_jf_functions import WaveFunction, FConst, HOFunc, GenerateHamV, GenerateSparseHamV, GenerateSparseHamAnharmV, VCISparseHamFromVSCF, HeatBath_Sort_FC, SpectralFrequencyPrune, SpectralFrequencyPruneFromVSCF, DoSpectralPT2, VCISparseHamNMode
from vstr.spectra.dipole import GetDipoleSurface, MakeDipoleList
from scipy import sparse
import matplotlib.pyplot as plt
import gc

def GetTransitionDipoleMatrix(mIR, xi = None, IncludeZeroth = False):
    if xi is None:
        mIR.D = []
        for x in range(3):
            Dx = GenerateSparseHamAnharmV(mIR.mVCI.Basis, list(mIR.mVCI.Frequencies), mIR.DipoleSurfaceList[x], mIR.DipoleSurface[x][3], mIR.DipoleSurface[x][4], mIR.DipoleSurface[x][5], mIR.DipoleSurface[x][6])
            if IncludeZeroth:
                D0 = Dx.diagonal()
                D0 += mIR.DipoleSurface[x][0][0]
                Dx.setdiag(D0)
            #else:
            #    Dx.setdiag(0)
            mIR.D.append(Dx)
    else:
        Dx = GenerateSparseHamAnharmV(mIR.mVCI.Basis, list(mIR.mVCI.Frequencies), mIR.DipoleSurfaceList[xi], mIR.DipoleSurface[xi][3], mIR.DipoleSurface[xi][4], mIR.DipoleSurface[xi][5], mIR.DipoleSurface[xi][6])
        if IncludeZeroth:
            D0 = Dx.diagonal()
            D0 += mIR.DipoleSurface[x][0][0]
            Dx.setdiag(D0)
        else:
            Dx.setdiag(0)
        mIR.D[xi] = Dx

def GetTransitionDipoleMatrixNMode(mIR, xi = None, IncludeZeroth = False):
    if xi is None:
        mIR.D = []
        for x in range(3):
            Dx = VCISparseHamNMode(mIR.mVCI.Basis, mIR.mVCI.Basis, list(np.zeros_like(mIR.mVCI.Frequencies)), mIR.mol.mu0[x], mIR.mol.dip_ints[0][x].tolist(), mIR.mol.dip_ints[1][x].tolist(), mIR.mol.dip_ints[2][x].tolist(), True)
            if not IncludeZeroth:
                Dx.setdiag(0)
            mIR.D.append(Dx)
    else:
        Dx = VCISparseHamNMode(mIR.mVCI.Basis, mIR.mVCI.Basis, list(np.zeros_like(mIR.mVCI.Frequencies)), mIR.mol.mu0[xi], mIR.mol.dip_ints[0][xi].tolist(), mIR.mol.dip_ints[1][xi].tolist(), mIR.mol.dip_ints[2][xi].tolist(), True)
        if not IncludeZeroth:
            Dx.setdiag(0)
        mIR.D[xi] = Dx

def GetTransitionDipoleMatrixFromVSCF(mIR, xi = None, IncludeZeroth = False):
    X0 = []
    for X in mIR.Xs:
        X0.append(np.zeros(X.shape))

    if xi is None:
        mIR.D = []
        for x in range(3):
            Dx = VCISparseHamFromVSCF(mIR.mVCI.Basis, mIR.mVCI.Basis, mIR.Frequencies, mIR.DipoleSurfaceList[x], mIR.Ys, X0, True)
            if IncludeZeroth:
                D0 = Dx.diagonal()
                D0 += mIR.DipoleSurface[x][0][0]
                Dx.setdiag(D0)
            mIR.D.append(Dx)
    else:
        Dx = VCISparseHamFromVSCF(mIR.mVCI.Basis, mIR.mVCI.Basis, mIR.Frequencies, mIR.DipoleSurfaceList[xi], mIR.Ys, X0, True)
        if IncludeZeroth:
            D0 = Dx.diagonal()
            D0 += mIR.DipoleSurface[xi][0][0]
            Dx.setdiag(D0)
        mIR.D[xi] = Dx


def GetAb(mIR, w, Basis = None, xi = None):
    if Basis is None:
        H = mIR.mVCI.H
        C = mIR.mVCI.C
        E = mIR.mVCI.E
    else:
        H = GenerateSparseHamV(Basis, mIR.mVCI.Frequencies, mIR.mVCI.PotentialList, mIR.mVCI.Potential[0], mIR.mVCI.Potential[1], mIR.mVCI.Potential[2], mIR.mVCI.Potential[3])
        E, C = sparse.linalg.eigsh(H, k = mIR.mVCI.NStates, which = 'SA')

    #A = np.eye(H.shape[0]) * (w + mIR.mVCI.E[0]) - H + np.eye(H.shape[0]) * mIR.eta * 1.j
    HDiag = H.diagonal()
    HDiag = (w + mIR.mVCI.E[0]) + mIR.eta * 1.j - HDiag
    A = -1 * H
    A.setdiag(HDiag)
    A = A.tocsr()
    b = [None] * 3
    if xi is None:
        for x in range(3):
            b[x] = np.asarray(mIR.D[x] @ C[:, 0]).T
    else:
        b[xi] = np.asarray(mIR.D[xi] @ C[:, 0]).T
    return A, b

def GetAbNMode(mIR, w, Basis = None, xi = None):
    if Basis is None:
        H = mIR.mVCI.H
        C = mIR.mVCI.C
        E = mIR.mVCI.E
    else:
        H = VCISparseHamNModeFromOM(Basis, Basis, mIR.mVCI.Frequencies, mIR.mVCI.mol.V0, mIR.mVCI.mol.onemode_eig, mIR.mVCI.mol.ints[1].tolist(), mIR.mVCI.mol.ints[2].tolist(), True)
        E, C = sparse.linalg.eigsh(H, k = mIR.mVCI.NStates, which = 'SA')

    #A = np.eye(H.shape[0]) * (w + mIR.mVCI.E[0]) - H + np.eye(H.shape[0]) * mIR.eta * 1.j
    HDiag = H.diagonal()
    HDiag = np.array((w + mIR.mVCI.E[0]) + mIR.eta * 1.j - HDiag, dtype = np.cdouble) 
    A = -1 * H
    A = A.astype(np.cdouble)
    A.setdiag(HDiag)
    A = A.tocsr()
    b = [None] * 3
    if xi is None:
        for x in range(3):
            b[x] = np.asarray(mIR.D[x] @ C[:, 0], dtype = np.cdouble).T
    else:
        b[xi] = np.asarray(mIR.D[xi] @ C[:, 0], dtype = np.cdouble).T
    return A, b


def GetAbFromVSCF(mIR, w, Basis = None, xi = None):
    if Basis is None:
        H = mIR.mVCI.H
        C = mIR.mVCI.C
        E = mIR.mVCI.E
    else:
        H = VCISparseHamFromVSCF(Basis, Basis, mIR.mVCI.Frequencies, mIR.mVCI.PotentialList, mIR.mVCI.Potential[0], mIR.mVCI.Potential[1], mIR.mVCI.Potential[2], mIR.mVCI.Potential[3], mIR.Ys, mIR.Xs, True)
        E, C = sparse.linalg.eigsh(H, k = mIR.mVCI.NStates, which = 'SA')

    #A = np.eye(H.shape[0]) * (w + mIR.mVCI.E[0]) - H + np.eye(H.shape[0]) * mIR.eta * 1.j
    HDiag = H.diagonal()
    HDiag = (w + mIR.mVCI.E[0]) + mIR.eta * 1.j - HDiag
    A = -1 * H
    A.setdiag(HDiag)
    A = A.tocsr()
    b = [None] * 3
    if xi is None:
        for x in range(3):
            b[x] = np.asarray(mIR.D[x] @ C[:, 0]).T
    else:
        b[xi] = np.asarray(mIR.D[xi] @ C[:, 0]).T

    return A, b

def SolveAxb(A, b):
    return sparse.linalg.spsolve(A, b)

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

def SpectralScreenBasis(mIR, Ws = None, C = None, eps = 0.01, InitState = None):
    if InitState is None:
        InitState = mIR.InitState
    if Ws is None:
        Ws = mIR.DipoleSurfaceList
    if C is None:
        C = abs(mIR.mVCI.C[:, InitState])

    return mIR.mVCI.ScreenBasis(Ws = Ws, C = C, eps = eps)

def SpectralHCIStep(mIR, w, xi, eps = 0.01, InitState = 0):
    if mIR.SpectralHBMethod == 1: # means perturbing the ground state
        # HB using dipole operator first
        mIR.mVCI.NewBasis, NAdded1 = mIR.SpectralScreenBasis(Ws = mIR.DipoleSurfaceList[xi], C = abs(mIR.mVCI.C[:, InitState]), eps = eps, InitState = InitState)
        mIR.mVCI.Basis += mIR.mVCI.NewBasis
        mIR.mVCI.SparseDiagonalize()
        mIR.mVCI.NewBasis = None
        del mIR.mVCI.NewBasis
        
        # HB using H and b/D
        mIR.GetTransitionDipoleMatrix(xi = xi)
        A, b = mIR.GetAb(w, xi = xi)
        b[xi] = b[xi].ravel()
        x = SolveAxb(A, b[xi])
        D = np.sqrt((w - mIR.mVCI.H.diagonal() + mIR.mVCI.E[0])**2.0 + mIR.eta**2.0)
        bD = abs(b[xi] / D)
        NewBasis, NAdded2 = mIR.SpectralScreenBasis(Ws = mIR.mVCI.PotentialListFull, C = bD, eps = eps, InitState = InitState)
        mIR.mVCI.Basis += NewBasis

    elif mIR.SpectralHBMethod == 2: # means perturbing x
        # HB using dipole operator first
        '''
        mIR.mVCI.NewBasis, NAdded1 = mIR.SpectralScreenBasis(Ws = mIR.DipoleSurfaceList[xi], C = abs(mIR.mVCI.C[:, InitState]), eps = eps, InitState = InitState)
        mIR.mVCI.Basis += mIR.mVCI.NewBasis
        mIR.mVCI.SparseDiagonalize()
        mIR.mVCI.NewBasis = None
        del mIR.mVCI.NewBasis
        '''
        NAdded1 = 0

        # HB using H and x
        mIR.GetTransitionDipoleMatrix(xi = xi)
        A, b = mIR.GetAb(w, xi = xi)
        b[xi] = b[xi].ravel()
        x = SolveAxb(A, b[xi])
        NewBasis, NAdded2 = mIR.SpectralScreenBasis(Ws = mIR.mVCI.PotentialListFull, C = abs(x.real), eps = eps, InitState = InitState)
        mIR.mVCI.Basis += NewBasis

    return NewBasis, NAdded1 + NAdded2

def SpectralHCI(mIR, w, xi):
    CartCoord = ['x', 'y', 'z']
    # First, get new basis based on the dipole operator
    NAdded = len(mIR.mVCI.Basis)
    it = 1
    while (float(NAdded) / float(len(mIR.mVCI.Basis))) > mIR.mVCI.tol:
        mIR.mVCI.NewBasis, NAdded = mIR.SpectralHCIStep(w, xi = xi, eps = mIR.eps1)
        #print("VHCI Iteration", it, "for w =", w, "complete with", NAdded, "new configurations and a total of", len(mIR.mVCI.Basis), flush = True)
        mIR.mVCI.SparseDiagonalize()
        it += 1
        if it > mIR.mVCI.MaxIter:
            raise RuntimeError("VHCI did not converge")
    print("VHCI for coordinate", CartCoord[xi], "converged for w =", w, "with a total of", len(mIR.mVCI.Basis), "configurations.", flush = True)
    mIR.mVCI.NewBasis = None
    del mIR.mVCI.NewBasis

def ResetVCI(mIR):
    del mIR.mVCI.Basis
    del mIR.mVCI.H
    del mIR.mVCI.C
    del mIR.mVCI.E

    mIR.mVCI.Basis = mIR.Basis0.copy()
    mIR.mVCI.H = mIR.H0.copy()
    mIR.mVCI.C = mIR.C0.copy()
    mIR.mVCI.E = mIR.E0.copy()
    
    gc.collect()

def Intensity(mIR, w, state_thr = 1e-6):
    I = np.zeros((3,3))
    # Should define new basis with HCI and then solve VHCI here, be sure to update mVCI object 
    for xi in range(3):
        mIR.SpectralHCI(w, xi = xi)
        mIR.GetTransitionDipoleMatrix(IncludeZeroth = False)
        # Solve for intensity using updated VCI object
        A, b = mIR.GetAb(w)
        x = SolveAxb(A, b[xi])
        x = np.asarray(x).ravel()
        mIR.XString[xi].append(mIR.mVCI.LCLine(0, thr = state_thr, C = np.reshape(abs(x), (x.shape[0], 1))))
        for xj in range(3):
            b[xj] = np.asarray(b[xj]).ravel()
            I[xj, xi] = (1.j * np.dot(b[xj], x)).real / np.pi

            # Include PT2 corrections
            if mIR.DoPT2:
                # Only do it for diagonal elements
                if xj == xi:
                    I[xj, xi] += DoSpectralPT2(x.reshape(x.shape[0], 1), mIR.mVCI.E, mIR.mVCI.C, mIR.mVCI.Basis, mIR.mVCI.PotentialListFull, mIR.mVCI.PotentialList, mIR.mVCI.Potential[0], mIR.mVCI.Potential[1], mIR.mVCI.Potential[2], mIR.mVCI.Potential[3], mIR.DipoleSurfaceList[xi], mIR.mVCI.Ys, mIR.eps2, 1, w, mIR.eta).real / np.pi
        # Reset VCI object
        mIR.ResetVCI()
    return I

def PlotSpectrum(mIR, PlotName, XLabel = "Frequency", YLabel = "Intensity", Title = "IR Spectrum"):
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

    StateFile = open(SaveName + "_states.txt", "w")
    for i in range(mIR.ws.shape[0]):
        StateFile.write('%.6f\t%.12f\tx: %s\ty: %s\tz: %s\n' % (mIR.ws[i], mIR.Is[i], mIR.XString[0][i], mIR.XString[1][i], mIR.XString[2][i]))
        

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
    SpectralScreenBasis = SpectralScreenBasis
    SpectralHCIStep = SpectralHCIStep
    SpectralHCI = SpectralHCI
    ResetVCI = ResetVCI
    Intensity = Intensity
    ApproximateAInv = ApproximateAInv
    PlotSpectrum = PlotSpectrum
    SaveSpectrum = SaveSpectrum

    TestPowerSeries = TestPowerSeries

    def __init__(self, mf, mVCI, FreqRange = [0, 5000], NPoints = 100, eta = 10, NormalModes = None, DipoleSurface = None, SpectralHBMethod = 2, **kwargs):
        self.mf = mf
        self.mVCI = mVCI
        #self.mVCI.HBMethod = 'coupling'
        
        self.SpectralHBMethod = SpectralHBMethod
        
        self.Basis0 = mVCI.Basis.copy()
        self.H0 = mVCI.H.copy()
        self.C0 = mVCI.C.copy()
        self.E0 = mVCI.E.copy()
        self.Frequencies = mVCI.Frequencies
        self.eps1 = mVCI.eps1 / 100
        self.eps2 = self.eps1 / 100
        self.NormalModes = NormalModes
        self.InitState = 0
        self.FreqRange = FreqRange
        self.NPoints = NPoints
        self.eta = eta
        self.Order = 1
        if DipoleSurface is not None:
            self.Order = len(DipoleSurface[0]) - 1
        self.DipoleSurface = DipoleSurface
        self.Normalize = False
        self.DoPT2 = False

        self.__dict__.update(kwargs)

    def kernel(self):
        # Make dipole surface
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
                DSListX += Dn.copy()
            self.DipoleSurfaceList.append(DSListX.copy())

        # Cubic/linear and quadratic/quartic terms must stack, and there must be something at the highest even order, so we will make sure that happens here
        for x in range(3):
            self.DipoleSurface[x][3] += self.DipoleSurface[x][1].copy()
            self.DipoleSurface[x][4] += self.DipoleSurface[x][2].copy()
            self.DipoleSurface[x][1] = []
            self.DipoleSurface[x][2] = []
            self.DipoleSurfaceList[x].append(FConst(0.0, [0] * 6, False))
            self.DipoleSurfaceList[x] = HeatBath_Sort_FC(self.DipoleSurfaceList[x])

        self.GetTransitionDipoleMatrix(IncludeZeroth = False)

        self.ws = np.linspace(self.FreqRange[0], self.FreqRange[1], num = self.NPoints)
        self.ITensors = []
        self.XString = [[], [], []]
        for w in self.ws:
            self.ITensors.append(self.Intensity(w))
        self.Is = []
        for I in self.ITensors:
            self.Is.append(I[0, 0] + I[1, 1] + I[2, 2])

        self.Is = np.asarray(self.Is)
        if self.Normalize:
            self.Is = np.asarray(self.Is) / max(self.Is)

class VSCFLinearResponseIR(LinearResponseIR):
    GetTransitionDipoleMatrix = GetTransitionDipoleMatrixFromVSCF
    GetAb = GetAbFromVSCF

    def __init__(self, mf, mVCI, FreqRange = [0, 5000], NPoints = 100, eta = 10, NormalModes = None, DipoleSurface = None, **kwargs):
        LinearResponseIR.__init__(self, mf, mVCI, FreqRange = FreqRange, NPoints = NPoints, eta = eta, NormalModes = NormalModes, DipoleSurface = DipoleSurface)

        self.Ys = mVCI.Ys
        self.Xs = mVCI.Xs

        self.__dict__.update(kwargs)

class LinearResponseIRNMode(LinearResponseIR):
    GetTransitionDipoleMatrix = GetTransitionDipoleMatrixNMode
    GetAb = GetAbNMode
    
    def __init__(self, mVCI, FreqRange = [0, 5000], NPoints = 100, eta = 10, SpectralHBMethod = 2, NormalModes = None, **kwargs):
        self.mVCI = mVCI
        self.mol = mVCI.mol
        
        self.SpectralHBMethod = SpectralHBMethod
        
        self.Basis0 = mVCI.Basis.copy()
        self.H0 = mVCI.H.copy()
        self.C0 = mVCI.C.copy()
        self.E0 = mVCI.E.copy()
        self.Frequencies = mVCI.Frequencies
        self.eps1 = mVCI.eps1 / 100
        self.eps2 = self.eps1 / 100
        self.NormalModes = NormalModes
        self.InitState = 0
        self.FreqRange = FreqRange
        self.NPoints = NPoints
        self.eta = eta
        self.Normalize = False
        self.DoPT2 = False

        self.__dict__.update(kwargs)

    def kernel(self):
        self.GetTransitionDipoleMatrix(IncludeZeroth = False)

        self.ws = np.linspace(self.FreqRange[0], self.FreqRange[1], num = self.NPoints)
        self.ITensors = []
        self.XString = [[], [], []]
        for w in self.ws:
            self.ITensors.append(self.Intensity(w))
        self.Is = []
        for I in self.ITensors:
            self.Is.append(I[0, 0] + I[1, 1] + I[2, 2])

        self.Is = np.asarray(self.Is)
        if self.Normalize:
            self.Is = np.asarray(self.Is) / max(self.Is)


if __name__ == "__main__":
    from vstr.ff.normal_modes import GetNormalModes
    from vstr.ff.force_field import GetFF
    from vstr.vhci.vhci import VHCI
    from pyscf import gto, scf

    mol = gto.M()
    mol.atom ='''
    H      0.5288      0.1610      0.9359
  C      0.0000      0.0000      0.0000
  H      0.2051      0.8240     -0.6786
  H      0.3345     -0.9314     -0.4496
  H     -1.0685     -0.0537      0.1921
  '''
    
    '''
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

    from vstr.spectra.ir_exact import IRSpectra

    mVHCI = VHCI(w, V, MaxQuanta = 10, MaxTotalQuanta = 3, eps1 = 10, eps2 = 0.01, eps3 = -1, NWalkers = 50, NSamples = 50, NStates = 10)
    mVHCI.kernel()

    mIR = IRSpectra(mf, mVHCI, NormalModes = NormalModes, Order = 2)
    mIR.kernel()
    mIR.PlotSpectrum("water_spectrum_lr.png", L = 100, XMin = 0, XMax = 5000)

    mVHCI = VHCI(w, V, MaxQuanta = 10, MaxTotalQuanta = 1, eps1 = 100, eps2 = 0.001, eps3 = -1, NWalkers = 50, NSamples = 50, NStates = 1)
    mVHCI.kernel()
    #mVHCI.E, mVHCI.C = np.linalg.eigh(mVHCI.H.todense())
    #mVHCI.E_HCI = mVHCI.E

    '''
    mIR2 = LinearResponseIR(mf, mVHCI, FreqRange = [0, 5000], NPoints = 100, eps1 = 1, eta = 100, DipoleSurface = mIR.DipoleSurface, Order = 2, SpectralHBMethod = 1)
    mIR2.kernel()
    mIR2.PlotSpectrum("water_spectrum_lr.png")
    '''

    #mIR2 = LinearResponseIR(mf, mVHCI, FreqRange = [0, 5000], NPoints = 100, eps1 = 1, eta = 100, DipoleSurface = mIR.DipoleSurface, Order = 2)
    #mIR2.kernel()
    #mIR2.PlotSpectrum("water_spectrum_lr.png")
    #mIR2.SaveSpectrum("water_spectrum")


    mIR2 = LinearResponseIR(mf, mVHCI, FreqRange = [0, 5000], NPoints = 100, eps1 = 1, eps2 = 0.001, DoPT2 = True, eta = 100, DipoleSurface = mIR.DipoleSurface, Order = 2)
    mIR2.kernel()
    mIR2.PlotSpectrum("water_spectrum_lr.png")
    mIR2.SaveSpectrum("water_spectrum")

    '''
    from vstr.ci.vci import VCI
    from vstr.mf.vscf import VSCF
    vmf = VSCF(w, V, MaxQuanta = 10, NStates = 1)
    vmf.kernel()
    mVCI = VCI(vmf, MaxTotalQuanta = 1, eps1 = 10, eps2 = 0.001, eps3 = -1, NWalkers = 50, NSamples = 50, NStates = 1)
    mVCI.kernel()

    mVSCFIR = VSCFLinearResponseIR(mf, mVCI, FreqRange = [0, 5000], NPoints = 1000, eps1 = 0.01, eta = 100, NormalModes = NormalModes, Order = 2)
    mVSCFIR.kernel()
    mVSCFIR.PlotSpectrum("water_spectrum_lr.png")
    '''
