import numpy as np
from vstr.cpp_wrappers.vhci_jf.vhci_jf_functions import FConst, ProdU, SetUij, SetUs
from vstr.mf.vscf import VSCF
from itertools import permutations

def ContractFC(mCO, U):
    FCs = []
    N = mCO.mf.NModes
    # First, make tensors
    VTensor3 = np.zeros((N, N, N))
    VTensor4 = np.zeros((N, N, N, N))
    VTensor5 = np.zeros((N, N, N, N, N))
    VTensor6 = np.zeros((N, N, N, N, N, N))
    VTensor = [VTensor3, VTensor4, VTensor5, VTensor6]

    # Fill them in
    for FC in mCO.PotentialList:
        Is = list(permutations(FC.QIndices))
        for I in Is:
            VTensor[FC.Order - 3][I] = FC.fc
    
    # Contract the anharmonic tensors
    VTensor[0] = np.einsum('ijk,ia,jb,kc->abc', VTensor[0], U, U, U, optimize = True)
    VTensor[1] = np.einsum('ijkl,ia,jb,kc,ld->abcd', VTensor[1], U, U, U, U, optimize = True)
    
    for i in range(N):
        for j in range(N):
            VTensor[2][i, j] = np.einsum('klm,kc,ld,me->cde', VTensor[2][i, j], U, U, U, optimize = True)
    for c in range(N):
        for d in range(N):
            for e in range(N):
                VTensor[2][:, :, c, d, e] = U.T @ VTensor[2][:, :, c, d, e] @ U
    #VTensor[2] = np.einsum('ijklm,ia,jb,kc,ld,me->abcde', VTensor[2], U, U, U, U, U,  optimize = True)
    for i in range(N):
        for j in range(N):
            for k in range(N):
                VTensor[3][i, j, k] = np.einsum('lmn,ld,me,nf->def', VTensor[3][i, j, k], U, U, U, optimize = True)
    for d in range(N):
        for e in range(N):
            for f in range(N):
                VTensor[3][:, :, :, d, e, f] = np.einsum('ijk,ia,jb,kc->abc', VTensor[3][:, :, :, d, e, f], U, U, U, optimize = True)
    #VTensor[3] = np.einsum('ijklmn,ia,jb,kc,ld,me,nf->abcdef', VTensor[3], U, U, U, U, U, U, optimize = True)

    # Make a new list of FCs
    for i in range(N):
        for j in range(i, N):
            for k in range(j, N):
                if abs(VTensor[0][i, j, k]) > 1e-12:
                    FCs.append(FConst(VTensor[0][i, j, k], [i, j, k], False))

    for i in range(N):
        for j in range(i, N):
            for k in range(j, N):
                for l in range(k, N):
                    if abs(VTensor[1][i, j, k, l]) > 1e-12:
                        FCs.append(FConst(VTensor[1][i, j, k, l], [i, j, k, l], False))
    for i in range(N):
        for j in range(i, N):
            for k in range(j, N):
                for l in range(k, N):
                    for m in range(l, N):
                        if abs(VTensor[2][i, j, k, l, m]) > 1e-12:
                            FCs.append(FConst(VTensor[2][i, j, k, l, m], [i, j, k, l, m], False))

    for i in range(N):
        for j in range(i, N):
            for k in range(j, N):
                for l in range(k, N):
                    for m in range(l, N):
                        for n in range(m, N):
                            if abs(VTensor[3][i, j, k, l, m, n]) > 1e-12:
                                FCs.append(FConst(VTensor[0][i, j, k, l, m, n], [i, j, k, l, m, n], False))

    return FCs

def E_SCF(mCO, U):
    mf_tmp = VSCF(mCO.mf.Frequencies, [], MaxQuanta = mCO.mf.MaxQuanta, verbose = 0)
    NewPotentialList = mCO.ContractFC(U)
    mf_tmp.UpdateFC(NewPotentialList)
    return mf_tmp.kernel()

def E_SCF_ij(mCO, ij, theta):
    NewUs = mCO.Us.copy()
    NewUs[ij] = SetUij(theta)
    U = ProdU(NewUs, mCO.NModes)
    return mCO.E_SCF(U)

def F(mCO, t, fn, tn):
    F = 0.0
    dF = 0.0
    ddF = 0.0

    # x0
    for f in fn:
        F += f

    # The rest
    for m in range(mCO.p):
        for n in range(len(fn)):
            F += fn[n] * np.cos(4*((m + 1) * t - (n - mCO.p) * tn[n]))
            dF -= fn[n] * 4 * (m + 1) * np.sin(4 * ((m + 1) * t - (n - mCO.p) * tn[n]))
            ddF -= fn[n] * (4 * (m + 1))**2 * np.cos(4 * ((m + 1) * t - (n - mCO.p) * tn[n]))

    F *= 1 / (2 * mCO.p + 1)
    dF *= 1 / (2 * mCO.p + 1)
    ddF *= 1 / (2 * mCO.p + 1)
    return F, dF, ddF

def OptF(mCO, fn, tn, thr = 1e-6):
    t = 0
    F, dF, ddF = mCO.F(t, fn, tn)
    it = 1
    while abs(dF) > thr:
        t = t - (dF / ddF)
        F, dF, ddF = mCO.F(t, fn, tn)
        print('\tFourier minimum NR iteration %d complete with value %.6f and error %.12f' % (it, t, abs(dF)), flush = True)
        it += 1
        if it > 20:
            print('\t!! Fourier minimum failed to converge !!', flush = True)
            t = 0
            break
    return t

def OptE(mCO, t, ij, fn, tn, dt = 1e-3, thr = 1e-6):
    E = mCO.E_SCF_ij(ij, t)
    Em = mCO.E_SCF_ij(ij, t - dt)
    Ep = mCO.E_SCF_ij(ij, t + dt) 
    dE = (Ep - Em) / (2 * dt)
    ddE = (Ep - 2 * E + Em) / (dt * dt)
    it = 1
    while abs(dE) > thr:
        t = t - dE / ddE
        E = mCO.E_SCF_ij(ij, t)
        Em = mCO.E_SCF_ij(ij, t - dt)
        Ep = mCO.E_SCF_ij(ij, t + dt) 
        dE = (Ep - Em) / (2 * dt)
        ddE = (Ep - 2 * E + Em) / (dt * dt)
        print('\tExact minimum NR iteration %d complete with energy %.6f and value %.6f and error %.12f' % (it, E, t, abs(dE)), flush = True)
        it += 1
    return t

def JacobiSweepIteration(mCO):
    # Optimize each Uij individually
    tn = []
    for n in range(-mCO.p, mCO.p + 1):
        tn.append(n * np.pi / (2 * (2 * mCO.p + 1)))

    OptUs = []
    thetas = []
    ij = 0
    # ij = (i - 2) (i - 1) / 2 + j
    for i in range(mCO.mf.NModes):
        for j in range(i + 1, mCO.mf.NModes):
            fn = []
            for tp in tn:
                fn.append(mCO.E_SCF_ij(ij, tp))
            t0 = mCO.OptF(fn, tn)
            t = mCO.OptE(t0, ij, fn, tn)
            U = SetUij(t)
            OptUs.append(U)
            thetas.append(t)
            ij += 1
    mCO.Us = OptUs
    mCO.thetas = thetas
    mCO.U = ProdU(OptUs, mCO.mf.NModes)

def JacobiSweep(mCO, thr = 1e-6):
    Ei = 0
    Ef = mCO.mf.ESCF
    it = 1
    while abs(Ef - Ei) > 1e-6:
        Ei = Ef
        mCO.JacobiSweepIteration()
        print(mCO.Us)
        print(mCO.U)
        Ef = mCO.E_SCF(mCO.U)
        print('== Jacobi Sweep iteration %d complete with SCF Energy %.6f and error %.12f' % (it, Ef, Ef - Ei), flush = True)
        it += 1

def InitU(mCO):
    mCO.Us = []
    for i in range(mCO.NModes):
        for j in range(i + 1, mCO.NModes):
            mCO.Us.append(np.eye(2))
    mCO.U = np.eye(mCO.NModes)

def MakeOCVSCF(mCO, U = None):
    if U is None:
        U = mCO.U
    mf_tmp = VSCF(mCO.mf.Frequencies, [], MaxQuanta = mCO.mf.MaxQuanta)
    NewPotentialList = mCO.ContractFC(U)
    mf_tmp.UpdateFC(NewPotentialList)
    #mf_tmp.HamHO = U.T @ mf_tmp.HamHO @ U
    mf_tmp.kernel()
    return mf_tmp

class CoordinateOptimizer:
    InitU = InitU
    JacobiSweepIteration = JacobiSweepIteration
    JacobiSweep = JacobiSweep
    OptF = OptF
    OptE = OptE
    F = F
    E_SCF_ij = E_SCF_ij
    E_SCF = E_SCF
    ContractFC = ContractFC

    MakeOCVSCF = MakeOCVSCF

    def __init__(self, mf):
        self.mf = mf
        self.PotentialList = []
        for Wp in self.mf.Potential:
            self.PotentialList += Wp
        self.p = 3
        self.NModes = mf.NModes
    
    def kernel(self):
        self.InitU()
        self.JacobiSweep()

if __name__ == "__main__":
    from vstr.utils.read_jf_input import Read
    w, MaxQuanta, MaxTotalQuanta, Vs, eps1, eps2, eps3, NWalkers, NSamples, NStates = Read('test.inp')
    mf = VSCF(w, Vs, MaxQuanta = MaxQuanta, NStates = NStates, SlowV = False)
    mf.kernel()
    mCO = CoordinateOptimizer(mf)
    mCO.kernel()

