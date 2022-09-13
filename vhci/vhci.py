import numpy as np
from functools import reduce

def FormW(mVHCI, V):
    Order = len(V.shape)
    Modes = mVHCI.w.shape[0]
    ws = []
    for i in range(Order):
        ws.append(mVHCI.w)
    wProd = reduce(np.multiply, np.ix_(*ws))
    wProd *= 2**Order
    wProd = wProd**-0.5

    return V * wProd

'''
Forms the explicit vibrational Hamiltonian in the basis given by Ds
'''
def HamV(mVHCI):
    N = len(mVHCI.Ds)
    mVHCI.H = np.zeros((N, N))
    for i in range(N):
        for j, n in enumerate(mVHCI.Ds[i]):
            mVHCI.H[i, i] += (n + 0.5) * mVHCI.w[j] # HO diagonal elements

    # Now we loop through each V
    for D in mVHCI.Ds:
        for W in mVHCI.Ws:
            WShape = W.shape
            WN = np.prod(W.shape)
            WFlat = W.reshape(WN)
            ConnectedD = []
            for i in range(WN):
                WElement = WFlat[i]
                QIndices = np.unravel_index(i, WShape)
                Conn0 = [D.copy()]
                Conn1 = []
                for q in QIndices:
                    for Con0 in Conn0:
                        Co0p = Con0.copy()
                        Co0p[q] += 1
                        Conn1.append(Co0p)
                        Co0m = Con0.copy()
                        Co0m[q] -= 1
                        if Co0m[q] < 0:
                            continue
                        Conn1.append(Co0m)
                    Conn0 = Conn1.copy()
                    Conn1 = []
                for Co0 in Conn0:
                    if Co0 not in Conn1:
                        Conn1.append(Co0)
                Conn0 = Conn1
                for DConn in Conn0:
                    try:
                        j = mVHCI.Ds.index(DConn)
                    except:
                        continue
                    Hij = 1
                    for Di, Dv in enumerate(D):
                        if D[Di] != DConn[Di]:
                            n = max(D[Di], DConn[Di])
                            Hij *= np.sqrt(n)
                    mVHCI.H[i, j] = Hij
                    mVHCI.H[j, i] = Hij
            
'''
Class that handles VHCI
'''
class VHCI:
    FormW = FormW
    HamV = HamV

    def __init__(self, Vs, w, Ds):
        self.Vs = Vs # PES Derivatives
        self.w = w # Harmonic Frequencies
        self.Ws = []
        self.Ds = Ds
        for V in self.Vs:
            self.Ws.append(self.FormW(V))

    def kernel(self):
        pass

if __name__ == "__main__":
    V2 = np.asarray([[1, 1], [1, 1]])
    Ds = [[0, 0], [1, 0], [2, 0], [0, 1], [1, 1], [1, 2], [2, 0], [2, 1], [2, 2]]
    w = np.asarray([10, 20])

    mVHCI = VHCI([V2], w, Ds)
    mVHCI.HamV()
    print(mVHCI.H)
