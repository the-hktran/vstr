'''
This is a quick script to translate Jonathan Fetherolf's input files for my code.
Details on the original file can be found here: https://github.com/berkelbach-group/VHCI
'''
import numpy as np

def Read(FilePath):
    f = open(FilePath, 'r')

    Comment = f.readline() # don't need this
    eps1 = float(f.readline().split()[1])
    NStates = int(f.readline().split()[1])
    doPT2 = bool(f.readline().split()[1])
    eps2 = float(f.readline().split()[1])
    MaxQuanta = int(f.readline().split()[1])
    NModes = int(f.readline().split()[1])

    ws = [None] * NModes
    MaxNs = [None] * NModes
    for i in range(NModes):
        a, ws[i], MaxNs[i] = f.readline().split()
    for i in range(NModes):
        ws[i] = float(ws[i])
        MaxNs[i] = int(MaxNs[i]) + 1
    ws = np.asarray(ws)

    NFC = int(f.readline().split()[1])
    Vs = [None] * 6 # will need to remove unused orders later
    for i in range(6):
        #shape = [NModes] * (i + 1)
        Vs[i] = []
    for i in range(NFC):
        FCLine = f.readline().split()
        order = int(FCLine[0])
        #I = ()
        #for j in range(order):
        #    I += (int(FCLine[j + 1]),)
        #Vs[order - 1][I] = FCLine[-1]
        I = []
        for j in range(order):
            I.append(int(FCLine[j + 1]))
        Vs[order].append((float(FCLine[-1]), I))
    VsFinal = []
    for i in range(6):
        #shape = [NModes] * (i + 1)
        #if not np.allclose(np.zeros(shape), Vs[i]):
        #    VsFinal.append(Vs[i])
        if len(Vs[i]) != 0:
            VsFinal.append(Vs[i])

    return ws, MaxNs, MaxQuanta, VsFinal, eps1, eps2, NStates

if __name__ == '__main__':
    Read('../../VHCI/input/C2H4.inp')
