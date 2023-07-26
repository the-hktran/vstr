import numpy as np
import sys

def ReadHessian(FilePath, NAtom):
    N = int(NAtom * 3)
    NElements = int(int(N * (N + 1)) / 2)
    Hessian = np.zeros((N, N))
    with open(FilePath, 'r') as f:
        I = 0
        J = 0
        for i in range(NElements):
            line = f.readline()
            #line = line.split()
            #I = int(line[0]) - 1
            #J = int(line[1]) - 1
            Hessian[I, J] = float(line)
            Hessian[J, I] = float(line)
            J += 1
            if J == N:
                I += 1
                J = I            
    return Hessian

def GetAtomList(N, psfFile):
    NStr = str(N)
    AtomList = []
    # read line of file one by one
    with open(psfFile, 'r') as f:
        for line in f:
            try:
                if line.split()[0] == NStr:
                    for i in range(N):
                        line = f.readline()
                        AtomList.append(line.split()[5][0])
                    break
            except:
                pass
    return AtomList

if __name__ == '__main__':
    '''
    FilePath = sys.argv[1]
    NAtom = int(sys.argv[2])
    Hessian = ReadHessian(FilePath, NAtom)
    E, V = np.linalg.eigh(Hessian)
    print(np.diag(Hessian))
    for e in E:
        print(e)
    print('evec')
    for v in V[:,-1]:
        print(v)
    '''

    psfFile = sys.argv[1]
    N = int(sys.argv[2])
    AtomList = GetAtomList(N, psfFile)
    print(AtomList)
    print(AtomList.count('C'))
    print(AtomList.count('H'))
    print(AtomList.count('N'))
    print(AtomList.count('O'))
    print(AtomList.count('S'))
