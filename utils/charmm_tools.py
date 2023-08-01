import numpy as np
import sys
from vstr.utils import constants

def GetCoordFromCRD(CRDFile, NAtom):
    Coord = np.zeros(3 * NAtom)
    with open(CRDFile, 'r') as f:
        for line in f:
            try:
                if line.split()[0] == str(NAtom):
                    for i in range(NAtom):
                        line = f.readline()
                        Coord[3 * i] = float(line.split()[4])
                        Coord[3 * i + 1] = float(line.split()[5])
                        Coord[3 * i + 2] = float(line.split()[6])
                    break
            except:
                pass
    return Coord * constants.ANGSTROM_TO_AU

def WriteNewCRD(OrigCRDFile, NewCRDFile, NAtom, Coord):
    with open(OrigCRDFile, 'r') as f:
        with open(NewCRDFile, 'w') as g:
            for line in f:
                if line.split()[0] == str(NAtom):
                    g.write(line)
                    for i in range(NAtom):
                        line = f.readline()
                        g.write('{}\t{}\t{}\t{}\t{:16.8f}{:16.8f}{:16.8f}\t{}\t{}\t{}\n'.format(line.split()[0], line.split()[1], line.split()[2], line.split()[3], Coord[3 * i], Coord[3 * i + 1], Coord[3 * i + 2], line.split()[7], line.split()[8], line.split()[9]))
                    break
                else:
                    g.write(line)

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
    # Hessian is in kcal/mol/A^2, convert to a.u.
    Hessian = Hessian * constants.KCAL_TO_AU / constants.ANGSTROM_TO_AU**2
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

def GetAtomMassList(AtomList):
    AtomMassList = []
    for atom in AtomList:
        if atom == 'C':
            AtomMassList.append(12.011)
        elif atom == 'H':
            AtomMassList.append(1.008)
        elif atom == 'N':
            AtomMassList.append(14.007)
        elif atom == 'O':
            AtomMassList.append(15.999)
        elif atom == 'S':
            AtomMassList.append(32.06)
        else:
            print('Unknown atom type')
            sys.exit()
    return np.asarray(AtomMassList)

def MassWeightHessian(Hessian, AtomMassList):
    sqrtmass = np.repeat(AtomMassList, 3)**(-0.5)
    return np.outer(sqrtmass, sqrtmass) * Hessian

def GetNormalModes(Hessian):
    E, V = np.linalg.eigh(Hessian)
    Frequencies = np.sqrt(E) / np.sqrt(constants.AMU_TO_ME) / (2 * np.pi * constants.C_AU * constants.BOHR_TO_CM)
    return Frequencies, V

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

    '''
    psfFile = sys.argv[1]
    N = int(sys.argv[2])
    HessFile = sys.argv[3]

    AtomList = GetAtomList(N, psfFile)
    print(AtomList)
    print(AtomList.count('C'))
    print(AtomList.count('H'))
    print(AtomList.count('N'))
    print(AtomList.count('O'))
    print(AtomList.count('S'))
    AtomMassList = GetAtomMassList(AtomList)
    MolarMass = np.sum(AtomMassList)
    Hessian = ReadHessian(HessFile, N)
    Hessian = MassWeightHessian(Hessian, AtomMassList)

    E, V = GetNormalModes(Hessian)
    for e in E:
        print(e)
    '''

    '''
    crdFile = sys.argv[1]
    newFile = sys.argv[2]
    NAtom = int(sys.argv[3])
    Coord = GetCoordFromCRD(crdFile, NAtom)
    Coord /= constants.ANGSTROM_TO_AU
    print(Coord)
    Coord = Coord + 10.0
    print(Coord)
    WriteNewCRD(crdFile, newFile, NAtom, Coord)
    '''
