#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <string>
#include <tuple>
#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <fstream>
#include <iostream>
#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/SparseSymMatProd.h>

class VDeriv
{
    public:
        double W; // Value of force constant
        std::vector<int> QIndices; // Index of derivatives
        int Order;
        std::vector<int> QUnique;
        std::vector<int> QPowers;
        bool doScaleW;

        VDeriv(double, std::vector<int>, bool);
        void ScaleW();
};

class BasisState
{
    public:
        std::vector<int> Quanta;
        double Coeff;

        BasisState(std::vector<int>);
        void Raise(int);
        void Lower(int);
};

class VHCI
{
    public:
        std::vector<std::vector<VDeriv>> Potential; // The scaled anharmonic potential terms
        std::vector<std::vector<VDeriv>> PotentialSD;
        std::vector<std::vector<VDeriv>> PotentialWithSD;
        std::vector<double> Frequencies; // Harmonic frequencies
        int NModes;
        Eigen::MatrixXd C;
        Eigen::VectorXd E;

        std::vector<int> MaxQuanta;
        int MaxTotalQuanta;
        std::vector<std::vector<int>> Basis;
        std::vector<std::vector<std::tuple<std::vector<int>, double>>> BasisConn;
        
        double Tolerance = 0.01; // Threshold for HCI termination condition
        double Epsilon1 = 0.01; // Threshold for basis screening
        double Epsilon2 = 0.0001; // Threshold for PT2 basis screening
        double Epsilon3 = 0.000001; // Threshold for Stochastic PT2
        int MaxIteration = 1000;
        int NStates;
        int NStatesOriginal = 0;
        int NWalkers;

        int NCPUs;
        std::fstream InputFile;
        std::fstream OutputFile;

        std::vector<std::vector<VDeriv>> FormW(std::vector<std::vector<VDeriv>>& V);
        void FormWSD();
        std::vector<std::vector<int>> ScreenBasis(std::vector<std::vector<VDeriv>>& Potential, Eigen::VectorXd& C, double Epsilon); 
        int HCIStep(double Tolerance);
        void HCI();
        void PT2();
        void Diagonalize();
        void SparseDiagonalize();
        void InitTruncatedBasis();
        void ReadInput();
        void RunVHCI();
        VHCI(int, char**);

        void ReadArgs(int, char**);
        void ReadInput(std::fstream&);
};


std::vector<std::tuple<std::vector<int>, double>> BasisConnectionsCPP(std::vector<int>& BasisFunction, std::vector<int>& Qs);
std::vector<std::vector<std::tuple<std::vector<int>, double>>> FormBasisConnectionsCPP(std::vector<std::vector<VDeriv>>& Ws, std::vector<std::vector<int>>& Basis);
std::string VectorToString(std::vector<int> V);
Eigen::MatrixXd HamVCPP(std::vector<std::vector<int>>& Basis, std::vector<std::vector<std::tuple<std::vector<int>, double>>>& BasisConn, std::vector<std::vector<int>>& BasisBras, std::vector<double>& Freq, std::vector<std::vector<VDeriv>>& Ws, bool DiagonalOnly, bool OffDiagonal);
Eigen::SparseMatrix<double, 0, ptrdiff_t> SpHamVCPP(std::vector<std::vector<int>>& Basis, std::vector<std::vector<std::tuple<std::vector<int>, double>>>& BasisConn, std::vector<std::vector<int>>& BasisBras, std::vector<double>& Freq, std::vector<std::vector<VDeriv>>& Ws, bool DiagonalOnly, bool OffDiagonal);

template <typename T>
std::vector<size_t> SortIndices(const std::vector<T>& V)
{
    std::vector<size_t> Indices(V.size());
    std::iota(Indices.begin(), Indices.end(), 0);

    // This sorts smallest to largest
    std::stable_sort(Indices.begin(), Indices.end(),
        [&V](size_t i1, size_t i2) {return V[i2] > V[i1];});

    return Indices;
}

template <typename T>
bool VectorContains(std::vector<T> V, const T& Elem)
{
    bool Result = false;
    if (std::find(V.begin(), V.end(), Elem) != V.end()) Result = true;
    return Result;
}

#include "input.cpp"
