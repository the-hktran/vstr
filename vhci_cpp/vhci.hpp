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
        std::vector<std::vector<Vderiv>> PotentialWithSD;
        std::vector<double> Frequencies; // Harmonic frequencies
        int NModes;

        std::vector<int> MaxQuanta;
        int MaxTotalQuanta;
        std::vector<std::vector<int>> Basis;
        std::vector<std::vector<std::vector<int>>> BasisConn;
        
        double Tolerance = 0.01; // Threshold for HCI termination condition
        double Epsilon1 = 0.01; // Threshold for basis screening
        double Epsilon2 = 0.0001; // Threshold for PT2 basis screening
        double Epsilon3 = 0.000001; // Threshold for Stochastic PT2
        int MaxIteration = 1000;
        int NStates;
        int NWalkers;

        int NCPUs;
        std::fstream InputFile;
        std::fstream OutputFile;

        std::vector<std::vector<VDeriv>> FormW(std::vector<std::vector<VDeriv>>& V);
        std::vector<std::vector<VDeriv>> FormWSD();
        std::vector<std::vector<int>> ScreenBasis(std::vector<std::vector<VDeriv>>& Potential, Eigen::VectorXd& C, double Epsilon); 
        int HCIStep(double Tolerance);
        void HCI();
        void PT2();
        void Diagonalize();
        void InitTruncatedBasis();
        void ReadInput();
        VHCI(std::string);

        void ReadArgs(int, char*, std::fstream&, std::fstream&, std::string&);
        void ReadInput(std::fstream&);
};


std::vector<std::tuple<std::vector<int>, double>> BasisConnectionsCPP(std::vector<int>& BasisFunction, std::vector<int>& Qs);
std::vector<std::vector<std::tuple<std::vector<int>, double>>> FormBasisConnectionsCPP(std::vector<std::vector<VDeriv>>& Ws, std::vector<std::vector<int>>& Basis);
std::string VectorToString(std::vector<int> V);
Eigen::MatrixXd HamVCPP(std::vector<std::vector<int>> Basis, std::vector<std::vector<std::tuple<std::vector<int>, double>>> BasisConn, std::vector<std::vector<int>> BasisBras, std::vector<double> Freq, std::vector<std::vector<VDeriv>> Ws, bool DiagonalOnly, bool OffDiagonal);
Eigen::SparseMatrix<double> SpHamVCPP(std::vector<std::vector<int>>& Basis, std::vector<std::vector<std::tuple<std::vector<int>, double>>>& BasisConn, std::vector<std::vector<int>>& BasisBras, std::vector<double>& Freq, std::vector<std::vector<VDeriv>>& Ws, bool DiagonalOnly, bool OffDiagonal);
