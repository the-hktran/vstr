#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <string>
#include <tuple>
#include <algorithm>
#include <cmath>
#include <unordered_map>

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


std::vector<std::tuple<std::vector<int>, double>> BasisConnectionsCPP(std::vector<int>& BasisFunction, std::vector<int>& Qs);
std::vector<std::vector<std::tuple<std::vector<int>, double>>> FormBasisConnectionsCPP(std::vector<std::vector<VDeriv>>& Ws, std::vector<std::vector<int>>& Basis);
std::string VectorToString(std::vector<int> V);
Eigen::MatrixXd HamVCPP(std::vector<std::vector<int>> Basis, std::vector<std::vector<std::tuple<std::vector<int>, double>>> BasisConn, std::vector<std::vector<int>> BasisBras, std::vector<double> Freq, std::vector<std::vector<VDeriv>> Ws, bool DiagonalOnly, bool OffDiagonal);
Eigen::SparseMatrix<double> SpHamVCPP(std::vector<std::vector<int>>& Basis, std::vector<std::vector<std::tuple<std::vector<int>, double>>>& BasisConn, std::vector<std::vector<int>>& BasisBras, std::vector<double>& Freq, std::vector<std::vector<VDeriv>>& Ws, bool DiagonalOnly, bool OffDiagonal);
