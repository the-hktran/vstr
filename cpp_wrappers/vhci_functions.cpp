// #include <Eigen/Dense>
// #include <Eigen/Sparse>
#include <vector>
#include <string>
#include <tuple>
#include <algorithm>

class VDeriv
{
    public:
        double W; // Value of force constant
        std::vector<int> QIndices; // Index of derivatives
        int Order;
        std::vector<int> QUnique;
        std::vector<int> QPowers;

        VDeriv(double, std::vector<int>);
}

VDeriv::VDeriv (double dV, std::vector<int> Qs)
{
    W = dV;
    QIndices = Qs;
    Order = QIndices.size();
    
    QUnique = QIndices;
    std::sort(QUnique.begin(), QUnique.end());
    auto Last = std::unique(QUnique.begin(), QUnique.end());
    QUnique.erase(Last, QUnique.end());
    for (int i : QUnique)
    {
        QPowers.push_back(std::count(QIndices.begin(), QIndices.edn(), i));
    }
}

std::vector<std::tuple<std::vector<int>, double>> BasisConnection(std::vector<std::vector<double, std::vector<int>>>> Ws, std::vector<int> BasisFunction)
{
    std::vector<std::tuple<std::vector<int>, double>> Conn
    Conn.push_back(std::tuple<std::vector<int>, double>(BasisFunction, 1.0))
    std::vector<std::tuple<std::vector<int>, double>> tmpC


    

}
std::vector<std::vector<std::tuple<std::vector<int>, double>>> FormBasisConnections(std::vector<std::vector<std::vector<double, std::vector<int>>>> Ws, std::vector<std::vector<int>> Basis)
{
}
