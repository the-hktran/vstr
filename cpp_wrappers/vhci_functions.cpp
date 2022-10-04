// #include <Eigen/Dense>
// #include <Eigen/Sparse>
#include <vector>
#include <string>
#include <tuple>
#include <algorithm>
#include <cmath>

// Class that holds values for 
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

std::vector<std::tuple<std::vector<int>, double>> BasisConnectionCPP(std::vector<int> BasisFunction, std::vector<int> Qs)
{
    std::vector<std::tuple<std::vector<int>, double>> Conn
    Conn.push_back(std::tuple<std::vector<int>, double>(BasisFunction, 1.0))
    std::vector<std::tuple<std::vector<int>, double>> Conn2
    for (int i = 0; i < Qs.size(); i++)
    {
        for (std::tuple<std::vector<int>, double> C : Conn)
        {
            std::vector<int> tmpQ = std::get<0>(C);
            double tmpCoeff = std::get<1>(C);
            tmpQ[i] += 1;
            tmpCoeff *= sqrt(tmpQ[i]);
            std::tuple<std::vector<int>, double> tmpC = std::make_tuple(tmpQ, tmpCoeff);
            Conn2.push_back(tmpC);

            std::vector<int> tmpQ = std::get<0>(C);
            double tmpCoeff = std::get<1>(C);
            tmpCoeff *= sqrt(tmpQ[i]);
            tmpQ[i] -= 1;
            if (tmpQ[i] < 0) continue;
            std::tuple<std::vector<int>, double> tmpC = std::make_tuple(tmpQ, tmpCoeff);
            Conn2.push_back(tmpC);
        }
        std::vector<std::tuple<std::vector<int>, double> Conn = Conn2;
        std::vector<std::tuple<std::vector<int>, double> Conn2;
    }
    return Conn;
}

std::vector<std::vector<std::tuple<std::vector<int>, double>>> FormBasisConnectionsCPP(std::vector<std::vector<std::vector<VDeriv>>> Ws, std::vector<std::vector<int>> Basis)
{
    std::vector<std::vector<std::tuple<std::vector<int>, double>>> BasisConnections;
    for (std::vector<int> B : Basis)
    {
        std::vector<std::tuple<std::vector<int>, double>> ConnB;
        for (int p = 0; p < Ws.size(); p++)
        {
            for (VDeriv W : Ws[p])
            {
                std::vector<std::tuple<std::vector<int>, double> ConnByW = BasisConnectionsCPP(B, W.QIndices);
                for (std::tuple<std::vector<int>, double>> C : ConnByW)
                {
                    std::get<0>(C) *= W.W;
                }
                ConnB.insert(ConnB.end(), C.begin(), C.end());
            }
        }
        BasisConnections.push_back(ConnB);
    }
    return BasisConnections; 
}
