#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <string>
#include <tuple>
#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <map>
#include "vhci.hpp"
#include <iostream>

VDeriv::VDeriv (double dV, std::vector<int> Qs, bool doScale)
{
    W = dV;
    QIndices = Qs;
    Order = QIndices.size();
    doScaleW = doScale;
    
    QUnique = QIndices;
    std::sort(QUnique.begin(), QUnique.end());
    auto Last = std::unique(QUnique.begin(), QUnique.end());
    QUnique.erase(Last, QUnique.end());
    for (int i : QUnique)
    {
        QPowers.push_back(std::count(QIndices.begin(), QIndices.end(), i));
    }
    if (doScaleW) ScaleW();
}

double Factorial(int n)
{
    double nFac = 1.0;
    for (int i = 0; i < n; i++)
    {
        nFac *= (double)(i + 1);
    }
    return nFac;
}

void VDeriv::ScaleW()
{
    W = W / sqrt(pow(2.0, Order));
    for (int i = 0; i < QPowers.size(); i++)
    {
        W /= Factorial(QPowers[i]);
    }
}

BasisState::BasisState (std::vector<int> Q)
{
    Quanta = Q;
    Coeff = 1.0;
}

void BasisState::Raise(int i)
{
    Quanta[i]++;
    Coeff *= sqrt(Quanta[i]);
}

void BasisState::Lower(int i)
{
    Coeff *= sqrt(Quanta[i]);
    Quanta[i]--;
}

void PrintVec(std::vector<int> V)
{
    for (int i = 0; i < V.size(); i++)
    {
        std::cout << V[i] << " ";
    }
    std::cout << std::endl;
}

std::vector<std::tuple<std::vector<int>, double>> BasisConnectionsCPP(std::vector<int>& BasisFunction, std::vector<int>& Qs)
{
    std::vector<std::tuple<std::vector<int>, double>> Conn;
    Conn.push_back(std::tuple<std::vector<int>, double>(BasisFunction, 1.0));
    std::vector<std::tuple<std::vector<int>, double>> Conn2;
    for (int i = 0; i < Qs.size(); i++)
    {
        for (std::tuple<std::vector<int>, double> C : Conn)
        {
            std::vector<int> tmpQR = std::get<0>(C);
            double tmpCoeffR = std::get<1>(C);
            tmpQR[Qs[i]] += 1;
            tmpCoeffR *= sqrt(tmpQR[Qs[i]]);
            std::tuple<std::vector<int>, double> tmpCR = std::make_tuple(tmpQR, tmpCoeffR);
            Conn2.push_back(tmpCR);

            std::vector<int> tmpQL = std::get<0>(C);
            double tmpCoeffL = std::get<1>(C);
            tmpCoeffL *= sqrt(tmpQL[Qs[i]]);
            tmpQL[Qs[i]] -= 1;
            if (tmpQL[Qs[i]] < 0) continue;
            std::tuple<std::vector<int>, double> tmpCL = std::make_tuple(tmpQL, tmpCoeffL);
            Conn2.push_back(tmpCL);
        }
        Conn.clear();
        for (std::tuple<std::vector<int>, double> C2 : Conn2)
        {
            Conn.push_back(C2);
        }
        Conn2.clear();
    }
    return Conn;
}

std::vector<std::vector<std::tuple<std::vector<int>, double>>> FormBasisConnectionsCPP(std::vector<std::vector<VDeriv>>& Ws, std::vector<std::vector<int>>& Basis)
{
    std::vector<std::vector<std::tuple<std::vector<int>, double>>> BasisConnections;
    for (std::vector<int> B : Basis)
    {
        std::vector<std::tuple<std::vector<int>, double>> ConnB;
        for (int p = 0; p < Ws.size(); p++)
        {
            for (VDeriv W : Ws[p])
            {
                std::vector<std::tuple<std::vector<int>, double>> ConnByW = BasisConnectionsCPP(B, W.QIndices);
                for (std::tuple<std::vector<int>, double>& C : ConnByW)
                {
                    std::get<1>(C) = std::get<1>(C) * W.W;
                }
                ConnB.insert(ConnB.end(), ConnByW.begin(), ConnByW.end());
            }
        }
        BasisConnections.push_back(ConnB);
    }
    return BasisConnections; 
}

std::string VectorToString(std::vector<int> V)
{
    std::string OutString = "";
    for (int i = 0; i < V.size(); i++)
    {
        OutString += std::to_string(V[i]) + " ";
    }
    return OutString;
}

Eigen::MatrixXd HamVCPP(std::vector<std::vector<int>> Basis, std::vector<std::vector<std::tuple<std::vector<int>, double>>> BasisConn, std::vector<std::vector<int>> BasisBras, std::vector<double> Freq, std::vector<std::vector<VDeriv>> Ws, bool DiagonalOnly = false, bool OffDiagonal = false)
{
    int N = Basis.size();
    int NL = BasisBras.size();
    // Make a map of the basis vectors and their positions
    std::unordered_map<std::string, int> BasisBrasDict;
    for (int i = 0; i < BasisBras.size(); i++)
    {
        std::string BBString = VectorToString(BasisBras[i]);
        BasisBrasDict.insert({BBString, i});
    }
    
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(NL, N);
    if (!OffDiagonal)
    {
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < Basis[i].size(); j++)
            {
                H(i, i) += (Basis[i][j] + 0.5) * Freq[j];
            }
        }
    }
    
    for (int i = 0; i < Basis.size(); i++)
    {
        for (int ic = 0; ic < BasisConn[i].size(); ic++)
        {
            std::string CStr = VectorToString(std::get<0>(BasisConn[i][ic]));
            if (BasisBrasDict.find(CStr) != BasisBrasDict.end())
            {
                double Val = std::get<1>(BasisConn[i][ic]);
                int j = BasisBrasDict[CStr];
                if (DiagonalOnly && j != i) continue;
                if (!OffDiagonal && j > i) continue;
                H(j, i) += Val;
                if (!OffDiagonal && i != j)
                {
                    H(i, j) += Val;
                }
            }
        }
    }
    return H;
}

typedef Eigen::Triplet<double> T;

Eigen::SparseMatrix<double> SpHamVCPP(std::vector<std::vector<int>>& Basis, std::vector<std::vector<std::tuple<std::vector<int>, double>>>& BasisConn, std::vector<std::vector<int>>& BasisBras, std::vector<double>& Freq, std::vector<std::vector<VDeriv>>& Ws, bool DiagonalOnly = false, bool OffDiagonal = false)
{
    int N = Basis.size();
    int NL = BasisBras.size();
    // Make a map of the basis vectors and their positions
    std::unordered_map<std::string, int> BasisBrasDict;
    for (int i = 0; i < BasisBras.size(); i++)
    {
        std::string BBString = VectorToString(BasisBras[i]);
        BasisBrasDict.insert({BBString, i});
    }

    std::map<std::tuple<int, int>, double> HijDict;
    
    Eigen::SparseMatrix<double> H(NL, N);
    if (!OffDiagonal)
    {
        for (int i = 0; i < N; i++)
        {
            double Val = 0.0;
            for (int j = 0; j < Basis[i].size(); j++)
            {
                Val += (Basis[i][j] + 0.5) * Freq[j];
            }
            HijDict.insert({std::make_tuple(i, i), Val});
        }
    }
    
    for (int i = 0; i < Basis.size(); i++)
    {
        for (int ic = 0; ic < BasisConn[i].size(); ic++)
        {
            std::string CStr = VectorToString(std::get<0>(BasisConn[i][ic]));
            if (BasisBrasDict.find(CStr) != BasisBrasDict.end())
            {
                double Val = std::get<1>(BasisConn[i][ic]);
                int j = BasisBrasDict[CStr];
                if (DiagonalOnly && j != i) continue;
                if (!OffDiagonal && j > i) continue;
                
                std::tuple<int, int> jiTuple = std::make_tuple(j, i);
                std::tuple<int, int> ijTuple = std::make_tuple(i, j);
                if (HijDict.count(jiTuple))
                {
                    HijDict[jiTuple] += Val;
                    if (!OffDiagonal && i != j) HijDict[ijTuple] += Val;
                }
                else
                {
                    HijDict.insert({jiTuple, Val});
                    if (!OffDiagonal && i != j) HijDict.insert({ijTuple, Val});
                }
            }
        }
    }
    
    std::vector<T> TripletList;
    for (const auto& it : HijDict)
    {
        TripletList.push_back(T(std::get<0>(it.first), std::get<1>(it.first), it.second));
    }
    H.setFromTriplets(TripletList.begin(), TripletList.end());
    return H;
}

void VHCI::InitTruncatedBasis()
{
    std::vector<std::vector<int>> Bs;
    std::vector<int> B0(NModes, 0);
    Bs.push_back(B0);
    Basis.push_back(B0);
    for (unsigned int m = 0; m < MaxTotalQuanta; m++)
    {
        std::vector<std::vector<int>> BNext;
        for (std::vector<int> B : Bs)
        {
            for (unsigned int i = 0; i < B.size(); i++)
            {
                NewB = B;
                NewB[i]++;
                if (NewB[i] < MaxQuanta[i]) BNext.push_back(NewB);
            }
        }
        Basis.insert(Basis.end(); BNext.begin(); BNext.end());
        Bs = BNext;
        BNext.clear();
    }
}

void VHCI::FormWSD()
{
    std::vector<VDeriv> VSingles;
    std::vector<VDeriv> VDoubles;

    for (std::vector<VDeriv> Wp : Potential)
    {
        if (Wp[0].Order == 3)
        {
            for (unsigned int i = 0; i < NModes; i++)
            {
                double Wi = 0.0;
                for (VDeriv W : Wp)
                {
                    // Two cases, Wiii and Wijj
                    if (std::count(W.QIndices.begin(), W.QIndices.end(), i) == 1) Wi += 2.0 * W.W;
                    else if (std::count(W.QIndices.begin(), W.QIndices.end(), i) == 3) Wi += 3.0 * W.W;
                }
                if (abs(Wi) > 1e-12)
                {
                    std::vector<int> tmpQ;
                    tmpQ.push_back(i);
                    mW = VDeriv(Wi, tmpQ, false);
                    VSingles.push_back(mW);
                }
            }
        }
        if (Wp[0].Order == 4)
        {
            for (unsigned int i = 0; i < NModes; i++)
            {
                for (unsigned int j = i; j < NModes, j++)
                {
                    double Wij = 0.0;
                    for (VDeriv W : Wp)
                    {
                        // Four cases, Wiiii, Wiikk, Wijjj, Wijkk
                        if (std::count(W.QIndices.begin(), W.QIndices.end(), i) == 1 && std::count(W.QIndices.begin(), W.QIndices.end(), j) == 1 && i != j) Wij += 2.0 * W.W;
                        else if (std::count(W.QIndices.begin(), W.QIndices.end(), i) == 2 && W.QUnique.size() == 2 && i == j) Wij += 2.0 * W.W;
                        else if ((std::count(W.QIndices.begin(), W.QIndices.end(), i) == 1 && std::count(W.QIndices.begin(), W.QIndices.end(), j) == 3) || (std::count(W.QIndices.begin(), W.QIndices.end(), i) == 3 && std::count(W.QIndices.begin(), W.QIndices.end(), j) == 1)) Wij += 3.0 * W.W;
                        else if (std::count(W.QIndices.begin(), W.QIndices.end(), i) == 4 && i == j) Wij += 4.0 * W.W;
                    }
                    if (abs(Wij) > 1e-12)
                    {
                        std::vector<int> tmpQ;
                        tmpQ.push_back(i);
                        tmpQ.push_back(j);
                        mW = VDeriv(Wi, tmpQ, false);
                        VDoubles.push_back(mW);
                    }
                }
            }
        }
    }
    PotentialSD.clear()
    PotentialSD.push_back(VSingles);
    PotentialSD.push_back(VDoubles);
}
