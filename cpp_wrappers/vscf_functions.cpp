#include "vscf_functions.hpp"

std::vector<Eigen::MatrixXd> GetVEffCPP(std::vector<Eigen::SparseMatrix<double>> &AnharmTensor, std::vector<WaveFunction> Basis, std::vector<Eigen::VectorXd> &CByModes, std::vector<std::vector<std::vector<int>>> &ModalSlices, std::vector<int> &MaxQuanta, std::vector<int> &ModeOcc)
{
    int NModes = CByModes.size();
    std::vector<Eigen::MatrixXd> VEff;

    for (unsigned int Mode = 0; Mode < NModes; Mode++)
    {
        Eigen::MatrixXd V = Eigen::Matrix::Zero(MaxQuanta[Mode], MaxQuanta[Mode]);
        for (unsigned int n = 0; n < MaxQuanta[Mode]; n++)
        {
            for (unsigned int m = n; m < MaxQuanta[Mode]; m++)
            {
                double Vnm = 0.0;
                for (auto HAnharm : AnharmTensor)
                {
                    for (unsigned int k = 0; k < HAnharm.outerSize(); k++)
                    {
                        for (Eigen::SparseMatrix<int, Eigen::ColMajor>::InnerIterator it(HAnharm, k); it; ++it)
                        {
                            int i = it.row();
                            int j = it.col();
                            if (std::count(ModalSlices[Mode][n].begin(), ModalSlices[Mode][n].end(), i)) // Means this element has n quanta in the mode we want.
                            {
                                std::vector<int> BraModeBasis; // List of quanta in all modes, which is the index of each mode's basis
                                for (auto HO : Basis[i].Modes)
                                {
                                    BraModeBasis.push_back(HO.Quanta);
                                }
                                if (std::count(ModalSlices[Mode][m].begin(), ModalSlices[Mode][m].end(), j))
                                {
                                    std::vector<int> KetModeBasis;
                                    for (auto HO: Basis[j].Modes)
                                    {
                                        KetModeBasis.push_back(HO.Quanta);
                                    }
                                    
                                    double CVC = HAnharm(i, j);
                                    for (unsigned int a = 0; a < BraModeBasis.size(); a++)
                                    {
                                        if (a != Mode)
                                        {
                                            CVC *= CByModes[a](BraModeBasis[a], ModeOcc[Mode]);
                                        }
                                    }
                                    for (unsigned int a = 0; a < KetModeBasis.size(); a++)
                                    {
                                        if (a != Mode)
                                        {
                                            CVC *= CByModes[a](KetModeBasis[a], ModeOcc[Mode]);
                                        }
                                    }
                                    Vnm += CVC;
                                }
                            }
                        }
                    }
                }
                V(n, m) = Vnm;
                V(m, n) = Vnm;
            }
        }
        VEff.push_back(V);
    }
    return VEff;
}
