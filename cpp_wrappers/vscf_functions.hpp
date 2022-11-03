#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <string>
#include <tuple>
#include <algorithm>
#include <cmath>
#include <unordered_map>
#include "vhci_jf/VCI_headers.h"

using namespace std;

std::vector<Eigen::MatrixXd> GetVEffCPP(std::vector<Eigen::SparseMatrix<double>> &AnharmTensor, std::vector<WaveFunction> Basis, std::vector<Eigen::VectorXd> &CByModes, std::vector<std::vector<std::vector<int>>> &ModalSlices, std::vector<int> &MaxQuanta, std::vector<int> &ModeOcc);
