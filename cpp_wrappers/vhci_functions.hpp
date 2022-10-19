#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <string>
#include <tuple>
#include <algorithm>
#include <cmath>
#include <unordered_map>

//Custom data structures
class HOFunc
{
    public:
        //Data structure for harmonic oscillator basis functions
        double Freq; //Frequency
        int Quanta; //Number of quanta in the mode
        // condition for equality
        bool operator == (const HOFunc &other) const {
            return (Freq == other.Freq) && (Quanta == other.Quanta);
        }
        // condition for inequality
        bool operator != (const HOFunc &other) const {
            return (Freq != other.Freq) || (Quanta != other.Quanta);
        }
};

class WaveFunction
{
    public:
        //Data structure for storing a VCI wavefunction
        int M; //Number of modes
        vector<HOFunc> Modes; // WaveFunctions
        bool operator == (const WaveFunction &other) const {
            bool same_wfn = 1;
            if(other.M != M){
                return 0;
            }
            for(int n=0; n<other.M; n++){
                if(other.Modes[n] != Modes[n]){
                    return 0;
                }
            }
            return same_wfn;
        }
        // condition for inequality
        bool operator != (const WaveFunction &other) const {
            bool diff_wfn = 0;
            if(other.M != M){
                return 1;
            }
            for(int n=0; n<other.M; n++){
                if(other.Modes[n] != Modes[n]){
                    return 1;
                }
            }
            return diff_wfn;
        }
};

struct FConst
{
    //Data structure for anharmonic force constants
    //Note: fc should include the permutation term
    double fc; //Value of the force constant
    vector<int> fcpow; //Modes and powers for the force constant
    vector<int> ShortModes; // Shortened list of modes only including those affected
    vector<int> ModePowers; // Power for each affected mode
};

struct WfnHasher
{
    size_t operator () (const WaveFunction& key) const 
    {
        // function to generate unique hash for a WaveFunction using Boost
        size_t seed = 0;
        for(int n=0; n<key.M; n++){
            hash_combine(seed, hash_value(key.Modes[n].Quanta));
            //hash_combine(seed, key.Modes[n].Freq);
        }
        return seed;
    }
};


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
