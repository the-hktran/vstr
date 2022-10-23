#include "VCI_headers.h"

WaveFunction::WaveFunction(std::vector<int> Quantas, std::vector<double> Frequencies)
{
    for (unsigned int i = 0; i < Quantas.size(); i++)
    {
        HOFunc myHO;
        myHO.Freq = Frequencies[i];
        myHO.Quanta = Quantas[i];
        Modes.push_back(myHO);
    }
}

FConst::FConst(double dV, std::vector<int> Qs, bool doScale)
{
    fc = dV;
    QIndices = Qs;
    Order = QIndices.size();
    QUnique = Qs;
    std::sort(QUnique.begin(), QUnique.end());
    auto Last = std::unique(QUnique.begin(), QUnique.end());
    QUnique.erase(Last, QUnique.end());
    for (int i : QUnique) QPowers.push_back(std::count(QIndices.begin(), QIndices.end(), i));
    for (unsigned int i = 0; i < Order; i++) fcpow.push_back(std::count(QIndices.begin(), QIndices.end(), i));
    if (doScale) ScaleW();
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

void FConst::ScaleW()
{
    fc = fc / sqrt(pow(2.0, Order));
    for (int i = 0; i < QPowers.size(); i++)
    {
        fc /= Factorial(QPowers[i]);
    }
}


inline void CreationLO(double& ci, int& ni)
{
    //Creation ladder operator
    ci *= sqrt(ni+1); //Update coefficient
    ni += 1; //Update state
    return;
};

inline void AnnihilationLO(double& ci, int& ni)
{
    //Annihilation ladder operator
    ci *= sqrt(ni); //Update coefficient
    ni -= 1; //Update state
    //Check for impossible states
    if (ni < 0)
    {
        ni = -1; //Used later to remove the state
        ci = 0; //Delete state
    }
    return;
};

inline void QDiffVec(WaveFunction &Bn, WaveFunction &Bm, int &qtot, int &mchange, vector<int> &DiffVec)
{
    for (unsigned int i = 0; i < Bn.M; i++)
    {
        int qdiff = 0;
        qdiff += Bn.Modes[i].Quanta;
        qdiff -= Bm.Modes[i].Quanta;
        DiffVec[i] = abs(qdiff);
        qtot += abs(qdiff);
        if (qdiff != 0) mchange += 1;
    }
    return;
}

inline bool ScreenState(int qdiff, int mchange, const vector<int>& QDiffVec, const FConst& fc)
{
    //Function for ignoring states that have no overlap
    //qdiff is the number of changed quanta, mchange is number modes with changed quanta, QDiffVec is a vector with number of changed quanta per mode, fc is force constant
    bool keepstate = 1; //Assume the state is good
    if(qdiff > fc.QPowers.size() || 
            mchange > fc.QPowers.size()
            || qdiff%2 != fc.QPowers.size()%2
            ){
        return 0;
    }
    //Check based on force constant powers (check that raising and lowering results in quanta match)
    for (unsigned int i=0;i<fc.QPowers.size();i++)
    {
        if ( QDiffVec[fc.QPowers[i]] > fc.QPowers[i] || 
                QDiffVec[fc.QPowers[i]] % 2 != fc.QPowers[i] % 2){
            //Impossible for the states to overlap if mode power is too small or wrong even/odd parity
            //Skip the rest of the checks
            return 0;
        }
    }
    //Check overlap of all other modes (modes not involved in FC)
    for (int i=0;i<QDiffVec.size();i++)
    {
        bool cont = 1; //Continue the check
        for (unsigned int j=0;j<fc.QPowers.size();j++)
        {
            if (fc.QPowers[j] == i)
            {
                //Ignore this mode since it is in the FC
                cont = 0;
                break; // No need to continue checking mode against fc.QPowers
            }
        }
        if (cont)
        {
            if ( QDiffVec[i] != 0)
            {
                //Remove state due to zero overlap in mode i
                return 0; // No need to check if the other modes match 
            }
        }
    }
    //Return decision
    return keepstate;
};


