#include "VCI_headers.h"

// General Functions

WaveFunction::WaveFunction(std::vector<int> Quantas, std::vector<double> Frequencies)
{
    M = Quantas.size(); 
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
    if(qdiff > fc.fcpow.size() || 
            mchange > fc.QUnique.size()
            || qdiff%2 != fc.fcpow.size()%2
            ){
        return 0;
    }
    //Check based on force constant powers (check that raising and lowering results in quanta match)
    for (unsigned int i=0;i<fc.QUnique.size();i++)
    {
        if ( QDiffVec[fc.QUnique[i]] > fc.QPowers[i] || 
                QDiffVec[fc.QUnique[i]] % 2 != fc.QPowers[i] % 2){
            //Impossible for the states to overlap if mode power is too small or wrong even/odd parity
            //Skip the rest of the checks
            return 0;
        }
    }
    //Check overlap of all other modes (modes not involved in FC)
    for (int i=0;i<QDiffVec.size();i++)
    {
        bool cont = 1; //Continue the check
        for (unsigned int j=0;j<fc.QUnique.size();j++)
        {
            if (fc.QUnique[j] == i)
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

// HB.cpp

bool sortByFC(const FConst &lhs, const FConst &rhs) { return abs(lhs.fc) > abs(rhs.fc); } 
// Sort force constants from large to small magnitude

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


std::vector<FConst> HeatBath_Sort_FC(std::vector<FConst> &AnharmHB){
    sort(AnharmHB.begin(), AnharmHB.end(), sortByFC); // Sort force constants from large to small magnitude
    return AnharmHB;
    
    /*
    std::vector<double> FCs;
    for (FConst FC : AnharmHB) FCs.push_back(abs(FC.fc));
    std::vector<long unsigned int> FCInd = SortIndices(FCs);
    std::vector<FConst> SortedFC;

    for (int i = FCInd.size() - 1; i >= 0; i--) {SortedFC.push_back(AnharmHB[i]); std::cout << AnharmHB[i].fc << std::endl;}

    return SortedFC;
    */
}

std::vector<WaveFunction> AddStatesHB(std::vector<WaveFunction> &BasisSet, std::vector<FConst> &AnharmHB, Eigen::Ref<Eigen::VectorXd> C, double eps){ // Expand basis via Heat Bath algorithm
    HashedStates HashedBasisInit; // hashed unordered_set containing BasisSet to check for duplicates
    HashedStates HashedNewStates; // hashed unordered_set of new states that only allows unique states to be inserted
    for( WaveFunction& wfn : BasisSet){
        HashedBasisInit.insert(wfn); // Populate hashed unordered_set with initial basis states
    }

    std::vector<double> CVec;
    for (unsigned int n = 0; n < C.rows(); n++) CVec.push_back(abs(C[n]));
    std::vector<long unsigned int> CSortedInd = SortIndices(CVec);
    for(unsigned int i=0; i<AnharmHB.size(); ++i){ // Loop over sorted force constants
        if (abs(AnharmHB[i].fc * C[CSortedInd[CSortedInd.size() - 1]]) < eps) break; // means that the largest Cn doesn't meet the criteria so we are done
        for (int nn = CSortedInd.size() - 1; nn >= 0; nn--)
        {
            unsigned int n = CSortedInd[nn];
            double Cn = C[n];
            if(abs(Cn*AnharmHB[i].fc) >= eps){ // States connected by fc will be added if |fc*Cn| >= eps
                if(AnharmHB[i].QUnique.size()==1){ // Stupid way to enumerate new states from F_ij
                    for( int a=-(AnharmHB[i].QPowers[0]); a<(AnharmHB[i].QPowers[0]+1); a+=2 ){
                        if( a != 0){// Skip if no change
                            WaveFunction tmp = BasisSet[n];
                            tmp.Modes[AnharmHB[i].QUnique[0]].Quanta += a;
                            if( tmp.Modes[AnharmHB[i].QUnique[0]].Quanta >=0 &&
                                                                HashedBasisInit.count(tmp) == 0
                                    ){ //make sure a|0> = 0 and tmp does not exist in original basis
                                HashedNewStates.insert(tmp); // add new state to set
                            }
                        }
                    }
                }
                if(AnharmHB[i].QUnique.size()==2){ // Stupid way to enumerate new states from F_ij
                    for( int a=-(AnharmHB[i].QPowers[0]); a<(AnharmHB[i].QPowers[0]+1); a+=2 ){
                        for( int b=-(AnharmHB[i].QPowers[1]); b<(AnharmHB[i].QPowers[1]+1); b+=2 ){
                            if(abs(a)+abs(b) != 0){// Skip if no change
                                WaveFunction tmp = BasisSet[n];
                                tmp.Modes[AnharmHB[i].QUnique[0]].Quanta += a;
                                tmp.Modes[AnharmHB[i].QUnique[1]].Quanta += b;
                                if( tmp.Modes[AnharmHB[i].QUnique[0]].Quanta >=0 &&
                                       tmp.Modes[AnharmHB[i].QUnique[1]].Quanta >=0 && 
                                                                       HashedBasisInit.count(tmp) == 0
                                       ){ //make sure a|0> = 0
                                    HashedNewStates.insert(tmp); // add new state to set
                                }
                            }
                        }
                    }
                }
                if(AnharmHB[i].QUnique.size()==3){ // Stupid way to enumerate new states from F_ijk
                    for( int a=-(AnharmHB[i].QPowers[0]); a<(AnharmHB[i].QPowers[0]+1); a+=2 ){
                        for( int b=-(AnharmHB[i].QPowers[1]); b<(AnharmHB[i].QPowers[1]+1); b+=2 ){
                            for( int c=-(AnharmHB[i].QPowers[2]); c<(AnharmHB[i].QPowers[2]+1); c+=2 ){
                                if(abs(a)+abs(b)+abs(c) != 0){// Skip if no change
                                    WaveFunction tmp = BasisSet[n];
                                    tmp.Modes[AnharmHB[i].QUnique[0]].Quanta += a;
                                    tmp.Modes[AnharmHB[i].QUnique[1]].Quanta += b;
                                    tmp.Modes[AnharmHB[i].QUnique[2]].Quanta += c;
                                    if( tmp.Modes[AnharmHB[i].QUnique[0]].Quanta >=0 &&
                                           tmp.Modes[AnharmHB[i].QUnique[1]].Quanta >=0  &&
                                           tmp.Modes[AnharmHB[i].QUnique[2]].Quanta >=0 &&
                                           HashedBasisInit.count(tmp) == 0
                                           ){ //make sure a|0> = 0
                                        HashedNewStates.insert(tmp); // add new state
                                    }
                                }      
                            }
                        }
                    }
                }
                if(AnharmHB[i].QUnique.size()==4){ // Stupid way to enumerate new states from F_ijkl
                    for( int a=-(AnharmHB[i].QPowers[0]); a<(AnharmHB[i].QPowers[0]+1); a+=2 ){
                        for( int b=-(AnharmHB[i].QPowers[1]); b<(AnharmHB[i].QPowers[1]+1); b+=2 ){
                            for( int c=-(AnharmHB[i].QPowers[2]); c<(AnharmHB[i].QPowers[2]+1); c+=2 ){
                                for( int d=-(AnharmHB[i].QPowers[3]); d<(AnharmHB[i].QPowers[3]+1); d+=2 ){
                                    if(abs(a)+abs(b)+abs(c)+abs(d) != 0){ // Skip if no change
                                        WaveFunction tmp = BasisSet[n];
                                        tmp.Modes[AnharmHB[i].QUnique[0]].Quanta += a;
                                        tmp.Modes[AnharmHB[i].QUnique[1]].Quanta += b;
                                        tmp.Modes[AnharmHB[i].QUnique[2]].Quanta += c;
                                        tmp.Modes[AnharmHB[i].QUnique[3]].Quanta += d;
                                        if( tmp.Modes[AnharmHB[i].QUnique[0]].Quanta >=0 &&
                                               tmp.Modes[AnharmHB[i].QUnique[1]].Quanta >=0  &&
                                               tmp.Modes[AnharmHB[i].QUnique[2]].Quanta >=0 &&
                                               tmp.Modes[AnharmHB[i].QUnique[3]].Quanta >=0 &&
                                               HashedBasisInit.count(tmp) == 0
                                               ){ //make sure a|0> = 0
                                            HashedNewStates.insert(tmp); // add new state
                                        }
                                    }
                                }
                            }
                        }
                    }
                }   
                if(AnharmHB[i].QUnique.size()==5){ // Stupid way to enumerate new states from F_ijklm 
                    for( int a=-(AnharmHB[i].QPowers[0]); a<(AnharmHB[i].QPowers[0]+1); a+=2 ){
                        for( int b=-(AnharmHB[i].QPowers[1]); b<(AnharmHB[i].QPowers[1]+1); b+=2 ){
                            for( int c=-(AnharmHB[i].QPowers[2]); c<(AnharmHB[i].QPowers[2]+1); c+=2 ){
                                for( int d=-(AnharmHB[i].QPowers[3]); d<(AnharmHB[i].QPowers[3]+1); d+=2 ){
                                    for( int f=-(AnharmHB[i].QPowers[4]); f<(AnharmHB[i].QPowers[4]+1); f+=2 ){
                                        if(abs(a)+abs(b)+abs(c)+abs(d)+abs(f) != 0){ // Skip if no change
                                            WaveFunction tmp = BasisSet[n];
                                            tmp.Modes[AnharmHB[i].QUnique[0]].Quanta += a;
                                            tmp.Modes[AnharmHB[i].QUnique[1]].Quanta += b;
                                            tmp.Modes[AnharmHB[i].QUnique[2]].Quanta += c;
                                            tmp.Modes[AnharmHB[i].QUnique[3]].Quanta += d;
                                            tmp.Modes[AnharmHB[i].QUnique[4]].Quanta += f;
                                            if( tmp.Modes[AnharmHB[i].QUnique[0]].Quanta >=0 &&
                                                   tmp.Modes[AnharmHB[i].QUnique[1]].Quanta >=0  &&
                                                   tmp.Modes[AnharmHB[i].QUnique[2]].Quanta >=0 &&
                                                   tmp.Modes[AnharmHB[i].QUnique[3]].Quanta >=0 &&
                                                   tmp.Modes[AnharmHB[i].QUnique[4]].Quanta >=0 && 
                                                   HashedBasisInit.count(tmp) == 0
                                                   ){ //make sure a|0> = 0
                                                HashedNewStates.insert(tmp); // add new state
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                if(AnharmHB[i].QUnique.size()==6){ // Stupid way to enumerate new states from F_ijklmn 
                    for( int a=-(AnharmHB[i].QPowers[0]); a<(AnharmHB[i].QPowers[0]+1); a+=2 ){
                        for( int b=-(AnharmHB[i].QPowers[1]); b<(AnharmHB[i].QPowers[1]+1); b+=2 ){
                            for( int c=-(AnharmHB[i].QPowers[2]); c<(AnharmHB[i].QPowers[2]+1); c+=2 ){
                                for( int d=-(AnharmHB[i].QPowers[3]); d<(AnharmHB[i].QPowers[3]+1); d+=2 ){
                                    for( int f=-(AnharmHB[i].QPowers[4]); f<(AnharmHB[i].QPowers[4]+1); f+=2 ){
                                        for( int g=-(AnharmHB[i].QPowers[5]); g<(AnharmHB[i].QPowers[5]+1); g+=2 ){
                                            if(abs(a)+abs(b)+abs(c)+abs(d)+abs(f)+abs(g) != 0){ // Skip if no change
                                                WaveFunction tmp = BasisSet[n];
                                                tmp.Modes[AnharmHB[i].QUnique[0]].Quanta += a;
                                                tmp.Modes[AnharmHB[i].QUnique[1]].Quanta += b;
                                                tmp.Modes[AnharmHB[i].QUnique[2]].Quanta += c;
                                                tmp.Modes[AnharmHB[i].QUnique[3]].Quanta += d;
                                                tmp.Modes[AnharmHB[i].QUnique[4]].Quanta += f;
                                                tmp.Modes[AnharmHB[i].QUnique[5]].Quanta += g;
                                                if( tmp.Modes[AnharmHB[i].QUnique[0]].Quanta >=0 &&
                                                       tmp.Modes[AnharmHB[i].QUnique[1]].Quanta >=0 &&
                                                       tmp.Modes[AnharmHB[i].QUnique[2]].Quanta >=0 &&
                                                       tmp.Modes[AnharmHB[i].QUnique[3]].Quanta >=0 &&
                                                       tmp.Modes[AnharmHB[i].QUnique[4]].Quanta >=0 &&
                                                       tmp.Modes[AnharmHB[i].QUnique[5]].Quanta >=0 &&
                                                                                                       HashedBasisInit.count(tmp) == 0
                                                       ){ //make sure a|0> = 0
                                                    HashedNewStates.insert(tmp); // add new state
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                if(AnharmHB[i].fcpow.size() == 0 || AnharmHB[i].fcpow.size() > 6 ){
                    //Print error message
                    cout << "Error: HCI only works with force constants up to 6th order." << endl;
                    cout.flush();
                    exit(0);
                }
            }
            else{break;}// break loop if you've reached element < eps, since all future elements will be smaller (HB sorting)
        }
    }
    std::vector<WaveFunction> NewBasis;
    for (const WaveFunction &WF : HashedNewStates) NewBasis.push_back(WF);
    return NewBasis;
}

std::vector<WaveFunction> AddStatesHB(std::vector<WaveFunction> &BasisSet, std::vector<FConst> &AnharmHB, int n, double Cn, double eps){ // Expand basis via Heat Bath algorithm
    HashedStates HashedBasisInit; // hashed unordered_set containing BasisSet to check for duplicates
    HashedStates HashedNewStates; // hashed unordered_set of new states that only allows unique states to be inserted
    for( WaveFunction& wfn : BasisSet){
        HashedBasisInit.insert(wfn); // Populate hashed unordered_set with initial basis states
    }

    for(unsigned int i=0; i<AnharmHB.size(); ++i){ // Loop over sorted force constants
            if(abs(Cn*AnharmHB[i].fc) >= eps){ // States connected by fc will be added if |fc*Cn| >= eps
                if(AnharmHB[i].QUnique.size()==1){ // Stupid way to enumerate new states from F_ij
                    for( int a=-(AnharmHB[i].QPowers[0]); a<(AnharmHB[i].QPowers[0]+1); a+=2 ){
                        if( a != 0){// Skip if no change
                            WaveFunction tmp = BasisSet[n];
                            tmp.Modes[AnharmHB[i].QUnique[0]].Quanta += a;
                            if( tmp.Modes[AnharmHB[i].QUnique[0]].Quanta >=0 &&
                                                                HashedBasisInit.count(tmp) == 0
                                    ){ //make sure a|0> = 0 and tmp does not exist in original basis
                                HashedNewStates.insert(tmp); // add new state to set
                            }
                        }
                    }
                }
                if(AnharmHB[i].QUnique.size()==2){ // Stupid way to enumerate new states from F_ij
                    for( int a=-(AnharmHB[i].QPowers[0]); a<(AnharmHB[i].QPowers[0]+1); a+=2 ){
                        for( int b=-(AnharmHB[i].QPowers[1]); b<(AnharmHB[i].QPowers[1]+1); b+=2 ){
                            if(abs(a)+abs(b) != 0){// Skip if no change
                                WaveFunction tmp = BasisSet[n];
                                tmp.Modes[AnharmHB[i].QUnique[0]].Quanta += a;
                                tmp.Modes[AnharmHB[i].QUnique[1]].Quanta += b;
                                if( tmp.Modes[AnharmHB[i].QUnique[0]].Quanta >=0 &&
                                       tmp.Modes[AnharmHB[i].QUnique[1]].Quanta >=0 && 
                                                                       HashedBasisInit.count(tmp) == 0
                                       ){ //make sure a|0> = 0
                                    HashedNewStates.insert(tmp); // add new state to set
                                }
                            }
                        }
                    }
                }
                if(AnharmHB[i].QUnique.size()==3){ // Stupid way to enumerate new states from F_ijk
                    for( int a=-(AnharmHB[i].QPowers[0]); a<(AnharmHB[i].QPowers[0]+1); a+=2 ){
                        for( int b=-(AnharmHB[i].QPowers[1]); b<(AnharmHB[i].QPowers[1]+1); b+=2 ){
                            for( int c=-(AnharmHB[i].QPowers[2]); c<(AnharmHB[i].QPowers[2]+1); c+=2 ){
                                if(abs(a)+abs(b)+abs(c) != 0){// Skip if no change
                                    WaveFunction tmp = BasisSet[n];
                                    tmp.Modes[AnharmHB[i].QUnique[0]].Quanta += a;
                                    tmp.Modes[AnharmHB[i].QUnique[1]].Quanta += b;
                                    tmp.Modes[AnharmHB[i].QUnique[2]].Quanta += c;
                                    if( tmp.Modes[AnharmHB[i].QUnique[0]].Quanta >=0 &&
                                           tmp.Modes[AnharmHB[i].QUnique[1]].Quanta >=0  &&
                                           tmp.Modes[AnharmHB[i].QUnique[2]].Quanta >=0 &&
                                           HashedBasisInit.count(tmp) == 0
                                           ){ //make sure a|0> = 0
                                        HashedNewStates.insert(tmp); // add new state
                                    }
                                }      
                            }
                        }
                    }
                }
                if(AnharmHB[i].QUnique.size()==4){ // Stupid way to enumerate new states from F_ijkl
                    for( int a=-(AnharmHB[i].QPowers[0]); a<(AnharmHB[i].QPowers[0]+1); a+=2 ){
                        for( int b=-(AnharmHB[i].QPowers[1]); b<(AnharmHB[i].QPowers[1]+1); b+=2 ){
                            for( int c=-(AnharmHB[i].QPowers[2]); c<(AnharmHB[i].QPowers[2]+1); c+=2 ){
                                for( int d=-(AnharmHB[i].QPowers[3]); d<(AnharmHB[i].QPowers[3]+1); d+=2 ){
                                    if(abs(a)+abs(b)+abs(c)+abs(d) != 0){ // Skip if no change
                                        WaveFunction tmp = BasisSet[n];
                                        tmp.Modes[AnharmHB[i].QUnique[0]].Quanta += a;
                                        tmp.Modes[AnharmHB[i].QUnique[1]].Quanta += b;
                                        tmp.Modes[AnharmHB[i].QUnique[2]].Quanta += c;
                                        tmp.Modes[AnharmHB[i].QUnique[3]].Quanta += d;
                                        if( tmp.Modes[AnharmHB[i].QUnique[0]].Quanta >=0 &&
                                               tmp.Modes[AnharmHB[i].QUnique[1]].Quanta >=0  &&
                                               tmp.Modes[AnharmHB[i].QUnique[2]].Quanta >=0 &&
                                               tmp.Modes[AnharmHB[i].QUnique[3]].Quanta >=0 &&
                                               HashedBasisInit.count(tmp) == 0
                                               ){ //make sure a|0> = 0
                                            HashedNewStates.insert(tmp); // add new state
                                        }
                                    }
                                }
                            }
                        }
                    }
                }   
                if(AnharmHB[i].QUnique.size()==5){ // Stupid way to enumerate new states from F_ijklm 
                    for( int a=-(AnharmHB[i].QPowers[0]); a<(AnharmHB[i].QPowers[0]+1); a+=2 ){
                        for( int b=-(AnharmHB[i].QPowers[1]); b<(AnharmHB[i].QPowers[1]+1); b+=2 ){
                            for( int c=-(AnharmHB[i].QPowers[2]); c<(AnharmHB[i].QPowers[2]+1); c+=2 ){
                                for( int d=-(AnharmHB[i].QPowers[3]); d<(AnharmHB[i].QPowers[3]+1); d+=2 ){
                                    for( int f=-(AnharmHB[i].QPowers[4]); f<(AnharmHB[i].QPowers[4]+1); f+=2 ){
                                        if(abs(a)+abs(b)+abs(c)+abs(d)+abs(f) != 0){ // Skip if no change
                                            WaveFunction tmp = BasisSet[n];
                                            tmp.Modes[AnharmHB[i].QUnique[0]].Quanta += a;
                                            tmp.Modes[AnharmHB[i].QUnique[1]].Quanta += b;
                                            tmp.Modes[AnharmHB[i].QUnique[2]].Quanta += c;
                                            tmp.Modes[AnharmHB[i].QUnique[3]].Quanta += d;
                                            tmp.Modes[AnharmHB[i].QUnique[4]].Quanta += f;
                                            if( tmp.Modes[AnharmHB[i].QUnique[0]].Quanta >=0 &&
                                                   tmp.Modes[AnharmHB[i].QUnique[1]].Quanta >=0  &&
                                                   tmp.Modes[AnharmHB[i].QUnique[2]].Quanta >=0 &&
                                                   tmp.Modes[AnharmHB[i].QUnique[3]].Quanta >=0 &&
                                                   tmp.Modes[AnharmHB[i].QUnique[4]].Quanta >=0 && 
                                                   HashedBasisInit.count(tmp) == 0
                                                   ){ //make sure a|0> = 0
                                                HashedNewStates.insert(tmp); // add new state
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                if(AnharmHB[i].QUnique.size()==6){ // Stupid way to enumerate new states from F_ijklmn 
                    for( int a=-(AnharmHB[i].QPowers[0]); a<(AnharmHB[i].QPowers[0]+1); a+=2 ){
                        for( int b=-(AnharmHB[i].QPowers[1]); b<(AnharmHB[i].QPowers[1]+1); b+=2 ){
                            for( int c=-(AnharmHB[i].QPowers[2]); c<(AnharmHB[i].QPowers[2]+1); c+=2 ){
                                for( int d=-(AnharmHB[i].QPowers[3]); d<(AnharmHB[i].QPowers[3]+1); d+=2 ){
                                    for( int f=-(AnharmHB[i].QPowers[4]); f<(AnharmHB[i].QPowers[4]+1); f+=2 ){
                                        for( int g=-(AnharmHB[i].QPowers[5]); g<(AnharmHB[i].QPowers[5]+1); g+=2 ){
                                            if(abs(a)+abs(b)+abs(c)+abs(d)+abs(f)+abs(g) != 0){ // Skip if no change
                                                WaveFunction tmp = BasisSet[n];
                                                tmp.Modes[AnharmHB[i].QUnique[0]].Quanta += a;
                                                tmp.Modes[AnharmHB[i].QUnique[1]].Quanta += b;
                                                tmp.Modes[AnharmHB[i].QUnique[2]].Quanta += c;
                                                tmp.Modes[AnharmHB[i].QUnique[3]].Quanta += d;
                                                tmp.Modes[AnharmHB[i].QUnique[4]].Quanta += f;
                                                tmp.Modes[AnharmHB[i].QUnique[5]].Quanta += g;
                                                if( tmp.Modes[AnharmHB[i].QUnique[0]].Quanta >=0 &&
                                                       tmp.Modes[AnharmHB[i].QUnique[1]].Quanta >=0 &&
                                                       tmp.Modes[AnharmHB[i].QUnique[2]].Quanta >=0 &&
                                                       tmp.Modes[AnharmHB[i].QUnique[3]].Quanta >=0 &&
                                                       tmp.Modes[AnharmHB[i].QUnique[4]].Quanta >=0 &&
                                                       tmp.Modes[AnharmHB[i].QUnique[5]].Quanta >=0 &&
                                                                                                       HashedBasisInit.count(tmp) == 0
                                                       ){ //make sure a|0> = 0
                                                    HashedNewStates.insert(tmp); // add new state
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                if(AnharmHB[i].fcpow.size() == 0 || AnharmHB[i].fcpow.size() > 6 ){
                    //Print error message
                    cout << "Error: HCI only works with force constants up to 6th order." << endl;
                    cout.flush();
                    exit(0);
                }
            }
            else{break;}// break loop if you've reached element < eps, since all future elements will be smaller (HB sorting)
        }
    std::vector<WaveFunction> NewBasis;
    for (const WaveFunction &WF : HashedNewStates) NewBasis.push_back(WF);
    return NewBasis;
}

void InternalAddStatesHB(std::vector<WaveFunction> &BasisSet, HashedStates &HashedNewStates, std::vector<FConst> &AnharmHB, int n, double Cn, double eps){ // Expand basis via Heat Bath algorithm
    HashedStates HashedBasisInit; // hashed unordered_set containing BasisSet to check for duplicates
    for( WaveFunction& wfn : BasisSet){
        HashedBasisInit.insert(wfn); // Populate hashed unordered_set with initial basis states
    }

    for(unsigned int i=0; i<AnharmHB.size(); ++i){ // Loop over sorted force constants
            if(abs(Cn*AnharmHB[i].fc) >= eps){ // States connected by fc will be added if |fc*Cn| >= eps
                if(AnharmHB[i].QUnique.size()==1){ // Stupid way to enumerate new states from F_ij
                    for( int a=-(AnharmHB[i].QPowers[0]); a<(AnharmHB[i].QPowers[0]+1); a+=2 ){
                        if( a != 0){// Skip if no change
                            WaveFunction tmp = BasisSet[n];
                            tmp.Modes[AnharmHB[i].QUnique[0]].Quanta += a;
                            if( tmp.Modes[AnharmHB[i].QUnique[0]].Quanta >=0 &&
                                                                HashedBasisInit.count(tmp) == 0
                                    ){ //make sure a|0> = 0 and tmp does not exist in original basis
                                HashedNewStates.insert(tmp); // add new state to set
                            }
                        }
                    }
                }
                if(AnharmHB[i].QUnique.size()==2){ // Stupid way to enumerate new states from F_ij
                    for( int a=-(AnharmHB[i].QPowers[0]); a<(AnharmHB[i].QPowers[0]+1); a+=2 ){
                        for( int b=-(AnharmHB[i].QPowers[1]); b<(AnharmHB[i].QPowers[1]+1); b+=2 ){
                            if(abs(a)+abs(b) != 0){// Skip if no change
                                WaveFunction tmp = BasisSet[n];
                                tmp.Modes[AnharmHB[i].QUnique[0]].Quanta += a;
                                tmp.Modes[AnharmHB[i].QUnique[1]].Quanta += b;
                                if( tmp.Modes[AnharmHB[i].QUnique[0]].Quanta >=0 &&
                                       tmp.Modes[AnharmHB[i].QUnique[1]].Quanta >=0 && 
                                                                       HashedBasisInit.count(tmp) == 0
                                       ){ //make sure a|0> = 0
                                    HashedNewStates.insert(tmp); // add new state to set
                                }
                            }
                        }
                    }
                }
                if(AnharmHB[i].QUnique.size()==3){ // Stupid way to enumerate new states from F_ijk
                    for( int a=-(AnharmHB[i].QPowers[0]); a<(AnharmHB[i].QPowers[0]+1); a+=2 ){
                        for( int b=-(AnharmHB[i].QPowers[1]); b<(AnharmHB[i].QPowers[1]+1); b+=2 ){
                            for( int c=-(AnharmHB[i].QPowers[2]); c<(AnharmHB[i].QPowers[2]+1); c+=2 ){
                                if(abs(a)+abs(b)+abs(c) != 0){// Skip if no change
                                    WaveFunction tmp = BasisSet[n];
                                    tmp.Modes[AnharmHB[i].QUnique[0]].Quanta += a;
                                    tmp.Modes[AnharmHB[i].QUnique[1]].Quanta += b;
                                    tmp.Modes[AnharmHB[i].QUnique[2]].Quanta += c;
                                    if( tmp.Modes[AnharmHB[i].QUnique[0]].Quanta >=0 &&
                                           tmp.Modes[AnharmHB[i].QUnique[1]].Quanta >=0  &&
                                           tmp.Modes[AnharmHB[i].QUnique[2]].Quanta >=0 &&
                                           HashedBasisInit.count(tmp) == 0
                                           ){ //make sure a|0> = 0
                                        HashedNewStates.insert(tmp); // add new state
                                    }
                                }      
                            }
                        }
                    }
                }
                if(AnharmHB[i].QUnique.size()==4){ // Stupid way to enumerate new states from F_ijkl
                    for( int a=-(AnharmHB[i].QPowers[0]); a<(AnharmHB[i].QPowers[0]+1); a+=2 ){
                        for( int b=-(AnharmHB[i].QPowers[1]); b<(AnharmHB[i].QPowers[1]+1); b+=2 ){
                            for( int c=-(AnharmHB[i].QPowers[2]); c<(AnharmHB[i].QPowers[2]+1); c+=2 ){
                                for( int d=-(AnharmHB[i].QPowers[3]); d<(AnharmHB[i].QPowers[3]+1); d+=2 ){
                                    if(abs(a)+abs(b)+abs(c)+abs(d) != 0){ // Skip if no change
                                        WaveFunction tmp = BasisSet[n];
                                        tmp.Modes[AnharmHB[i].QUnique[0]].Quanta += a;
                                        tmp.Modes[AnharmHB[i].QUnique[1]].Quanta += b;
                                        tmp.Modes[AnharmHB[i].QUnique[2]].Quanta += c;
                                        tmp.Modes[AnharmHB[i].QUnique[3]].Quanta += d;
                                        if( tmp.Modes[AnharmHB[i].QUnique[0]].Quanta >=0 &&
                                               tmp.Modes[AnharmHB[i].QUnique[1]].Quanta >=0  &&
                                               tmp.Modes[AnharmHB[i].QUnique[2]].Quanta >=0 &&
                                               tmp.Modes[AnharmHB[i].QUnique[3]].Quanta >=0 &&
                                               HashedBasisInit.count(tmp) == 0
                                               ){ //make sure a|0> = 0
                                            HashedNewStates.insert(tmp); // add new state
                                        }
                                    }
                                }
                            }
                        }
                    }
                }   
                if(AnharmHB[i].QUnique.size()==5){ // Stupid way to enumerate new states from F_ijklm 
                    for( int a=-(AnharmHB[i].QPowers[0]); a<(AnharmHB[i].QPowers[0]+1); a+=2 ){
                        for( int b=-(AnharmHB[i].QPowers[1]); b<(AnharmHB[i].QPowers[1]+1); b+=2 ){
                            for( int c=-(AnharmHB[i].QPowers[2]); c<(AnharmHB[i].QPowers[2]+1); c+=2 ){
                                for( int d=-(AnharmHB[i].QPowers[3]); d<(AnharmHB[i].QPowers[3]+1); d+=2 ){
                                    for( int f=-(AnharmHB[i].QPowers[4]); f<(AnharmHB[i].QPowers[4]+1); f+=2 ){
                                        if(abs(a)+abs(b)+abs(c)+abs(d)+abs(f) != 0){ // Skip if no change
                                            WaveFunction tmp = BasisSet[n];
                                            tmp.Modes[AnharmHB[i].QUnique[0]].Quanta += a;
                                            tmp.Modes[AnharmHB[i].QUnique[1]].Quanta += b;
                                            tmp.Modes[AnharmHB[i].QUnique[2]].Quanta += c;
                                            tmp.Modes[AnharmHB[i].QUnique[3]].Quanta += d;
                                            tmp.Modes[AnharmHB[i].QUnique[4]].Quanta += f;
                                            if( tmp.Modes[AnharmHB[i].QUnique[0]].Quanta >=0 &&
                                                   tmp.Modes[AnharmHB[i].QUnique[1]].Quanta >=0  &&
                                                   tmp.Modes[AnharmHB[i].QUnique[2]].Quanta >=0 &&
                                                   tmp.Modes[AnharmHB[i].QUnique[3]].Quanta >=0 &&
                                                   tmp.Modes[AnharmHB[i].QUnique[4]].Quanta >=0 && 
                                                   HashedBasisInit.count(tmp) == 0
                                                   ){ //make sure a|0> = 0
                                                HashedNewStates.insert(tmp); // add new state
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                if(AnharmHB[i].QUnique.size()==6){ // Stupid way to enumerate new states from F_ijklmn 
                    for( int a=-(AnharmHB[i].QPowers[0]); a<(AnharmHB[i].QPowers[0]+1); a+=2 ){
                        for( int b=-(AnharmHB[i].QPowers[1]); b<(AnharmHB[i].QPowers[1]+1); b+=2 ){
                            for( int c=-(AnharmHB[i].QPowers[2]); c<(AnharmHB[i].QPowers[2]+1); c+=2 ){
                                for( int d=-(AnharmHB[i].QPowers[3]); d<(AnharmHB[i].QPowers[3]+1); d+=2 ){
                                    for( int f=-(AnharmHB[i].QPowers[4]); f<(AnharmHB[i].QPowers[4]+1); f+=2 ){
                                        for( int g=-(AnharmHB[i].QPowers[5]); g<(AnharmHB[i].QPowers[5]+1); g+=2 ){
                                            if(abs(a)+abs(b)+abs(c)+abs(d)+abs(f)+abs(g) != 0){ // Skip if no change
                                                WaveFunction tmp = BasisSet[n];
                                                tmp.Modes[AnharmHB[i].QUnique[0]].Quanta += a;
                                                tmp.Modes[AnharmHB[i].QUnique[1]].Quanta += b;
                                                tmp.Modes[AnharmHB[i].QUnique[2]].Quanta += c;
                                                tmp.Modes[AnharmHB[i].QUnique[3]].Quanta += d;
                                                tmp.Modes[AnharmHB[i].QUnique[4]].Quanta += f;
                                                tmp.Modes[AnharmHB[i].QUnique[5]].Quanta += g;
                                                if( tmp.Modes[AnharmHB[i].QUnique[0]].Quanta >=0 &&
                                                       tmp.Modes[AnharmHB[i].QUnique[1]].Quanta >=0 &&
                                                       tmp.Modes[AnharmHB[i].QUnique[2]].Quanta >=0 &&
                                                       tmp.Modes[AnharmHB[i].QUnique[3]].Quanta >=0 &&
                                                       tmp.Modes[AnharmHB[i].QUnique[4]].Quanta >=0 &&
                                                       tmp.Modes[AnharmHB[i].QUnique[5]].Quanta >=0 &&
                                                                                                       HashedBasisInit.count(tmp) == 0
                                                       ){ //make sure a|0> = 0
                                                    HashedNewStates.insert(tmp); // add new state
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                if(AnharmHB[i].fcpow.size() == 0 || AnharmHB[i].fcpow.size() > 6 ){
                    //Print error message
                    cout << "Error: HCI only works with force constants up to 6th order." << endl;
                    cout.flush();
                    exit(0);
                }
            }
            else{break;}// break loop if you've reached element < eps, since all future elements will be smaller (HB sorting)
        }
}

// Ham.cpp

double AnharmPot(WaveFunction &Bn, WaveFunction &Bm, const FConst& fc)
{
    //Calculate anharmonic matrix elements for <m|H|n>
    double Vnm = 0;
    // Initialize new states
    vector<vector<int> > NewStates(fc.QUnique.size());
    vector<vector<double> > StateCoeffs(fc.QUnique.size());
    //Create new states
    for (unsigned int i=0;i<fc.QPowers.size();i++)
    {
        //Apply operators
        NewStates[i].push_back(Bn.Modes[fc.QUnique[i]].Quanta);
        StateCoeffs[i].push_back(1.0);
        for (int j=0;j<fc.QPowers[i];j++)
        {
            //Create new state for mode i
            vector<int> stateupdate;
            vector<double> coeffupdate;
            for (unsigned int k=0;k<NewStates[i].size();k++)
            {
                int quant;
                double coeff;
                //Creation
                quant = NewStates[i][k];
                coeff = StateCoeffs[i][k];
                CreationLO(coeff,quant);
                stateupdate.push_back(quant);
                coeffupdate.push_back(coeff);
                //Annihilation
                quant = NewStates[i][k];
                coeff = StateCoeffs[i][k];
                AnnihilationLO(coeff,quant);
                if (quant >= 0)
                {
                    stateupdate.push_back(quant);
                    coeffupdate.push_back(coeff);
                }
            }
            //Save states
            NewStates[i] = stateupdate;
            StateCoeffs[i] = coeffupdate; // Accounting for permutations
        }
    }
    //Sum energies
    vector<double> CoeffSum;
    for (unsigned int i=0;i<fc.QUnique.size();i++)
    {
        CoeffSum.push_back(0.0);
        for (unsigned int j=0;j<NewStates[i].size();j++)
        {
            int quantn = NewStates[i][j];
            int quantm = Bm.Modes[fc.QUnique[i]].Quanta;
            if (quantn == quantm)
            {
                CoeffSum[i] += StateCoeffs[i][j];
            }
        }
    }
    //Scale by the force constant
    Vnm = fc.fc;
    //Combine coeffcients
    for (unsigned int i=0;i<CoeffSum.size();i++)
    {
        Vnm *= CoeffSum[i];
    }
    return Vnm;
};


//Hamiltonian operators
void ZerothHam(Eigen::MatrixXd &H, std::vector<WaveFunction> &BasisSet)
{
    //Calculate the harmonic Hamiltonian matrix elements
    //Harmonic matrix elements
    #pragma omp parallel for
    for (unsigned int i=0;i<BasisSet.size();i++)
    {
        //Loop over all modes
        double Ei = 0.; //Hii matrix element
        for (int j=0;j<BasisSet[i].M;j++)
        {
            //Calculate partial energies
            double Ej = 0.5;
            Ej += BasisSet[i].Modes[j].Quanta;
            Ej *= BasisSet[i].Modes[j].Freq;
            //Update matrix element
            Ei += Ej;
        }
        //Update Hamiltonian
        H(i,i) += Ei;
    }
    return;
};

void AnharmHam(Eigen::MatrixXd &H, std::vector<WaveFunction> &BasisSet, std::vector<double> &Frequencies, std::vector<FConst> &AnharmFC, std::vector<FConst> &CubicFC, std::vector<FConst> &QuarticFC, std::vector<FConst> &QuinticFC, std::vector<FConst> &SexticFC)
{
    //Add anharmonic terms to the Hamiltonian
    int fcmax=0;
    for (unsigned int k=0;k<AnharmFC.size();k++){ //Order of maximum anharmonic term
        if(AnharmFC[k].fcpow.size()>fcmax){
            fcmax = AnharmFC[k].fcpow.size();
        }
    }
    #pragma omp parallel for
    for (unsigned int i=0;i<BasisSet.size();i++)
    {
        vector<int> qdiffvec(BasisSet[0].M,0);
        for (unsigned int j=i;j<BasisSet.size();j++) // starts with j=i to exploit Hermiticity (Hij = Hji)
        {
            double Vij = 0;
            int qdiff = 0; // total number of quanta difference between states
            int mchange = 0; // number of modes with nonzero change in quanta
            QDiffVec(BasisSet[i], BasisSet[j], qdiff,mchange,qdiffvec);
            if(qdiff <= fcmax && mchange <= fcmax && qdiff%2==0){ 
                // States cannot differ by more than fcmax quanta 
                for (unsigned int k=0;k<QuarticFC.size();k++)
                {
                    if (ScreenState(qdiff,mchange,qdiffvec,QuarticFC[k]))
                    {// Screen force constants that cannot connect basis states i and j
                        //Add anharmonic matrix elements
                        Vij += AnharmPot(BasisSet[i], BasisSet[j], QuarticFC[k]);
                    }
                }
                for (unsigned int k=0;k<SexticFC.size();k++)
                {
                    if (ScreenState(qdiff,mchange,qdiffvec,SexticFC[k]))
                    {// Screen force constants that cannot connect basis states i and j
                        //Add anharmonic matrix elements
                        Vij += AnharmPot(BasisSet[i], BasisSet[j], SexticFC[k]);
                    }
                }
                H(i,j) += Vij;
                if(i!=j){
                    H(j,i) += Vij; // Exploiting Hermiticity (Hij=Hji)
                }
            }
            if(qdiff <= fcmax-1 && mchange <= fcmax-1 && qdiff%2==1){ 
                // fcmax-1 assumes max order is even (4th or 6th)
                // States cannot differ by more than fcmax quanta 
                for (unsigned int k=0;k<CubicFC.size();k++)
                {
                    if (ScreenState(qdiff,mchange,qdiffvec,CubicFC[k]))
                    {// Screen force constants that cannot connect basis states i and j
                        //Add anharmonic matrix elements
                        Vij += AnharmPot(BasisSet[i], BasisSet[j], CubicFC[k]);
                    }
                }
                for (unsigned int k=0;k<QuinticFC.size();k++)
                {
                    if (ScreenState(qdiff,mchange,qdiffvec,QuinticFC[k]))
                    {// Screen force constants that cannot connect basis states i and j
                        //Add anharmonic matrix elements
                        Vij += AnharmPot(BasisSet[i], BasisSet[j], QuinticFC[k]);
                    }
                }
                H(i,j) += Vij;
                if(i!=j){
                    H(j,i) += Vij; // Exploiting Hermiticity (Hij=Hji)
                }
            }
        }
    }
    return;
};

//Hamiltonian operators
void ZerothHamSparse(vector<Trip>& HTrip, std::vector<WaveFunction> &BasisSet, std::vector<double> &Frequencies)
{
    //Calculate the harmonic Hamiltonian matrix elements in sparse form of Eigen::Triplet
    //Harmonic matrix elements
    #pragma omp parallel for 
    for (unsigned int i=0;i<BasisSet.size();i++)
    {
        //Loop over all modes
        double Ei = 0.; //Hii matrix element
        for (int j=0;j<BasisSet[i].M;j++)
        {
            //Calculate partial energies
            double Ej = 0.5;
            Ej += BasisSet[i].Modes[j].Quanta;
            Ej *= BasisSet[i].Modes[j].Freq;
            //Update matrix element
            Ei += Ej;
        }
        //Update Hamiltonian
        #pragma omp critical
        HTrip.push_back(Trip(i,i,Ei/2.)); // Dividing by 2 so I can do H=H+H*
    }
    return;
};

void AnharmHamSparse(vector<Trip>& HTrip, std::vector<WaveFunction> &BasisSet, std::vector<double> &Frequencies, std::vector<FConst> &AnharmFC, std::vector<FConst> &CubicFC, std::vector<FConst> &QuarticFC, std::vector<FConst> &QuinticFC, std::vector<FConst> &SexticFC)
{
    //Add anharmonic terms to the Sparse Hamiltonian in sparse form of Eigen::Triplet
    int fcmax=0;
    for (unsigned int k=0;k<AnharmFC.size();k++)
    { //Order of maximum anharmonic term
        if(AnharmFC[k].fcpow.size()>fcmax)
        {
            fcmax = AnharmFC[k].fcpow.size();
        }
    }
    #pragma omp parallel for
    for (unsigned int i=0;i<BasisSet.size();i++)
    {
        vector<int> qdiffvec(BasisSet[0].M,0);
        for (unsigned int j=i;j<BasisSet.size();j++) // starts with j=i to exploit Hermiticity (Hij = Hji)
        {
            double Vij = 0;
            int mchange = 0; // number of modes with nonzero change in quanta
            int qdiff = 0; // total number of quanta difference between states
            QDiffVec(BasisSet[i], BasisSet[j],qdiff,mchange,qdiffvec);
            if(qdiff <= fcmax && mchange <= fcmax && qdiff%2==0){ 
                // States cannot differ by more than fcmax quanta
                double W = 0.;
                double V = 0.;
                for (unsigned int k=0;k<QuarticFC.size();k++){
                    if (ScreenState(qdiff,mchange,qdiffvec,QuarticFC[k])){
                        // Screen force constants that cannot connect basis states i and j
                        //Add anharmonic matrix elements
                        double val = AnharmPot(BasisSet[i], BasisSet[j], QuarticFC[k]);
                        Vij += val; 
                        if(abs(val) > abs(W)){
                            W = val;
                            V = QuarticFC[k].fc;
                        }
                    }
                }
                for (unsigned int k=0;k<SexticFC.size();k++){
                    if (ScreenState(qdiff,mchange,qdiffvec,SexticFC[k])){
                        // Screen force constants that cannot connect basis states i and j
                        //Add anharmonic matrix elements
                        double val = AnharmPot(BasisSet[i], BasisSet[j], SexticFC[k]); 
                        Vij += val;
                        if(abs(val) > abs(W)){
                            W = val;
                            V = SexticFC[k].fc;
                        }
                    }
                }
                #pragma omp critical
                if (abs(Vij) > 1e-12)
                {
                    if(i==j){
                        HTrip.push_back(Trip(i,j,Vij/2.)); // Dividing by 2 so I can do H=H+H*

                    }else{
                        HTrip.push_back(Trip(i,j,Vij)); // Exploiting Hermiticity (Hij=Hji)
                    }
                }
            }
            if(qdiff <= fcmax-1 && mchange <= fcmax-1 && qdiff%2==1){
                // fcmax-1 assumes max order is even (4th or 6th)
                // States cannot differ by more than fcmax quanta 
                double W = 0.;
                double V = 0.;
                for (unsigned int k=0;k<CubicFC.size();k++){
                    if (ScreenState(qdiff,mchange,qdiffvec,CubicFC[k])){
                        // Screen force constants that cannot connect basis states i and j
                        //Add anharmonic matrix elements
                        
                        double val = AnharmPot(BasisSet[i], BasisSet[j], CubicFC[k]);
                        Vij += val; 
                        if(abs(val) > abs(W)){
                            W = val;
                            V = CubicFC[k].fc;
                        }
                        //Vij += AnharmPot(i,j,CubicFC[k]);
                    }
                }
                for (unsigned int k=0;k<QuinticFC.size();k++){
                    if (ScreenState(qdiff,mchange,qdiffvec,QuinticFC[k])){
                        // Screen force constants that cannot connect basis states i and j
                        //Add anharmonic matrix elements
                        double val = AnharmPot(BasisSet[i], BasisSet[j], QuinticFC[k]);
                        Vij += val; 
                        if(abs(val) > abs(W)){
                            W = val;
                            V = QuinticFC[k].fc;
                        }
                        //Vij += AnharmPot(i,j,QuinticFC[k]);
                    }
                }
                #pragma omp critical
                if (abs(Vij) > 1e-12)
                {
                    if(i==j){
                        HTrip.push_back(Trip(i,j,Vij/2.)); // Dividing by 2 so I can do H=H+H*

                    }else{
                        HTrip.push_back(Trip(i,j,Vij)); // Exploiting Hermiticity (Hij=Hji)
                    }
                }
            }
        }
    }
    return;
};

inline void MakeHamSparse(SpMat &HSp, std::vector<WaveFunction> &BasisSet, std::vector<double> &Frequencies, std::vector<FConst> &AnharmFC, std::vector<FConst> &CubicFC, std::vector<FConst> &QuarticFC, std::vector<FConst> &QuinticFC, std::vector<FConst> &SexticFC, bool MakeZeroth, bool MakeAnharm)
{
    //Build the sparse CI Hamiltonian
    vector< Trip > HTrip;
    if (MakeZeroth) ZerothHamSparse(HTrip, BasisSet, Frequencies);
    if (MakeAnharm) AnharmHamSparse(HTrip, BasisSet, Frequencies, AnharmFC, CubicFC, QuarticFC, QuinticFC, SexticFC);
    HSp.setFromTriplets(HTrip.begin(),HTrip.end());
    HSp.makeCompressed();
    HTrip = vector< Trip >(); // Free memory
    SpMat HSpT = HSp.transpose();
    HSpT.makeCompressed();
    HSp += HSpT; // Complete symmetric matrix
    HSpT = SpMat(1,1); // Free memory
    HSp.makeCompressed();
    cout << "The Hamiltonian is " << fixed << setprecision(2) << 
        100*(1.-(double)HSp.nonZeros()/(double)HSp.size()) << "% sparse." << endl;
    return; 
}

inline void MakeHamDense(Eigen::MatrixXd &H, std::vector<WaveFunction> &BasisSet, std::vector<double> &Frequencies, std::vector<FConst> &AnharmFC, std::vector<FConst> &CubicFC, std::vector<FConst> &QuarticFC, std::vector<FConst> &QuinticFC, std::vector<FConst> &SexticFC)
{
    ZerothHam(H, BasisSet);
    AnharmHam(H, BasisSet, Frequencies, AnharmFC, CubicFC, QuarticFC, QuinticFC, SexticFC);
    return;
}

inline void MakeHamAnharmDense(Eigen::MatrixXd &H, std::vector<WaveFunction> &BasisSet, std::vector<double> &Frequencies, std::vector<FConst> &AnharmFC, std::vector<FConst> &CubicFC, std::vector<FConst> &QuarticFC, std::vector<FConst> &QuinticFC, std::vector<FConst> &SexticFC)
{
    AnharmHam(H, BasisSet, Frequencies, AnharmFC, CubicFC, QuarticFC, QuinticFC, SexticFC);
    return;
}

//Utility functions
std::tuple<Eigen::VectorXd, Eigen::MatrixXd> SparseDiagonalizeCPP(std::vector<WaveFunction> &BasisSet, std::vector<double> &Frequencies, std::vector<FConst> &AnharmFC, std::vector<FConst> &CubicFC, std::vector<FConst> &QuarticFC, std::vector<FConst> &QuinticFC, std::vector<FConst> &SexticFC, int NEig)
{
    SpMat H(BasisSet.size(), BasisSet.size());
    MakeHamSparse(H, BasisSet, Frequencies, AnharmFC, CubicFC, QuarticFC, QuinticFC, SexticFC, true, true);
    typedef SparseSymMatProd<double,Eigen::Lower,0,ptrdiff_t> SparseMVProd;
    SparseMVProd op(H);
    // Construct eigen solver object, requesting the largest three eigenvalues
    int NCV = 0;
//    int NState = 0;
//    if(NEig > BasisSet.size()){
//        NState = BasisSet.size()-1;
//    }else{
//        NState = NEig;
//    }
    NCV = max(2*NEig+1,20); // Default from Scipy's Lanczos/Arnoldi implementation
    Spectra::SymEigsSolver<double, Spectra::SMALLEST_ALGE, SparseMVProd> eigs(&op, NEig, NCV);
    // Initialize and compute
    eigs.init();
    int nconv = eigs.compute(1000,1e-10,Spectra::SMALLEST_ALGE);
    Eigen::VectorXd E;
    Eigen::MatrixXd Psi;
    if(eigs.info() == SUCCESSFUL){
        E = eigs.eigenvalues().real();
        Psi = eigs.eigenvectors().real();
    }else{
        cout << "Error: Eigenvalues did not converge." << endl; exit(0);}
    return std::make_tuple(E, Psi);

}

std::tuple<Eigen::VectorXd, Eigen::MatrixXd> DenseDiagonalizeCPP(std::vector<WaveFunction> &BasisSet, std::vector<double> &Frequencies, std::vector<FConst> &AnharmFC, std::vector<FConst> &CubicFC, std::vector<FConst> &QuarticFC, std::vector<FConst> &QuinticFC, std::vector<FConst> &SexticFC)
{
    //Wrapper for the Eigen diagonalization
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(BasisSet.size(), BasisSet.size());
    MakeHamDense(H, BasisSet, Frequencies, AnharmFC, CubicFC, QuarticFC, QuinticFC, SexticFC);
    SelfAdjointEigenSolver<MatrixXd> SE; //Schrodinger equation
    SE.compute(H); //Diagonalize the matrix
    Eigen::VectorXd E = SE.eigenvalues().real(); //Extract frequencies
    Eigen::MatrixXd Psi = SE.eigenvectors().real(); //Extract CI vectors
    return std::make_tuple(E, Psi);
}

Eigen::MatrixXd GenerateHamV(std::vector<WaveFunction> &BasisSet, std::vector<double> &Frequencies, std::vector<FConst> &AnharmFC, std::vector<FConst> &CubicFC, std::vector<FConst> &QuarticFC, std::vector<FConst> &QuinticFC, std::vector<FConst> &SexticFC)
{
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(BasisSet.size(), BasisSet.size());
    MakeHamDense(H, BasisSet, Frequencies, AnharmFC, CubicFC, QuarticFC, QuinticFC, SexticFC);
    return H;
}

Eigen::MatrixXd GenerateHam0V(std::vector<WaveFunction> &BasisSet, std::vector<double> &Frequencies)
{
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(BasisSet.size(), BasisSet.size());
    ZerothHam(H, BasisSet);
    return H;
}

Eigen::MatrixXd GenerateHamAnharmV(std::vector<WaveFunction> &BasisSet, std::vector<double> &Frequencies, std::vector<FConst> &AnharmFC, std::vector<FConst> &CubicFC, std::vector<FConst> &QuarticFC, std::vector<FConst> &QuinticFC, std::vector<FConst> &SexticFC)
{
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(BasisSet.size(), BasisSet.size());
    AnharmHam(H, BasisSet, Frequencies, AnharmFC, CubicFC, QuarticFC, QuinticFC, SexticFC);
    return H;
}

SpMat GenerateSparseHamV(std::vector<WaveFunction> &BasisSet, std::vector<double> &Frequencies, std::vector<FConst> &AnharmFC, std::vector<FConst> &CubicFC, std::vector<FConst> &QuarticFC, std::vector<FConst> &QuinticFC, std::vector<FConst> &SexticFC)
{
    SpMat H(BasisSet.size(), BasisSet.size());
    MakeHamSparse(H, BasisSet, Frequencies, AnharmFC, CubicFC, QuarticFC, QuinticFC, SexticFC, true, true);
    return H;
}

SpMat GenerateSparseHamAnharmV(std::vector<WaveFunction> &BasisSet, std::vector<double> &Frequencies, std::vector<FConst> &AnharmFC, std::vector<FConst> &CubicFC, std::vector<FConst> &QuarticFC, std::vector<FConst> &QuinticFC, std::vector<FConst> &SexticFC)
{
    SpMat H(BasisSet.size(), BasisSet.size());
    MakeHamSparse(H, BasisSet, Frequencies, AnharmFC, CubicFC, QuarticFC, QuinticFC, SexticFC, false, true);
    return H;
}

// PT2.cpp

std::vector<double> StateProbability(std::vector<double>& CMax)
{
    double CTotal = 0.0;
    for (unsigned int i = 0; i < CMax.size(); i++)
    {
        CTotal += abs(CMax[i]);
    }

    std::vector<double> CProbability;
    for (int i = 0; i < CMax.size(); i++)
    {
        CProbability.push_back(abs(CMax[i]) / CTotal);
    }
    return CProbability;
}

void FillWalkers(std::map<int, int>& WalkerPopulation, std::vector<double>& C, int Nd)
{
    std::random_device RD;
    std::mt19937 Gen(RD());

    std::discrete_distribution<> Distribution(C.begin(), C.end());

    for (unsigned int i = 0; i < Nd; i++) WalkerPopulation[Distribution(Gen)]++;
}

std::vector<double> DoPT2(MatrixXd& Evecs, VectorXd& Evals, std::vector<WaveFunction> &BasisSet, std::vector<FConst> &AnharmHB, std::vector<FConst> &AnharmFC, std::vector<FConst> &CubicFC, std::vector<FConst> &QuarticFC, std::vector<FConst> &QuinticFC, std::vector<FConst> &SexticFC, double PT2_Eps, int NEig)
{
    int N_opt;
    if(NEig > BasisSet.size()){ // If we don't have enough states to optimize for yet
        N_opt = BasisSet.size();
    }else{ // If number of states exceeds the number selected to optimize
        N_opt = NEig;
    }

    vector<double> DeltaE(N_opt,0.);  // Vector will contain the PT correction for each eigenvalue
    int fcmax=0;
    for (unsigned int k=0;k<AnharmFC.size();k++){ //Order of maximum anharmonic term
        if(AnharmFC[k].fcpow.size()>fcmax){
            fcmax = AnharmFC[k].fcpow.size();
        }
    }

    for (unsigned int n = 0; n < N_opt; n++)
    {
        std::vector<WaveFunction> PTBasisSet = AddStatesHB(BasisSet, AnharmHB, Evecs.col(n), PT2_Eps);
        std::cout << " Perturbative space for state " << n << " contains " << PTBasisSet.size() << " basis states." << std::endl;

        #pragma omp parallel for
        for (unsigned int a = 0; a < PTBasisSet.size(); a++)
        {
            double HaiCi = 0.0;
            vector<int> qdiffvec(BasisSet[0].M,0);
            for (unsigned int i = 0; i < BasisSet.size(); i++)
            {
                double Hai = 0;
                int mchange = 0; // number of modes with nonzero change in quanta
                int qdiff = 0; // total number of quanta difference 
                QDiffVec(PTBasisSet[a], BasisSet[i], qdiff, mchange, qdiffvec);
                if(qdiff <= fcmax && mchange <= fcmax && qdiff%2==0)
                { 
                    // States cannot differ by more than fcmax quanta
                    for (unsigned int k=0;k<QuarticFC.size();k++)
                    {
                        if ( ScreenState(qdiff,mchange,qdiffvec,QuarticFC[k]) ){
                            //Screen force constants for connection                      
                            //Add anharmonic matrix elements
                            Hai += AnharmPot(PTBasisSet[a], BasisSet[i], QuarticFC[k]);
                        }
                    }
                    for (unsigned int k=0;k<SexticFC.size();k++)
                    {
                        if ( ScreenState(qdiff,mchange,qdiffvec,SexticFC[k]) ){
                            //Screen force constants for connection                      
                            //Add anharmonic matrix elements
                            Hai += AnharmPot(PTBasisSet[a], BasisSet[i], SexticFC[k]);
                        }
                    }
                    HaiCi += Hai * Evecs(i, n); // C_i Hai for each eigenvalue of interest
                }
                if(qdiff <= fcmax-1 && mchange <= fcmax-1 && qdiff%2==1)
                { 
                    // fcmax-1 assumes 4th or 6th max order 
                    // States cannot differ by more than fcmax quanta
                    for (unsigned int k=0;k<CubicFC.size();k++)
                    {
                        if ( ScreenState(qdiff,mchange,qdiffvec,CubicFC[k]) ){
                            //Screen force constants for connection                      
                            //Add anharmonic matrix elements
                            Hai += AnharmPot(PTBasisSet[a], BasisSet[i], CubicFC[k]);
                        }
                    }
                    for (unsigned int k=0;k<QuinticFC.size();k++)
                    {
                        if ( ScreenState(qdiff,mchange,qdiffvec,QuinticFC[k]) ){
                            //Screen force constants for connection                      
                            //Add anharmonic matrix elements
                            Hai += AnharmPot(PTBasisSet[a], BasisSet[i], QuinticFC[k]);
                        }
                    }
                    HaiCi += Hai * Evecs(i, n); // C_i Hai for each eigenvalue of interest
                }
            }
            double Ea = 0.; //Hii matrix element
            for (unsigned int j = 0; j < PTBasisSet[a].M; j++)
            {
              //Calculate partial energies
              double Ej = 0.5;
              Ej += PTBasisSet[a].Modes[j].Quanta;
              Ej *= PTBasisSet[a].Modes[j].Freq;
              //Update matrix element
              Ea += Ej;
            }
            vector<int> zerodiffvec(BasisSet[0].M,0);
            int qdiff=0;
            int mchange=0;
            for (unsigned int k = 0; k < QuarticFC.size(); k++) // Only even-ordered fc can affect this
            {
                if (ScreenState(qdiff, mchange, zerodiffvec, QuarticFC[k]))
                {
                    // Screen force constants that cannot connect basis states a and a
                    //Add anharmonic matrix elements
                    Ea += AnharmPot(PTBasisSet[a], PTBasisSet[a], QuarticFC[k]);
                }
            }
            for (unsigned int k = 0; k < SexticFC.size(); k++) // Only even-ordered fc can affect this
            {    
                if (ScreenState(qdiff, mchange, zerodiffvec, SexticFC[k]))
                {
                    // Screen force constants that cannot connect basis states a and a
                    //Add anharmonic matrix elements
                    Ea += AnharmPot(PTBasisSet[a], PTBasisSet[a], SexticFC[k]);
                }
            }
            #pragma omp atomic // Will cause floating point error if blindly done in parallel
            DeltaE[n] += pow(HaiCi, 2) / (Evals(n) - Ea);
        }
    }
    return DeltaE;    
}

std::tuple<std::vector<double>, std::vector<double>> DoSPT2(MatrixXd& Evecs, VectorXd& Evals, std::vector<WaveFunction> &BasisSet, std::vector<FConst> &AnharmHB, std::vector<FConst> &AnharmFC, std::vector<FConst> &CubicFC, std::vector<FConst> &QuarticFC, std::vector<FConst> &QuinticFC, std::vector<FConst> &SexticFC, double PT2_Eps, int NEig, int Nd, int Ns, bool SemiStochastic = false, double PT2_Eps2 = 0.0)
{
    int N_opt;
    if(NEig > BasisSet.size()){ // If we don't have enough states to optimize for yet
        N_opt = BasisSet.size();
    }else{ // If number of states exceeds the number selected to optimize
        N_opt = NEig;
    }

    vector<double> DeltaE(N_opt,0.);  // Vector will contain the PT correction for each eigenvalue
    std::vector<std::vector<double>> DeltaESample;
    std::vector<double> DeltaEDet(N_opt, 0.0);
    int fcmax=0;
    for (unsigned int k=0;k<AnharmFC.size();k++){ //Order of maximum anharmonic term
        if(AnharmFC[k].fcpow.size()>fcmax){
            fcmax = AnharmFC[k].fcpow.size();
        }
    }

    for (unsigned int n = 0; n < N_opt; n++)
    {
        // For each state, we determine the perturbative basis and populate the state with walkers.
        std::vector<WaveFunction> DetPTBasisSet;
        if (SemiStochastic)
        {
            DetPTBasisSet = AddStatesHB(BasisSet, AnharmHB, Evecs.col(n), PT2_Eps);
            std::cout << "Perturbative space for state " << n << " contains " << DetPTBasisSet.size() << " deterministic basis states." << std::endl;
        }

        std::vector<double> Cn;
        for (unsigned int i = 0; i < Evecs.rows(); i++) Cn.push_back(Evecs(i ,n));

        // If we are doing semi-stochastic, I need to do a deterministic calculation first. That will be executed here with the DetPTBasisSet.
        if (SemiStochastic)
        {
            #pragma omp parallel for
            for (unsigned int a = 0; a < DetPTBasisSet.size(); a++)
            {
                double HaiCi = 0.0;
                vector<int> qdiffvec(BasisSet[0].M,0);
                for (unsigned int i = 0; i < BasisSet.size(); i++)
                {
                    double Hai = 0;
                    int mchange = 0; // number of modes with nonzero change in quanta
                    int qdiff = 0; // total number of quanta difference 
                    QDiffVec(DetPTBasisSet[a], BasisSet[i], qdiff, mchange, qdiffvec);
                    if(qdiff <= fcmax && mchange <= fcmax && qdiff%2==0)
                    { 
                        // States cannot differ by more than fcmax quanta
                        for (unsigned int k=0;k<QuarticFC.size();k++)
                        {
                            if ( ScreenState(qdiff,mchange,qdiffvec,QuarticFC[k]) ){
                                //Screen force constants for connection                      
                                //Add anharmonic matrix elements
                                Hai += AnharmPot(DetPTBasisSet[a], BasisSet[i], QuarticFC[k]);
                            }
                        }
                        for (unsigned int k=0;k<SexticFC.size();k++)
                        {
                            if ( ScreenState(qdiff,mchange,qdiffvec,SexticFC[k]) ){
                                //Screen force constants for connection                      
                                //Add anharmonic matrix elements
                                Hai += AnharmPot(DetPTBasisSet[a], BasisSet[i], SexticFC[k]);
                            }
                        }
                        HaiCi += Hai * Evecs(i, n); // C_i Hai for each eigenvalue of interest
                    }
                    if(qdiff <= fcmax-1 && mchange <= fcmax-1 && qdiff%2==1)
                    { 
                        // fcmax-1 assumes 4th or 6th max order 
                        // States cannot differ by more than fcmax quanta
                        for (unsigned int k=0;k<CubicFC.size();k++)
                        {
                            if ( ScreenState(qdiff,mchange,qdiffvec,CubicFC[k]) ){
                                //Screen force constants for connection                      
                                //Add anharmonic matrix elements
                                Hai += AnharmPot(DetPTBasisSet[a], BasisSet[i], CubicFC[k]);
                            }
                        }
                        for (unsigned int k=0;k<QuinticFC.size();k++)
                        {
                            if ( ScreenState(qdiff,mchange,qdiffvec,QuinticFC[k]) ){
                                //Screen force constants for connection                      
                                //Add anharmonic matrix elements
                                Hai += AnharmPot(DetPTBasisSet[a], BasisSet[i], QuinticFC[k]);
                            }
                        }
                        HaiCi += Hai * Evecs(i, n); // C_i Hai for each eigenvalue of interest
                    }
                }
                double Ea = 0.; //Hii matrix element
                for (unsigned int j = 0; j < DetPTBasisSet[a].M; j++)
                {
                  //Calculate partial energies
                  double Ej = 0.5;
                  Ej += DetPTBasisSet[a].Modes[j].Quanta;
                  Ej *= DetPTBasisSet[a].Modes[j].Freq;
                  //Update matrix element
                  Ea += Ej;
                }
                vector<int> zerodiffvec(BasisSet[0].M,0);
                int qdiff=0;
                int mchange=0;
                for (unsigned int k = 0; k < QuarticFC.size(); k++) // Only even-ordered fc can affect this
                {
                    if (ScreenState(qdiff, mchange, zerodiffvec, QuarticFC[k]))
                    {
                        // Screen force constants that cannot connect basis states a and a
                        //Add anharmonic matrix elements
                        Ea += AnharmPot(DetPTBasisSet[a], DetPTBasisSet[a], QuarticFC[k]);
                    }
                }
                for (unsigned int k = 0; k < SexticFC.size(); k++) // Only even-ordered fc can affect this
                {    
                    if (ScreenState(qdiff, mchange, zerodiffvec, SexticFC[k]))
                    {
                        // Screen force constants that cannot connect basis states a and a
                        //Add anharmonic matrix elements
                        Ea += AnharmPot(DetPTBasisSet[a], DetPTBasisSet[a], SexticFC[k]);
                    }
                }
                #pragma omp atomic // Will cause floating point error if blindly done in parallel
                DeltaEDet[n] += pow(HaiCi, 2) / (Evals(n) - Ea);
            }
        } // end deterministic calculation for semi-stochastic

        std::vector<double> DeltaEs(Ns, 0.0);
        int PTBasisSize = 0;
        #pragma omp parallel for
        for (unsigned int s = 0; s < Ns; s++)
        {
            std::vector<double> WalkerProbability = StateProbability(Cn);
            std::map<int, int> WalkerPopulation;
            FillWalkers(WalkerPopulation, WalkerProbability, Nd);

            std::vector<WaveFunction> PTBasisSet;
            HashedStates PTBasisSetHashed;
            if (SemiStochastic)
            {
                std::vector<WaveFunction> VarDetBasisSet; // Actually I think this destructs outside of this if statement.
                for (const WaveFunction &WF : BasisSet) VarDetBasisSet.push_back(WF);
                for (const WaveFunction &WF : DetPTBasisSet) VarDetBasisSet.push_back(WF);
                for (std::map<int, int>::iterator it = WalkerPopulation.begin(); it != WalkerPopulation.end(); ++it)
                {
                    int i = it->first;
                    InternalAddStatesHB(VarDetBasisSet, PTBasisSetHashed, AnharmHB, i, Evecs(i, n), PT2_Eps2);
                }
            }
            else
            {
                for (std::map<int, int>::iterator it = WalkerPopulation.begin(); it != WalkerPopulation.end(); ++it)
                {
                    int i = it->first;
                    InternalAddStatesHB(BasisSet, PTBasisSetHashed, AnharmHB, i, Evecs(i, n), PT2_Eps);
                }
            }
            for (const WaveFunction &WF : PTBasisSetHashed) PTBasisSet.push_back(WF);
            PTBasisSetHashed = HashedStates();
            PTBasisSize += PTBasisSet.size();
            //std::cout << " Perturbative space for state " << n << " and sample " << s << " contains " << PTBasisSet.size() << " stochastic basis states." << std::endl;

            for (unsigned int a = 0; a < PTBasisSet.size(); a++)
            {
                double HaiCi = 0.0;
                double Hai2Ci2 = 0.0;
                vector<int> qdiffvec(BasisSet[0].M,0);
                for (std::map<int, int>::iterator it = WalkerPopulation.begin(); it != WalkerPopulation.end(); ++it)
                {
                    int i = it->first;
                    double Hai = 0;
                    int mchange = 0; // number of modes with nonzero change in quanta
                    int qdiff = 0; // total number of quanta difference 
                    QDiffVec(PTBasisSet[a], BasisSet[i], qdiff, mchange, qdiffvec);
                    if(qdiff <= fcmax && mchange <= fcmax && qdiff%2==0)
                    { 
                        // States cannot differ by more than fcmax quanta
                        for (unsigned int k=0;k<QuarticFC.size();k++)
                        {
                            if ( ScreenState(qdiff,mchange,qdiffvec,QuarticFC[k]) ){
                                //Screen force constants for connection                      
                                //Add anharmonic matrix elements
                                Hai += AnharmPot(PTBasisSet[a], BasisSet[i], QuarticFC[k]);
                            }
                        }
                        for (unsigned int k=0;k<SexticFC.size();k++)
                        {
                            if ( ScreenState(qdiff,mchange,qdiffvec,SexticFC[k]) ){
                                //Screen force constants for connection                      
                                //Add anharmonic matrix elements
                                Hai += AnharmPot(PTBasisSet[a], BasisSet[i], SexticFC[k]);
                            }
                        }
                        HaiCi += (Hai * Evecs(i, n) * WalkerPopulation[i]) / WalkerProbability[i]; // C_i Hai for each eigenvalue of interest
                        Hai2Ci2 += (pow(Hai, 2) * pow(Evecs(i, n), 2)) * (WalkerPopulation[i] * (Nd - 1) / WalkerProbability[i] - pow(WalkerPopulation[i], 2) / pow(WalkerProbability[i], 2));
                    }
                    if(qdiff <= fcmax-1 && mchange <= fcmax-1 && qdiff%2==1)
                    { 
                        // fcmax-1 assumes 4th or 6th max order 
                        // States cannot differ by more than fcmax quanta
                        for (unsigned int k=0;k<CubicFC.size();k++)
                        {
                            if ( ScreenState(qdiff,mchange,qdiffvec,CubicFC[k]) ){
                                //Screen force constants for connection                      
                                //Add anharmonic matrix elements
                                Hai += AnharmPot(PTBasisSet[a], BasisSet[i], CubicFC[k]);
                            }
                        }
                        for (unsigned int k=0;k<QuinticFC.size();k++)
                        {
                            if ( ScreenState(qdiff,mchange,qdiffvec,QuinticFC[k]) ){
                                //Screen force constants for connection                      
                                //Add anharmonic matrix elements
                                Hai += AnharmPot(PTBasisSet[a], BasisSet[i], QuinticFC[k]);
                            }
                        }
                        HaiCi += (Hai * Evecs(i, n) * WalkerPopulation[i]) / WalkerProbability[i]; // C_i Hai for each eigenvalue of interest
                        Hai2Ci2 += (pow(Hai, 2) * pow(Evecs(i, n), 2)) * (WalkerPopulation[i] * (Nd - 1) / WalkerProbability[i] - pow(WalkerPopulation[i], 2) / pow(WalkerProbability[i], 2));
                    }
                }
                double Ea = 0.; //Hii matrix element
                for (unsigned int j = 0; j < PTBasisSet[a].M; j++)
                {
                  //Calculate partial energies
                  double Ej = 0.5;
                  Ej += PTBasisSet[a].Modes[j].Quanta;
                  Ej *= PTBasisSet[a].Modes[j].Freq;
                  //Update matrix element
                  Ea += Ej;
                }
                vector<int> zerodiffvec(BasisSet[0].M,0);
                int qdiff=0;
                int mchange=0;
                for (unsigned int k = 0; k < QuarticFC.size(); k++) // Only even-ordered fc can affect this
                {
                    if (ScreenState(qdiff, mchange, zerodiffvec, QuarticFC[k]))
                    {
                        // Screen force constants that cannot connect basis states a and a
                        //Add anharmonic matrix elements
                        Ea += AnharmPot(PTBasisSet[a], PTBasisSet[a], QuarticFC[k]);
                    }
                }
                for (unsigned int k = 0; k < SexticFC.size(); k++) // Only even-ordered fc can affect this
                {    
                    if (ScreenState(qdiff, mchange, zerodiffvec, SexticFC[k]))
                    {
                        // Screen force constants that cannot connect basis states a and a
                        //Add anharmonic matrix elements
                        Ea += AnharmPot(PTBasisSet[a], PTBasisSet[a], SexticFC[k]);
                    }
                }
                #pragma omp atomic // Will cause floating point error if blindly done in parallel
                DeltaEs[s] += (pow(HaiCi, 2) + Hai2Ci2) / ((Evals(n) - Ea) * Nd * (Nd - 1));
            }
        }
        DeltaESample.push_back(DeltaEs);
        std::cout << " Perturbative space for state " << n << " contains " << (float)PTBasisSize / (float)Ns << " stochastic basis states on average." << std::endl;
    }

    std::vector<double> SigmaDeltaE(N_opt, 0.0);
    for (unsigned int n = 0; n < N_opt; n++)
    {
        for (unsigned int s = 0; s < Ns; s++)
        {
            DeltaE[n] += DeltaESample[n][s];
            //SigmaDeltaE[n] += pow(DeltaESample[n][s], 2);
        }
        //SigmaDeltaE[n] -= (pow(DeltaE[n], 2) / Ns);
        DeltaE[n] /= Ns;
        std::cout << setprecision(12);
        for (unsigned int s = 0; s < Ns; s++)
        {
            SigmaDeltaE[n] += pow(DeltaESample[n][s] - DeltaE[n], 2);
        }
        if (SemiStochastic) DeltaE[n] += DeltaEDet[n];
        SigmaDeltaE[n] /= (Ns - 1);
        SigmaDeltaE[n] = sqrt(SigmaDeltaE[n]); // This is sample standard deviation
        SigmaDeltaE[n] /= sqrt(Ns); // This is standard error
        //cout << DeltaE[n] << " " << SigmaDeltaE[n] << std::endl;
    }

    return std::make_tuple(DeltaE, SigmaDeltaE);
}

/*************************************************************************************
************************************* VSCF Functions *********************************
*************************************************************************************/

std::vector<Eigen::MatrixXd> GetVEffSLOW1CPP(std::vector<Eigen::SparseMatrix<double>> &AnharmTensor, std::vector<WaveFunction> Basis, std::vector<Eigen::MatrixXd> &CByModes, std::vector<std::vector<std::vector<int>>> &ModalSlices, std::vector<int> &MaxQuanta, std::vector<int> &ModeOcc, bool FirstV = false)
{
    int NModes = CByModes.size();
    std::vector<Eigen::MatrixXd> VEff;
    if (FirstV) NModes = 1;

    for (unsigned int Mode = 0; Mode < NModes; Mode++)
    {
        Eigen::MatrixXd V = Eigen::MatrixXd::Zero(MaxQuanta[Mode], MaxQuanta[Mode]);
        #pragma omp parallel for
        for (unsigned int n = 0; n < MaxQuanta[Mode]; n++)
        {
            for (unsigned int m = n; m < MaxQuanta[Mode]; m++)
            {
                double Vnm = 0.0;
                for (auto HAnharm : AnharmTensor)
                {
                    for (unsigned int k = 0; k < HAnharm.outerSize(); k++)
                    {
                        for (Eigen::SparseMatrix<double, Eigen::ColMajor>::InnerIterator it(HAnharm, k); it; ++it)
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
                                    
                                    double CVC = it.value();
                                    for (unsigned int a = 0; a < BraModeBasis.size(); a++)
                                    {
                                        if (a != Mode)
                                        {
                                            CVC *= CByModes[a](BraModeBasis[a], ModeOcc[a]);
                                        }
                                    }
                                    for (unsigned int a = 0; a < KetModeBasis.size(); a++)
                                    {
                                        if (a != Mode)
                                        {
                                            CVC *= CByModes[a](KetModeBasis[a], ModeOcc[a]);
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

std::vector<Eigen::MatrixXd> MakeCTensorsCPP(std::vector<Eigen::MatrixXd> &Cs, std::vector<std::vector<int>> QUniques, std::vector<int> ModeOcc, std::vector<std::vector<WaveFunction>> &RestrictedBases)
{
    std::vector<Eigen::MatrixXd> CTensors;
    for (unsigned int q = 0; q < RestrictedBases.size(); q++)
    {
        Eigen::MatrixXd CTensor = Eigen::MatrixXd::Ones(RestrictedBases[q].size(), 1);
        #pragma omp parallel for
        for (unsigned int i = 0; i < RestrictedBases[q].size(); i++)
        {
            for (unsigned int m = 0; m < QUniques[q].size(); m++)
            {
                CTensor(i, 0) *= Cs[QUniques[q][m]](RestrictedBases[q][i].Modes[QUniques[q][m]].Quanta, ModeOcc[QUniques[q][m]]);
            }
        }
        CTensors.push_back(CTensor);
    }
    return CTensors;
}

Eigen::MatrixXd MakeCTensorCPP(std::vector<Eigen::MatrixXd> &Cs, std::vector<int> QUnique, std::vector<std::vector<int>> ModeOccs, std::vector<WaveFunction> &RestrictedBasis)
{
    Eigen::MatrixXd CTensor = Eigen::MatrixXd::Ones(RestrictedBasis.size(), ModeOccs.size());
    for (unsigned int k = 0; k < ModeOccs.size(); k++)
    {
        #pragma omp parallel for
        for (unsigned int i = 0; i < RestrictedBasis.size(); i++)
        {
            for (unsigned int m = 0; m < QUnique.size(); m++)
            {
                CTensor(i, k) *= Cs[QUnique[m]](RestrictedBasis[i].Modes[QUnique[m]].Quanta, ModeOccs[k][QUnique[m]]);
            }
        }
    }
    return CTensor;
}

std::tuple<std::vector<int>, std::vector<WaveFunction>> GetRestrictedSlice(std::vector<WaveFunction> &RestrictedBasis, unsigned int n, unsigned int Mode)
{
    std::vector<int> SlicedIndex;
    std::vector<WaveFunction> SlicedBasis;
    for (unsigned int i = 0; i < RestrictedBasis.size(); i++)
    {
        if (RestrictedBasis[i].Modes[Mode].Quanta == n)
        {
            SlicedBasis.push_back(RestrictedBasis[i]);
            SlicedIndex.push_back(i);
        }
    }
    return std::make_tuple(SlicedIndex, SlicedBasis);
}

Eigen::MatrixXd SliceMatrix(Eigen::SparseMatrix<double> M, std::vector<int> Ind, bool Row)
{
    Eigen::MatrixXd SlicedM;
    if (Row)
    {
        SlicedM = Eigen::MatrixXd::Zero(Ind.size(), M.cols());
        for (unsigned int i = 0; i < Ind.size(); i++)
        {
            SlicedM.row(i) = M.row(Ind[i]);
        }
    }
    else
    {
        SlicedM = Eigen::MatrixXd::Zero(M.rows(), Ind.size());
        for (unsigned int i = 0; i < Ind.size(); i++)
        {
            SlicedM.col(i) = M.col(Ind[i]);
        }
    }
    return SlicedM;
}

Eigen::MatrixXd SliceMatrix(Eigen::MatrixXd M, std::vector<int> Ind, bool Row)
{
    Eigen::MatrixXd SlicedM;
    if (Row)
    {
        SlicedM = Eigen::MatrixXd::Zero(Ind.size(), M.cols());
        for (unsigned int i = 0; i < Ind.size(); i++)
        {
            SlicedM.row(i) = M.row(Ind[i]);
        }
    }
    else
    {
        SlicedM = Eigen::MatrixXd::Zero(M.rows(), Ind.size());
        for (unsigned int i = 0; i < Ind.size(); i++)
        {
            SlicedM.col(i) = M.col(Ind[i]);
        }
    }
    return SlicedM;
}

std::vector<Eigen::MatrixXd> GetVEffSLOW2CPP(std::vector<Eigen::SparseMatrix<double>> &AnharmTensor, std::vector<std::vector<WaveFunction>> &RestrictedBases, std::vector<std::vector<int>> &QUniques, std::vector<Eigen::MatrixXd> &Cs, std::vector<int> MaxQuanta, std::vector<int> ModeOcc, bool FirstV)
{
    std::vector<Eigen::MatrixXd> VEff;
    int NModes = Cs.size();
    if (FirstV) NModes = 1;
    std::vector<std::vector<int>> ModeOccVec;
    ModeOccVec.push_back(ModeOcc);

    for (unsigned int Mode = 0; Mode < NModes; Mode++)
    {
        Eigen::MatrixXd Vm = Eigen::MatrixXd::Zero(MaxQuanta[Mode], MaxQuanta[Mode]);
        #pragma omp parallel for 
        for (unsigned int q = 0; q < AnharmTensor.size(); q++)
        {
            Eigen::MatrixXd Vmq = Eigen::MatrixXd::Zero(MaxQuanta[Mode], MaxQuanta[Mode]);
            if (!std::count(QUniques[q].begin(), QUniques[q].end(), Mode))
            {
                Eigen::MatrixXd CTensor = MakeCTensorCPP(Cs, QUniques[q], ModeOccVec, RestrictedBases[q]);
                Vmq = Eigen::MatrixXd::Identity(MaxQuanta[Mode], MaxQuanta[Mode]) * (CTensor.transpose() * AnharmTensor[q] * CTensor)(0, 0);
            }
            else
            {
                std::vector<int> Qs = QUniques[q];
                std::vector<int>::iterator iMode = std::find(Qs.begin(), Qs.end(), Mode);
                Qs.erase(iMode);
                for (unsigned int i = 0; i < MaxQuanta[Mode]; i++)
                {
                    std::tuple<std::vector<int>, std::vector<WaveFunction>> SliceI = GetRestrictedSlice(RestrictedBases[q], i, Mode);
                    Eigen::MatrixXd CTensorI = MakeCTensorCPP(Cs, Qs, ModeOccVec, std::get<1>(SliceI));
                    Eigen::MatrixXd HI = SliceMatrix(AnharmTensor[q], std::get<0>(SliceI), true);
                    for (unsigned int j = i; j < MaxQuanta[Mode]; j++)
                    {
                        std::tuple<std::vector<int>, std::vector<WaveFunction>> SliceJ;
                        Eigen::MatrixXd CTensorJ;
                        if (j == i)
                        {
                            SliceJ = SliceI;
                            CTensorJ = CTensorI;
                        }
                        else
                        {
                            SliceJ = GetRestrictedSlice(RestrictedBases[q], j, Mode);
                            CTensorJ = MakeCTensorCPP(Cs, Qs, ModeOccVec, std::get<1>(SliceJ));
                        }
                        Eigen::MatrixXd HIJ = SliceMatrix(HI, std::get<0>(SliceJ), false);
                        Vmq(i, j) = (CTensorI.transpose() * HIJ * CTensorJ)(0, 0);
                        Vmq(j, i) = Vmq(i, j);
                    }
                }
            }
            #pragma omp critical
            Vm += Vmq;
        }
        VEff.push_back(Vm);
    }
    return VEff;
}

/*************************************************************************************
************************************ HO2MO Functions *********************************
*************************************************************************************/

// Computers C.T @ V @ C for each power and mode. First index is Mode, second index is power.
std::vector<std::vector<double>> ContractedAnharmonicPotential(std::vector<Eigen::MatrixXd> &Cs, std::vector<Eigen::SparseMatrix<double>> &Vs, std::vector<int> &ModeOcc1, std::vector<int> &ModeOcc2)
{
    unsigned int NumP = Vs.size();
    unsigned int NumM = Cs.size();
    
    std::vector<std::vector<double>> Ys;
    for (unsigned int m = 0; m < NumM; m++)
    {
        std::vector<double> Ym;
        Eigen::VectorXd C1 = Cs[m].col(ModeOcc1[m]);
        Eigen::VectorXd C2 = Cs[m].col(ModeOcc2[m]);
        for (int p = 0; p < NumP; p++)
        {
            double Y = C1.transpose() * Vs[p] * C2;
            /*
            double Y = 0.0;
            for (int i = 0; i < C1.rows(); i++)
            {
                for (int pi = i - p; pi <= i + p; pi += 2)
                {
                    if (pi >= 0 && pi < C2.rows()) Y += C1[i] * C2[pi] * Vs[p].coeff(i, pi);
                }
            }
            */
            Ym.push_back(Y);
        }
        Ys.push_back(Ym);
    }
    return Ys;
}

std::vector<Eigen::MatrixXd> GetVEffCPP(std::vector<Eigen::SparseMatrix<double>> &AnharmTensor, std::vector<double> &FCs, std::vector<std::vector<int>> &QUniques, std::vector<std::vector<int>> &QPowers, std::vector<Eigen::MatrixXd> &Cs, std::vector<int> &MaxQuanta, std::vector<int> &ModeOcc1, std::vector<int> &ModeOcc2, bool FirstV)
{
    std::vector<Eigen::MatrixXd> VEff;
    int NModes = Cs.size();
    if (FirstV) NModes = 1;

    std::vector<std::vector<double>> Ys = ContractedAnharmonicPotential(Cs, AnharmTensor, ModeOcc1, ModeOcc2);

    for (unsigned int Mode = 0; Mode < NModes; Mode++)
    {
        Eigen::MatrixXd Vm = Eigen::MatrixXd::Zero(MaxQuanta[Mode], MaxQuanta[Mode]);
        #pragma omp parallel for 
        for (unsigned int q = 0; q < FCs.size(); q++)
        {
            Eigen::MatrixXd Vmq; // = Eigen::MatrixXd::Zero(MaxQuanta[Mode], MaxQuanta[Mode]);
            if (!std::count(QUniques[q].begin(), QUniques[q].end(), Mode))
            {
                Vmq = Eigen::MatrixXd::Identity(MaxQuanta[Mode], MaxQuanta[Mode]) * FCs[q];
                for (unsigned int Q = 0; Q < QUniques[q].size(); Q++)
                {
                    Vmq *= Ys[QUniques[q][Q]][QPowers[q][Q]];
                }
            }
            else
            {
                std::vector<int>::iterator itMode = std::find(QUniques[q].begin(), QUniques[q].end(), Mode);
                int iMode = itMode - QUniques[q].begin();
                Vmq = AnharmTensor[QPowers[q][iMode]] * FCs[q]; 
                for (unsigned int i = 0; i < Vmq.rows(); i++)
                {
                    for (unsigned int j = i; j < Vmq.cols(); j++)
                    {
                        if (abs(Vmq(i, j)) < 1e-12) continue;
                        else
                        {
                            for (unsigned int Q = 0; Q < QUniques[q].size(); Q++)
                            {
                                if (QUniques[q][Q] == Mode) continue;
                                Vmq(i, j) *= Ys[QUniques[q][Q]][QPowers[q][Q]];
                            }
                            Vmq(j, i) = Vmq(i, j);
                        }
                    }
                }
            }
            #pragma omp critical
            Vm += Vmq;
        }
        VEff.push_back(Vm);
    }
    return VEff;
}

/*************************************************************************************
************************************* VCI Functions **********************************
*************************************************************************************/

std::vector<int> CalcDiffModes(WaveFunction &Bi, WaveFunction &Bj)
{
    std::vector<int> DiffModes;
    for (unsigned int i = 0; i < Bi.Modes.size(); i++)
    {
        if (Bi.Modes[i].Quanta - Bj.Modes[i].Quanta != 0) DiffModes.push_back(i);
    }
    return DiffModes;
}

bool VectorContainedIn(std::vector<int> &b, std::vector<int> &B)
{
    for (auto &a : b)
    {
        if (std::find(B.begin(), B.end(), a) == B.end()) return false;
    }
    return true;
}

std::vector<double> ContractedHOTerms(std::vector<Eigen::MatrixXd> &Cs, std::vector<double> &Frequencies, std::vector<int> ModeOccI, std::vector<int> ModeOccJ)
{
    std::vector<double> Xs;
    for (unsigned int m = 0; m < Cs.size(); m++)
    {
        double Xm = 0.0;
        for (unsigned int i = 0; i < Cs[m].rows(); i++)
        {
            Xm += Cs[m](i, ModeOccI[m]) * Cs[m](i, ModeOccJ[m]) * (i + 0.5);
        }
        Xm *= Frequencies[m];
        Xs.push_back(Xm);
    }
    return Xs;
}

Eigen::MatrixXd VCIHamFromVSCF(std::vector<WaveFunction> &BasisSet, std::vector<double> &Frequencies, std::vector<FConst> &FCs, std::vector<Eigen::MatrixXd> &Cs, std::vector<Eigen::SparseMatrix<double>> &GenericV)
{
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(BasisSet.size(), BasisSet.size());

    std::vector<int> MaxQuanta;
    for (Eigen::MatrixXd &C : Cs) MaxQuanta.push_back(C.rows());

    #pragma omp parallel for
    for (unsigned int i = 0; i < BasisSet.size(); i++)
    {
        std::vector<int> ModeOccI;
        for (unsigned int m = 0; m < BasisSet[i].Modes.size(); m++) ModeOccI.push_back(BasisSet[i].Modes[m].Quanta);
        for (unsigned int j = i; j < BasisSet.size(); j++) // starts with j=i to exploit Hermiticity (Hij = Hji)
        {
            double Vij = 0;
            std::vector<int> ModeOccJ;
            for (unsigned int m = 0; m < BasisSet[j].Modes.size(); m++) ModeOccJ.push_back(BasisSet[j].Modes[m].Quanta);
            std::vector<int> DiffModes = CalcDiffModes(BasisSet[i], BasisSet[j]);
           
            // Harmonic Part
            if (DiffModes.size() < 2)
            {
                if (DiffModes.size() == 0)
                {
                    std::vector<double> Xs = ContractedHOTerms(Cs, Frequencies, ModeOccI, ModeOccJ);
                    double HOij = 0.0;
                    for (unsigned int m = 0; m < Xs.size(); m++)
                    {
                        HOij += Xs[m];
                    }
                    Vij += HOij;
                }
                if (DiffModes.size() == 1)
                {
                    Eigen::MatrixXd Cm = Cs[DiffModes[0]];
                    double Xm = 0.0;
                    for (unsigned int im = 0; im < Cs[DiffModes[0]].rows(); im++)
                    {
                        Xm += Cs[DiffModes[0]](im, ModeOccI[DiffModes[0]]) * Cs[DiffModes[0]](im, ModeOccJ[DiffModes[0]]) * (im + 0.5);
                    }
                    Xm *= Frequencies[DiffModes[0]];
                    Vij += Xm;
                }
            }

            // Anharmonic Part
            std::vector<std::vector<double>> Ys = ContractedAnharmonicPotential(Cs, GenericV, ModeOccI, ModeOccJ);
            for (FConst &FC : FCs)
            {
                double Vijq = 0.0;
                if (VectorContainedIn(DiffModes, FC.QUnique))
                {      
                    Vijq = FC.fc;
                    for (unsigned int m = 0; m < FC.QUnique.size(); m++)
                    {
                        Vijq *= Ys[FC.QUnique[m]][FC.QPowers[m]];
                    }
                }
                Vij += Vijq;
            }
            H(i, j) += Vij;
            if (i != j) H(j, i) += Vij;
        }
    }
    return H;
};

bool FitsMaxQuanta(WaveFunction &B, std::vector<int> &MaxQuanta)
{
    for (unsigned int m = 0; m < B.Modes.size(); m++)
    {
        if (B.Modes[m].Quanta >= MaxQuanta[m]) return false;
    }
    return true;
}

std::vector<WaveFunction> AddStatesHBWithMax(std::vector<WaveFunction> &BasisSet, std::vector<FConst> &AnharmHB, Eigen::Ref<Eigen::VectorXd> C, double eps, std::vector<int> &MaxQuanta){ // Expand basis via Heat Bath algorithm
    HashedStates HashedBasisInit; // hashed unordered_set containing BasisSet to check for duplicates
    HashedStates HashedNewStates; // hashed unordered_set of new states that only allows unique states to be inserted
    for( WaveFunction& wfn : BasisSet){
        HashedBasisInit.insert(wfn); // Populate hashed unordered_set with initial basis states
    }

    std::vector<double> CVec;
    for (unsigned int n = 0; n < C.rows(); n++) CVec.push_back(abs(C[n]));
    std::vector<long unsigned int> CSortedInd = SortIndices(CVec);
    for(unsigned int i=0; i<AnharmHB.size(); ++i){ // Loop over sorted force constants
        if (abs(AnharmHB[i].fc * C[CSortedInd[CSortedInd.size() - 1]]) < eps) break; // means that the largest Cn doesn't meet the criteria so we are done
        for (int nn = CSortedInd.size() - 1; nn >= 0; nn--)
        {
            unsigned int n = CSortedInd[nn];
            double Cn = C[n];
            if(abs(Cn*AnharmHB[i].fc) >= eps){ // States connected by fc will be added if |fc*Cn| >= eps
                if(AnharmHB[i].QUnique.size()==1){ // Stupid way to enumerate new states from F_ij
                    for( int a=-(AnharmHB[i].QPowers[0]); a<(AnharmHB[i].QPowers[0]+1); a+=2 ){
                        if( a != 0){// Skip if no change
                            WaveFunction tmp = BasisSet[n];
                            tmp.Modes[AnharmHB[i].QUnique[0]].Quanta += a;
                            if( tmp.Modes[AnharmHB[i].QUnique[0]].Quanta >=0 &&
                                                                HashedBasisInit.count(tmp) == 0
                                    ){ //make sure a|0> = 0 and tmp does not exist in original basis
                                if (FitsMaxQuanta(tmp, MaxQuanta))
                                {
                                    HashedNewStates.insert(tmp); // add new state to set
                                }
                            }
                        }
                    }
                }
                if(AnharmHB[i].QUnique.size()==2){ // Stupid way to enumerate new states from F_ij
                    for( int a=-(AnharmHB[i].QPowers[0]); a<(AnharmHB[i].QPowers[0]+1); a+=2 ){
                        for( int b=-(AnharmHB[i].QPowers[1]); b<(AnharmHB[i].QPowers[1]+1); b+=2 ){
                            if(abs(a)+abs(b) != 0){// Skip if no change
                                WaveFunction tmp = BasisSet[n];
                                tmp.Modes[AnharmHB[i].QUnique[0]].Quanta += a;
                                tmp.Modes[AnharmHB[i].QUnique[1]].Quanta += b;
                                if( tmp.Modes[AnharmHB[i].QUnique[0]].Quanta >=0 &&
                                       tmp.Modes[AnharmHB[i].QUnique[1]].Quanta >=0 && 
                                                                       HashedBasisInit.count(tmp) == 0
                                       ){ //make sure a|0> = 0
                                        if (FitsMaxQuanta(tmp, MaxQuanta))
                                        {
                                            HashedNewStates.insert(tmp); // add new state to set
                                        }
                                }
                            }
                        }
                    }
                }
                if(AnharmHB[i].QUnique.size()==3){ // Stupid way to enumerate new states from F_ijk
                    for( int a=-(AnharmHB[i].QPowers[0]); a<(AnharmHB[i].QPowers[0]+1); a+=2 ){
                        for( int b=-(AnharmHB[i].QPowers[1]); b<(AnharmHB[i].QPowers[1]+1); b+=2 ){
                            for( int c=-(AnharmHB[i].QPowers[2]); c<(AnharmHB[i].QPowers[2]+1); c+=2 ){
                                if(abs(a)+abs(b)+abs(c) != 0){// Skip if no change
                                    WaveFunction tmp = BasisSet[n];
                                    tmp.Modes[AnharmHB[i].QUnique[0]].Quanta += a;
                                    tmp.Modes[AnharmHB[i].QUnique[1]].Quanta += b;
                                    tmp.Modes[AnharmHB[i].QUnique[2]].Quanta += c;
                                    if( tmp.Modes[AnharmHB[i].QUnique[0]].Quanta >=0 &&
                                           tmp.Modes[AnharmHB[i].QUnique[1]].Quanta >=0  &&
                                           tmp.Modes[AnharmHB[i].QUnique[2]].Quanta >=0 &&
                                           HashedBasisInit.count(tmp) == 0
                                           ){ //make sure a|0> = 0
                                        if (FitsMaxQuanta(tmp, MaxQuanta))
                                        {
                                            HashedNewStates.insert(tmp); // add new state
                                        }
                                    }
                                }      
                            }
                        }
                    }
                }
                if(AnharmHB[i].QUnique.size()==4){ // Stupid way to enumerate new states from F_ijkl
                    for( int a=-(AnharmHB[i].QPowers[0]); a<(AnharmHB[i].QPowers[0]+1); a+=2 ){
                        for( int b=-(AnharmHB[i].QPowers[1]); b<(AnharmHB[i].QPowers[1]+1); b+=2 ){
                            for( int c=-(AnharmHB[i].QPowers[2]); c<(AnharmHB[i].QPowers[2]+1); c+=2 ){
                                for( int d=-(AnharmHB[i].QPowers[3]); d<(AnharmHB[i].QPowers[3]+1); d+=2 ){
                                    if(abs(a)+abs(b)+abs(c)+abs(d) != 0){ // Skip if no change
                                        WaveFunction tmp = BasisSet[n];
                                        tmp.Modes[AnharmHB[i].QUnique[0]].Quanta += a;
                                        tmp.Modes[AnharmHB[i].QUnique[1]].Quanta += b;
                                        tmp.Modes[AnharmHB[i].QUnique[2]].Quanta += c;
                                        tmp.Modes[AnharmHB[i].QUnique[3]].Quanta += d;
                                        if( tmp.Modes[AnharmHB[i].QUnique[0]].Quanta >=0 &&
                                               tmp.Modes[AnharmHB[i].QUnique[1]].Quanta >=0  &&
                                               tmp.Modes[AnharmHB[i].QUnique[2]].Quanta >=0 &&
                                               tmp.Modes[AnharmHB[i].QUnique[3]].Quanta >=0 &&
                                               HashedBasisInit.count(tmp) == 0
                                               ){ //make sure a|0> = 0
                                                if (FitsMaxQuanta(tmp, MaxQuanta))
                                                {
                                                    HashedNewStates.insert(tmp); // add new state
                                                }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }   
                if(AnharmHB[i].QUnique.size()==5){ // Stupid way to enumerate new states from F_ijklm 
                    for( int a=-(AnharmHB[i].QPowers[0]); a<(AnharmHB[i].QPowers[0]+1); a+=2 ){
                        for( int b=-(AnharmHB[i].QPowers[1]); b<(AnharmHB[i].QPowers[1]+1); b+=2 ){
                            for( int c=-(AnharmHB[i].QPowers[2]); c<(AnharmHB[i].QPowers[2]+1); c+=2 ){
                                for( int d=-(AnharmHB[i].QPowers[3]); d<(AnharmHB[i].QPowers[3]+1); d+=2 ){
                                    for( int f=-(AnharmHB[i].QPowers[4]); f<(AnharmHB[i].QPowers[4]+1); f+=2 ){
                                        if(abs(a)+abs(b)+abs(c)+abs(d)+abs(f) != 0){ // Skip if no change
                                            WaveFunction tmp = BasisSet[n];
                                            tmp.Modes[AnharmHB[i].QUnique[0]].Quanta += a;
                                            tmp.Modes[AnharmHB[i].QUnique[1]].Quanta += b;
                                            tmp.Modes[AnharmHB[i].QUnique[2]].Quanta += c;
                                            tmp.Modes[AnharmHB[i].QUnique[3]].Quanta += d;
                                            tmp.Modes[AnharmHB[i].QUnique[4]].Quanta += f;
                                            if( tmp.Modes[AnharmHB[i].QUnique[0]].Quanta >=0 &&
                                                   tmp.Modes[AnharmHB[i].QUnique[1]].Quanta >=0  &&
                                                   tmp.Modes[AnharmHB[i].QUnique[2]].Quanta >=0 &&
                                                   tmp.Modes[AnharmHB[i].QUnique[3]].Quanta >=0 &&
                                                   tmp.Modes[AnharmHB[i].QUnique[4]].Quanta >=0 && 
                                                   HashedBasisInit.count(tmp) == 0
                                                   ){ //make sure a|0> = 0
                                                    if (FitsMaxQuanta(tmp, MaxQuanta))
                                                    {
                                                        HashedNewStates.insert(tmp); // add new state
                                                    }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                if(AnharmHB[i].QUnique.size()==6){ // Stupid way to enumerate new states from F_ijklmn 
                    for( int a=-(AnharmHB[i].QPowers[0]); a<(AnharmHB[i].QPowers[0]+1); a+=2 ){
                        for( int b=-(AnharmHB[i].QPowers[1]); b<(AnharmHB[i].QPowers[1]+1); b+=2 ){
                            for( int c=-(AnharmHB[i].QPowers[2]); c<(AnharmHB[i].QPowers[2]+1); c+=2 ){
                                for( int d=-(AnharmHB[i].QPowers[3]); d<(AnharmHB[i].QPowers[3]+1); d+=2 ){
                                    for( int f=-(AnharmHB[i].QPowers[4]); f<(AnharmHB[i].QPowers[4]+1); f+=2 ){
                                        for( int g=-(AnharmHB[i].QPowers[5]); g<(AnharmHB[i].QPowers[5]+1); g+=2 ){
                                            if(abs(a)+abs(b)+abs(c)+abs(d)+abs(f)+abs(g) != 0){ // Skip if no change
                                                WaveFunction tmp = BasisSet[n];
                                                tmp.Modes[AnharmHB[i].QUnique[0]].Quanta += a;
                                                tmp.Modes[AnharmHB[i].QUnique[1]].Quanta += b;
                                                tmp.Modes[AnharmHB[i].QUnique[2]].Quanta += c;
                                                tmp.Modes[AnharmHB[i].QUnique[3]].Quanta += d;
                                                tmp.Modes[AnharmHB[i].QUnique[4]].Quanta += f;
                                                tmp.Modes[AnharmHB[i].QUnique[5]].Quanta += g;
                                                if( tmp.Modes[AnharmHB[i].QUnique[0]].Quanta >=0 &&
                                                       tmp.Modes[AnharmHB[i].QUnique[1]].Quanta >=0 &&
                                                       tmp.Modes[AnharmHB[i].QUnique[2]].Quanta >=0 &&
                                                       tmp.Modes[AnharmHB[i].QUnique[3]].Quanta >=0 &&
                                                       tmp.Modes[AnharmHB[i].QUnique[4]].Quanta >=0 &&
                                                       tmp.Modes[AnharmHB[i].QUnique[5]].Quanta >=0 &&
                                                                                                       HashedBasisInit.count(tmp) == 0
                                                       ){ //make sure a|0> = 0
                                                        if (FitsMaxQuanta(tmp, MaxQuanta))
                                                        {
                                                            HashedNewStates.insert(tmp); // add new state
                                                        }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                if(AnharmHB[i].fcpow.size() == 0 || AnharmHB[i].fcpow.size() > 6 ){
                    //Print error message
                    cout << "Error: HCI only works with force constants up to 6th order." << endl;
                    cout.flush();
                    exit(0);
                }
            }
            else{break;}// break loop if you've reached element < eps, since all future elements will be smaller (HB sorting)
        }
    }
    std::vector<WaveFunction> NewBasis;
    for (const WaveFunction &WF : HashedNewStates) NewBasis.push_back(WF);
    return NewBasis;
}

SpMat VCISparseHamFromVSCF(std::vector<WaveFunction> &BasisSet, std::vector<double> &Frequencies, std::vector<FConst> &FCs, std::vector<Eigen::MatrixXd> &Cs, std::vector<Eigen::SparseMatrix<double>> &GenericV)
{
    SpMat H(BasisSet.size(), BasisSet.size());
    std::vector<Trip> HTrip;

    std::vector<int> MaxQuanta;
    for (Eigen::MatrixXd &C : Cs) MaxQuanta.push_back(C.rows());

    #pragma omp parallel for
    for (unsigned int i = 0; i < BasisSet.size(); i++)
    {
        std::vector<int> ModeOccI;
        for (unsigned int m = 0; m < BasisSet[i].Modes.size(); m++) ModeOccI.push_back(BasisSet[i].Modes[m].Quanta);
        for (unsigned int j = i; j < BasisSet.size(); j++) // starts with j=i to exploit Hermiticity (Hij = Hji)
        {
            double Vij = 0;
            std::vector<int> ModeOccJ;
            for (unsigned int m = 0; m < BasisSet[j].Modes.size(); m++) ModeOccJ.push_back(BasisSet[j].Modes[m].Quanta);
            std::vector<int> DiffModes = CalcDiffModes(BasisSet[i], BasisSet[j]);
           
            // Harmonic Part
            if (DiffModes.size() < 2)
            {
                if (DiffModes.size() == 0)
                {
                    std::vector<double> Xs = ContractedHOTerms(Cs, Frequencies, ModeOccI, ModeOccJ);
                    double HOij = 0.0;
                    for (unsigned int m = 0; m < Xs.size(); m++)
                    {
                        HOij += Xs[m];
                    }
                    Vij += HOij;
                }
                if (DiffModes.size() == 1)
                {
                    Eigen::MatrixXd Cm = Cs[DiffModes[0]];
                    double Xm = 0.0;
                    for (unsigned int im = 0; im < Cs[DiffModes[0]].rows(); im++)
                    {
                        Xm += Cs[DiffModes[0]](im, ModeOccI[DiffModes[0]]) * Cs[DiffModes[0]](im, ModeOccJ[DiffModes[0]]) * (im + 0.5);
                    }
                    Xm *= Frequencies[DiffModes[0]];
                    Vij += Xm;
                }
            }

            // Anharmonic Part
            std::vector<std::vector<double>> Ys = ContractedAnharmonicPotential(Cs, GenericV, ModeOccI, ModeOccJ);
            for (FConst &FC : FCs)
            {
                double Vijq = 0.0;
                if (VectorContainedIn(DiffModes, FC.QUnique))
                {      
                    Vijq = FC.fc;
                    for (unsigned int m = 0; m < FC.QUnique.size(); m++)
                    {
                        Vijq *= Ys[FC.QUnique[m]][FC.QPowers[m]];
                    }
                }
                Vij += Vijq;
            }
            #pragma omp critical
            if (abs(Vij) > 1e-12)
            {
                if (i == j) HTrip.push_back(Trip(i, j, Vij / 2.0));
                else HTrip.push_back(Trip(i, j, Vij));
            }
        }
    }
    
    H.setFromTriplets(HTrip.begin(), HTrip.end());
    H.makeCompressed();
    HTrip = std::vector<Trip>(); // Free memory
    SpMat HT = H.transpose();
    HT.makeCompressed();
    H += HT; // Complete symmetric matrix
    HT = SpMat(1,1); // Free memory
    H.makeCompressed();
    cout << "The Hamiltonian is " << fixed << setprecision(2) << 
        100*(1.-(double)H.nonZeros()/(double)H.size()) << "% sparse." << endl;
    return H;
};


