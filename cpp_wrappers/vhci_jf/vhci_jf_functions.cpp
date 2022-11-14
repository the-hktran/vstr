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

std::vector<Eigen::MatrixXd> GetVEffCPP(std::vector<Eigen::SparseMatrix<double>> &AnharmTensor, std::vector<WaveFunction> Basis, std::vector<Eigen::MatrixXd> &CByModes, std::vector<std::vector<std::vector<int>>> &ModalSlices, std::vector<int> &MaxQuanta, std::vector<int> &ModeOcc, bool FirstV = false)
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
