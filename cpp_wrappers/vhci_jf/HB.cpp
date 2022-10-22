/*

#############################################################################
#                                                                           #
#              Vibrational Heat-Bath Configuration Interaction              # 
#                         By: Jonathan H. Fetherolf                         #
#                                                                           #
#                      Based on LOVCI by Eric G. Kratz                      #
#                                                                           #
#############################################################################

Functions for generating the heat bath product basis

*/
#include "VCI_headers.h"

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


void HeatBath_Sort_FC(std::vector<FConst> &AnharmHB){
    sort(AnharmHB.begin(), AnharmHB.end(), sortByFC); // Sort force constants from large to small magnitude
}

std::vector<WaveFunction> AddStatesHB(std::vector<WaveFunction> &BasisSet, std::vector<FConst> &AnharmHB, Eigen::VectorXd &C, double eps){ // Expand basis via Heat Bath algorithm
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
