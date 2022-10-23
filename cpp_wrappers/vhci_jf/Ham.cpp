/*

#############################################################################
#                                                                           #
#              Vibrational Heat-Bath Configuration Interaction              # 
#                         By: Jonathan H. Fetherolf                         #
#                                                                           #
#                      Based on LOVCI by Eric G. Kratz                      #
#                                                                           #
#############################################################################

Functions for building and diagonalizing the Hamiltonian matrix

*/
#include "VCI_headers.h"

double AnharmPot(WaveFunction &Bn, WaveFunction &Bm, const FConst& fc)
{
    //Calculate anharmonic matrix elements for <m|H|n>
    double Vnm = 0;
    // Initialize new states
    vector<vector<int> > NewStates(fc.QPowers.size());
    vector<vector<double> > StateCoeffs(fc.QPowers.size());
    //Create new states
    for (unsigned int i=0;i<fc.QPowers.size();i++)
    {
        //Apply operators
        NewStates[i].push_back(Bn.Modes[fc.QPowers[i]].Quanta);
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
    for (unsigned int i=0;i<fc.QPowers.size();i++)
    {
        CoeffSum.push_back(0.0);
        for (unsigned int j=0;j<NewStates[i].size();j++)
        {
            int quantn = NewStates[i][j];
            int quantm = Bm.Modes[fc.QPowers[i]].Quanta;
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
        if(AnharmFC[k].QPowers.size()>fcmax){
            fcmax = AnharmFC[k].QPowers.size();
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
        if(AnharmFC[k].QPowers.size()>fcmax)
        {
            fcmax = AnharmFC[k].QPowers.size();
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
                if(i==j){
                    HTrip.push_back(Trip(i,j,Vij/2.)); // Dividing by 2 so I can do H=H+H*

                }else{
                    HTrip.push_back(Trip(i,j,Vij)); // Exploiting Hermiticity (Hij=Hji)
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
                if(i==j){
                    HTrip.push_back(Trip(i,j,Vij/2.)); // Dividing by 2 so I can do H=H+H*

                }else{
                    HTrip.push_back(Trip(i,j,Vij)); // Exploiting Hermiticity (Hij=Hji)
                }
            }
        }
    }
    return;
};

inline void MakeHamSparse(SpMat &HSp, std::vector<WaveFunction> &BasisSet, std::vector<double> &Frequencies, std::vector<FConst> &AnharmFC, std::vector<FConst> &CubicFC, std::vector<FConst> &QuarticFC, std::vector<FConst> &QuinticFC, std::vector<FConst> &SexticFC)
{
    //Build the sparse CI Hamiltonian
    vector< Trip > HTrip;
    ZerothHamSparse(HTrip, BasisSet, Frequencies);
    AnharmHamSparse(HTrip, BasisSet, Frequencies, AnharmFC, CubicFC, QuarticFC, QuinticFC, SexticFC);
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
    MakeHamSparse(H, BasisSet, Frequencies, AnharmFC, CubicFC, QuarticFC, QuinticFC, SexticFC);
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
    SymEigsSolver< double, SMALLEST_ALGE, SparseMVProd > eigs(&op, NEig, NCV);
    // Initialize and compute
    eigs.init();
    int nconv = eigs.compute(1000,1e-10,SMALLEST_ALGE);
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
