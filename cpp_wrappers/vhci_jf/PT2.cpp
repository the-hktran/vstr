/*

#############################################################################
#                                                                           #
#              Vibrational Heat-Bath Configuration Interaction              # 
#                         By: Jonathan H. Fetherolf                         #
#                                                                           #
#                      Based on LOVCI by Eric G. Kratz                      #
#                                                                           #
#############################################################################

 Implementation of Epstein-Nesbet PT2 with Heat Bath sorting

*/
#include "VCI_headers.h"

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
    cout << " Finding connected states..." << endl;
    int N_opt;
    if(NEig > BasisSet.size()){ // If we don't have enough states to optimize for yet
        N_opt = BasisSet.size();
    }else{ // If number of states exceeds the number selected to optimize
        N_opt = NEig;
    }

    if(NEig==1){
        cout << " Calculating the 2nd-order perturbative correction on the ground state energy" << endl;
    }else{
        cout << " Calculating the 2nd-order perturbative correction for on first " << N_opt << " eigenvalues." << endl;
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

std::tuple<std::vector<double>, std::vector<double>> DoSPT2(MatrixXd& Evecs, VectorXd& Evals, std::vector<WaveFunction> &BasisSet, std::vector<WaveFunction> &PTBasisSet, std::vector<FConst> &AnharmFC, std::vector<FConst> &CubicFC, std::vector<FConst> &QuarticFC, std::vector<FConst> &QuinticFC, std::vector<FConst> &SexticFC, double PT2_Eps, int NEig, int Nd, int Ns, bool SemiStochastic = false, double PT2_Eps2 = 0.0)
{
    cout << " Starting Stocastic PT2 corrections." << endl;
    int N_opt;
    if(NEig > BasisSet.size()){ // If we don't have enough states to optimize for yet
        N_opt = BasisSet.size();
    }else{ // If number of states exceeds the number selected to optimize
        N_opt = NEig;
    }

    if(NEig==1){
        cout << " Calculating the 2nd-order perturbative correction on the ground state energy" << endl;
    }else{
        cout << " Calculating the 2nd-order perturbative correction for on first " << N_opt << " eigenvalues." << endl;
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
        if (SemiStochastic)
        {
            std::vector<WaveFunction> VarDetBasisSet;
            std::vector<WaveFunction> DetPTBasisSet = AddStatesHB(BasisSet, AnharmHB, Evecs.col(n), PT2_Eps);
            for (const WaveFunction &WF : BasisSet) VarDetBasisSet.push_back(WF);
            for (const WaveFunction &WF : DetPTBasisSet) VarDetBasisSet.push_back(WF);
            std::vector<WaveFunction> PTBasisSet = AddStatesHB(DetPTBasisSet, AnharmHB, Evecs.col(n), PT2_Eps2);
            VarDetBasisSet = HashedStates(); // Clear memory, this entity doesn't really need to exist if we are careful with slicing, but I don't feel like doing that now.
            // If I were feeling like it, I would store the original basis set size and just use the variational and deterministic parts together, like how JF used to do it, 
            // and reset after each iteration.
        }
        else
        {
            std::vector<WaveFunction> PTBasisSet = AddStatesHB(BasisSet, AnharmHB, Evecs.col(n), PT2_Eps);
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
        } // end variational calculation for semi-stochastic

        std::vector<double> DeltaEs(Ns, 0.0);
        #pragma omp parallel for
        for (unsigned int s = 0; s < Ns; s++)
        {
            std::vector<double> WalkerProbability = StateProbability(Cn);
            std::map<int, int> WalkerPopulation;
            FillWalkers(WalkerPopulation, WalkerProbability, Nd);

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
    }

    std::vector<double> SigmaDeltaE(N_opt, 0.0);
    for (unsigned int n = 0; n < N_opt; n++)
    {
        for (unsigned int s = 0; s < Ns; s++)
        {
            DeltaE[n] += DeltaESample[n][s];
            SigmaDeltaE[n] += pow(DeltaESample[n][s], 2);
        }
        SigmaDeltaE[n] -= (pow(DeltaE[n], 2) / Ns);
        DeltaE[n] /= Ns;
        if (SemiStochastic) DeltaE[n] += DeltaEDet[n];
        SigmaDeltaE[n] /= (Ns - 1);
        SigmaDeltaE[n] = sqrt(SigmaDeltaE[n]);
        cout << DeltaE[n] << " " << SigmaDeltaE[n] << std::endl;
    }

    return std::make_tuple<std::vector<double, std::vector<double>>(DeltaE, SigmaDeltaE);
}
