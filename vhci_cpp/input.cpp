#include "vhci.hpp"
using namespace std;

void VHCI::ReadArgs(int argc, char* argv[], fstream& vcidata, fstream& freqfile, string& cpfilename)
{
    //Function to read the command line arguments
    cpfilename = "none";
    bool DoQuit = 0; //Exit if an error is detected
    string dummy; //Generic string
    stringstream call; //Stream for system calls and reading/writing files
    //Read command line arguments
    if (argc == 1)
    {
        //Escape if there are no arguments
        cout << '\n';
        cout << "Missing arguments...";
        cout << '\n' << '\n';
        cout << "Usage: vhci -n Ncpus -i Modes.inp -o Freqs.dat -c Checkpoint.txt";
        cout << '\n' << '\n';
        cout << "Use -h or --help for detailed instructions.";
        cout << '\n' << '\n';
        cout.flush();
        exit(0);
    }
    if ((argc % 2) != 1)
    {
        dummy = string(argv[1]);
        if ((dummy != "-h") and (dummy != "--help"))
        {
            //Escape if there are missing arguments
            cout << '\n';
            cout << "Wrong number of arguments...";
            cout << '\n' << '\n';
            cout << "Usage: vhci -n Ncpus -i Modes.inp -o Freqs.dat -c Checkpoint.txt";
            cout << '\n' << '\n';
            cout << "Use -h or --help for detailed instructions.";
            cout << '\n' << '\n';
            cout.flush();
            exit(0);
        }
    }
    for (int i=0;i<argc;i++)
    {
        //Read file names and CPUs
        dummy = string(argv[i]);
        if ((dummy == "-h") or (dummy == "--help"))
        {
            //Print helpful information and exit
            cout << '\n';
            cout << "Usage: vhci -n Ncpus -i Modes.inp -o Freqs.txt";
            cout << '\n' << '\n';
            cout << "Command line arguments:";
            cout << '\n' << '\n';
            cout << "  -n    Number of CPUs to use for the calculation.";
            cout << '\n' << '\n';
            cout << "  -i    Input file.";
            cout << '\n' << '\n';
            cout << "  -o    Output file for the first NEig VCI frequencies.";
            cout << '\n' << '\n';
            cout << "  -c    Checkpoint file for HCI iterations on large systems.";
            cout << '\n' << '\n';
            cout.flush();
            exit(0);
        }
        if (dummy == "-n")
        {
            Ncpus = atoi(argv[i+1]);
        }
        if (dummy == "-i")
        {
            InputFile.open(argv[i+1],ios_base::in);
        }
        if (dummy == "-o")
        {
            OutputFile.open(argv[i+1],ios_base::out);
        }
        if (dummy == "-c"){
            cpfilename = argv[i+1];
        }
    }
    //Check for argument errors
    if (Ncpus < 1)
    {
        //Checks the number of threads and continue
        cout << " Warning: Calculations cannot run with ";
        cout << Ncpus << " CPUs.";
        cout << '\n';
        cout << "  Do you know how computers work?";
        cout << " Ncpus set to 1";
        cout << '\n';
        Ncpus = 1;
        cout.flush(); //Print warning
    }
    if (Ncpus > FindMaxThreads())
    {
        cout << " Warning: Too many threads requested.";
        cout << '\n';
        cout << "  Requested: " << Ncpus << ",";
        cout << " Available: " << FindMaxThreads();
        cout << '\n';
        cout << "  Ncpus set to " << FindMaxThreads();
        cout << '\n';
        Ncpus = FindMaxThreads();
        cout.flush(); //Print warning
    }
    if (!InputFile.good())
    {
        //Check input file
        cout << " Error: Could not open input file.";
        cout << '\n';
        DoQuit = 1;
    }
    if (!OutputFile.good())
    {
        //Check output file
        cout << " Error: Could not output file.";
        cout << '\n';
        DoQuit = 1;
    }

    if (DoQuit)
    {
        //Quit if there is an error
        cout << '\n';
        cout.flush();
        exit(0);
    }
    else
    {
        //Sarcastically continue
        cout << '\n';
        cout << "No fatal errors detected.";
        cout << '\n';
        cout << " And there was much rejoicing. Yay...";
        cout << '\n' << '\n';
        cout.flush();
    }
    if(cpfilename=="none"){
        cout << "No checkpoint filename specified. Checkpoints will not be saved." << '\n' << '\n'; 
    }
    else{
        cout << "Locating checkpoint file..." << '\n';
        IsRestart(cpfilename);
        if(restart==1){
            cout << "Checkpoint file found. Calculation will be resumed." << '\n' << '\n';
        }else{
            cout << "No checkpoint file found. Calculation will start from the beginning." << '\n' << '\n';
        }
        cout.flush();
    }
    //Set threads
    omp_set_num_threads(Ncpus);
    Eigen::setNbThreads(Ncpus);
    return;
};

void VHCI::ReadInput(fstream& vcidata)
{
    //Function to read the input files
    string dummy; //Generic sting
    //Count basis functions and read modes
    int Nmodes = 0; //Number of different modes
    int Ntot = 0; //Maximum quanta in a single product state
    vector<HOFunc> BasisCount; //Temp. storage of modes
    int Nfc = 0; //Number of force constants
    //Heat bath parameters 
    getline(vcidata,dummy); //Junk (String with comment for input file)
    vcidata >> dummy; //Junk
    vcidata >> Epsilon1; // Variational HCI cutoff energy
    vcidata >> dummy;
    vcidata >> NStates; // Number of eigenstates to include in Heat Bath optimization
    // Set PT2 parameters
    vcidata >> dummy;
    vcidata >> perturb;
    vcidata >> dummy;
    vcidata >> Epsilon2;
    vcidata >> dummy;
    vcidata >> Epsilon3;
    vcidata >> NWalkers;
    if (perturb != 0 && perturb != 1)
    {
        //Print an error message
        cout << "Error: perturb must be either 0 or 1 ";
        cout << '\n' << '\n';
        cout.flush(); //Print message
        //Quit
        exit(0);
    }
    if(perturb==0){
        Epsilon2=0;
    }
    //Set maximum simultaneous excitations 
    vcidata >> dummy; //Junk
    vcidata >> MaxTotalQuanta; //Maximum quanta in a single product state
    //Read active modes
    vcidata >> dummy; //Clear junk
    vcidata >> NModes; //Read modes
    for (int i=0;i<Nmodes;i++)
    {
        //Actual vibrational modes
        double Freq;
        int QMax;
        int modeid;
        vcidata >> modeid; //Read mode ID
        //Check mode order
        if (modeid != i)
        {
            //Print an error message
            cout << "Error: Expected mode " << i;
            cout << " but read data for mode " << modeid;
            cout << '\n' << '\n';
            cout.flush(); //Print message
            //Quit
            exit(0);
        }
        vcidata >> Freq; //Frequency
        vcidata >> QMax; //Max number of quanta
        Frequencies.push_back(Freq);
        MaxQuanta.push_back(QMax);
    }
    //Read anharmonic force constants
    vcidata >> dummy; //Clear junk
    int Nfc;
    vcidata >> Nfc; //Read number of anharmonic force constants
    int oldfcpower = 0;
    std::vector<VDeriv> PotentialByOrder;
    for (int i=0;i<Nfc;i++)
    {
        //Save force constant data
        std::vector<int> Qs;
        int fcpower = 0;
        vcidata >> fcpower;
        if (oldfcpower == 0) oldfcpower = fcpower;
        for (int j=0;j<fcpower;j++)
        {
            int modej = 0;
            vcidata >> modej;
            Qs.push_back(modej);
        }
        double fc;
        vcidata >> fc; //Read force constant value
        VTerm = VDeriv(fc, Qs, true);
        if (oldfcpower != fcpower)
        {
            Potential.push_back(PotentialByOrder);
            PotentialByOrder.clear();
        }
        PotentialByOrder.push_back(VTerm);
        oldfcpower = fcpower;
    }
    Potential.push_back(PotentialByOrder);

    //Print settings
    cout << "General settings:" << '\n';
    cout << "  CPU threads: " << NCPUs << '\n' << '\n'; 
    cout << "  Normal modes: " << Frequencies.size() << '\n';
    cout << "  Max quanta per state: " << MaxTotalQuanta << '\n';
    if(Epsilon1 > 1e-12){
        cout.precision(3);
        cout << "  Using HCI with energy cutoff: " << Epsilon1 << '\n';
        if(NStates == 1){
            cout << "  Optimizing just the ground states (ZPE)." << '\n';
        }
        else{
            cout << "  Optimizing the first " << NStates << " eigenstates." << '\n';
        }
    }
    cout << '\n';
    cout.precision(12); //Replace settings
    cout.flush();
    return;
};


