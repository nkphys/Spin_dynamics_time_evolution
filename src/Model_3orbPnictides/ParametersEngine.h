#ifndef Parameters_class
#define Parameters_class

class Parameters{

public:
    int lx, ly, ns, orbs, IterMax, MCNorm, RandomSeed;
    int TBC_mx, TBC_my;
    int TBC_cellsX, TBC_cellsY;
    int lx_cluster, ly_cluster;
    double mus, mus_Cluster, Fill,pi,J_NN,J_NNN,J_HUND,k_const,lamda_12,lamda_66;

    bool Cooling_;
    bool ED_;

    bool Metropolis_Algorithm;
    bool Heat_Bath_Algorithm;

    double temp_max, beta_min;
    double temp_min, beta_max;
    double d_Temp;

    int Last_n_sweeps_for_measurement;
    int Measurement_after_each_m_sweeps;

    double temp,beta,Eav,maxmoment;
    double WindowSize, AccCount[2];
    char Dflag;

    void Initialize(string inputfile_);
    double matchstring(string file,string match);

};


void Parameters::Initialize(string inputfile_){

    maxmoment=1.0;

    cout << "____________________________________" << endl;
    cout << "Reading the inputfile: " << inputfile_ << endl;
    cout << "____________________________________" << endl;
    lx = int(matchstring(inputfile_,"Xsite"));
    ly = int(matchstring(inputfile_,"Ysite"));


    ns = lx*ly;
    cout << "TotalNumberOfSites = "<< ns << endl;
    orbs = int(matchstring(inputfile_,"Orbitals"));
    Fill = matchstring(inputfile_,"Fill");
    cout << "TotalNumberOfParticles = "<< ns*Fill*orbs*2.0 << endl;


    MCNorm = 0.0; ; //matchstring(inputfile,"MCNorm")
    RandomSeed = matchstring(inputfile_,"RandomSeed");
    Dflag = 'N';
    J_NN=double(matchstring(inputfile_,"J_NN"));
    J_NNN=double(matchstring(inputfile_,"J_NNN"));
    J_HUND=double(matchstring(inputfile_,"J_HUND"));


        Cooling_=false;

        temp = double(matchstring(inputfile_,"Temperature"));   // temperature in kelvin
        beta=double(11604.0/ temp);    //Beta which is (T*k_b)^-1


    pi=4.00*atan(double(1.0));
    mus=0.25;
    cout << "____________________________________" << endl;
}


double Parameters::matchstring(string file,string match) {
    string test;
    string line;
    ifstream readFile(file);
    double amount;
    bool pass=false;
    while (std::getline(readFile, line)) {
        std::istringstream iss(line);
        if (std::getline(iss, test, '=') && pass==false) {
            // ---------------------------------
            if (iss >> amount && test==match) {
                // cout << amount << endl;
                pass=true;
            }
            else {
                pass=false;
            }
            // ---------------------------------
            if(pass) break;
        }
    }
    if (pass==false) {
        string errorout=match;
        errorout+="= argument is missing in the input file!";
        throw std::invalid_argument(errorout);
    }
    cout << match << " = " << amount << endl;
    return amount;
}

#endif



