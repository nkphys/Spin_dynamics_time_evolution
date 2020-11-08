#include "ParametersEngine_MultiOrbSF.h"

void Parameters_MultiOrbSF::Initialize(string inputfile_)
{


    maxmoment = 10.0;

    cout << "____________________________________" << endl;
    cout << "Reading the inputfile: " << inputfile_ << endl;
    cout << "____________________________________" << endl;


    lx = int(matchstring(inputfile_, "Xsite"));
    ly = int(matchstring(inputfile_, "Ysite"));

    TBC_mx = int(matchstring(inputfile_, "TwistedBoundaryCond_mx"));
    n_orbs = int(matchstring(inputfile_, "N_Orbs"));
    J_Hund.resize(n_orbs);
    OnSiteE.resize(n_orbs);

    TBC_my = int(matchstring(inputfile_, "TwistedBoundaryCond_my"));
    TBC_cellsX = int(matchstring(inputfile_, "TBC_cellsX"));
    TBC_cellsY = int(matchstring(inputfile_, "TBC_cellsY"));
    fix_mu = matchstring(inputfile_, "Fix_mu");
    fixed_mu_value = double(matchstring(inputfile_, "fixed_mu_value")) * 1.0;
    BoundaryConnection = double(matchstring(inputfile_, "PBC"));


    ns = lx * ly;
    cout << "TotalNumberOfSites = " << ns << endl;

    Fill = matchstring(inputfile_, "Filling_per_site");
    Total_Particles = ns * Fill;
    cout << "TotalNumberOfParticles = " << Total_Particles << endl;


    MCNorm = 0.0; //matchstring(inputfile,"MCNorm")
    RandomSeed = matchstring(inputfile_, "RandomSeed");
    RandomDisorderSeed = matchstring(inputfile_, "RandomDisorderSeed");
    Disorder_Strength = matchstring(inputfile_, "Disorder_Strength");
    Boltzman_constant = matchstring(inputfile_, "Boltzman_constant");


    string J_Hund_str = "J_Hund";
    string OnSiteE_str = "OnSiteE";
    string J_Hund_str_temp, OnSiteE_str_temp;

    for(int n=0;n<n_orbs;n++){
        J_Hund_str_temp = J_Hund_str + to_string(n);
        OnSiteE_str_temp = OnSiteE_str + to_string(n);

        J_Hund[n] = matchstring(inputfile_, J_Hund_str_temp);
        OnSiteE[n] = matchstring(inputfile_, OnSiteE_str_temp);
    }

    //Hopping matrices -------------------


    hopping_NN_X.resize(n_orbs,n_orbs);
    hopping_NN_Y.resize(n_orbs,n_orbs);
    string Nearest_Neigh_Hopping_basis_X;
    string Nearest_Neigh_Hopping_basis_Y;

    string NN_X_str, NN_Y_str;
    for (int m=0;m<n_orbs;m++){

        NN_X_str = "Nearest_Neigh_Hopping_X_basis_row" + to_string(m);
        NN_Y_str = "Nearest_Neigh_Hopping_Y_basis_row" + to_string(m);
        Nearest_Neigh_Hopping_basis_X=matchstring2(inputfile_, NN_X_str);
        Nearest_Neigh_Hopping_basis_Y=matchstring2(inputfile_, NN_Y_str);

        stringstream stream_X(Nearest_Neigh_Hopping_basis_X);
        stringstream stream_Y(Nearest_Neigh_Hopping_basis_Y);

        for(int n=0;n<n_orbs;n++){
            stream_X >> hopping_NN_X(m,n);
            stream_Y >> hopping_NN_Y(m,n);

        }
    }



    //Next Nearest hopping------------
    hopping_NNN_PXPY.resize(n_orbs,n_orbs);
    hopping_NNN_PXMY.resize(n_orbs,n_orbs);
    //If needed read from input file

    //Hopping matrices done---------------




    lambda_lattice = matchstring (inputfile_, "lambda_lattice");
    K1x = matchstring(inputfile_, "K");
    K1y = K1x;
    cout << "K1x= " << K1x << endl;

    Dflag = 'N';

    Cooling_ = false;

    temp = double(matchstring(inputfile_, "Temperature")); // temperature in kelvin
    beta = double(1.0/(Boltzman_constant*temp));                         //Beta which is (T*k_b)^-1

    temp_min = temp;
    temp_max = temp;
    d_Temp = 10.0; //arbitrary positive number

    Temp_values.resize(1);
    Temp_values[0]=temp_min;

    pi = 4.00 * atan(double(1.0));
    Eav = 0.0;
    AccCount[0] = 0;
    AccCount[1] = 0;

    WindowSize = double(0.01);
    mus = 0.25;
    cout << "____________________________________" << endl;
}

double Parameters_MultiOrbSF::matchstring(string file, string match)
{
    string test;
    string line;
    ifstream readFile(file);
    double amount;
    bool pass = false;
    while (std::getline(readFile, line))
    {
        std::istringstream iss(line);
        if (std::getline(iss, test, '=') && pass == false)
        {
            // ---------------------------------
            if (iss >> amount && test == match)
            {
                // cout << amount << endl;
                pass = true;
            }
            else
            {
                pass = false;
            }
            // ---------------------------------
            if (pass)
                break;
        }
    }
    if (pass == false)
    {
        string errorout = match;
        errorout += "= argument is missing in the input file!";
        throw std::invalid_argument(errorout);
    }
    cout << match << " = " << amount << endl;
    return amount;
}

string Parameters_MultiOrbSF::matchstring2(string file, string match)
{

    string line;
    ifstream readFile(file);
    string amount;
    int offset;

    if (readFile.is_open())
    {
        while (!readFile.eof())
        {
            getline(readFile, line);

            if ((offset = line.find(match, 0)) != string::npos)
            {
                amount = line.substr(offset + match.length() + 1);
            }
        }
        readFile.close();
    }
    else
    {
        cout << "Unable to open input file while in the Parameters class." << endl;
    }

    cout << match << " = " << amount << endl;
    return amount;
}


