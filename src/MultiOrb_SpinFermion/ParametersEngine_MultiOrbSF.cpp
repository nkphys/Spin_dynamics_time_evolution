#include "ParametersEngine_MultiOrbSF.h"

void Parameters_MultiOrbSF::Initialize(string inputfile_)
{


    maxmoment = 10.0;

    cout << "____________________________________" << endl;
    cout << "Reading the inputfile: '" << inputfile_ <<"'"<< endl;
    cout << "____________________________________" << endl;


    lx = int(matchstring(inputfile_, "Xsite"));
    ly = int(matchstring(inputfile_, "Ysite"));

    TBC_mx = int(matchstring(inputfile_, "TwistedBoundaryCond_mx"));
    n_orbs = int(matchstring(inputfile_, "N_Orbs"));
    n_Spins = int(matchstring(inputfile_, "No_of_classical_spins_per_site"));
    J_Hund.resize(n_orbs);
    OnSiteE.resize(n_orbs);

    TBC_my = int(matchstring(inputfile_, "TwistedBoundaryCond_my"));
    TBC_cellsX = int(matchstring(inputfile_, "TBC_cellsX"));
    TBC_cellsY = int(matchstring(inputfile_, "TBC_cellsY"));
    fix_mu = matchstring(inputfile_, "Fix_mu");
    fixed_mu_value = double(matchstring(inputfile_, "fixed_mu_value")) * 1.0;
    hz_mag=double(matchstring(inputfile_, "hz_mag")) * 1.0;
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
    hopping_0X_0Y.resize(n_orbs,n_orbs);
    hopping_1X_0Y.resize(n_orbs,n_orbs);
    hopping_0X_1Y.resize(n_orbs,n_orbs);
    hopping_m1X_1Y.resize(n_orbs,n_orbs);

    string Hopping_0X_0Y, Hopping_1X_0Y, Hopping_0X_1Y, Hopping_m1X_1Y;
    string Hop_0X_0Y_str, Hop_1X_0Y_str, Hop_0X_1Y_str, Hop_m1X_1Y_str;
    for (int m=0;m<n_orbs;m++){

        Hop_0X_0Y_str = "Hopping_0X_0Y_row" + to_string(m);
        Hop_1X_0Y_str = "Hopping_1X_0Y_row" + to_string(m);
        Hop_0X_1Y_str = "Hopping_0X_1Y_row" + to_string(m);
        Hop_m1X_1Y_str = "Hopping_m1X_1Y_row" + to_string(m);

        Hopping_0X_0Y=matchstring2(inputfile_, Hop_0X_0Y_str);
        Hopping_1X_0Y=matchstring2(inputfile_, Hop_1X_0Y_str);
        Hopping_0X_1Y=matchstring2(inputfile_, Hop_0X_1Y_str);
        Hopping_m1X_1Y=matchstring2(inputfile_, Hop_m1X_1Y_str);

        stringstream stream_Hopping_0X_0Y(Hopping_0X_0Y);
        stringstream stream_Hopping_1X_0Y(Hopping_1X_0Y);
        stringstream stream_Hopping_0X_1Y(Hopping_0X_1Y);
        stringstream stream_Hopping_m1X_1Y(Hopping_m1X_1Y);

        for(int n=0;n<n_orbs;n++){
            stream_Hopping_0X_0Y >> hopping_0X_0Y(m,n);
            stream_Hopping_1X_0Y >> hopping_1X_0Y(m,n);
            stream_Hopping_0X_1Y >> hopping_0X_1Y(m,n);
            stream_Hopping_m1X_1Y >> hopping_m1X_1Y(m,n);
        }
    }

    //Hopping matrices done---------------



    //Superexchange matrices -------------------
    K_0X_0Y.resize(n_Spins,n_Spins);
    K_1X_0Y.resize(n_Spins,n_Spins);
    K_0X_1Y.resize(n_Spins,n_Spins);
    K_m1X_1Y.resize(n_Spins,n_Spins);

    string K_0X_0Y_str, K_1X_0Y_str, K_0X_1Y_str, K_m1X_1Y_str;

    string KExc_0X_0Y_str, KExc_1X_0Y_str, KExc_0X_1Y_str, KExc_m1X_1Y_str ;
    for (int m=0;m<n_Spins;m++){

        KExc_0X_0Y_str = "K_0X_0Y_row" + to_string(m);
        KExc_1X_0Y_str = "K_1X_0Y_row" + to_string(m);
        KExc_0X_1Y_str = "K_0X_1Y_row" + to_string(m);
        KExc_m1X_1Y_str = "K_m1X_1Y_row" + to_string(m);

        K_0X_0Y_str=matchstring2(inputfile_, KExc_0X_0Y_str);
        K_1X_0Y_str=matchstring2(inputfile_, KExc_1X_0Y_str);
        K_0X_1Y_str=matchstring2(inputfile_, KExc_0X_1Y_str);
        K_m1X_1Y_str=matchstring2(inputfile_, KExc_m1X_1Y_str);

        stringstream stream_K_0X_0Y(K_0X_0Y_str);
        stringstream stream_K_1X_0Y(K_1X_0Y_str);
        stringstream stream_K_0X_1Y(K_0X_1Y_str);
        stringstream stream_K_m1X_1Y(K_m1X_1Y_str);

        for(int n=0;n<n_Spins;n++){
            stream_K_0X_0Y >> K_0X_0Y(m,n);
            stream_K_1X_0Y >> K_1X_0Y(m,n);
            stream_K_0X_1Y >> K_0X_1Y(m,n);
            stream_K_m1X_1Y >> K_m1X_1Y(m,n);
        }
    }




     //Right now only b/w Classical Spin_i=0 to Spin_j=0--------------------------------------
      
      assert(n_Spins==1);
	
      J_px.resize(3,3);J_py.resize(3,3);J_mxpy.resize(3,3);
      string J_px_file_str = matchstring2(inputfile_, "Nearest_neighbour_Exc_px_Matrix_file_name");
      string J_py_file_str = matchstring2(inputfile_, "Nearest_neighbour_Exc_py_Matrix_file_name");
      string J_mxpy_file_str = matchstring2(inputfile_, "Nearest_neighbour_Exc_mxpy_Matrix_file_name");

        ifstream J_px_stream(J_px_file_str.c_str());
        ifstream J_py_stream(J_py_file_str.c_str());
        ifstream J_mxpy_stream(J_mxpy_file_str.c_str());


        for(int alpha=0;alpha<3;alpha++){
        for(int beta=0;beta<3;beta++){
        J_px_stream>>J_px(alpha,beta);
        J_py_stream>>J_py(alpha,beta);
        J_mxpy_stream>>J_mxpy(alpha,beta);
        }
        }





    //Superexchange matrices done---------------





    lambda_lattice = matchstring (inputfile_, "lambda_lattice");
//    K1x = matchstring(inputfile_, "K");
//    K1y = K1x;
//    cout << "K1x= " << K1x << endl;

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

    lx_cluster=lx;
    ly_cluster=ly;

    assert( (n_Spins==n_orbs) || (n_Spins==1));
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


