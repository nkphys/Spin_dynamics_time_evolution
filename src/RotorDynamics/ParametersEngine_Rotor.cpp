#include "ParametersEngine_Rotor.h"

void Parameters_Rotor::Initialize(string inputfile_)
{



    cout << "____________________________________" << endl;
    cout << "Reading the inputfile: '" << inputfile_ <<"'"<< endl;
    cout << "____________________________________" << endl;


    PrintingNoOfTimeSlices=int(matchstring(inputfile_,"PrintingNoOfTimeSlices"));

    DampingConst = double(matchstring(inputfile_, "DampingConst"));
    hz_mag=double(matchstring(inputfile_, "hz_mag")) * 1.0;
    hx_mag=double(matchstring(inputfile_, "hx_mag")) * 1.0;

    ns = int(matchstring(inputfile_,"Total_Sites"));
    cout << "TotalNumberOfSites = " << ns << endl;


    RandomDisorderSeed = matchstring(inputfile_, "RandomDisorderSeed");
    Disorder_Strength = matchstring(inputfile_, "Disorder_Strength");
    Boltzman_constant = matchstring(inputfile_, "Boltzman_constant");


    string RandomNoiseSeed_str = matchstring2(inputfile_, "RandomNoiseSeed");

    stringstream RandomNoiseSeed_stream(RandomNoiseSeed_str);

    int No_NoiseSeeds;
    RandomNoiseSeed_stream>>No_NoiseSeeds;
    RandomNoiseSeed_array.clear();
    int temp_int_min, temp_int_max;
    RandomNoiseSeed_stream>>temp_int_min;
    RandomNoiseSeed_stream>>temp_int_max;
    for(int i=temp_int_min;i<=temp_int_max;i++){
    RandomNoiseSeed_array.push_back(i);
    }




     //Right now only b/w Classical Spin_i=0 to Spin_j=0--------------------------------------
      
	
      Jzz_longrange.resize(ns,ns);
      string Jzz_longrange_file = matchstring2(inputfile_, "Jzz_longrange_file_name");
      ifstream J_zz_stream(Jzz_longrange_file.c_str());

        for(int i=0;i<ns;i++){
        for(int j=0;j<ns;j++){
        J_zz_stream>>Jzz_longrange(i,j);
        }}





    //Superexchange matrices done---------------


    Dflag = 'N';


    temp = double(matchstring(inputfile_, "Temperature")); // temperature in kelvin
    beta = double(1.0/(Boltzman_constant*temp));                         //Beta which is (T*k_b)^-1

    cout << "____________________________________" << endl;
}

double Parameters_Rotor::matchstring(string file, string match)
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

string Parameters_Rotor::matchstring2(string file, string match)
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


