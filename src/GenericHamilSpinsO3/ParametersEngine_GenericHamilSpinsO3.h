#include <assert.h>
#include <complex>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <sstream>
#include <cmath>
#include "../../Matrix.h"
#include "../../tensor_type.h"

using namespace std;

#ifndef Parameters_GenericHamilSpinsO3_class
#define Parameters_GenericHamilSpinsO3_class

class Parameters_GenericHamilSpinsO3
{

public:
    int IterMax, MCNorm, RandomSeed, RandomNoiseSeed;
    Mat_1_int RandomNoiseSeed_array;
    int ns;

    int lx,ly;
    double DampingConst;
    int n_Spins;
    string ModelType;

    double Disorder_Strength, RandomDisorderSeed;
    double Boltzman_constant;

    bool Read_Seed_from_file_;
    string Seed_file_name_;

    int PrintingNoOfTimeSlices;

    int N_ConnectionsFiles;
    Mat_1_string ConnectionFiles;
    Mat_2_string Connections;

    /*Matrix<double> hopping_NN_X;
    Matrix<double> hopping_NN_Y;
    Matrix<double> hopping_NNN_PXPY;
    Matrix<double> hopping_NNN_PXMY;*/


    /*
SavingMicroscopicStates=1
NoOfMicroscopicStates=50
      */


    void Initialize(string inputfile_);
    double matchstring(string file, string match);
    string matchstring2(string file, string match);
};

#endif
