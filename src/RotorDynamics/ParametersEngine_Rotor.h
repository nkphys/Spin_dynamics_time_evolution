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

#ifndef Parameters_Rotor_class
#define Parameters_Rotor_class

class Parameters_Rotor
{

    //ROTOR is in XZ plane
public:
    int ns;
    int RandomNoiseSeed;
    Mat_1_int RandomNoiseSeed_array;

    double DampingConst;

    string ModelType;

    double hz_mag, hx_mag;


    Matrix<double> Jzz_longrange;

    string DisorderType;
    double Disorder_Strength, RandomDisorderSeed;

    double Boltzman_constant;

    int PrintingNoOfTimeSlices;


    double temp, beta;
    char Dflag;

    void Initialize(string inputfile_);
    double matchstring(string file, string match);
    string matchstring2(string file, string match);
};

#endif
