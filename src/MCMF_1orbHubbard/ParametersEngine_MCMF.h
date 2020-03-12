#include <assert.h>
#include <complex>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <sstream>
#include <cmath>
using namespace std;

#ifndef Parameters_class_MCMF
#define Parameters_class_MCMF

class Parameters_MCMF{

public:
    int lx, ly, ns, orbs, IterMax, MCNorm, RandomSeed;
    int TBC_mx, TBC_my;
    int TBC_cellsX, TBC_cellsY;
    int lx_cluster, ly_cluster;
    double mus, mus_Cluster, Fill,pi,U_COUL,k_const,lamda_12,lamda_66;
    double t_hopping;

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

#endif



