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

#ifndef Parameters_MultiOrbSF_class
#define Parameters_MultiOrbSF_class

class Parameters_MultiOrbSF
{

public:
    int lx, ly, ns, IterMax, MCNorm, RandomSeed;
    int n_orbs;
    string ModelType;
    int TBC_mx, TBC_my;
    int TBC_cellsX, TBC_cellsY;
    int lx_cluster, ly_cluster;
    double mus, mus_Cluster, Fill, pi;
    double Total_Particles;
    double K1x, K1y;
    Mat_1_doub J_Hund, OnSiteE;
    double lambda_lattice;
    double k_const;
    double Disorder_Strength, RandomDisorderSeed;
    double Boltzman_constant;
    double BoundaryConnection;
    Mat_1_doub Temp_values;
    bool Read_Seed_from_file_;
    string Seed_file_name_;

    Matrix<double> hopping_NN_X;
    Matrix<double> hopping_NN_Y;
    Matrix<double> hopping_NNN_PXPY;
    Matrix<double> hopping_NNN_PXMY;

    Mat_1_string MC_DOF;

    bool Cooling_;
    bool ED_;

    bool Metropolis_Algorithm;
    bool Heat_Bath_Algorithm;

    bool MC_on_theta, MC_on_phi, MC_on_theta_and_phi, MC_on_theta_and_phi_and_u, MC_on_moment_size, MC_on_local_density;

    bool fix_mu;

    double fixed_mu_value;

    /*
SavingMicroscopicStates=1
NoOfMicroscopicStates=50
      */
    bool Saving_Microscopic_States;
    int No_Of_Microscopic_States;

    double temp_max, beta_min;
    double temp_min, beta_max;
    double d_Temp;

    int Last_n_sweeps_for_measurement;
    int Measurement_after_each_m_sweeps;

    double temp, beta, Eav, maxmoment;
    double WindowSize, AccCount[2];
    char Dflag;

    void Initialize(string inputfile_);
    double matchstring(string file, string match);
    string matchstring2(string file, string match);
};

#endif
