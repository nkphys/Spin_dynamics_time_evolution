#include <iostream>
#include <math.h>  //fabs(double x) =|x|
#include <algorithm>
#include <stdlib.h>  //for div(q,n).rem(quot),rand
#include <time.h>
#include <fstream>
#include <sstream>
#include <limits>
#include <iomanip>
#include <stdio.h>
#include <random>
#include "ParametersEngine_Rotor.h"
#include "../../tensor_type.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#define PI acos(-1.0)

using namespace std;

#ifndef SC_SW_engine_Rotor
#define SC_SW_engine_Rotor



class SC_SW_ENGINE_Rotor{

public:
    SC_SW_ENGINE_Rotor(Parameters_Rotor& Parameters__)
        : Parameters_(Parameters__)
    {

    }



    int ns;
    double Rotor_mass;
    double time_max;
    double dt_;
    int Runge_Kutta_order;
    int no_of_processors;

    string spins_r_t_out;


    int time_steps;
    double StdDev;
    int TIME_STEP_GLOBAL;

    Mat_1_doub Theta;   //[site]
    Mat_1_doub Theta_eq;   //[site]

    Mat_1_doub Momentum;   //[site]
    Mat_1_doub Momentum_eq;   //[site]

    Mat_1_doub YVec0, YVec1;

    Mat_2_Complex_doub Pauli_x,Pauli_y,Pauli_z;

    Mat_2_int Jzz_connections_neigh;
    Mat_2_doub Jzz_connections_vals;

    bool Use_Scheduler;
    string Scheduler_File;


    Mat_1_doub Gamma_bare;//Field along x component; Gamma_[]*Sx
    Mat_1_doub Js_bare;//Exchange constant in z components Js_[]*Sz*Sz
    Mat_1_doub Time_bare;

    Mat_1_doub Gamma_;//Field along x component; Gamma_[]*Sx
    Mat_1_doub Js_;//Exchange constant in z components Js_[]*Sz*Sz
    Mat_1_doub Time_;


    Mat_1_int index_to_theta, index_to_momentum;
    Mat_1_int theta_to_index, momentum_to_index;

    double Rotor_KinkDen_type1, Rotor_KinkDen_type2;

    mt19937_64 GeneratorNoise_;
    normal_distribution<double> GaussianDistribution;

    mt19937_64 Generator_;
    uniform_real_distribution<double> dis1_;

    double EPSILON;

    Parameters_Rotor& Parameters_;



    void Read_parameters(string filename);
    double random1();
    void Initialize_engine();
    void IndexMapping_bw_Y_and_Variables();
    void Map_Variables_to_Y(Mat_1_doub & Theta_, Mat_1_doub & Momentum_,  Mat_1_doub & Y_);
    void Map_Y_to_Variables(Mat_1_doub & Y_, Mat_1_doub & Theta_, Mat_1_doub & Momentum_);
    void Initial_condition();
    double Get_Classical_Energy(Mat_1_doub & Y_);
    void Evolve_classical_spins_Runge_Kutta(int ts);
    void RungeKuttaOne(Mat_1_doub & Yn, Mat_1_doub & Ynp1);
    void RungeKuttaFour(Mat_1_doub & Yn, Mat_1_doub & Ynp1);
    void RungeKuttaSix(Mat_1_doub & Yn, Mat_1_doub & Ynp1);
    void RungeKuttaEight(Mat_1_doub & Yn, Mat_1_doub & Ynp1);
    void Create_Scheduler();
    void Derivative(Mat_1_doub & Y_, Mat_1_doub & dYbydt);
    void Start_Engine();
    double Get_Kink_Density_1d();
    double Get_Kink_Density_1d_type2();
    double sign(double x);

};



#endif
