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
#include "Coordinates_MultiOrbSF.h"
#include "MFParams_MultiOrbSF.h"
#include "Hamiltonian_MultiOrbSF.h"
#include "Observables_MultiOrbSF.h"
#include "ParametersEngine_MultiOrbSF.h"
#include "../../tensor_type.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#define PI acos(-1.0)

using namespace std;

#ifndef SC_SW_engine_VNE_MultiOrbSF
#define SC_SW_engine_VNE_MultiOrbSF



class SC_SW_ENGINE_VNE_MultiOrbSF{

public:
    SC_SW_ENGINE_VNE_MultiOrbSF(Parameters_MultiOrbSF& Parameters__, Coordinates_MultiOrbSF& Coordinates__,
                                   MFParams_MultiOrbSF& MFParams__, Hamiltonian_MultiOrbSF& Hamiltonian__, Observables_MultiOrbSF& Observables__)
        : Parameters_(Parameters__), Coordinates_(Coordinates__), MFParams_(MFParams__),
          Hamiltonian_(Hamiltonian__), Observables_(Observables__)
    {

    }




    string Dynamic_Spin_Type;

    double Rotor_mass;
    int RANDOM_NO_SEED;
    double w_min,w_max,dw;
    double time_max;
    double dt_;
    int Runge_Kutta_order;
    int no_of_processors;

    string conf_input;
    string conf_initialize;
    int conf_seed;
    string spins_r_t_out;
    string Skw_out;
    string Skw_out_full;
    string RK4_type;
    string final_conf_out;
    double w_conv;
    bool IgnoreFermions;


    int n_wpoints;
    int time_steps;

    int n_Spins_;
    Mat_3_doub Theta;   //[Spin_no][lx][ly]
    Mat_3_doub Phi;
    Mat_3_doub Moment_Size;   //[Spin_no][lx][ly]

    Mat_3_doub Theta_eq;   //[Spin_no][lx][ly]
    Mat_3_doub Phi_eq;

    Mat_2_doub quant_s_x;  //[time_step][pos]
    Mat_2_doub quant_s_y;
    Mat_2_doub quant_s_z;
    Mat_1_doub quant_s_x_eq;
    Mat_1_doub quant_s_y_eq;
    Mat_1_doub quant_s_z_eq;

    Mat_2_doub Aux_S_x;  //[time_step][pos]
    Mat_2_doub Aux_S_y;
    Mat_2_doub Aux_S_z;
    Mat_1_doub Aux_S_x_eq;
    Mat_1_doub Aux_S_y_eq;
    Mat_1_doub Aux_S_z_eq;

    Mat_4_Complex_doub Red_Den_mat; //c^dagger[i][orb + n_orb*spin]c[j][orb + n_orb*spin2]

    Mat_1_Complex_doub YVec0, YVec1;


    Mat_5_Complex_doub Red_Den_mat_time; //[ts]c^dagger[i][orb + n_orb*spin]c[j][orb + n_orb*spin2]
    Mat_3_doub Theta_time;   //[time_step][lx][ly]
    Mat_3_doub Phi_time;     //[time_step][lx][ly]

    Mat_1_doub Jval_array;
    Mat_1_doub Jval_array_in_H;
    double J_Classical;
    Mat_1_doub Bval_array;
    int center;
    double S_mag;

    Mat_3_Complex_doub S_rw;
    Mat_2_Complex_doub S_kw;

    Mat_3_Complex_doub s_quantum_rw;
    Mat_2_Complex_doub s_quantum_kw;

    Mat_3_Complex_doub T_rw;
    Mat_2_Complex_doub T_kw;


    Mat_2_int AuxSx_to_index, AuxSy_to_index, AuxSz_to_index;
    Mat_4_int RedDen_to_index;
    Mat_1_intpair index_to_AuxSx, index_to_AuxSy, index_to_AuxSz;
    Mat_1_fourint index_to_RedDen;

    Mat_2_Complex_doub Pauli_x,Pauli_y,Pauli_z;


    Mat_2_doub SE_connections_vals;
    vector<Mat_1_inttrio> SE_connections;
    vector<Mat_1_intpair> Hopp_connections;
    Mat_2_doub Hopp_connections_vals;


    bool Use_Scheduler;
    string Scheduler_File;



    Mat_1_doub Gamma_bare;//Field along y component; Gamma_[]*Sy
    Mat_1_doub Js_bare;//Exchange constant in x components Js_[]*Sx*Sx
    Mat_1_doub Time_bare;

    Mat_1_doub Gamma_;//Field along y component; Gamma_[]*Sy
    Mat_1_doub Js_;//Exchange constant in x components Js_[]*Sx*Sx
    Mat_1_doub Time_;

    double AnnealingTime;


    bool RESTART;
    double Restart_Time;
    string Restart_DM_file;
    string Restart_Classical_file;



    bool SAVE;
    string Save_DM_file;
    string Save_Classical_file;


    bool Insitu_SpaceTimeFourier;


    bool Predictor_Corrector;
    bool Use_FFT;

    bool SelfConsistentEvolution;

    mt19937_64 Generator_;
    normal_distribution<double> GaussianDistribution;

    mt19937_64 Generator1_;
    uniform_real_distribution<double> dis1_;


    Parameters_MultiOrbSF& Parameters_;
    Coordinates_MultiOrbSF& Coordinates_;
    MFParams_MultiOrbSF& MFParams_;
    Hamiltonian_MultiOrbSF& Hamiltonian_;
    Observables_MultiOrbSF& Observables_;


    double random1();
    void Create_Scheduler();
    void Initialize_engine();
    void Read_parameters(string filename);
    void Read_equilibrium_configuration();
    void Start_Engine();
    void Read_Restart_Data();
    void Evolve_classical_spins_Runge_Kutta(int ts);
    void Evolve_classical_spins_Predictor_Corrector();
    void Map_Variables_to_Y(Mat_4_Complex_doub & Red_Den_mat_temp, Mat_2_doub & AuxSx, Mat_2_doub & AuxSy, Mat_2_doub & AuxSz, Mat_1_Complex_doub & Y_);
    void Map_Y_to_Variables(Mat_1_Complex_doub & Y_, Mat_4_Complex_doub & Red_Den_mat_temp, Mat_2_doub & AuxSx, Mat_2_doub & AuxSy, Mat_2_doub & AuxSz);
    void RungeKuttaOne(Mat_1_Complex_doub & Yn, Mat_1_Complex_doub & Ynp1);
    void RungeKuttaFour(Mat_1_Complex_doub & Yn, Mat_1_Complex_doub & Ynp1);
    void RungeKuttaSix(Mat_1_Complex_doub & Yn, Mat_1_Complex_doub & Ynp1);
    void RungeKuttaEight(Mat_1_Complex_doub & Yn, Mat_1_Complex_doub & Ynp1);
    void Derivative(Mat_1_Complex_doub & Y_, Mat_1_Complex_doub & dYbydt);
    void Derivative_O3(Mat_1_Complex_doub & Y_, Mat_1_Complex_doub & dYbydt);
    void Derivative_Rotor(Mat_1_Complex_doub & Y_, Mat_1_Complex_doub & dYbydt);
    void Write_final_time_result();
    void IndexMapping_bw_Y_and_Variables();
    double Get_Classical_Energy(Mat_1_Complex_doub & Y_);
    void Set_Initial_configuration();


};



#endif
