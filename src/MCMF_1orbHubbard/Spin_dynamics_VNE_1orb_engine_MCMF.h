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
#include "Coordinates_MCMF.h"
#include "MFParams_MCMF.h"
#include "Hamiltonian_MCMF.h"
#include "Observables_MCMF.h"
#include "ParametersEngine_MCMF.h"
#include "../../tensor_type.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#define PI acos(-1.0)

using namespace std;

#ifndef SC_SW_engine_VNE_1orb_MCMF
#define SC_SW_engine_VNE_1orb_MCMF



class SC_SW_ENGINE_VNE_1orb_MCMF{

public:
    SC_SW_ENGINE_VNE_1orb_MCMF(Parameters_MCMF& Parameters__, Coordinates_MCMF& Coordinates__,
                                   MFParams_MCMF& MFParams__, Hamiltonian_MCMF& Hamiltonian__, Observables_MCMF& Observables__)
        : Parameters_(Parameters__), Coordinates_(Coordinates__), MFParams_(MFParams__),
          Hamiltonian_(Hamiltonian__), Observables_(Observables__)
    {

    }


    int RANDOM_NO_SEED;
    double w_min,w_max,dw;
    double time_max;
    double dt_;
    int Runge_Kutta_order;
    int no_of_processors;

    string conf_input;
    string spins_r_t_out;
    string Skw_out;
    string Skw_out_full;
    string RK4_type;
    string final_conf_out;
    double w_conv;


    int n_wpoints;
    int time_steps;

    Mat_3_doub Theta;   //[time_step][lx][ly]
    Mat_3_doub Phi;
    Mat_3_doub Moment_Size;   //[time_step][lx][ly]

    Mat_2_doub Theta_eq;   //[time_step][lx][ly]
    Mat_2_doub Phi_eq;

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

    Mat_4_Complex_doub Red_Den_mat; //c^dagger[i][spin]c[j][spin2]


    Mat_5_Complex_doub Red_Den_mat_time; //[ts]c^dagger[i][spin]c[j][spin2]
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


    Parameters_MCMF& Parameters_;
    Coordinates_MCMF& Coordinates_;
    MFParams_MCMF& MFParams_;
    Hamiltonian_MCMF& Hamiltonian_;
    Observables_MCMF& Observables_;


    void Initialize_engine();
    void Read_parameters(string filename);
    void Read_equilibrium_configuration();
    void Start_Engine();
    void Read_Restart_Data();
    void Evolve_classical_spins_Runge_Kutta(int ts);
    void Evolve_classical_spins_Predictor_Corrector();

    void Write_final_time_result();



};



#endif
