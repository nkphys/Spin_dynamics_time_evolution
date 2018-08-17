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
#include "Coordinates.h"
#include "MFParams.h"
#include "Hamiltonian.h"
#include "Observables.h"
#include "ParametersEngine.h"
#include "../../tensor_type.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#define PI acos(-1.0)

using namespace std;

#ifndef SC_SW_engine_VNE_3orbPnictides
#define SC_SW_engine_VNE_3orbPnictides



class SC_SW_ENGINE_VNE_3orbPnictides{

public:
    SC_SW_ENGINE_VNE_3orbPnictides(Parameters& Parameters__, Coordinates& Coordinates__,
                                   MFParams& MFParams__, Hamiltonian& Hamiltonian__, Observables& Observables__)
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
    Mat_2_doub Theta_eq;   //[time_step][lx][ly]
    Mat_2_doub Phi_eq;

    Mat_3_doub quant_s_x;  //[time_step][pos][orb]
    Mat_3_doub quant_s_y;
    Mat_3_doub quant_s_z;
    Mat_2_doub quant_s_x_eq;
    Mat_2_doub quant_s_y_eq;
    Mat_2_doub quant_s_z_eq;

    Mat_6_Complex_doub Red_Den_mat; //c^dagger[i][orb][spin]c[j][orb2][spin2]

    Mat_1_doub Jval_array;
    Mat_1_doub Jval_array_in_H;
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



    Parameters& Parameters_;
    Coordinates& Coordinates_;
    MFParams& MFParams_;
    Hamiltonian& Hamiltonian_;
    Observables& Observables_;


    void Initialize_engine();
    void Read_parameters(string filename);
    void Read_equilibrium_configuration();
    void Start_Engine();
    void Read_Restart_Data();
    void Evolve_classical_spins(int ts);

    void Write_final_time_result();



};



#endif
