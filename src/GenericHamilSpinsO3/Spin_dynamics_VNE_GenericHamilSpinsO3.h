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
#include "ParametersEngine_GenericHamilSpinsO3.h"
#include "../../tensor_type.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#define PI acos(-1.0)

using namespace std;

#ifndef SC_SW_engine_VNE_GenericHamilSpinsO3
#define SC_SW_engine_VNE_GenericHamilSpinsO3



class SC_SW_ENGINE_VNE_GenericHamilSpinsO3{

public:
    SC_SW_ENGINE_VNE_GenericHamilSpinsO3(Parameters_GenericHamilSpinsO3& Parameters__)
        : Parameters_(Parameters__)
    {

    }



    Matrix<complex<double>> Mat_B;

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


    int n_wpoints;
    int time_steps;

    int TIME_STEP_GLOBAL;

    int n_Spins_;
    Mat_1_doub Theta;   //[ns]
    Mat_1_doub Phi;
    Mat_1_doub Moment_Size;

    Mat_1_doub Theta_eq;   //[ns]
    Mat_1_doub Phi_eq;


    Mat_1_doub Aux_S_x;  //[time_step][pos]
    Mat_1_doub Aux_S_y;
    Mat_1_doub Aux_S_z;

    Mat_1_doub Aux_S_x_eq;
    Mat_1_doub Aux_S_y_eq;
    Mat_1_doub Aux_S_z_eq;


    Matrix<complex<double>> Sp0_Smt;

    Mat_1_Complex_doub YVec0, YVec1;

    Mat_3_doub Theta_time;   //[time_step][lx][ly]
    Mat_3_doub Phi_time;     //[time_step][lx][ly]

    int center;
    double S_mag;

    Mat_3_Complex_doub S_rw;
    Mat_2_Complex_doub S_kw;

    Mat_3_Complex_doub T_rw;
    Mat_2_Complex_doub T_kw;


    Mat_1_int AuxSx_to_index, AuxSy_to_index, AuxSz_to_index;
    Mat_1_int index_to_AuxSx, index_to_AuxSy, index_to_AuxSz;

    Mat_2_Complex_doub Pauli_x,Pauli_y,Pauli_z;


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


    Parameters_GenericHamilSpinsO3& Parameters_;


    Mat_2_int LocalConnectionSize;
    Mat_3_int LocalConnectionSites;
    Mat_3_string LocalConnectionOprs;
    Mat_2_doub LocalConnectionValue;

    string AnsatzType;

    int ns;


    void Read_parameters(string filename);
    void Create_Connections_wrt_site();
    void Initialize_engine();
    double random1();
    void Set_Initial_configuration();
    void Set_Initial_configuration_using_Ansatz();
    void Read_equilibrium_configuration();
    void Map_Variables_to_Y(Mat_1_doub & AuxSx, Mat_1_doub & AuxSy, Mat_1_doub & AuxSz, Mat_1_Complex_doub & Y_);
    void Map_Y_to_Variables(Mat_1_Complex_doub & Y_, Mat_1_doub & AuxSx, Mat_1_doub & AuxSy, Mat_1_doub & AuxSz);
    void IndexMapping_bw_Y_and_Variables();
    complex<double> GetLocalOprExp(string opr_str, int opr_site, Mat_1_Complex_doub & Y_);
    double Calculate_TotalE(Mat_1_Complex_doub & Y_);
    void Start_Engine();

    void Evolve_classical_spins_Runge_Kutta(int ts);
    void RungeKuttaOne(Mat_1_Complex_doub & Yn, Mat_1_Complex_doub & Ynp1);
    void RungeKuttaFour(Mat_1_Complex_doub & Yn, Mat_1_Complex_doub & Ynp1);
    void RungeKuttaSix(Mat_1_Complex_doub & Yn, Mat_1_Complex_doub & Ynp1);
    void RungeKuttaEight(Mat_1_Complex_doub & Yn, Mat_1_Complex_doub & Ynp1);
    void Derivative(Mat_1_Complex_doub & Y_, Mat_1_Complex_doub & dYbydt);
    void Write_final_time_result();
    void Read_Restart_Data();





    /*
    void Create_Scheduler();
    void Read_parameters(string filename);
    void Evolve_classical_spins_Predictor_Corrector();


    double Get_Classical_Energy(Mat_1_Complex_doub & Y_);
    void Update_Hamiltonian_Classical_dof();
    void Update_matrix_B(int ts);
    void Update_UnequalTimeRedDenMatrix();
    void Get_UnequalTimeSpinSpinCorr();
    double Get_Kink_Density_1d();
    double Get_Kink_Density_1d_type2();
*/

};



#endif
