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
#include "../../tensor_type.h"
#include "Model_1orb_SF.h"
#include "Basis_1orb_SF.h"
#define PI acos(-1.0)

using namespace std;

#ifndef SC_SW_engine_VNE_1orbHubbard
#define SC_SW_engine_VNE_1orbHubbard



class SC_SW_ENGINE_VNE_1orbHubbard{

public:

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

void Initialize_engine(MODEL_1_orb_SF & model, BASIS_1_orb_SF &basis);
void Read_equilibrium_configuration(MODEL_1_orb_SF & model, BASIS_1_orb_SF &basis);
void Start_Engine(MODEL_1_orb_SF & model, BASIS_1_orb_SF &basis);
void Evolve_classical_spins(int ts, MODEL_1_orb_SF & model, BASIS_1_orb_SF &basis);
void Read_parameters(string filename, MODEL_1_orb_SF &model, BASIS_1_orb_SF &basis);
void Write_final_time_result(MODEL_1_orb_SF & model, BASIS_1_orb_SF &basis);
void Read_Restart_Data(MODEL_1_orb_SF & model, BASIS_1_orb_SF &basis);


};
#endif
