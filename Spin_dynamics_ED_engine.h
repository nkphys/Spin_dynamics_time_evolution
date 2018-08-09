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
#include "tensor_type.h"
#include "Model_1orb_SF.h"
#include "Basis_1orb_SF.h"
#define PI acos(-1.0)

using namespace std;

#ifndef SC_SW_engine_ED
#define SC_SW_engine_ED

class SC_SW_ENGINE_ED{

public:


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

void Initialize_engine(BASIS_1_orb_SF & Basis,MODEL_1_orb_SF & Model);
void Read_equilibrium_configuration(BASIS_1_orb_SF & Basis,MODEL_1_orb_SF & Model);
void Start_Engine(BASIS_1_orb_SF & Basis,MODEL_1_orb_SF & Model);
void Evolve_classical_spins(int ts,MODEL_1_orb_SF &Model,BASIS_1_orb_SF & Basis);
void Read_parameters(string filename, MODEL_1_orb_SF &Model,BASIS_1_orb_SF & Basis);

};
#endif
