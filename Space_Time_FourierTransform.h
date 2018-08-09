#include <iostream>
#include <math.h>  //fabs(double x) =|x|
#include <algorithm>
#include <stdlib.h>  //for div(q,n).rem(quot),rand
#include <assert.h>
#include <time.h>
#include <fstream>
#include <sstream>
#include <limits>
#include <iomanip>
#include <stdio.h>
#include "tensor_type.h"
#include "src/Model_1orbHubbard/Model_1orb_SF.h"
#include "src/Model_1orbHubbard/Basis_1orb_SF.h"
#define PI acos(-1.0)

using namespace std;

#ifndef ST_fourier
#define ST_fourier

class ST_Fourier{

public:

double w_min,w_max,dw;
double w_conv;
int n_wpoints;

double time_max;
double dt_;
int time_steps;

int no_of_processors;

Mat_1_string conf_inputs;
int No_Of_Inputs;


string Skw_out;
string Skw_out_full;


/*
Mat_3_doub SS_rt;
Mat_2_doub S_rt;


Mat_3_doub sqsq_rt;
Mat_2_doub sq_rt;

Mat_3_doub TT_rt;
Mat_2_doub T_rt;
*/

Mat_2_doub S_tr;

Mat_3_doub S_rw; //
Mat_3_doub S_kw; //

Mat_3_doub s_quantum_rw; //
Mat_3_doub s_quantum_kw; //

Mat_3_doub T_rw; //
Mat_3_doub T_kw; //


Mat_1_doub Sz_eq,Sx_eq,Sy_eq;
Mat_1_doub sz_eq,sx_eq,sy_eq;
Mat_1_doub Tz_eq,Tx_eq,Ty_eq;


Mat_1_doub Sz_t,Sx_t,Sy_t;
Mat_1_doub sz_t,sx_t,sy_t;
Mat_1_doub Tz_t,Tx_t,Ty_t;




void Initialize_engine(BASIS_1_orb_SF & Basis);
//void Read_equilibrium_configuration(BASIS_1_orb_SF & Basis,MODEL_1_orb_SF & Model);
//void Start_Engine(BASIS_1_orb_SF & Basis,MODEL_1_orb_SF & Model);
//void Evolve_classical_spins(int ts,MODEL_1_orb_SF &Model,BASIS_1_orb_SF & Basis);

//void Write_final_time_result(MODEL_1_orb_SF &Model);
//void Read_Restart_Data(MODEL_1_orb_SF &Model);



void Read_parameters(string filename, BASIS_1_orb_SF & Basis);

void Perform_Averaging_on_one_point(BASIS_1_orb_SF & Basis);
void Perform_Smarter_Averaging_on_one_point(BASIS_1_orb_SF & Basis);


void Calculate_Skw_from_Srt_file(BASIS_1_orb_SF & Basis, string filename, string fileout);

};
#endif

