#include "../../tensor_type.h"
#include "../../functions.h"
#include "Basis_1orb_SF.h"

using namespace std;

#ifndef Model_1_orb_SF
#define Model_1_orb_SF


class MODEL_1_orb_SF{

public:
Mat_2_Complex_doub Hamil;
Mat_1_doub Evals;
Mat_2_Complex_doub Eigvecs;
Mat_2_doub Theta;   //[time_step][pos]
Mat_2_doub Phi;
Mat_1_doub Theta_eq;   //[time_step][pos]
Mat_1_doub Phi_eq;
Mat_2_doub quant_s_x;
Mat_2_doub quant_s_y;
Mat_2_doub quant_s_z;
Mat_1_doub quant_s_x_eq;
Mat_1_doub quant_s_y_eq;
Mat_1_doub quant_s_z_eq;

Mat_1_doub Jval_array;
Mat_1_doub Jval_array_in_H;
Mat_1_doub Bval_array;

Mat_4_Complex_doub Red_Den_mat;


double S_mag, Jval_p, Jval, B_mag;
int N_total;
double t_hop;
double Beta;
double Mu_;


//******************To reproduce New J. Phys. 17 (2015) 113058********//
bool Single_classical_S_in_B;
int center;
//--------------------------------------------------------------------//

bool Decide_initial_cond_1D;

void Create_Hamil(BASIS_1_orb_SF &Basis, int time_step_);
void Create_Hamil(BASIS_1_orb_SF &Basis, int time_step_, string str);
void Diagonalize_Hamil();
double Calculate_mu(Mat_1_doub &Evals_Temp);
void Get_quantum_Spins(int ts);
void Get_red_den_mat(int ts);

double fermi(Mat_1_doub &Evals_Temp,int n, double mu);
void Create_Hamil(BASIS_1_orb_SF &Basis, Mat_1_doub &Phi_vec,
                  Mat_1_doub &delta_phi_vec, Mat_1_doub &Theta_vec,
                  Mat_1_doub &delta_theta_vec, Mat_2_Complex_doub &Hamil_Temp, double factor);
void Diagonalize_Hamil(Mat_2_Complex_doub &Hamil_Temp, Mat_1_doub &Evals_Temp, Mat_2_Complex_doub &Eigvecs_Temp);
void Get_quantum_Spins(Mat_1_doub &quant_s_y_new,Mat_1_doub &quant_s_x_new,Mat_1_doub &quant_s_z_new,
                       Mat_2_Complex_doub &Hamil_Temp, Mat_2_Complex_doub &Eigvecs_Temp,Mat_1_doub &Evals_Temp, double mu);

};
#endif
