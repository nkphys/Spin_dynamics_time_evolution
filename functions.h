#include "tensor_type.h"
#include <assert.h>
#include <complex>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <sstream>
#include <cmath>
//#include "mkl_spblas.h"
//#include <mkl_types.h>
//#include <mkl_cblas.h>
//#include <mkl_lapacke.h>
//#include <mkl.h>
#ifdef _OPENMP
#include <omp.h>
#endif
//#include <bits/stdc++.h>

int Find_int_in_intarray(int num, Mat_1_int &array);
void Print_Matrix_COO(Matrix_COO &A);
void Print_vector_in_file(Mat_1_doub vec, string filename);
void Diagonalize(Matrix_COO &X, double & EG, Mat_1_doub & vecG);
void Matrix_COO_vector_multiplication(string COO_type, Matrix_COO &A,Mat_1_doub &u,Mat_1_doub &v);
double dot_product(Mat_1_doub &vec1, Mat_1_doub &vec2);
void Subtract( Mat_1_doub &temp1, double x, Mat_1_doub &temp2);
void Diagonalize(Mat_2_Complex_doub &Matrix, Mat_1_doub & Evals, Mat_2_Complex_doub &Eigvecs);
void Sum(Matrix_COO A, Matrix_COO B, Matrix_COO & C, double value1, double value2);
void Calculate_recursive_GF(Mat_1_doub A, Mat_1_doub B2, complex<double> &Recursive_GF, double omega,
                            double eta, double GS_energy);
bool comp_greater(double i, double j);
bool comp_greater_pair_double_int(pair_double_int i, pair_double_int j);
static bool sort_using_greater_than(double u, double v);
string NumberToString ( int Number );
bool present_before(Mat_1_int nup_2, Mat_1_int ndn_2, Mat_2_int nup_2_group, Mat_2_int ndn_2_group, int &pos);
void Get_NI_Skw(string filename_full , string filename_specific_k, int dim, double thop);
double E_NI(double kx, double ky, double thop);
double Lorentzian(double x);


