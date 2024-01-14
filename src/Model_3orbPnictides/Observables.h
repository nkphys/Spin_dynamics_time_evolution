#include "ParametersEngine.h"
#include "Coordinates.h"
#include "MFParams.h"
#include "Hamiltonian.h"
#include <iomanip>
#include "../../tensor_type.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#ifndef OBSERVABLES_H
#define OBSERVABLES_H

#define PI acos(-1.0)

class Observables{
public:

    Observables(Parameters& Parameters__, Coordinates& Coordinates__,
                MFParams& MFParams__, Hamiltonian& Hamiltonian__)
        : Parameters_(Parameters__), Coordinates_(Coordinates__), MFParams_(MFParams__),
          Hamiltonian_(Hamiltonian__),
          lx_(Parameters_.lx), ly_(Parameters_.ly), ns_(Parameters_.ns)
    {
        Initialize();
    }

    void Initialize();

    void OccDensity();
    void Calculate_Akw();
    void Calculate_Skw(double mu);
    void Calculate_Nw();
    void Get_red_den_mat(Mat_6_Complex_doub &Red_Den_mat, double mu);
    void Get_Non_Interacting_dispersion();
    double Lorentzian(double x, double brd);
    void TotalOccDensity();
    void DensityOfStates();
    void SiSjFULL();


    void OccDensity(int tlabel);
    void DOSprint(int tlabel);
    complex<double> SiSjQ(int i,int j);
    double SiSj(int i,int j);
    double Omega(int i);

    complex<double> SiSjQ_Mean(int i, int j);
    complex<double> SiSjQ_square_Mean(int i, int j);

    double SiSj_Mean(int i, int j);
    double SiSj_square_Mean(int i, int j);

    double BandWidth;
    double nia_t,nib_t,nic_t,n_t;
    Matrix<complex<double>> SiSjQ_, SiSjQ_Mean_, SiSjQ_square_Mean_;
    Matrix<double> SiSj_Mean_, SiSj_square_Mean_;
    double Nematic_order_mean_, Nematic_order_square_mean_;
    Parameters& Parameters_;
    Coordinates& Coordinates_;
    MFParams& MFParams_;
    Hamiltonian& Hamiltonian_;
    int lx_,ly_,ns_;
    double dosincr_,tpi_;
    vector<double> nia_,nib_,nic_;
    Matrix<double> SiSj_,dos;
    vector<double> sx_,sy_,sz_;
    double AVG_Total_Energy, AVG_Total_Energy_sqr;
};

#endif // OBSERVABLES_H
