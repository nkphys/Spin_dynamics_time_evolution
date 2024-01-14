#include "ParametersEngine_MultiOrbSF.h"
#include "Coordinates_MultiOrbSF.h"
#include "MFParams_MultiOrbSF.h"
#include "Hamiltonian_MultiOrbSF.h"
#include <iomanip>
#include "../../tensor_type.h"

#ifndef OBSERVABLES_MultiOrbSF_H
#define OBSERVABLES_MultiOrbSF_H

#define PI acos(-1.0)

class Observables_MultiOrbSF
{
public:
    Observables_MultiOrbSF(Parameters_MultiOrbSF &Parameters__, Coordinates_MultiOrbSF &Coordinates__,
                MFParams_MultiOrbSF &MFParams__, Hamiltonian_MultiOrbSF &Hamiltonian__)
        : Parameters_(Parameters__), Coordinates_(Coordinates__), MFParams_(MFParams__),
          Hamiltonian_(Hamiltonian__),
          lx_(Parameters_.lx), ly_(Parameters_.ly), ncells_(Coordinates_.ncells_), nbasis_(Coordinates_.nbasis_), n_orbs_(Coordinates_.n_orbs_)
    {
        Initialize();
    }

    void Initialize();
    void Get_red_den_mat(Mat_4_Complex_doub &Red_Den_mat, double mu);
    void OccDensity();
    void Get_Non_Interacting_dispersion();
    double Lorentzian(double x, double brd);
    void TotalOccDensity();
    void DensityOfStates();
    void SiSjFULL();
    double fermi_function(int n);

    void calculate_quantum_SiSj();
    void quantum_SiSjQ_Average();
    void quantum_SiSj_Average();

    void SiSjQ_Average();
    void SiSj_Average();
    void Total_Energy_Average(double Curr_QuantE, double CurrE);

    void OccDensity(int tlabel);
    void DOSprint(int tlabel);
    complex<double> SiSjQ(int i, int j);
    double SiSj(int i, int j);
    double Omega(int i);

    complex<double> SiSjQ_Mean(int i, int j);
    complex<double> SiSjQ_square_Mean(int i, int j);

    double SiSj_Mean(int i, int j);
    double SiSj_square_Mean(int i, int j);

    double BandWidth;
    Matrix<complex<double>> SiSjQ_, SiSjQ_Mean_, SiSjQ_square_Mean_;
    Matrix<double> SiSj_Mean_, SiSj_square_Mean_;
    double Nematic_order_mean_, Nematic_order_square_mean_;
    Parameters_MultiOrbSF &Parameters_;
    Coordinates_MultiOrbSF &Coordinates_;
    MFParams_MultiOrbSF &MFParams_;
    Hamiltonian_MultiOrbSF &Hamiltonian_;
    int lx_, ly_, ncells_, nbasis_;
    int n_orbs_;
    double dosincr_, tpi_;
    Matrix<double> SiSj_, dos;
    vector<double> sx_, sy_, sz_;
    double AVG_Total_Energy, AVG_Total_Energy_sqr;

    Mat_2_doub local_density;
    Mat_2_doub local_density_Mean;
    Mat_2_doub local_density_square_Mean;

    Matrix<complex<double>> quantum_SiSjQ_, quantum_SiSjQ_Mean_, quantum_SiSjQ_square_Mean_;
    Matrix<complex<double>> quantum_SiSj_, quantum_SiSj_Mean_, quantum_SiSj_square_Mean_;

    void calculate_local_density();
    void local_density_average();
};

#endif // OBSERVABLES_H
