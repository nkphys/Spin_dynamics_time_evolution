#include <algorithm>
#include <functional>
#include <math.h>
#include "../../tensor_type.h"
#include "ParametersEngine_MCMF.h"
#include "Coordinates_MCMF.h"
#include "MFParams_MCMF.h"
#define PI acos(-1.0)

#ifndef Hamiltonian_class_MCMF
#define Hamiltonian_class_MCMF

extern "C" void   zheev_(char *,char *,int *,std::complex<double> *, int *, double *,
                         std::complex<double> *,int *, double *, int *);


class Hamiltonian_MCMF {
public:

    Hamiltonian_MCMF(Parameters_MCMF& Parameters__, Coordinates_MCMF&  Coordinates__, MFParams_MCMF& MFParams__ )
        :Parameters_(Parameters__),Coordinates_(Coordinates__),MFParams_(MFParams__)

    {
        Initialize();
        Hoppings();
        HTBCreate();
    }


    void Initialize();
    void Hoppings();
    double GetCLEnergy();
    void InteractionsCreate();
    void Check_Hermiticity();
    void Check_up_down_symmetry();
    void HTBCreate();
    double chemicalpotential(double muin,double filling);

    double TotalDensity();
    double ClusterDensity();
    double E_QM();

    void Diagonalize(char option);
    void copy_eigs(int i);

    Parameters_MCMF &Parameters_;
    Coordinates_MCMF &Coordinates_;
    MFParams_MCMF &MFParams_;
    int lx_, ly_, ns_, orbs_;
    Matrix<complex<double>> HTB_;
    Matrix<complex<double>> Ham_;
    vector<double> potential_local;

    Matrix<double> Tx,Ty,Tpxpy,Tpxmy;
    vector<double> eigs_,eigs_saved_,sx_,sy_,sz_;

};

#endif
