#include <algorithm>
#include <functional>
#include <math.h>
#include "../../tensor_type.h"
#include "ParametersEngine_MultiOrbSF.h"
#include "Coordinates_MultiOrbSF.h"
#include "MFParams_MultiOrbSF.h"
#define PI acos(-1.0)

#ifndef Hamiltonian_MultiOrbSF_class
#define Hamiltonian_MultiOrbSF_class

extern "C" void zheev_(char *, char *, int *, std::complex<double> *, int *, double *,
                       std::complex<double> *, int *, double *, int *);

class Hamiltonian_MultiOrbSF
{
public:
    Hamiltonian_MultiOrbSF(Parameters_MultiOrbSF &Parameters__, Coordinates_MultiOrbSF &Coordinates__, Coordinates_MultiOrbSF &CoordinatesCluster__, MFParams_MultiOrbSF &MFParams__)
        : Parameters_(Parameters__), Coordinates_(Coordinates__), CoordinatesCluster_(CoordinatesCluster__), MFParams_(MFParams__)

    {
        Initialize();
        Hoppings();
        HTBCreate();
        HTBClusterCreate();
    }

    void Initialize();                                     //::DONE
    void Hoppings();                                       //::DONE
    double GetCLEnergy();                                  //::DONE
    void InteractionsCreate();                             //::DONE
    void InteractionsCreateType2();
    void InteractionsClusterCreate(int Center_site);       //::DONE
    void Check_Hermiticity();                              //::DONE
    void Check_up_down_symmetry();                         //::DONE
    void HTBCreate();                                      //::DONE
    void HTBClusterCreate();                               //::DONE
    double chemicalpotential(double muin, double filling); //::DONE

    double chemicalpotentialCluster(double muin, double filling); //::DONE

    double TotalDensity();   //::DONE
    double ClusterDensity(); //::DONE
    double E_QM();           //::DONE

    double E_QMCluster();                 //::DONE
    void Diagonalize(char option);        //::DONE
    void DiagonalizeCluster(char option); //::DONE
    void copy_eigs(int i);                //::DONE
    void copy_eigs_Cluster(int i);        //::DONE

    void Get_matrix_A();

    Parameters_MultiOrbSF &Parameters_;
    Coordinates_MultiOrbSF &Coordinates_;
    Coordinates_MultiOrbSF &CoordinatesCluster_;
    MFParams_MultiOrbSF &MFParams_;
    int lx_, ly_, ncells_, n_orbs_, n_Spins_;
    int lx_cluster_, ly_cluster_, ncells_cluster;
    Matrix<complex<double>> HTB_;
    Matrix<complex<double>> HTBCluster_;
    Matrix<complex<double>> Ham_;
    Matrix<complex<double>> HamCluster_;
    Matrix<double> Tx, Ty, Tpxpy, Tpxmy;
    vector<double> eigs_, eigsCluster_, eigsCluster_saved_, eigs_saved_;
    Mat_2_doub sx_, sy_, sz_;
    Mat_2_doub Classical_Sx_, Classical_Sy_, Classical_Sz_;
    Matrix<double> IntraCell_Hopp, InterCell_px_Hopp, InterCell_py_Hopp, InterCell_pxmy_Hopp ;

    Matrix<complex<double>> Mat_A;

    double dt_;
    double HS_factor;
};

#endif
