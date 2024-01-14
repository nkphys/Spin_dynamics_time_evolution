#include <math.h>
#include "../../tensor_type.h"
#include "ParametersEngine_MultiOrbSF.h"
#include "Coordinates_MultiOrbSF.h"
#include "random"
#include <stdlib.h>
#include <iomanip>
#define PI acos(-1.0)

#ifndef MFParams_MultiOrbSF_class
#define MFParams_MultiOrbSF_class

class MFParams_MultiOrbSF
{
public:
    // Define Fields
    vector<Matrix<double>> etheta, ephi;
    vector<Matrix<double>> Sz, Sx, Sy;
    vector<Matrix<double>> etheta_avg, ephi_avg;
    vector<Matrix<double>> Moment_Size;
    Matrix<double> Disorder;

    // Constructor
    MFParams_MultiOrbSF(Parameters_MultiOrbSF &Parameters__, Coordinates_MultiOrbSF &Coordinates__, mt19937_64 &Generator1__, mt19937_64 &Generator2__)
        : Parameters_(Parameters__), Coordinates_(Coordinates__), Generator1_(Generator1__), Generator2_(Generator2__)
    {
        //setupm_arr();
        initialize();
    }

    double random1();
    double random2();
    void FieldThrow(int site,  int Spin_no, string mc_dof_type);
    void initialize();
    void Adjust_MCWindow();
    void Calculate_Fields_Avg();
    void Read_classical_DOFs(string filename);

    Parameters_MultiOrbSF &Parameters_;
    Coordinates_MultiOrbSF &Coordinates_;
    mt19937_64 &Generator1_; //for random fields
    mt19937_64 &Generator2_; //for random disorder
    int lx_, ly_, ns_;

    uniform_real_distribution<double> dis1_; //for random fields
    uniform_real_distribution<double> dis2_; //for random disorder

    //mt19937_64 mt_rand(Parameters_.RandomSeed);
};

#endif
