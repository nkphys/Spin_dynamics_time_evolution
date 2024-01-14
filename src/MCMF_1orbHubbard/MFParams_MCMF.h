#include <math.h>
#include "../../tensor_type.h"
#include "ParametersEngine_MCMF.h"
#include "Coordinates_MCMF.h"
#include "random"
#include <stdlib.h>
#define PI acos(-1.0)

#ifndef MFParams_class_MCMF
#define MFParams_class_MCMF

class MFParams_MCMF{
public:

    // Define Fields
    Matrix<double> etheta, ephi;
    Matrix<double> etheta_avg, ephi_avg;
    Matrix<double> Moment_Size;
    Matrix<double> Local_density;

    // Constructor
    MFParams_MCMF(Parameters_MCMF& Parameters__, Coordinates_MCMF&  Coordinates__, mt19937_64& Generator__)
        :Parameters_(Parameters__),Coordinates_(Coordinates__), Generator_(Generator__)
    {
        //setupm_arr();
        initialize();
    }


    double random();
    void initialize();
    void Read_classical_DOFs(string filename);


    Parameters_MCMF &Parameters_;
    Coordinates_MCMF &Coordinates_;
    mt19937_64 &Generator_;
    int lx_,ly_,ns_;
    uniform_real_distribution<double> dis_;
    //mt19937_64 mt_rand(Parameters_.RandomSeed);


};




#endif
