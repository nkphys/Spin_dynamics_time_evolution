#include <math.h>
#include "../../tensor_type.h"
#include "ParametersEngine.h"
#include "Coordinates.h"
#include "random"
#include <stdlib.h>
#define PI acos(-1.0)

#ifndef MFParams_class
#define MFParams_class

class MFParams{ 
public:

    // Define Fields
    Matrix<double> etheta, ephi;
    Matrix<double> etheta_avg, ephi_avg;

    // Constructor
    MFParams(Parameters& Parameters__, Coordinates&  Coordinates__, mt19937_64& Generator__)
        :Parameters_(Parameters__),Coordinates_(Coordinates__), Generator_(Generator__)
    {
        //setupm_arr();
        initialize();
    }


    double random();
    void initialize();
    void Read_classical_DOFs(string filename);


    Parameters &Parameters_;
    Coordinates &Coordinates_;
    mt19937_64 &Generator_;
    int lx_,ly_,ns_;
    uniform_real_distribution<double> dis_;
    //mt19937_64 mt_rand(Parameters_.RandomSeed);


};




#endif
