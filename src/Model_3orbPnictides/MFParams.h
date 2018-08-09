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




double MFParams::random(){

    return dis_(Generator_);

}

void MFParams::initialize(){


    lx_=Coordinates_.lx_;
    ly_=Coordinates_.ly_;

   // srand(Parameters_.RandomSeed);


    etheta_avg.resize(lx_,ly_);
    ephi_avg.resize(lx_,ly_);

    etheta.resize(lx_,ly_);
    ephi.resize(lx_,ly_);

    for(int j=0;j<ly_;j++){
        for(int i=0;i<lx_;i++){
            //ephi(i,j)=(0.5+0.5*pow(-1.0f,i))*Parameters_.pi + grnd()*0.2;
            //etheta(i,j)=0.5*Parameters_.pi + grnd()*0.2;

            //q=(pi,pi)
           // ephi(i,j)=0.0; //(0.5+0.5*pow(-1.0f,i))*Parameters_.pi + grnd()*0.2;
           // etheta(i,j)=0.5*(pow(-1.0,j+i)  + 1.0 )*PI ;//+ grnd()*0.2;

            //q=(0,pi)
            //ephi(i,j)=0.0; //(0.5+0.5*pow(-1.0f,i))*Parameters_.pi + grnd()*0.2;
            //etheta(i,j)=0.5*(pow(-1.0,j)  + 1.0 )*PI; //+ grnd()*0.2;

            //q=(0,0)
           // ephi(i,j)=0.0; //(0.5+0.5*pow(-1.0f,i))*Parameters_.pi + grnd()*0.2;
           // etheta(i,j)=0.0; //+ grnd()*0.2;

            //RANDOM
            ephi(i,j)=2.0*random()*PI;
            etheta(i,j)=random()*PI;
        }
    }


} // ----------



void MFParams::Read_classical_DOFs(string filename){

    string tmp_str;
    double tmp_double;
    ifstream fl_in(filename.c_str());
    fl_in >> tmp_str;

    for(int i=0;i<lx_;i++){
        for(int j=0;j<ly_;j++){
            fl_in>>tmp_double>>tmp_double>>etheta(i,j)>>ephi(i,j);
        }
    }


} // ----------


#endif
