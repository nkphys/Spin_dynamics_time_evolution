
#include "MFParams_MCMF.h"

double MFParams_MCMF::random(){

    return dis_(Generator_);

}

void MFParams_MCMF::initialize(){


    lx_=Coordinates_.lx_;
    ly_=Coordinates_.ly_;

   // srand(Parameters_.RandomSeed);


    etheta_avg.resize(lx_,ly_);
    ephi_avg.resize(lx_,ly_);

    etheta.resize(lx_,ly_);
    ephi.resize(lx_,ly_);

    Moment_Size.resize(lx_,ly_);
    Local_density.resize(lx_, ly_);

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

            Moment_Size(i, j) = random();

            Local_density(i, j) = Parameters_.Fill * 2.0;
        }
    }


} // ----------



void MFParams_MCMF::Read_classical_DOFs(string filename){

    string tmp_str;
    double tmp_double;
    ifstream fl_in(filename.c_str());
    getline(fl_in, tmp_str);

    for(int i=0;i<lx_;i++){
        for(int j=0;j<ly_;j++){
            fl_in>>tmp_double>>tmp_double;

            fl_in>>tmp_double;
            etheta(i,j)=tmp_double;

            fl_in>>tmp_double;
            ephi(i,j)=tmp_double;

            fl_in>>tmp_double;
            Moment_Size(i,j)=tmp_double;

            fl_in>>tmp_double;
            Local_density(i,j)=tmp_double;

        }
    }


} // ----------

