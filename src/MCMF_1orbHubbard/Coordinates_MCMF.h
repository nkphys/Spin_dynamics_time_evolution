#include "../../tensor_type.h"
#include "ParametersEngine_MCMF.h"
#include "../../Matrix.h"


#ifndef Coordinates_class_MCMF
#define Coordinates_class_MCMF
class Coordinates_MCMF {

public:
    Coordinates_MCMF(int lx, int ly)
        :lx_(lx),ly_(ly)
    {
        Numbering();
    }

    enum {PX=0,MX,PY,MY,PXPY,MXPY,MXMY,PXMY};	// not needed - but keep as a "key"
    void Numbering();
    int indx(int i);
    int indy(int i);
    int Nc(int x, int y);
    int NNc(int x, int y);
    int neigh(int site, int wneigh);
    int getneigh(int site,int wneigh);

    int lx_,ly_,ns_;
    Mat_1_int indx_,indy_;
    Matrix<int> Nc_,neigh_;
};

#endif
