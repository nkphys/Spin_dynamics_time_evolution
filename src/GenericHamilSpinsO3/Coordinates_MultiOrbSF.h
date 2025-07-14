#include "ParametersEngine_MultiOrbSF.h"



#ifndef Coordinates_MultiOrbSF_class
#define Coordinates_MultiOrbSF_class
class Coordinates_MultiOrbSF {

public:
    Coordinates_MultiOrbSF(int lx, int ly, int n_orbs)
        :lx_(lx),ly_(ly),n_orbs_(n_orbs)
    {
        Numbering();
    }

    enum {PX=0,MX,PY,MY,PXPY,MXPY,MXMY,PXMY};	// not needed - but keep as a "key"
    void Numbering();
    int indx_basiswise(int i);
    int indy_basiswise(int i);
    int indorb_basiswise(int i);
    int Nbasis(int x, int y, int orb);

    int indx_cellwise(int i);
    int indy_cellwise(int i);
    int Ncell(int x, int y);

    int neigh(int cell, int wneigh);
    int getneigh(int cell,int wneigh);

    int lx_,ly_,n_orbs_,nbasis_, ncells_;
    Mat_1_int indx_basiswise_,indy_basiswise_, indorb_basiswise_;
    Mat_3_int Nbasis_;
    Matrix<int> Ncell_, neigh_;
    Mat_1_int indx_cellwise_, indy_cellwise_;
};

#endif
