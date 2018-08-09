#include "../../tensor_type.h"
#include "../../functions.h"
using namespace std;

#ifndef Basis_1_orb_SF
#define Basis_1_orb_SF

class BASIS_1_orb_SF{

public:
    string Geometry;
    bool PBC;

    Mat_1_int x_pos,y_pos;

    int Lx,Ly;
    int Max_pos;
    int No_neighs;
    Mat_2_doub Neighbour;



void Initialize_Basis();

};
#endif
