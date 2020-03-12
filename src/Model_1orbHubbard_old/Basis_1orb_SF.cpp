#include <iostream>
#include <math.h>
#include "Basis_1orb_SF.h"
using namespace std;


void BASIS_1_orb_SF::Initialize_Basis(){



    //all above read from inputfile

    if(Geometry=="2D"){


        Max_pos = Lx*Ly - 1;

        No_neighs=4;
        Neighbour.clear();
        Neighbour.resize(Lx*Ly);

        for(int pos=0;pos<=Max_pos;pos++){
            Neighbour[pos].clear();
        }


        int pos;
        int pos_neigh, y_neigh, x_neigh;

        for(int x=0;x<Lx;x++){
            for(int y=0;y<Ly;y++){

                pos = y*Lx + x;

                for(int neigh=0;neigh<4;neigh++){
                if(neigh==0){
                    // ----> +x
                    x_neigh = x + 1;
                    y_neigh = y;
                }
                if(neigh==1){
                    // ----> +y
                    x_neigh = x;
                    y_neigh = y + 1;
                }
                if(neigh==2){
                    // ----> -x
                    x_neigh = x - 1;
                    y_neigh = y;
                }
                if(neigh==3){
                    // ----> -y
                    x_neigh = x;
                    y_neigh = y - 1;
                }



                if( (x_neigh == Lx) && PBC){
                    x_neigh = 0;
                }
                if( (x_neigh == -1) && PBC){
                    x_neigh = Lx-1;
                }

                if( (y_neigh == Ly) && PBC){
                    y_neigh = 0;
                }

                if( (y_neigh == -1) && PBC){
                    y_neigh = Ly-1;
                }

                //neighboring sites have to be on lattice not on moon (unless lattice is on moon)
                if( (x_neigh < Lx) && (x_neigh > -1) && (y_neigh < Ly) && (y_neigh > -1) )
                {
                    pos_neigh = y_neigh*Lx + x_neigh;
                    Neighbour[pos].push_back(pos_neigh);
                }

            }

            }

        }

    }


}
