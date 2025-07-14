#include <iostream>  //for cin and cout
#include <math.h>  // for pow
#include <stdlib.h>  //for div(q,n).rem(quot),abs(int n)
#include <time.h>
#include <fstream>
#include <sstream>
#include <string>
#include "tensor_type.h"

int main(){

    double JS=1.0;
    double J1_ST=0.5;
    double J2_ST=0.5;
    double 
    double K=0.0;//Ring Exchange
    double JChi=0.0;

    //int L1=2; //X
    //int L2=2; //Y

    int N; 

    //sqrt(8)Xsqrt(8)
    /*
    N=8;
    Mat_1_int Bonds_Xtype_site1 = {0, 1, 2, 3, 4, 5, 6, 7};
    Mat_1_int Bonds_Xtype_site2 = {6, 7, 1, 0, 2, 3, 5, 4};

    Mat_1_int Bonds_Ytype_site1 = {0, 1, 2, 3, 4, 5, 6, 7};
    Mat_1_int Bonds_Ytype_site2 = {2, 3, 5, 4, 6, 7, 1, 0};
    */
    

   //single bond +X
    /*
    N=2;
    Mat_1_int Bonds_Xtype_site1 = {0};
    Mat_1_int Bonds_Xtype_site2 = {1};

    Mat_1_int Bonds_Ytype_site1;
    Mat_1_int Bonds_Ytype_site2;
    Bonds_Ytype_site1.clear();
    Bonds_Ytype_site2.clear();
    */

  //2x2 system
    /*
    N=4;
    Mat_1_int Bonds_Xtype_site1 = {0, 2};
    Mat_1_int Bonds_Xtype_site2 = {1, 3};
    
    Mat_1_int Bonds_Ytype_site1 = {0, 1};
    Mat_1_int Bonds_Ytype_site2 = {2, 3};
   */

   //4x4 OBC(x) PBC(y)
   /*    
    N=16;
    Mat_1_int Bonds_Ytype_site1 = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9,  10,  11,  12,  13,  14,  15};
    Mat_1_int Bonds_Ytype_site2 = {1, 2, 3, 0, 5, 6, 7, 4, 9, 10, 11,   8,  13,  14,  15,  12};

    Mat_1_int Bonds_Xtype_site1 = {0, 1, 2, 3, 4, 5, 6,  7,  8,  9,  10, 11};
    Mat_1_int Bonds_Xtype_site2 = {4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
   */


  //8x8 PBCxPBC--------------------------
 
   int Lx=8;
   int Ly=8;
   N=Lx*Ly;
   
   //index = x +y*Lx
   Mat_1_int Bonds_Ytype_site1, Bonds_Ytype_site2;
   Mat_1_int Bonds_Xtype_site1, Bonds_Xtype_site2;

   int jx,jy,site_j,site_i;
   for(int ix=0;ix<Lx;ix++){
   for(int iy=0;iy<Ly;iy++){
   site_i=ix+iy*Lx;

   //+X neigh
   jx=(ix+1)%Lx;
   jy=iy;
   site_j = jx + jy*Lx;
   Bonds_Xtype_site1.push_back(site_i);Bonds_Xtype_site2.push_back(site_j);

   //+Y neigh
   jy=(iy+1)%Ly;
   jx=ix;
   site_j = jx + jy*Lx;
   Bonds_Ytype_site1.push_back(site_i);Bonds_Ytype_site2.push_back(site_j);

   }
   }


   //Plaquettes
   Mat_1_int Plaquette_site1, Plaquette_site2, Plaquette_site3, Plaquette_site4;
   int ix2, iy2, ix3, iy3, ix4, iy4;
   int site_1, site_2, site_3, site_4;
   for(int ix1=0;ix1<Lx;ix1++){
   for(int iy1=0;iy1<Ly;iy1++){
   
   ix2=(ix1+1)%Lx;iy2=(iy1+0)%Ly;
   ix3=(ix1+0)%Lx;iy3=(iy1+1)%Ly;
   ix4=(ix1+1)%Lx;iy4=(iy1+1)%Ly;

   site_1=ix1+iy1*Lx;	
   site_2=ix2+iy2*Lx;
   site_3=ix3+iy3*Lx;
   site_4=ix4+iy4*Lx;
   
   Plaquette_site1.push_back(site_1);Plaquette_site2.push_back(site_2);
   Plaquette_site4.push_back(site_4);Plaquette_site3.push_back(site_3);

   }
   }

  //-----------------------------------------
   

    //6x4 OBC(x) PBC(y)------------------------
/*    
    N=24;
    Mat_1_int Bonds_Ytype_site1 = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9,  10,  11,  12,  13,  14,  15, 16, 17, 18 ,19, 20, 21, 22, 23};
    Mat_1_int Bonds_Ytype_site2 = {1, 2, 3, 0, 5, 6, 7, 4, 9, 10, 11,   8,  13,  14,  15,  12, 17, 18, 19, 16, 21, 22, 23, 20};

    Mat_1_int Bonds_Xtype_site1 = {0, 1, 2, 3, 4, 5, 6,  7,  8,  9,  10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
    Mat_1_int Bonds_Xtype_site2 = {4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23};
  
 */  
   //--------------------------------------


    string H_SSConnfile_str = "H_SS.dat";
    ofstream H_SSConn_file_stream(H_SSConnfile_str.c_str());
 

    int site_neigh;
    int site;

    Mat_1_string Neigh_type;
    Mat_1_int x_offset, y_offset;
    Mat_1_doub hop_x_or_y;
    Mat_1_string opr_str, opr_str2;
    Neigh_type.push_back("plusX");x_offset.push_back(1);y_offset.push_back(0);
    Neigh_type.push_back("plusY");x_offset.push_back(0);y_offset.push_back(1);
    opr_str.push_back("y");opr_str.push_back("x");
    opr_str2.push_back("x");opr_str2.push_back("y");


    string str_temp;


    //H_SS
    for(int neigh_type=0;neigh_type<Neigh_type.size();neigh_type++){

        Mat_1_int Bonds_site1, Bonds_site2;
        if(Neigh_type[neigh_type]=="plusX"){
            Bonds_site1 = Bonds_Xtype_site1;
            Bonds_site2 = Bonds_Xtype_site2;
        }
        if(Neigh_type[neigh_type]=="plusY"){
            Bonds_site1 = Bonds_Ytype_site1;
            Bonds_site2 = Bonds_Ytype_site2;
        }

        for(int bond_no=0;bond_no<Bonds_site1.size();bond_no++){

    	    site=Bonds_site1[bond_no];
            site_neigh=Bonds_site2[bond_no];

            Mat_1_string SS_type;
            Mat_1_doub fac;
            SS_type.push_back("Sz Sz");fac.push_back(1.0);
            SS_type.push_back("Sx Sx");fac.push_back(1.0);
            SS_type.push_back("Sy Sy");fac.push_back(1.0);

            for(int type=0;type<3;type++){
               H_SSConn_file_stream<<"2 "<<SS_type[type]<<" ";
               H_SSConn_file_stream<<site<<" "<<site_neigh<<" ";
               H_SSConn_file_stream<<fac[type]*J<<endl;
            }

        }
    }//neigh

//--------------------------------------------------



//--------------------------------------------------
//----------------Ring Exchange-----------------------------

    string H_RExcConnfile_str = "H_RExc.dat";
    ofstream H_RExcConn_file_stream(H_RExcConnfile_str.c_str());

    int sitei_,sitej_,sitek_,sitel_;
    double fac_;
   for(int plqt_n=0;plqt_n<Plaquette_site1.size();plqt_n++){

	  //(i.j)(k.l)
   for(int type=0;type<3;type++){
	if(type==0){
	sitei_=Plaquette_site1[plqt_n];
	sitej_=Plaquette_site2[plqt_n];
	sitek_=Plaquette_site3[plqt_n];
	sitel_=Plaquette_site4[plqt_n];
	fac_=1.0;
	}	
	if(type==1){
        sitei_=Plaquette_site1[plqt_n];
        sitej_=Plaquette_site3[plqt_n];
        sitek_=Plaquette_site2[plqt_n];
        sitel_=Plaquette_site4[plqt_n];
        fac_=1.0;
	}
	if(type==2){
        sitei_=Plaquette_site1[plqt_n];
        sitej_=Plaquette_site4[plqt_n];
        sitek_=Plaquette_site2[plqt_n];
        sitel_=Plaquette_site3[plqt_n];
        fac_=-1.0;
	}


	  Mat_1_string SS_type;
            SS_type.push_back("Sz Sz");
            SS_type.push_back("Sx Sx");
            SS_type.push_back("Sy Sy");


   for(int bond_type_ij=0;bond_type_ij<3;bond_type_ij++){
      for(int bond_type_kl=0;bond_type_kl<3;bond_type_kl++){
  
      H_RExcConn_file_stream<<"4 "<<SS_type[bond_type_ij]<<" "<<SS_type[bond_type_kl]<<" ";
               H_RExcConn_file_stream<<sitei_<<" "<<sitej_<<" "<<sitek_<<" "<<sitel_<<" ";
               H_RExcConn_file_stream<<fac_*K<<endl;

   }
   }

   }//type


  }//plqt_n



//-----------------------------------------------------------------







    return 0;
}
