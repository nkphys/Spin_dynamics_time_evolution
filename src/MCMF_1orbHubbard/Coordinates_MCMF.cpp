#include "Coordinates_MCMF.h"
/*
 * ***********
 *  Functions in Class Coordinates ------
 *  ***********
*/

int Coordinates_MCMF::indx(int i){
    if(i>ns_-1){perror("Lattice.h:x-coordinate of lattice excede limit");}
    return indx_[i];
 // ----------
}


int Coordinates_MCMF::indy(int i){
    if(i>ns_-1){perror("Lattice.h:y-coordinate of lattice excede limit");}
    return indy_[i];
 // ----------
}

int Coordinates_MCMF::Nc(int x, int y){
    if(x>ly_&&y<ly_){perror("Lattice.h:ith-sitelabel of lattice excede limit");}
    return Nc_(x,y);
 // ----------
}


int Coordinates_MCMF::neigh(int site, int wneigh){
    if(site>ns_-1 || wneigh>=8){perror("Lattice.h:getneigh -> ifstatement-1");}
    return neigh_(site,wneigh);
} // ----------


void Coordinates_MCMF::Numbering(){

    ns_=lx_*ly_;

    indx_.clear(); 	indx_.resize(ns_);
    indy_.clear();	indy_.resize(ns_);
    Nc_.resize(lx_,ly_);
    neigh_.resize(ns_,8);

    // Site labeling
    int icount=0;
    for(int j=0;j<ly_;j++){
        for(int i=0;i<lx_;i++){
            indx_[icount]=i;
            indy_[icount]=j;
            Nc_(i,j)=icount;
            icount++;
        }}

    // Neighbors for each site
    for(int i=0;i<ns_;i++){ 	// ith site
        for(int j=0;j<8;j++) {		// jth neighbor
            neigh_(i,j)=getneigh(i,j);
        }
    }

} // ----------


int Coordinates_MCMF::getneigh(int site,int wneigh){
    if(site>ns_-1 || wneigh>=8){perror("Lattice.h:getneigh -> ifstatement-1");}
    int nx=indx(site);
    int ny=indy(site);
    int mx=0;
    int my=0;

    // Nearest!
    if(wneigh==0){ //PX
        mx=(nx+1)%(lx_);
        my=ny;
    }
    if(wneigh==1){ //MX
        mx=(nx+lx_-1)%(lx_);
        my=ny;
    }
    if(wneigh==2){ //PY
        mx=nx;
        my=(ny+1)%(ly_);
    }
    if(wneigh==3){ //MY
        mx=nx;
        my=(ny+ly_-1)%(ly_);
    }

    // Next-Nearest!
    if(wneigh==4){ //PXPY
        mx=(nx+1)%(lx_);
        my=(ny+1)%(ly_);
    }
    if(wneigh==5){ //MXPY
        mx=(nx+lx_-1)%(lx_);
        my=(ny+1)%(ly_);
    }
    if(wneigh==6){ //MXMY
        mx=(nx+lx_-1)%(lx_);
        my=(ny+ly_-1)%(ly_);
    }
    if(wneigh==7){ //PXMY
        mx=(nx+1)%(lx_);
        my=(ny+ly_-1)%(ly_);
    }
    return Nc_(mx,my); //Nc(mx,my);
} // ----------



