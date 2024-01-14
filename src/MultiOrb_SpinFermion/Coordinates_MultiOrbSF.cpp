#include "Coordinates_MultiOrbSF.h"


/*  
 * ***********
 *  Functions in Class Coordinates ------
 *  ***********
*/  

int Coordinates_MultiOrbSF::indx_basiswise(int i){
    if(i>nbasis_-1){perror("Coordinates.h:x-coordinate of lattice excede limit");}
    return indx_basiswise_[i];
    // ----------
}


int Coordinates_MultiOrbSF::indy_basiswise(int i){
    if(i>nbasis_-1){perror("Coordinates.h:y-coordinate of lattice excede limit");}
    return indy_basiswise_[i];
    // ----------
}

int Coordinates_MultiOrbSF::indx_cellwise(int i){
    if(i>ncells_-1){perror("Coordinates.h:x-coordinate of lattice excede limit");}
    return indx_cellwise_[i];
    // ----------
}


int Coordinates_MultiOrbSF::indy_cellwise(int i){
    if(i>ncells_-1){perror("Coordinates.h:y-coordinate of lattice excede limit");}
    return indy_cellwise_[i];
    // ----------
}

int Coordinates_MultiOrbSF::Nbasis(int x, int y, int orb){
    if(!( (x<lx_&& y<ly_) &&  orb<n_orbs_)){perror("Coordinates.h:ith-sitelabel of lattice excede limit");}
    return Nbasis_[x][y][orb];
    // ----------
}


int Coordinates_MultiOrbSF::Ncell(int x, int y){
    if(!(x<lx_&& y<ly_)){perror("Coordinates.h:ith-sitelabel of lattice excede limit");}
    return Ncell_(x,y);
    // ----------
}


int Coordinates_MultiOrbSF::neigh(int cell, int wneigh){
    if(cell> (lx_*ly_)-1 || wneigh>=8){perror("Coordinates.h:getneigh -> ifstatement-1");}
    return neigh_(cell,wneigh);
} // ----------


void Coordinates_MultiOrbSF::Numbering(){

    nbasis_=lx_*ly_*n_orbs_;
    ncells_= lx_*ly_;

    indx_basiswise_.clear(); 	indx_basiswise_.resize(nbasis_);
    indy_basiswise_.clear();	indy_basiswise_.resize(nbasis_);
    indorb_basiswise_.clear();  indorb_basiswise_.resize(nbasis_);

    indx_cellwise_.clear();indx_cellwise_.resize(ncells_);
    indy_cellwise_.clear();indy_cellwise_.resize(ncells_);

    Nbasis_.resize(lx_);
    for(int ix=0;ix<lx_;ix++){
        Nbasis_[ix].resize(ly_);
        for(int iy=0;iy<ly_;iy++){
            Nbasis_[ix][iy].resize(n_orbs_);
        }
    }

    neigh_.resize(ncells_,8);


    //basis labeling
    int icount=0;
    for(int j=0;j<ly_;j++){
        for(int i=0;i<lx_;i++){
            for(int orb=0;orb<n_orbs_;orb++){
                indx_basiswise_[icount]=i;
                indy_basiswise_[icount]=j;
                indorb_basiswise_[icount]=orb;
                Nbasis_[i][j][orb]=icount;
                icount++;
            }
        }}



    Ncell_.resize(lx_,ly_);
    //cell labeling
    icount=0;
    for(int j=0;j<ly_;j++){
        for(int i=0;i<lx_;i++){
            indx_cellwise_[icount]=i;
            indy_cellwise_[icount]=j;
            Ncell_(i,j)=icount;
            icount++;
        }
    }


    // Neighbors for each unit cell
    for(int i=0;i<ncells_;i++){ 	// ith site
        for(int j=0;j<8;j++) {		// jth neighbor
            neigh_(i,j)=getneigh(i,j);
        }
    }


} // ----------


int Coordinates_MultiOrbSF::getneigh(int site,int wneigh){
    if(site>ncells_-1 || wneigh>7){perror("Coordinates.h:getneigh -> ifstatement-1");}
    int nx=indx_cellwise(site);
    int ny=indy_cellwise(site);
    int mx=0;
    int my=0;

    // Nearest Neighbours
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


    return Ncell(mx,my); //Nc(mx,my);
} // ----------

