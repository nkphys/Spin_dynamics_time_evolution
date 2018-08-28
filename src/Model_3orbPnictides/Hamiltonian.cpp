
#include "Hamiltonian.h"

double Hamiltonian::chemicalpotential(double muin,double filling){
    double mu_out;
    double n1,N;
    double dMubydN;
    double nstate = eigs_.size();
    dMubydN = 0.05*(eigs_[nstate-1] - eigs_[0])/nstate;
    N=filling*double(eigs_.size());
    //temp=Parameters_.temp;
    mu_out = muin;
    bool converged=false;


    if(1==1){
    for(int i=0;i<50000;i++){
        n1=0.0;
        for(int j=0;j<nstate;j++){
            n1+=double(1.0/( exp( (eigs_[j]-mu_out)*Parameters_.beta ) + 1.0));
        }
        //cout <<"i  "<< i << "  n1  " << n1 << "  mu  " << mu_out<< endl;
        if(abs(N-n1)<double(0.0001)){
            //cout<<abs(N-n1)<<endl;
            converged=true;
            break;
        }
        else {
            mu_out += (N-n1)*dMubydN;
            //cout<<i<<"    "<<n1<<"    "<<N-n1<<endl;

        }
    }

    if(!converged){
        //cout<<"mu_not_converged, N = "<<n1<<endl;
    }
    else{
        //cout<<"mu converged, N = "<<n1<<endl;
    }

    }


    double mu1, mu2;
    double mu_temp = muin;
    //cout<<"mu_input = "<<mu_temp<<endl;
    if(1==2){
        mu1=eigs_[0];
        mu2=eigs_[nstate-1];
    for(int i=0;i<40000;i++){
        n1=0.0;
        for(int j=0;j<nstate;j++){
            n1+=double(1.0/( exp( (eigs_[j]-mu_temp)*Parameters_.beta ) + 1.0));
        }
        //cout <<"i  "<< i << "  n1  " << n1 << "  mu  " << mu_out<< endl;
        if(abs(N-n1)<double(0.0001)){
            //cout<<abs(N-n1)<<endl;
            converged=true;
            break;
        }
        else {
           if(n1 >N){
               mu2=mu_temp;
               mu_temp=0.5*(mu1 + mu_temp);
           }
           else{
               mu1=mu_temp;
               mu_temp=0.5*(mu2 + mu_temp);
           }

        }
        //cout<<"mu_temp = "<<mu_temp<<"   "<<n1<<endl;
    }

    if(!converged){
        //cout<<"mu_not_converged, N = "<<n1<<endl;
    }
    else{
        //cout<<"mu converged, N = "<<n1<<endl;
    }

    mu_out = mu_temp;
    }

    return mu_out;
} // ----------


void Hamiltonian::Initialize(){


    double DeltaXY=0.4;
    ly_=Parameters_.ly;
    lx_=Parameters_.lx;
    ns_=Parameters_.ns;
    orbs_=Parameters_.orbs;

    int space=2*orbs_*ns_;

    Tx.resize(orbs_,orbs_);
    Ty.resize(orbs_,orbs_);
    Tpxpy.resize(orbs_,orbs_);
    Tpxmy.resize(orbs_,orbs_);
    HTB_.resize(space,space);
    Ham_.resize(space,space);
    eigs_.resize(space);
    sx_.resize(space);
    sy_.resize(space);
    sz_.resize(space);
    eigs_saved_.resize(space);

    potential_local.resize(3);
    potential_local[0]=0.0;
    potential_local[1]=0.0;
    potential_local[2]=DeltaXY;


} // ----------

double Hamiltonian::TotalDensity(){
    double n1=0.0;
    for(int j=0;j<eigs_.size();j++){
        n1 +=  1.0f/( exp(Parameters_.beta*(eigs_[j]-Parameters_.mus) ) + 1.0);
    }
    return n1;
} // ----------


double Hamiltonian::E_QM(){
    double E=0.0;
    for(int j=0;j<eigs_.size();j++){
        //E +=  (eigs_[j]-Parameters_.mus)/( exp(Parameters_.beta*(eigs_[j]-Parameters_.mus) ) + 1.0);
        E +=  (eigs_[j])/( exp(Parameters_.beta*(eigs_[j]-Parameters_.mus) ) + 1.0);
    }
    return E;
} // ----------


double Hamiltonian::GetCLEnergy(){
    double EClassical;
    double ei, ai;
    int x,y,site;

    // Spin Components
    for(int i=0;i<lx_;i++){
        for(int j=0;j<ly_;j++){
            site = Coordinates_.Nc(i,j); //+x
            ei=MFParams_.etheta(i,j);
            ai=MFParams_.ephi(i,j);
            sx_[site] = cos(ai) * sin(ei);
            sy_[site] = sin(ai) * sin(ei);
            sz_[site] = cos(ei);
        }
    }

    // Classical Energy
    EClassical=double(0.0);

    // NN TERMS
    for(int i=0;i<ns_;i++){
        site = Coordinates_.neigh(i,0); //+x
        EClassical += 1.14*Parameters_.J_NN * (  sx_[i]*sx_[site]  + sy_[i]*sy_[site] + 1.0*sz_[i]*sz_[site] );
        site = Coordinates_.neigh(i,2); //+y
        EClassical += Parameters_.J_NN * (  sx_[i]*sx_[site]  + sy_[i]*sy_[site] + 1.0*sz_[i]*sz_[site] );

        // NNN TERMS
        site = Coordinates_.neigh(i,4); //+x+y
        EClassical += Parameters_.J_NNN * (  sx_[i]*sx_[site]  + sy_[i]*sy_[site] + 1.0*sz_[i]*sz_[site] );
        site = Coordinates_.neigh(i,7); //+x-y
        EClassical += Parameters_.J_NNN * (  sx_[i]*sx_[site]  + sy_[i]*sy_[site] + 1.0*sz_[i]*sz_[site] );

    }

    return EClassical;
} // ----------



void Hamiltonian::InteractionsCreate(){



    int a,b;
    int space=2*orbs_*ns_;
    double ei, ai,J_Hund=Parameters_.J_HUND;

    for(int i=0;i<space;i++) {
        for(int j=0;j<space;j++) {
            Ham_(i,j)=HTB_(i,j);
        }
    }
    //HUND COUPLING
    for(int i=0;i<ns_;i++) {  // For each site
        ei=MFParams_.etheta(Coordinates_.indx(i),Coordinates_.indy(i));
        ai=MFParams_.ephi(Coordinates_.indx(i),Coordinates_.indy(i));
        for(int k=0;k<orbs_;k++) {  // For each orb
            Ham_(i+k*ns_,i+k*ns_) -=  J_Hund*( cos(ei))*0.5;
            Ham_(i+k*ns_+orbs_*ns_,i+k*ns_+orbs_*ns_) -=  J_Hund*(-cos(ei))*0.5;
            Ham_(i+k*ns_,i+k*ns_+orbs_*ns_) -=  J_Hund* sin(ei)*complex<double>( cos(ai),-sin(ai) )*0.5 ; //S-
            Ham_(i+k*ns_+orbs_*ns_,i+k*ns_) -=  J_Hund* sin(ei)*complex<double>( cos(ai), sin(ai) )*0.5;  //S+
        }
    }


    /*
    Ham_.resize(3,3);
    Ham_(0,0)=one_complex*2.0;Ham_(0,1)=one_complex*3.0;Ham_(0,2)=one_complex*4.0;
    Ham_(1,0)=one_complex*3.0;Ham_(1,1)=one_complex*7.0;Ham_(1,2)=one_complex*5.0;
    Ham_(2,0)=one_complex*4.0;Ham_(2,1)=one_complex*5.0;Ham_(2,2)=one_complex*8.0;
*/



} // ----------



void Hamiltonian::Check_up_down_symmetry()

{
    complex<double> temp(0,0);
    complex<double> temp2;

    for(int i=0;i<orbs_*ns_;i++) {
        for(int j=0;j<orbs_*ns_;j++) {
            temp2 = Ham_(i,j) - Ham_(i+orbs_*ns_,j+orbs_*ns_); //+ Ham_(i+orbs_*ns_,j) + Ham_(i,j+orbs_*ns_);
            temp +=temp2*conj(temp2);
        }
    }

    cout<<"Assymetry in up-down sector: "<<temp<<endl;
}




void Hamiltonian::Check_Hermiticity()

{
    complex<double> temp(0,0);
    complex<double>temp2;

    for(int i=0;i<2*orbs_*ns_;i++) {
        for(int j=0;j<2*orbs_*ns_;j++) {
            assert(Ham_(i,j)==conj(Ham_(j,i))); //+ Ham_(i+orbs_*ns_,j) + Ham_(i,j+orbs_*ns_);
            //temp +=temp2*conj(temp2);
        }
    }

    // cout<<"Hermiticity: "<<temp<<endl;
}





void Hamiltonian::Diagonalize(char option){

    //extern "C" void   zheev_(char *,char *,int *,std::complex<double> *, int *, double *,
     //                      std::complex<double> *,int *, double *, int *);


    char jobz=option;
    char uplo='L'; //WHY ONLY 'L' WORKS?
    int n=Ham_.n_row();
    int lda=Ham_.n_col();
    vector<complex<double>> work(3);
    vector<double> rwork(3*n -2);
    int info;
    int lwork= -1;

    eigs_.resize(Ham_.n_row());
    fill(eigs_.begin(),eigs_.end(),0);
    // query:
    zheev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&(rwork[0]),&info);
    //lwork = int(real(work[0]))+1;
    lwork = int((work[0].real())) ;
    work.resize(lwork);
    // real work:
    zheev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&(rwork[0]),&info);
    if (info!=0) {
        std::cerr<<"info="<<info<<"\n";
        perror("diag: zheev: failed with info!=0.\n");
    }

    // Ham_.print();

    //  for(int i=0;i<eigs_.size();i++){
    //    cout<<eigs_[i]<<endl;
    //}


}



void Hamiltonian::HTBCreate(){

    int space=2*orbs_*ns_;
    int mx=0;
    int my=0;
    complex<double> phasex, phasey;
    int l,m,a,b;
    double l_i, DeltaXY=0.4;
    HTB_.fill(0.0);


    for(l=0;l<ns_;l++) {

        // Phase from As positions
        l_i=pow (-1.00, Coordinates_.indx(l) + Coordinates_.indy(l) );
        //  l_i=1.0;

        // * +x direction Neighbor
        if(Coordinates_.indx(l)==(Coordinates_.lx_ -1)){
            phasex=one_complex;
            phasey=one_complex;
        }
        else{
            phasex=one_complex;
            phasey=one_complex;
        }
        m = Coordinates_.neigh(l,0);
        for(int spin=0;spin<2;spin++) {
            for(int or1=0;or1<orbs_;or1++) {
                for(int or2=0;or2<orbs_;or2++) {
                    a = l + or1*ns_ + ns_*orbs_*spin;
                    b = m + or2*ns_ + ns_*orbs_*spin;
                    assert (a!=b);
                    if(a!=b){
                        if ( (or1==2) ^ (or2==2) ) {
                            HTB_(a,b)=complex<double>(1.0*l_i*Tx(or1,or2),0.0)*phasex;
                            //if(Tx(or1,or2) != 0.0){
                            //cout<<or1<<"   "<<or2<<"  "<<l<<"  "<<m<<"   "<<HTB_(a,b)<<"   "<<endl;
                            // }
                        }
                        else {
                            HTB_(a,b)=complex<double>(1.0*Tx(or1,or2), 0.0)*phasex;
                            //                            if(Tx(or1,or2) != 0.0){
                            //                                cout<<"here"<<endl;
                            //                            }
                        }
                        HTB_(b,a)=conj(HTB_(a,b));
                    }
                }
            }
        }


        // * +y direction Neighbor
        if(Coordinates_.indy(l)==(Coordinates_.ly_ -1)){
            phasex=one_complex;
            phasey=one_complex;
        }
        else{
            phasex=one_complex;
            phasey=one_complex;
        }
        m = Coordinates_.neigh(l,2);
        for(int spin=0;spin<2;spin++) {
            for(int or1=0;or1<orbs_;or1++) {
                for(int or2=0;or2<orbs_;or2++) {
                    a = l + or1*ns_ + ns_*orbs_*spin;
                    b = m + or2*ns_ + ns_*orbs_*spin;
                    assert (a!=b);
                    if(a!=b){
                        if ( (or1==2) ^ (or2==2) ) {
                            HTB_(a,b)=complex<double>(1.0*l_i*Ty(or1,or2),0.0)*phasey;
                        }
                        else {
                            HTB_(a,b)=complex<double>(1.0*Ty(or1,or2),0.0)*phasey;
                            //                        if(Ty(or1,or2) != 0.0){
                            //                            cout<<"here"<<endl;
                            //                        }
                        }
                        HTB_(b,a)=conj(HTB_(a,b));
                    }
                }
            }
        }


        // * +x+y direction Neighbor
        phasex=one_complex;
        phasey=one_complex;
        if( Coordinates_.indy(l)==(Coordinates_.ly_ -1)  ){
            phasey=one_complex;
        }
        if(Coordinates_.indx(l)==(Coordinates_.lx_ -1)){
            phasex=one_complex;
        }

        m = Coordinates_.neigh(l,4);
        for(int spin=0;spin<2;spin++) {
            for(int or1=0;or1<orbs_;or1++) {
                for(int or2=0;or2<orbs_;or2++) {
                    a = l + or1*ns_ + ns_*orbs_*spin;
                    b = m + or2*ns_ + ns_*orbs_*spin;
                    assert (a!=b);
                    if(a!=b){
                        if ( (or1==2) ^ (or2==2) ) {
                            HTB_(a,b)=complex<double>(1.0*l_i*Tpxpy(or1,or2),0.0)*phasex*phasey;
                        }
                        else {
                            HTB_(a,b)=complex<double>(1.0*Tpxpy(or1,or2),0.0)*phasex*phasey;
                        }

                        HTB_(b,a)=conj(HTB_(a,b));
                    }
                }
            }
        }


        // * +x-y direction Neighbor
        phasex=one_complex;
        phasey=one_complex;
        if( Coordinates_.indy(l)==0  ){
            phasey=one_complex;
        }
        if(Coordinates_.indx(l)==(Coordinates_.lx_ -1)){
            phasex=one_complex;
        }
        m = Coordinates_.neigh(l,7);

        for(int spin=0;spin<2;spin++) {
            for(int or1=0;or1<orbs_;or1++) {
                for(int or2=0;or2<orbs_;or2++) {
                    a = l + or1*ns_ + ns_*orbs_*spin;
                    b = m + or2*ns_ + ns_*orbs_*spin;
                    assert (a!=b);
                    if(a!=b){
                        if ( (or1==2) ^ (or2==2) ) {
                            HTB_(a,b)=complex<double>(1.0*l_i*Tpxmy(or1,or2),0.0)*phasex*phasey;
                        }
                        else {
                            HTB_(a,b)=complex<double>(1.0*Tpxmy(or1,or2),0.0)*phasex*phasey;
                        }
                        HTB_(b,a)=conj(HTB_(a,b));
                    }
                }
            }
        }

        // On-Site potential for orb = 2 (xy)
        for(int spin=0;spin<2;spin++) {
            a = l + 2*ns_ + ns_*orbs_*spin;
            HTB_(a,a)=complex<double>(1.0,0.0)*DeltaXY;
        }
    }



} // ----------


void Hamiltonian::Hoppings(){

    double t1,t2,t3,t4,t5,t6,t7,t8;

    //EXACT HAMILTONIAN FROM EQ-1,2,3 from PRB81, 014511 (2010) is implemented
    // Butfollowing have to be used to reproduce the bands shown in Fig-1 of same paper.


    t1 = -0.02;         t2 = -0.06;
    t3 = -0.03;         t4 = 0.01;
    t5 = -0.2;          t6 = -0.3;
    t7 = 0.2;          t8 = 0.1;


    Tx(0,0)=-t2;
    Ty(1,1)=-t2;

    Ty(0,0)=-t1;
    Tx(1,1)=-t1;

    Tpxpy(0,0)=-t3;
    Tpxmy(0,0)=-t3;
    Tpxpy(1,1)=-t3;
    Tpxmy(1,1)=-t3;

    Tpxpy(0,1)=t4;
    Tpxmy(1,0)=-t4;
    Tpxmy(0,1)=-t4;
    Tpxpy(1,0)=t4;

    Tx(2,2)=t5;
    Ty(2,2)=t5;

    Tpxmy(2,2)=-t6;
    Tpxpy(2,2)=-t6;

    /*
    Tx(0,2)=t7;
    Tx(2,0)=t7;
    Ty(1,2)=t7;
    Ty(2,1)=t7;


    Tpxpy(0,2)=-t8;
    Tpxpy(2,0)=t8;
    Tpxmy(0,2)=-t8;
    Tpxmy(2,0)=t8;

    Tpxpy(1,2)=t8;
    Tpxpy(2,1)=-t8;
    Tpxmy(1,2)=-t8;
    Tpxmy(2,1)=t8;
    */

    Tx(0,2)=-t7;
    Tx(2,0)=-t7;
    Ty(1,2)=-t7;
    Ty(2,1)=-t7;


    Tpxpy(0,2)=-t8;
    Tpxpy(2,0)=t8;
    Tpxmy(0,2)=-t8;
    Tpxmy(2,0)=t8;

    Tpxpy(1,2)=-t8;
    Tpxpy(2,1)=t8;
    Tpxmy(1,2)=t8;
    Tpxmy(2,1)=-t8;



} // ----------

void Hamiltonian::copy_eigs(int i){

    int space=2*orbs_*ns_;

    if (i == 0) {
        for(int j=0;j<space;j++) {
            eigs_[j] = eigs_saved_[j];
        }
    }
    else {
        for(int j=0;j<space;j++) {
            eigs_saved_[j] = eigs_[j];
        }
    }

}


