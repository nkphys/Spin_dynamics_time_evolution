
#include "Hamiltonian_MCMF.h"

double Hamiltonian_MCMF::chemicalpotential(double muin,double filling){
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


void Hamiltonian_MCMF::Initialize(){


    double DeltaXY=0.4;
    ly_=Parameters_.ly;
    lx_=Parameters_.lx;
    ns_=Parameters_.ns;
    orbs_=Parameters_.orbs;
    assert(orbs_==1);

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

    potential_local.resize(1);
    potential_local[0]=0.0;


} // ----------

double Hamiltonian_MCMF::TotalDensity(){
    double n1=0.0;
    for(int j=0;j<eigs_.size();j++){
        n1 +=  1.0f/( exp(Parameters_.beta*(eigs_[j]-Parameters_.mus) ) + 1.0);
    }
    return n1;
} // ----------


double Hamiltonian_MCMF::E_QM(){
    double E=0.0;
    for(int j=0;j<eigs_.size();j++){
        //E +=  (eigs_[j]-Parameters_.mus)/( exp(Parameters_.beta*(eigs_[j]-Parameters_.mus) ) + 1.0);
        E +=  (eigs_[j])/( exp(Parameters_.beta*(eigs_[j]-Parameters_.mus) ) + 1.0);
    }
    return E;
} // ----------


double Hamiltonian_MCMF::GetCLEnergy(){
    double EClassical;

    // Classical Energy
    EClassical=double(0.0);

    return EClassical;
} // ----------



void Hamiltonian_MCMF::InteractionsCreate(){


    double ei, ai;
    double den;

    Ham_ = HTB_;
    // Ham_.print();

    for (int i = 0; i < ns_; i++)
    { // For each site
        ei = MFParams_.etheta(Coordinates_.indx(i), Coordinates_.indy(i));
        ai = MFParams_.ephi(Coordinates_.indx(i), Coordinates_.indy(i));
        den = MFParams_.Local_density(Coordinates_.indx(i), Coordinates_.indy(i));

        Ham_(i, i) +=  (0.5) * Parameters_.U_COUL * (den);
        Ham_(i + ns_, i + ns_) += (0.5) * Parameters_.U_COUL * (den);
        Ham_(i, i) += (-1.0)*Parameters_.U_COUL * (cos(ei)) * MFParams_.Moment_Size(Coordinates_.indx(i), Coordinates_.indy(i));
        Ham_(i + ns_, i + ns_) += (-1.0)*Parameters_.U_COUL * (-cos(ei)) * MFParams_.Moment_Size(Coordinates_.indx(i), Coordinates_.indy(i));
        Ham_(i, i + ns_) += (-1.0)*Parameters_.U_COUL * sin(ei) * complex<double>(cos(ai), -sin(ai)) * MFParams_.Moment_Size(Coordinates_.indx(i), Coordinates_.indy(i)); //S-
        Ham_(i + ns_, i) += (-1.0)*Parameters_.U_COUL * sin(ei) * complex<double>(cos(ai), sin(ai)) * MFParams_.Moment_Size(Coordinates_.indx(i), Coordinates_.indy(i));  //S+

    }



} // ----------



void Hamiltonian_MCMF::Check_up_down_symmetry()

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




void Hamiltonian_MCMF::Check_Hermiticity()

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





void Hamiltonian_MCMF::Diagonalize(char option){

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



void Hamiltonian_MCMF::HTBCreate(){

    int l, m, a, b;

    HTB_.fill(0.0);

    for (l = 0; l < ns_; l++)
    {

        // * +x direction Neighbor

        m = Coordinates_.neigh(l, 0);
        for (int spin = 0; spin < 2; spin++)
        {

            a = l + ns_ * spin;
            b = m + ns_ * spin;
            assert(a != b);
            if (a != b)
            {
                HTB_(b, a) = complex<double>(1.0 * Parameters_.t_hopping, 0.0);
                HTB_(a, b) = conj(HTB_(b, a));
            }
        }

        // * +y direction Neighbor
        m = Coordinates_.neigh(l, 2);
        for (int spin = 0; spin < 2; spin++)
        {

            a = l + ns_ * spin;
            b = m + ns_ * spin;
            assert(a != b);
            if (a != b)
            {

                HTB_(b, a) = complex<double>(1.0 * Parameters_.t_hopping, 0.0);
                HTB_(a, b) = conj(HTB_(b, a));
            }
        }
    }



} // ----------


void Hamiltonian_MCMF::Hoppings(){


} // ----------

void Hamiltonian_MCMF::copy_eigs(int i){

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


