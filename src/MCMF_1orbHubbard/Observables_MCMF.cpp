#include "Observables_MCMF.h"
/*
 * ***********
 *  Functions in Class Observables ------
 *  ***********
*/




void Observables_MCMF::Calculate_Skw(double mu){


    int  no_of_processors=4;
#ifdef _OPENMP
    omp_set_num_threads(no_of_processors);
    int N_p = omp_get_max_threads();
    cout<<N_p<<" processors are used parallely"<<endl;
#endif

    double Beta=Parameters_.beta;
    complex<double> iota(0.0,1.0);
    Mat_2_Complex_doub Pauli_x,Pauli_y,Pauli_z;
    Pauli_x.resize(2);Pauli_y.resize(2);Pauli_z.resize(2);

    for(int i=0;i<2;i++){
        Pauli_x[i].resize(2); Pauli_y[i].resize(2);Pauli_z[i].resize(2);
        for(int j=0;j<2;j++){
            Pauli_x[i][j]=0;Pauli_y[i][j]=0;Pauli_z[i][j]=0;
        }
    }

    Pauli_x[0][1]=1.0;Pauli_x[1][0]=1.0;
    Pauli_y[0][1]=-1.0*iota;Pauli_y[1][0]=1.0*iota;
    Pauli_z[0][0]=1.0;Pauli_z[1][1]=-1.0;




    //---------Read from input file-----------------------//
    string fileout="Skw.txt";
    double omega_min, omega_max, d_omega;
    double eta = 0.05;
    omega_min=0.0;omega_max=4.0;d_omega=0.025;
    //---------------------------------------------------//


    int omega_index_max = int( (omega_max - omega_min)/(d_omega) );

    ofstream file_Skw_out(fileout.c_str());



    Mat_4_Complex_doub S_jl_nnp;

    /*
    S_lj_nnp.resize(Parameters_.ns);
    for(int l=0;l<Parameters_.ns;l++){
        S_lj_nnp[l].resize(Parameters_.ns) ;
        for(int j=0;j<Parameters_.ns;j++){
            S_lj_nnp[l][j].resize(Hamiltonian_.Ham_.n_row());
            for(int n=0;n<Hamiltonian_.Ham_.n_row();n++){
                S_lj_nnp[l][j][n].resize(Hamiltonian_.Ham_.n_row());
                for(int np=0;np<Hamiltonian_.Ham_.n_row();np++){
                    S_lj_nnp[l][j][n][np]=zero_complex;

                }
            }
        }
    }
    */




    S_jl_nnp.resize(Parameters_.ns);
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
    for(int j=0;j<Parameters_.ns;j++){
        S_jl_nnp[j].resize(Parameters_.ns) ;
        for(int l=0;l<Parameters_.ns;l++){
            S_jl_nnp[j][l].resize(Hamiltonian_.Ham_.n_row());
            for(int n=0;n<Hamiltonian_.Ham_.n_row();n++){
                S_jl_nnp[j][l][n].resize(Hamiltonian_.Ham_.n_row());
                for(int np=0;np<Hamiltonian_.Ham_.n_row();np++){
                    S_jl_nnp[j][l][n][np]=zero_complex;

                    for(int orb1=0;orb1<Parameters_.orbs;orb1++){
                        for(int orb2=0;orb2<Parameters_.orbs;orb2++){

                            for(int a_=0;a_<2;a_++){
                                for(int b_=0;b_<2;b_++){

                                    for(int ap_=0;ap_<2;ap_++){
                                        for(int bp_=0;bp_<2;bp_++){

                                            int c1,c2,c3,c4;
                                            c1 = l + orb1*ns_ + ns_*Parameters_.orbs*a_;
                                            c2 = l + orb1*ns_ + ns_*Parameters_.orbs*b_;

                                            c3 = j + orb2*ns_ + ns_*Parameters_.orbs*ap_;
                                            c4 = j + orb2*ns_ + ns_*Parameters_.orbs*bp_;

                                            S_jl_nnp[j][l][n][np] +=  (Hamiltonian_.Ham_(c1,n)*conj(Hamiltonian_.Ham_(c2,np))*
                                                                       (Hamiltonian_.Ham_(c3,np))*conj(Hamiltonian_.Ham_(c4,n)))*
                                                    ((Pauli_x[a_][b_]*Pauli_x[ap_][bp_]) +
                                                     (Pauli_y[a_][b_]*Pauli_y[ap_][bp_]) +
                                                     (Pauli_z[a_][b_]*Pauli_z[ap_][bp_])
                                                     );

                                        }
                                    }
                                }
                            }
                        }
                    }

                }
            }
        }
        cout <<"j = "<<j<<endl;
    }





    Mat_3_Complex_doub S_rw;
    S_rw.resize(Parameters_.ns);

    for (int i=0;i<Parameters_.ns;i++){
        S_rw[i].resize(Parameters_.ns);
        for(int j=0;j<Parameters_.ns;j++){
            S_rw[i][j].resize(omega_index_max);
        }
    }


    complex<double> Local_Quantum_moment(0,0);



#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
    for (int j=0;j<Parameters_.ns;j++){
        for (int l=0;l<Parameters_.ns;l++){
            for(int omega_ind=0;omega_ind<omega_index_max;omega_ind++){

                S_rw[j][l][omega_ind]=zero_complex;
                for(int n=0;n<Hamiltonian_.Ham_.n_row();n++){
                    for(int np=0;np<Hamiltonian_.Ham_.n_row();np++){
                        if(n!=np){



                            //c= l + orb*ns_ + ns_*Parameters_.orbs*spin;
                            //Hamiltonian_.Ham_(c2,n) is nth eigenvector and c2th component [checked];



                            S_rw[j][l][omega_ind] += S_jl_nnp[j][l][n][np]*
                                    Lorentzian( omega_min + (omega_ind*d_omega) -
                                                Hamiltonian_.eigs_[np] + Hamiltonian_.eigs_[n], eta)*
                                    (1.0/(1.0 + exp(Beta*(Hamiltonian_.eigs_[n] - mu))))*
                                    (1.0/(1.0 + exp(-1.0*Beta*(Hamiltonian_.eigs_[np] - mu))))
                                    ;



                        }

                    }
                }

                cout<<"j,l,omega_ind = "<<j<<"  "<<l<<"  "<<omega_ind<<"   "<<endl;
            }
        }

    }


    complex<double> temp;
    double kx,ky;
    int kx_i,ky_i;

    Mat_1_intpair k_path;
    k_path.clear();
    pair_int temp_pair;

    //--------\Gamma to X-----------------
    ky_i=0;
    for(kx_i=0;kx_i<=(Parameters_.lx/2);kx_i++){
        temp_pair.first = kx_i;
        temp_pair.second = ky_i;
        k_path.push_back(temp_pair);
    }
    //----------------------------------

    //--------X to M-----------------
    kx_i=(Parameters_.lx/2);
    for(ky_i=1;ky_i<=(Parameters_.lx/2);ky_i++){
        temp_pair.first = kx_i;
        temp_pair.second = ky_i;
        k_path.push_back(temp_pair);
    }
    //----------------------------------

    //--------M to \Gamma[with one extra point,
    //                  because in gnuplot use "set pm3d corners2color c1"
    //                  ]-----------------
    kx_i=(Parameters_.lx/2) - 1;
    ky_i=(Parameters_.lx/2) - 1;
    for(kx_i=(Parameters_.lx/2) - 1;kx_i>=-1;kx_i--){
        temp_pair.first = kx_i;
        temp_pair.second = kx_i;
        k_path.push_back(temp_pair);
    }
    //----------------------------------



    for(int k_point=0;k_point<k_path.size();k_point++){

        kx_i=k_path[k_point].first;
        ky_i=k_path[k_point].second;
        kx=(2.0*PI*kx_i)/(1.0*Parameters_.lx);
        ky=(2.0*PI*ky_i)/(1.0*Parameters_.ly);

        for(int omega_ind=0;omega_ind<omega_index_max;omega_ind++){
            temp=zero_complex;

            for(int j=0;j<ns_;j++){
                for(int l=0;l<ns_;l++){
                    temp+= one_complex*
                            exp(iota_complex*(kx*(Coordinates_.indx(j) - Coordinates_.indx(l)) +
                                              ky*(Coordinates_.indy(j) - Coordinates_.indy(l))))*
                            S_rw[j][l][omega_ind];
                }
            }
            //Use 1:6:7----for gnuplot
            file_Skw_out<< k_point<<"   "<<kx_i<<"   "<<ky_i<<"   "<<(ky_i*Parameters_.lx) + kx_i<<"    "<<
                           omega_min + (d_omega*omega_ind)<<"   "<<omega_ind<<"    "
                        <<temp.real()<<"    "<<temp.imag()<<"    "<<endl;

        }
        file_Skw_out<<endl;

        cout<<"k index ="<<k_point<<" done"<<endl;
    }



}

void Observables_MCMF::Calculate_Akw(){


    //---------Read from input file-----------------------//
    string fileout="Akw.txt";
    double omega_min, omega_max, d_omega;
    double eta = 0.001;
    omega_min=-1.6;omega_max=2.6;d_omega=0.03;
    //---------------------------------------------------//


    int omega_index_max = int( (omega_max - omega_min)/(d_omega) );

    ofstream file_Akw_out(fileout.c_str());

    int c1,c2;

    Mat_3_Complex_doub A_up_00, A_up_11, A_up_22;
    Mat_3_Complex_doub A_dn_00, A_dn_11, A_dn_22;
    A_up_00.resize(Parameters_.ns);
    A_up_11.resize(Parameters_.ns);
    A_up_22.resize(Parameters_.ns);
    A_dn_00.resize(Parameters_.ns);
    A_dn_11.resize(Parameters_.ns);
    A_dn_22.resize(Parameters_.ns);
    for (int i=0;i<Parameters_.ns;i++){
        A_up_00[i].resize(Parameters_.ns);
        A_up_11[i].resize(Parameters_.ns);
        A_up_22[i].resize(Parameters_.ns);
        A_dn_00[i].resize(Parameters_.ns);
        A_dn_11[i].resize(Parameters_.ns);
        A_dn_22[i].resize(Parameters_.ns);
        for(int j=0;j<Parameters_.ns;j++){
            A_up_00[i][j].resize(omega_index_max);
            A_up_11[i][j].resize(omega_index_max);
            A_up_22[i][j].resize(omega_index_max);
            A_dn_00[i][j].resize(omega_index_max);
            A_dn_11[i][j].resize(omega_index_max);
            A_dn_22[i][j].resize(omega_index_max);
        }
    }


    complex<double> Nup_check(0,0);
    complex<double> Ndn_check(0,0);

    for (int j=0;j<Parameters_.ns;j++){
        for (int l=0;l<Parameters_.ns;l++){
            for(int omega_ind=0;omega_ind<omega_index_max;omega_ind++){
                A_up_00[j][l][omega_ind]=zero_complex;
                A_up_11[j][l][omega_ind]=zero_complex;
                A_up_22[j][l][omega_ind]=zero_complex;
                A_dn_00[j][l][omega_ind]=zero_complex;
                A_dn_11[j][l][omega_ind]=zero_complex;
                A_dn_22[j][l][omega_ind]=zero_complex;

                for(int n=0;n<Hamiltonian_.Ham_.n_row();n++){

                    //c= l + or1*ns_ + ns_*orbs_*spin;

                    //Hamiltonian_.Ham_(c2,n) is nth eigenvector and c2th component [checked];
                    c1 = l + ns_*Parameters_.orbs; c2 = j+ ns_*Parameters_.orbs;
                    A_dn_00[j][l][omega_ind] +=  conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c2,n)*
                            Lorentzian( omega_min + (omega_ind*d_omega) - Hamiltonian_.eigs_[n], eta);

                    c1 = l+ ns_+ ns_*Parameters_.orbs; c2 = j+ ns_+ ns_*Parameters_.orbs;
                    A_dn_11[j][l][omega_ind] +=  conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c2,n)*
                            Lorentzian( omega_min + (omega_ind*d_omega) - Hamiltonian_.eigs_[n], eta);

                    c1 = l + 2*ns_+ ns_*Parameters_.orbs; c2 = j + 2*ns_+ ns_*Parameters_.orbs;
                    A_dn_22[j][l][omega_ind] +=  conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c2,n)*
                            Lorentzian( omega_min + (omega_ind*d_omega) - Hamiltonian_.eigs_[n], eta);



                    c1 = l; c2 = j;
                    A_up_00[j][l][omega_ind] +=  conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c2,n)*
                            Lorentzian( omega_min + (omega_ind*d_omega) - Hamiltonian_.eigs_[n], eta);

                    c1 = l+ ns_; c2 = j+ ns_;
                    A_up_11[j][l][omega_ind] +=  conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c2,n)*
                            Lorentzian( omega_min + (omega_ind*d_omega) - Hamiltonian_.eigs_[n], eta);

                    c1 = l + 2*ns_; c2 = j + 2*ns_;
                    A_up_22[j][l][omega_ind] +=  conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c2,n)*
                            Lorentzian( omega_min + (omega_ind*d_omega) - Hamiltonian_.eigs_[n], eta);






                }

                if(j==l){
                    Nup_check += (A_up_00[j][l][omega_ind] + A_up_11[j][l][omega_ind] + A_up_22[j][l][omega_ind])*d_omega;
                    Ndn_check += (A_dn_00[j][l][omega_ind] + A_dn_11[j][l][omega_ind] + A_dn_22[j][l][omega_ind])*d_omega;
                }
            }
        }
    }

    cout << "Nup_check = "<<Nup_check<<endl;
    cout << "Ndn_check = "<<Ndn_check<<endl;

    complex<double> temp_up_00,temp_up_11,temp_up_22;
    complex<double> temp_dn_00,temp_dn_11,temp_dn_22;
    double kx,ky;
    int kx_i,ky_i;

    Mat_1_intpair k_path;
    k_path.clear();
    pair_int temp_pair;

    //--------\Gamma to X-----------------
    ky_i=0;
    for(kx_i=0;kx_i<=(Parameters_.lx/2);kx_i++){
        temp_pair.first = kx_i;
        temp_pair.second = ky_i;
        k_path.push_back(temp_pair);
    }
    //----------------------------------

    //--------X to M-----------------
    kx_i=(Parameters_.lx/2);
    for(ky_i=1;ky_i<=(Parameters_.lx/2);ky_i++){
        temp_pair.first = kx_i;
        temp_pair.second = ky_i;
        k_path.push_back(temp_pair);
    }
    //----------------------------------

    //--------M to \Gamma[with one extra point,
    //                  because in gnuplor use "set pm3d corners2color c1"
    //                  ]-----------------
    kx_i=(Parameters_.lx/2) - 1;
    ky_i=(Parameters_.lx/2) - 1;
    for(kx_i=(Parameters_.lx/2) - 1;kx_i>=-1;kx_i--){
        temp_pair.first = kx_i;
        temp_pair.second = kx_i;
        k_path.push_back(temp_pair);
    }
    //----------------------------------


    double k22_offset=PI;
    for(int k_point=0;k_point<k_path.size();k_point++){

        kx_i=k_path[k_point].first;
        ky_i=k_path[k_point].second;
        kx=(2.0*PI*kx_i)/(1.0*Parameters_.lx);
        ky=(2.0*PI*ky_i)/(1.0*Parameters_.ly);

        for(int omega_ind=0;omega_ind<omega_index_max;omega_ind++){
            temp_up_00=zero_complex;temp_up_11=zero_complex;temp_up_22=zero_complex;
            temp_dn_00=zero_complex;temp_dn_11=zero_complex;temp_dn_22=zero_complex;

            for(int j=0;j<ns_;j++){
                for(int l=0;l<ns_;l++){
                    temp_up_00 += one_complex*
                            exp(iota_complex*(kx*(Coordinates_.indx(j) - Coordinates_.indx(l)) +
                                              ky*(Coordinates_.indy(j) - Coordinates_.indy(l))))*
                            A_up_00[j][l][omega_ind];
                    temp_up_11 += one_complex*
                            exp(iota_complex*(kx*(Coordinates_.indx(j) - Coordinates_.indx(l)) +
                                              ky*(Coordinates_.indy(j) - Coordinates_.indy(l))))*
                            A_up_11[j][l][omega_ind];
                    temp_up_22 += one_complex*
                            exp(iota_complex*((kx+k22_offset)*(Coordinates_.indx(j) - Coordinates_.indx(l)) +
                                              (ky+k22_offset)*(Coordinates_.indy(j) - Coordinates_.indy(l))))*
                            A_up_22[j][l][omega_ind];
                    temp_dn_00 += one_complex*
                            exp(iota_complex*(kx*(Coordinates_.indx(j) - Coordinates_.indx(l)) +
                                              ky*(Coordinates_.indy(j) - Coordinates_.indy(l))))*
                            A_dn_00[j][l][omega_ind];
                    temp_dn_11 += one_complex*
                            exp(iota_complex*(kx*(Coordinates_.indx(j) - Coordinates_.indx(l)) +
                                              ky*(Coordinates_.indy(j) - Coordinates_.indy(l))))*
                            A_dn_11[j][l][omega_ind];
                    temp_dn_22 += one_complex*
                            exp(iota_complex*((kx+k22_offset)*(Coordinates_.indx(j) - Coordinates_.indx(l)) +
                                              (ky+k22_offset)*(Coordinates_.indy(j) - Coordinates_.indy(l))))*
                            A_dn_22[j][l][omega_ind];

                }
            }
            //Use 1:6:7----for gnuplot
            file_Akw_out<< k_point<<"   "<<kx_i<<"   "<<ky_i<<"   "<<(ky_i*Parameters_.lx) + kx_i<<"    "<<
                           omega_min + (d_omega*omega_ind)<<"   "<<omega_ind<<"    "<<temp_up_00.real()<<"    "<<temp_up_11.real()
                        <<"    "<<
                          temp_up_22.real()<<"    "<<temp_dn_00.real()<<"    "<<temp_dn_11.real()<<"    "<<temp_dn_22.real()<<"    "<<temp_up_00.imag()<<"    "<<temp_up_11.imag()
                       <<"    "<<
                         temp_up_22.imag()<<"    "<<temp_dn_00.imag()<<"    "<<temp_dn_11.imag()<<"    "<<temp_dn_22.imag()<<"    "<<endl;

        }
        file_Akw_out<<endl;
    }



}

void Observables_MCMF::Get_red_den_mat(Mat_4_Complex_doub &Red_Den_mat, double mu){

    double Beta=Parameters_.beta;
    //Red_Den_mat; //c^dagger[j][spin]c[l][spin2]
    int c1, c2;
    for (int j=0;j<Parameters_.ns;j++){
        for (int l=0;l<Parameters_.ns;l++){
            for (int s=0;s<2;s++){
                for (int s2=0;s2<2;s2++){
                    Red_Den_mat[j][s][l][s2]=zero_complex;
                    for (int n=0;n<2*Parameters_.ns;n++){
                        c1 = j + s*ns_;
                        c2 = l + s2*ns_;
                        //remember psi goes with anhilation
                        //Red_Den_mat[j][s][l][s2] = <C_{l,s2}^{dagger} C_{j,s}>
                        Red_Den_mat[j][s][l][s2] += conj(Hamiltonian_.Ham_(c2,n))*(Hamiltonian_.Ham_(c1,n))*(1.0/(1.0 + exp(Beta*(Hamiltonian_.eigs_[n] - mu))));

                    }

                }
            }
        }
    }


}

void Observables_MCMF::Calculate_Nw(){

    //---------Read from input file-----------------------//
    string fileout="Nw.txt";
    double omega_min, omega_max, d_omega;
    double eta = 0.005;
    omega_min=1.0;omega_max=3.0;d_omega=0.00025;
    //---------------------------------------------------//

    int omega_index_max = int( (omega_max - omega_min)/(d_omega) );

    ofstream file_Nw_out(fileout.c_str());

    complex<double> temp_val11, temp_val22, temp_val00 ;
    int c1;


    for(int omega_ind=0;omega_ind<omega_index_max;omega_ind++){

        temp_val11=zero_complex;
        temp_val22=zero_complex;
        temp_val00=zero_complex;
        //l + 2*ns_ + ns_*orbs_*spin

        for (int j=0;j<Parameters_.ns;j++){
            for(int n=0;n<Hamiltonian_.Ham_.n_row();n++){

                c1 = j + ns_*Parameters_.orbs;
                temp_val00 +=  conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c1,n)*
                        Lorentzian( omega_min + (omega_ind*d_omega) - Hamiltonian_.eigs_[n], eta);

                c1 = j+ ns_+ ns_*Parameters_.orbs;
                temp_val11 +=  conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c1,n)*
                        Lorentzian( omega_min + (omega_ind*d_omega) - Hamiltonian_.eigs_[n], eta);

                c1 = j + 2*ns_+ ns_*Parameters_.orbs;
                temp_val22 +=  conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c1,n)*
                        Lorentzian( omega_min + (omega_ind*d_omega) - Hamiltonian_.eigs_[n], eta);



                c1 = j;
                temp_val00 +=  conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c1,n)*
                        Lorentzian( omega_min + (omega_ind*d_omega) - Hamiltonian_.eigs_[n], eta);

                c1 = j+ ns_;
                temp_val11 +=  conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c1,n)*
                        Lorentzian( omega_min + (omega_ind*d_omega) - Hamiltonian_.eigs_[n], eta);

                c1 = j + 2*ns_;
                temp_val22 +=  conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c1,n)*
                        Lorentzian( omega_min + (omega_ind*d_omega) - Hamiltonian_.eigs_[n], eta);


            }
        }

        file_Nw_out<<omega_min + (omega_ind*d_omega)<<"     "<<temp_val00.real()<<"     "
                  <<temp_val11.real()<<"     "<<temp_val22.real()<<endl;

    }


}

void Observables_MCMF::Get_Non_Interacting_dispersion(){

    complex<double> one(1,0);
    complex<double> iota(0,1);
    complex<double> t1_,t2_,t3_,t4_,t5_,t6_,t7_,t8_;
    complex<double> Delta_xy_;
    Matrix<complex<double>> H;
    H.resize(3,3);

    char option = 'V';

    /*
    t1_ = 0.02*one;  t2_ = 0.06*one;
    t3_ = 0.03*one;  t4_ = -0.01*one;
    t5_ = 0.2*one;   t6_ = 0.3*one;
    */

    t1_ = -0.02*one;  t2_ = -0.06*one;
    t3_ = -0.03*one;  t4_ = 0.01*one;
    t5_ = -0.2*one;   t6_ = -0.3*one;
    t7_ = 0.2*one;  t8_ = 0.1*one;






    Delta_xy_=0.4*one;


    string fileout="k_vs_E_NI.txt";

    int counter_k=0;
    double dk_ =0.01;
    ofstream file_out(fileout.c_str());



    double kx_=0.0;
    double ky_=0.0;

    while(kx_<PI){

        H(0,0)=one*(-2.0*t2_*cos(kx_)  -2.0*t1_*cos(ky_)  -4.0*t3_*cos(kx_)*cos(ky_));
        H(1,1)=one*(-2.0*t1_*cos(kx_) -2.0*t2_*cos(ky_) - 4.0*t3_*cos(kx_)*cos(ky_));
        H(2,2)=one*(-2.0*t5_*(cos(kx_) + cos(ky_)) - 4.0*t6_*cos(kx_)*cos(ky_) + Delta_xy_);
        H(0,1)=one*(-4.0*t4_*sin(kx_)*sin(ky_));

        H(0,2)=iota*(-2.0*t7_*sin(kx_) + 4.0*t8_*sin(kx_)*cos(ky_) );
        H(1,2)=iota*(-2.0*t7_*sin(ky_) + 4.0*t8_*sin(ky_)*cos(kx_) );

        /*
        H(0,2)=(2.0*t7_*cos(kx_) + 4.0*iota*t8_*sin(kx_)*cos(ky_) );
       H(1,2)=(2.0*t7_*cos(ky_) + 4.0*iota*t8_*sin(ky_)*cos(kx_) );
       */






        char jobz=option;
        char uplo='U';
        int n=H.n_row();
        int lda=H.n_col();
        vector<complex<double>> work(3);
        vector<double> rwork(3*n);
        int info,lwork= -1;

        vector<double> E_;
        E_.resize(3);
        fill(E_.begin(), E_.end(),0);
        // query:
        zheev_(&jobz,&uplo,&n,&(H(0,0)),&lda,&(E_[0]),&(work[0]),&lwork,&(rwork[0]),&info);
        lwork = int(real(work[0]))+1;
        work.resize(lwork+1);
        // real work:
        zheev_(&jobz,&uplo,&n,&(H(0,0)),&lda,&(E_[0]),&(work[0]),&lwork,&(rwork[0]),&info);
        if (info!=0) {
            std::cerr<<"info="<<info<<"\n";
            perror("diag: zheev: failed with info!=0.\n");
        }

        file_out<<counter_k<<"    "<<kx_<<"   "<<ky_<<"    "<<E_[0]<<"     "<<E_[1]<<"     "<<E_[2]<<endl;
        kx_= kx_ + dk_;
        counter_k++;
    }


    kx_=PI;
    ky_=0;
    while(ky_<PI){

        H(0,0)=one*(-2.0*t2_*cos(kx_)  -2.0*t1_*cos(ky_)  -4.0*t3_*cos(kx_)*cos(ky_));
        H(1,1)=one*(-2.0*t1_*cos(kx_) -2.0*t2_*cos(ky_) - 4.0*t3_*cos(kx_)*cos(ky_));
        H(2,2)=one*(-2.0*t5_*(cos(kx_) + cos(ky_)) - 4.0*t6_*cos(kx_)*cos(ky_) + Delta_xy_);
        H(0,1)=one*(-4.0*t4_*sin(kx_)*sin(ky_));

        H(0,2)=iota*(-2.0*t7_*sin(kx_) + 4.0*t8_*sin(kx_)*cos(ky_) );
        H(1,2)=iota*(-2.0*t7_*sin(ky_) + 4.0*t8_*sin(ky_)*cos(kx_) );

        char jobz=option;
        char uplo='U';
        int n=H.n_row();
        int lda=H.n_col();
        vector<complex<double>> work(3);
        vector<double> rwork(3*n);
        int info,lwork= -1;

        vector<double> E_;
        E_.resize(3);
        fill(E_.begin(), E_.end(),0);
        // query:
        zheev_(&jobz,&uplo,&n,&(H(0,0)),&lda,&(E_[0]),&(work[0]),&lwork,&(rwork[0]),&info);
        lwork = int(real(work[0]))+1;
        work.resize(lwork+1);
        // real work:
        zheev_(&jobz,&uplo,&n,&(H(0,0)),&lda,&(E_[0]),&(work[0]),&lwork,&(rwork[0]),&info);
        if (info!=0) {
            std::cerr<<"info="<<info<<"\n";
            perror("diag: zheev: failed with info!=0.\n");
        }

        file_out<<counter_k<<"    "<<kx_<<"   "<<ky_<<"    "<<E_[0]<<"     "<<E_[1]<<"     "<<E_[2]<<endl;
        ky_= ky_ + dk_;
        counter_k++;
    }


    kx_=PI;
    ky_=PI;
    while(ky_>=0){

        kx_=ky_;

        H(0,0)=one*(-2.0*t2_*cos(kx_)  -2.0*t1_*cos(ky_)  -4.0*t3_*cos(kx_)*cos(ky_));
        H(1,1)=one*(-2.0*t1_*cos(kx_) -2.0*t2_*cos(ky_) - 4.0*t3_*cos(kx_)*cos(ky_));
        H(2,2)=one*(-2.0*t5_*(cos(kx_) + cos(ky_)) - 4.0*t6_*cos(kx_)*cos(ky_) + Delta_xy_);
        H(0,1)=one*(-4.0*t4_*sin(kx_)*sin(ky_));

        H(0,2)=iota*(-2.0*t7_*sin(kx_) + 4.0*t8_*sin(kx_)*cos(ky_) );
        H(1,2)=iota*(-2.0*t7_*sin(ky_) + 4.0*t8_*sin(ky_)*cos(kx_) );

        char jobz=option;
        char uplo='U';
        int n=H.n_row();
        int lda=H.n_col();
        vector<complex<double>> work(3);
        vector<double> rwork(3*n);
        int info,lwork= -1;

        vector<double> E_;
        E_.resize(3);
        fill(E_.begin(), E_.end(),0);
        // query:
        zheev_(&jobz,&uplo,&n,&(H(0,0)),&lda,&(E_[0]),&(work[0]),&lwork,&(rwork[0]),&info);
        lwork = int(real(work[0]))+1;
        work.resize(lwork+1);
        // real work:
        zheev_(&jobz,&uplo,&n,&(H(0,0)),&lda,&(E_[0]),&(work[0]),&lwork,&(rwork[0]),&info);
        if (info!=0) {
            std::cerr<<"info="<<info<<"\n";
            perror("diag: zheev: failed with info!=0.\n");
        }

        file_out<<counter_k<<"    "<<kx_<<"   "<<ky_<<"    "<<E_[0]<<"     "<<E_[1]<<"     "<<E_[2]<<endl;
        ky_= ky_ - dk_;
        counter_k++;
    }


}


double Observables_MCMF::Lorentzian(double x, double brd){
    double temp;

    temp = (1.0/PI)*( (brd/2.0)/ ( (x*x) + ((brd*brd)/4.0) ) );

    return temp;

}

void Observables_MCMF::DensityOfStates(){
    //-----------Calculate Bandwidth------//
    BandWidth=2.0;
    //-----------------------------------//

} // ----------


void Observables_MCMF::OccDensity(){

} // ----------


void Observables_MCMF::TotalOccDensity(){

} // ----------



complex<double> Observables_MCMF::SiSjQ(int i,int j){return SiSjQ_(i,j);}

double Observables_MCMF::SiSj(int i,int j){return SiSj_(i,j);}


complex<double> Observables_MCMF::SiSjQ_Mean(int i,int j){return SiSjQ_Mean_(i,j);}

complex<double> Observables_MCMF::SiSjQ_square_Mean(int i,int j){return SiSjQ_square_Mean_(i,j);}


double Observables_MCMF::SiSj_Mean(int i,int j){return SiSj_Mean_(i,j);}

double Observables_MCMF::SiSj_square_Mean(int i,int j){return SiSj_square_Mean_(i,j);}


void Observables_MCMF::SiSjFULL(){

    double Cos_ij,Sin_ij,ei,ai,phase;
    int site_,site_p,ax,ay;


    for(int i=0;i<lx_;i++){
        for(int j=0;j<ly_;j++){
            site_ = Coordinates_.Nc(i,j);
            ei=MFParams_.etheta(i,j);
            ai=MFParams_.ephi(i,j);
            sx_[site_] = cos(ai) * sin(ei);
            sy_[site_] = sin(ai) * sin(ei);
            sz_[site_] = cos(ei);
        }
    }

    for(int xr=0;xr<lx_;xr++){
        for(int yr=0;yr<ly_;yr++){
            SiSj_(xr,yr)=double(0.0);
            for(int i=0;i<lx_;i++){
                for(int j=0;j<ly_;j++){
                    site_ = Coordinates_.Nc(i,j);
                    ax = (i+xr)%lx_;
                    ay = (j+yr)%ly_;
                    site_p = Coordinates_.Nc(ax,ay);
                    SiSj_(xr,yr) += sx_[site_]*sx_[site_p];
                    SiSj_(xr,yr) += sy_[site_]*sy_[site_p];
                    SiSj_(xr,yr) += sz_[site_]*sz_[site_p];
                }
            }
            SiSj_(xr,yr)*= double(1.0/(lx_*ly_));
            //cout << xr << " "<< yr<< " "<<  SiSj_(xr,yr) << endl;
        }
    }

    for(int qx=0; qx<lx_; qx++) {
        for(int qy=0; qy<ly_; qy++) {
            SiSjQ_(qx,qy)=complex<double>(0.0,0.0);
            for(int xr=0;xr<lx_;xr++){
                for(int yr=0;yr<ly_;yr++){
                    phase=2.0*Parameters_.pi*(double(qx*xr)/double(lx_)+double(qy*yr)/double(ly_));
                    Cos_ij = cos(phase);
                    Sin_ij = sin(phase);
                    SiSjQ_(qx,qy) += SiSj_(xr,yr)*complex<double>(Cos_ij,Sin_ij);
                }
            }
            SiSjQ_(qx,qy)*= double(1.0/(lx_*ly_));
            //cout << qx << " "<< qy<< " "<<  SiSjQ_(qx,qy) << endl;
        }
    }

    //     cout << 0 << " "<< 1 << " "<<  SiSj_(0,1) << endl;
    //     cout << 1 << " "<< 0 << " "<<  SiSj_(1,0) << endl;
    //     cout << 0 << " "<< 4 << " "<<  SiSjQ_(0,4) << endl;
    //     cout << 4 << " "<< 0 << " "<<  SiSjQ_(4,0) << endl;
    //     cout << 2 << " "<< 6 << " "<<  SiSjQ_(2,6) << endl;
    //     cout << 6 << " "<< 2 << " "<<  SiSjQ_(6,2) << endl;

} // ----------


void Observables_MCMF::OccDensity(int tlabel){

} // ----------


void Observables_MCMF::DOSprint(int tlabel){
    double omega;
    //create name
    std::string name="Output/OrbDOS_" + to_string(tlabel) + ".dat";
    ofstream myfile;
    myfile.open(name);
    for(int ll=0;ll<=800;ll++) {
        omega=Omega(ll);
        myfile << omega-Parameters_.mus << "\t"
               << setw(12) << dos(0,ll)/(ns_*Parameters_.MCNorm) << "\t"
               << setw(12) << dos(1,ll)/(ns_*Parameters_.MCNorm) << "\t"
               << setw(12) << dos(2,ll)/(ns_*Parameters_.MCNorm) << "\t"
               << setw(12) << dos(3,ll)/(ns_*Parameters_.MCNorm) << "\t"
               << endl;
    }
    myfile.close();
} // ----------



void Observables_MCMF::Initialize(){

    int space=2*Parameters_.orbs*ns_;
    sx_.resize(space);
    sy_.resize(space);
    sz_.resize(space);

    dos.resize(4,801);   dos.fill(0.0);

    SiSj_.resize(lx_,ly_);

    SiSjQ_.resize(lx_,ly_);

    dosincr_=0.05;
    tpi_=4.0f*atan(1.0f);
} // ----------


double Observables_MCMF::Omega(int i){
    return -20.0+double(i)*dosincr_;
} // ----------









