#include "Spin_dynamics_ED_engine.h"
#include <stdlib.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;



void SC_SW_ENGINE_ED::Initialize_engine(BASIS_1_orb_SF & Basis, MODEL_1_orb_SF & Model){



    RK4_type="type_1";   //Not used
    final_conf_out="final_config.txt"; //Not used

    //all above read from input file

    n_wpoints=(int) ((w_max-w_min)/dw + 0.5);

    time_steps=(int) (time_max/dt_ + 0.5);


    Model.Theta.resize(Runge_Kutta_order+1);
    Model.Phi.resize(Runge_Kutta_order+1);


    Model.Theta_eq.resize(Basis.Max_pos + 1);
    Model.Phi_eq.resize(Basis.Max_pos + 1);

    for(int i=0;i<Model.Theta.size();i++){

        Model.Theta[i].resize(Basis.Max_pos + 1);
        Model.Phi[i].resize(Basis.Max_pos + 1);
    }

    S_kw.resize(Basis.Max_pos + 1 );
    for(int i=0;i<=Basis.Max_pos;i++){
        S_kw[i].resize(n_wpoints);

    }
}

void SC_SW_ENGINE_ED::Read_equilibrium_configuration(BASIS_1_orb_SF & Basis, MODEL_1_orb_SF & Model){





    ifstream file_in(conf_input.c_str());

    double temp_pos,temp_theta,temp_phi;
    //Later on read from file
    for(int x=0;x<Basis.Lx;x++){
        for(int y=0;y<Basis.Ly;y++){

            //int pos= y*Basis.Lx + x;

            file_in>>temp_pos>>temp_theta>>temp_phi;
            Model.Theta[0][temp_pos-1]=temp_theta*PI;
            Model.Phi[0][temp_pos-1]=temp_phi*PI;

            //cout<<temp_pos<<"  "<<temp_theta<<"   "<<temp_phi<<endl;

            Model.Theta_eq[temp_pos-1]=temp_theta*PI;
            Model.Phi_eq[temp_pos-1]=temp_phi*PI;

        }

    }







}


void SC_SW_ENGINE_ED::Start_Engine(BASIS_1_orb_SF & Basis,MODEL_1_orb_SF & Model){
#ifdef _OPENMP
    omp_set_num_threads(no_of_processors);
    int N_p = omp_get_max_threads();
    cout<<N_p<<" processors are used parallely"<<endl;
#endif

    bool use_only_allowed_k = true;
    string mu_output = "mu.txt";

    complex<double> zero(0,0);
    complex<double> one(1,0);
    complex<double> iota(0,1);
    complex<double> temp;

    double kx, ky;
    double dk_ = 1.0/64.0;

    S_rw.resize(Basis.Max_pos +1);
    for(int pos_i=0;pos_i<Basis.Max_pos +1;pos_i++){
        S_rw[pos_i].resize(Basis.Max_pos +1);
        for(int pos_j=0;pos_j<Basis.Max_pos +1;pos_j++){
            S_rw[pos_i][pos_j].resize(n_wpoints);
            for(int wi=0;wi<n_wpoints;wi++){
                S_rw[pos_i][pos_j][wi]=zero;
            }
        }
    }


    s_quantum_rw.resize(Basis.Max_pos +1);
    for(int pos_i=0;pos_i<Basis.Max_pos +1;pos_i++){
        s_quantum_rw[pos_i].resize(Basis.Max_pos +1);
        for(int pos_j=0;pos_j<Basis.Max_pos +1;pos_j++){
            s_quantum_rw[pos_i][pos_j].resize(n_wpoints);
            for(int wi=0;wi<n_wpoints;wi++){
                s_quantum_rw[pos_i][pos_j][wi]=zero;
            }
        }
    }


    T_rw.resize(Basis.Max_pos +1);
    for(int pos_i=0;pos_i<Basis.Max_pos +1;pos_i++){
        T_rw[pos_i].resize(Basis.Max_pos +1);
        for(int pos_j=0;pos_j<Basis.Max_pos +1;pos_j++){
            T_rw[pos_i][pos_j].resize(n_wpoints);
            for(int wi=0;wi<n_wpoints;wi++){
                T_rw[pos_i][pos_j][wi]=zero;
            }
        }
    }


    ofstream file_out(spins_r_t_out.c_str());
    file_out<<"# Sz------, Sx-----,Sy-----"<<endl;

    ofstream file_mu_out(mu_output.c_str());



    for(int ts=0;ts<time_steps -1;ts++){

        //create Hamiltonian for time_step=ts
        //if(ts!=0){
        Model.Create_Hamil(Basis, 0);
        //}
        //else{
        //Model.Create_Hamil(Basis, 0, "biased_quantum_spins");}

        Model.Diagonalize_Hamil();

        Model.Mu_=Model.Calculate_mu(Model.Evals);

        file_mu_out<<ts<<"     "<<Model.Mu_<<endl;

        Model.Get_quantum_Spins(0);

        if(ts==0){
            Model.quant_s_x_eq = Model.quant_s_x[0];
            Model.quant_s_y_eq = Model.quant_s_y[0];
            Model.quant_s_z_eq = Model.quant_s_z[0];
        }

        Evolve_classical_spins(0,Model,Basis);


        file_out<<(ts*dt_)<<"  ";
        for(int pos=0;pos<=Basis.Max_pos;pos++){
            file_out<<cos(Model.Theta[0][pos])<<"  ";
        }
        for(int pos=0;pos<=Basis.Max_pos;pos++){
            file_out<<sin(Model.Theta[0][pos])*cos(Model.Phi[0][pos])<<"  ";
        }
        for(int pos=0;pos<=Basis.Max_pos;pos++){
            file_out<<sin(Model.Theta[0][pos])*sin(Model.Phi[0][pos])<<"  ";
        }

        for(int pos=0;pos<=Basis.Max_pos;pos++){
            file_out<<Model.quant_s_z[0][pos]<<"  ";
        }
        for(int pos=0;pos<=Basis.Max_pos;pos++){
            file_out<<Model.quant_s_x[0][pos]<<"  ";
        }
        for(int pos=0;pos<=Basis.Max_pos;pos++){
            file_out<<Model.quant_s_y[0][pos]<<"  ";
        }

        file_out<<endl;



#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
        for(int wi=0;wi<n_wpoints;wi++){

            for(int pos_i=0;pos_i<Basis.Max_pos +1;pos_i++){
                for(int pos_j=0;pos_j<Basis.Max_pos +1;pos_j++){



                    /*               S_rw[pos_i][pos_j][wi] += one*cos((wi * dw) * (ts* dt_))*exp(-0.5*(ts*dt_*w_conv*ts*dt_*w_conv))*dt_*Model.S_mag*Model.S_mag*(  (cos(Model.Theta_eq[pos_i])*cos(Model.Theta[0][pos_j]))
                            +
                            (sin(Model.Theta_eq[pos_i])*cos(Model.Phi_eq[pos_i])*sin(Model.Theta[0][pos_j])*cos(Model.Phi[0][pos_j]) )
                            +
                            (sin(Model.Theta_eq[pos_i])*sin(Model.Phi_eq[pos_i])*sin(Model.Theta[0][pos_j])*sin(Model.Phi[0][pos_j]) )
                            );
*/

                    S_rw[pos_i][pos_j][wi] += exp(iota * (-wi * dw) * (ts* dt_))*exp(-0.5*(ts*dt_*w_conv*ts*dt_*w_conv))*dt_*Model.S_mag*Model.S_mag*(  (cos(Model.Theta_eq[pos_i])*cos(Model.Theta[0][pos_j]))
                            +
                            (sin(Model.Theta_eq[pos_i])*cos(Model.Phi_eq[pos_i])*sin(Model.Theta[0][pos_j])*cos(Model.Phi[0][pos_j]) )
                            +
                            (sin(Model.Theta_eq[pos_i])*sin(Model.Phi_eq[pos_i])*sin(Model.Theta[0][pos_j])*sin(Model.Phi[0][pos_j]) )
                            );



                    s_quantum_rw[pos_i][pos_j][wi] += exp(iota * (-wi * dw) * (ts* dt_))*exp(-0.5*(ts*dt_*w_conv*ts*dt_*w_conv))*dt_*

                            ( (Model.quant_s_x_eq[pos_i]*Model.quant_s_x[0][pos_j])
                            +
                            (Model.quant_s_y_eq[pos_i]*Model.quant_s_y[0][pos_j])
                            +
                            (Model.quant_s_z_eq[pos_i]*Model.quant_s_z[0][pos_j])


                            );


                    T_rw[pos_i][pos_j][wi] +=exp(iota * (-wi * dw) * (ts* dt_))*exp(-0.5*(ts*dt_*w_conv*ts*dt_*w_conv))*dt_*  ( Model.S_mag*(
                                                                                                                                    (
                                                                                                                                        (cos(Model.Theta_eq[pos_i])*Model.quant_s_z[0][pos_j])
                                                                                                                                    +
                                                                                                                                    (sin(Model.Theta_eq[pos_i])*cos(Model.Phi_eq[pos_i])*Model.quant_s_x[0][pos_j])
                                                                                                                                +
                                                                                                                                (sin(Model.Theta_eq[pos_i])*sin(Model.Phi_eq[pos_i])*Model.quant_s_y[0][pos_j])

                            )
                            +
                            (
                                (Model.quant_s_z_eq[pos_i]*cos(Model.Theta[0][pos_j]))
                            +
                            (Model.quant_s_x_eq[pos_i]*sin(Model.Theta[0][pos_j])*cos(Model.Phi[0][pos_j]) )
                            +
                            (Model.quant_s_y_eq[pos_i]*sin(Model.Theta[0][pos_j])*sin(Model.Phi[0][pos_j]) )

                            )


                            )

                            +


                            (
                            (Model.quant_s_x_eq[pos_i]*Model.quant_s_x[0][pos_j])
                            +
                            (Model.quant_s_y_eq[pos_i]*Model.quant_s_y[0][pos_j])
                            +
                            (Model.quant_s_z_eq[pos_i]*Model.quant_s_z[0][pos_j])


                            )

                            +

                            (
                            Model.S_mag*Model.S_mag*(  (cos(Model.Theta_eq[pos_i])*cos(Model.Theta[0][pos_j]))
                            +
                            (sin(Model.Theta_eq[pos_i])*cos(Model.Phi_eq[pos_i])*sin(Model.Theta[0][pos_j])*cos(Model.Phi[0][pos_j]) )
                            +
                            (sin(Model.Theta_eq[pos_i])*sin(Model.Phi_eq[pos_i])*sin(Model.Theta[0][pos_j])*sin(Model.Phi[0][pos_j]) )
                            )

                            )



                            );



                    /*  S_rw[pos_i][pos_j][wi] += exp(iota * (wi * dw) * (ts* dt_))*dt_*(  cos(Model.Theta[0][pos_i])*cos(Model.Theta[0][pos_j])
                                +
                                sin(Model.Theta[0][pos_i])*cos(Model.Phi[0][pos_i])*sin(Model.Theta[0][pos_j])*cos(Model.Phi[0][pos_j])
                                +
                                sin(Model.Theta[0][pos_i])*sin(Model.Phi[0][pos_i])*sin(Model.Theta[0][pos_j])*sin(Model.Phi[0][pos_j])
                                );

                        */




                }
            }
        }




        Model.Theta[0]=Model.Theta[1];
        Model.Phi[0]=Model.Phi[1];
    }







    ofstream file_out_full(Skw_out_full.c_str());


    ofstream file_out2(Skw_out.c_str());
    file_out2<<"# kx (kx=2*PI*nx/Lx), ky, n(ny*Lx + nx),  w_value, w_no, Skw.real, Skw.imag , skw.real, skw.imag, Tkw.real, Tkw.imag"<<endl;
    file_out_full<<"# kx (kx=2*PI*nx/Lx), ky, n(ny*Lx + nx),  w_value, w_no, Skw.real, Skw.imag , skw.real, skw.imag, Tkw.real, Tkw.imag"<<endl;

    complex<double> temp2, temp3;
    int pos_i, pos_j;
    int nx, ny;
    int k_ind;
    int Lxby2 = (int) ( ((1.0*Basis.Lx)/2.0) + 0.5);
    int Lyby2 = (int) ( ((1.0*Basis.Ly)/2.0) + 0.5);
    int Lxby4 = (int) ( ((1.0*Basis.Lx)/4.0) + 0.5);
    int Lyby4 = (int) ( ((1.0*Basis.Ly)/4.0) + 0.5);

    k_ind=0;

    Mat_3_Complex_doub Skw_Mat;
    Skw_Mat.resize(2*Lxby2+1);

    Mat_3_Complex_doub squant_kw_Mat;
    squant_kw_Mat.resize(2*Lxby2+1);

    Mat_3_Complex_doub Tkw_Mat;
    Tkw_Mat.resize(2*Lxby2+1);

    for(int i=0;i<=2*Lxby2;i++){
        Skw_Mat[i].resize(2*Lyby2+1);
        squant_kw_Mat[i].resize(2*Lyby2+1);
        Tkw_Mat[i].resize(2*Lyby2+1);
        for(int j=0;j<=2*Lyby2;j++){
            Skw_Mat[i][j].resize(n_wpoints);
            squant_kw_Mat[i][j].resize(n_wpoints);
            Tkw_Mat[i][j].resize(n_wpoints);
        }
    }

    if(use_only_allowed_k != true){

        kx=-(PI*dk_);

        while(kx <  (2*PI + (PI*dk_))  ){
            ky=-(PI*dk_);
            while(ky < (2*PI + (PI*dk_)) ){

                // kx = (2*nx*PI)/(1.0*Basis.Lx);
                // ky = (2*ny*PI)/(1.0*Basis.Ly);


#ifdef _OPENMP
#pragma omp parallel for default(shared)  private(pos_i,pos_j,temp,temp2,temp3)
#endif

                for(int wi=0;wi<n_wpoints;wi++){
                    temp3=zero;
                    temp2=zero;
                    temp=zero;
                    for(int x_i=0;x_i<Basis.Lx;x_i++){
                        for(int y_i=0;y_i<Basis.Ly;y_i++){
                            pos_i = y_i*Basis.Lx + x_i;
                            for(int x_j=0;x_j<Basis.Lx;x_j++){
                                for(int y_j=0;y_j<Basis.Ly;y_j++){
                                    pos_j = y_j*Basis.Lx + x_j;

                                    temp += one*S_rw[pos_i][pos_j][wi]*cos(( (x_j - x_i)*kx +  (y_j - y_i)*ky ) );

                                    //                                   temp += S_rw[pos_i][pos_j][wi]*exp(iota*( (x_j - x_i)*kx +  (y_j - y_i)*ky ) );
                                    temp2 += one*s_quantum_rw[pos_i][pos_j][wi]*cos(( (x_j - x_i)*kx +  (y_j - y_i)*ky ) );
                                    temp3 += one*T_rw[pos_i][pos_j][wi]*cos(( (x_j - x_i)*kx +  (y_j - y_i)*ky ) ); ;

                                }
                            }
                        }
                    }
                    file_out_full<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<wi<<"    "<<temp.real()<<"   "<<temp.imag()<<"    "<<temp2.real()<<"   "<<temp2.imag()<<"    "<<temp3.real()<<"   "<<temp3.imag()<<endl;
                    //file_out_full<<kx<<"   "<<ky<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<wi<<"   "<<temp.real()<<"   "<<temp.imag()<<"    "<<"none"<<"   "<<"none"<<"    "<<"none"<<"   "<<"none"<<endl;
                }
                k_ind +=1;
                ky=ky+(PI*dk_);

                cout<<kx<<"   "<<ky<<endl;
            }

            kx=kx+(PI*dk_);
        }

    }

    else{
        for(nx=0;nx<=2*Lxby2;nx++){
            for(ny=0;ny<=2*Lyby2;ny++){

                kx = (2*nx*PI)/(1.0*Basis.Lx);
                ky = (2*ny*PI)/(1.0*Basis.Ly);

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(pos_i,pos_j,temp,temp2,temp3)
#endif

                for(int wi=0;wi<n_wpoints;wi++){
                    temp3=zero;
                    temp2=zero;
                    temp=zero;
                    for(int x_i=0;x_i<Basis.Lx;x_i++){
                        for(int y_i=0;y_i<Basis.Ly;y_i++){
                            pos_i = y_i*Basis.Lx + x_i;
                            for(int x_j=0;x_j<Basis.Lx;x_j++){
                                for(int y_j=0;y_j<Basis.Ly;y_j++){
                                    pos_j = y_j*Basis.Lx + x_j;

                                    temp += one*S_rw[pos_i][pos_j][wi]*cos(( (x_j - x_i)*kx +  (y_j - y_i)*ky ) );
                                    //temp += S_rw[pos_i][pos_j][wi]*exp(iota*( (x_j - x_i)*kx +  (y_j - y_i)*ky ) );
                                    temp2 += one*s_quantum_rw[pos_i][pos_j][wi]*cos(( (x_j - x_i)*kx +  (y_j - y_i)*ky ) );
                                    temp3 += one*T_rw[pos_i][pos_j][wi]*cos(( (x_j - x_i)*kx +  (y_j - y_i)*ky ) ) ;

                                }
                            }
                        }
                    }
                    //file_out2<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<temp.real()<<"   "<<temp.imag()<<"    "<<temp2.real()<<"   "<<temp2.imag()<<"    "<<temp3.real()<<"   "<<temp3.imag()<<endl;
                    //file_out_full<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<wi<<"   "<<temp.real()<<"   "<<temp.imag()<<"    "<<"none"<<"   "<<"none"<<"    "<<"none"<<"   "<<"none"<<endl;
                    Skw_Mat[nx][ny][wi]=temp;
                    squant_kw_Mat[nx][ny][wi]=temp2;
                    Tkw_Mat[nx][ny][wi]=temp3;
                }



            }

        }


        k_ind=0;
        for(nx=0;nx<=2*Lxby2;nx++){
            for(ny=0;ny<=2*Lyby2;ny++){

                for(int wi=0;wi<n_wpoints;wi++){

                    //file_out2<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<temp.real()<<"   "<<temp.imag()<<"    "<<temp2.real()<<"   "<<temp2.imag()<<"    "<<temp3.real()<<"   "<<temp3.imag()<<endl;
                    file_out_full<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<wi<<"   "<<Skw_Mat[nx][ny][wi].real()<<"   "<<Skw_Mat[nx][ny][wi].imag()<<"    "<<squant_kw_Mat[nx][ny][wi].real()<<"   "<<squant_kw_Mat[nx][ny][wi].imag()<<
                                   "     "<<Tkw_Mat[nx][ny][wi].real()<<"   "<<Tkw_Mat[nx][ny][wi].imag()<<endl;
                    //Skw_Mat[nx][ny][wi]=temp;
                }

                k_ind +=1;

            }

        }



    }



    if(use_only_allowed_k == true){

        k_ind=0;
        //From [(pi/2,pi/2)----to----(pi,pi) )
        for(nx=Lxby4;nx<Lxby2;nx++){
            ny=nx;
            for(int wi=0;wi<n_wpoints;wi++){

                //file_out2<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<temp.real()<<"   "<<temp.imag()<<"    "<<temp2.real()<<"   "<<temp2.imag()<<"    "<<temp3.real()<<"   "<<temp3.imag()<<endl;
                file_out2<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<wi<<"   "<<Skw_Mat[nx][ny][wi].real()<<"   "<<Skw_Mat[nx][ny][wi].imag()<<"    "<<squant_kw_Mat[nx][ny][wi].real()<<"   "<<squant_kw_Mat[nx][ny][wi].imag()<<
                           "     "<<Tkw_Mat[nx][ny][wi].real()<<"   "<<Tkw_Mat[nx][ny][wi].imag()<<endl;

            }
            file_out2<<endl;
            k_ind +=1;
        }




        //From [(pi,pi)----to----(pi,0) )
        nx = Lxby2;
        for(ny=Lyby2;ny>0;ny--){

            for(int wi=0;wi<n_wpoints;wi++){

                // file_out2<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<temp.real()<<"   "<<temp.imag()<<"   "<<temp2.real()<<"   "<<temp2.imag()   <<"    "<<temp3.real()<<"   "<<temp3.imag()<<   endl;
                file_out2<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<wi<<"   "<<Skw_Mat[nx][ny][wi].real()<<"   "<<Skw_Mat[nx][ny][wi].imag()<<"    "<<squant_kw_Mat[nx][ny][wi].real()<<"   "<<squant_kw_Mat[nx][ny][wi].imag()<<
                           "     "<<Tkw_Mat[nx][ny][wi].real()<<"   "<<Tkw_Mat[nx][ny][wi].imag()<<endl;
            }
            file_out2<<endl;
            k_ind +=1;
        }



        //From [(pi,0)----to----(pi/2,pi/2) )
        nx=Lxby2;
        ny=0;
        for(int i_n=0;i_n<Lxby4;i_n++){
            nx=nx - i_n;
            ny=ny + i_n;

            for(int wi=0;wi<n_wpoints;wi++){

                // file_out2<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<temp.real()<<"   "<<temp.imag()<<"   "<<temp2.real()<<"   "<<temp2.imag()   <<"    "<<temp3.real()<<"   "<<temp3.imag()<<endl;
                file_out2<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<wi<<"   "<<Skw_Mat[nx][ny][wi].real()<<"   "<<Skw_Mat[nx][ny][wi].imag()<<"    "<<squant_kw_Mat[nx][ny][wi].real()<<"   "<<squant_kw_Mat[nx][ny][wi].imag()<<
                           "     "<<Tkw_Mat[nx][ny][wi].real()<<"   "<<Tkw_Mat[nx][ny][wi].imag()<<endl;
            }
            file_out2<<endl;
            k_ind +=1;
        }

        //From [(pi/2,pi/2)----to----(0,0) )

        for(ny=Lyby4;ny>=0;ny--){

            nx=ny;

            for(int wi=0;wi<n_wpoints;wi++){

                // file_out2<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<temp.real()<<"   "<<temp.imag()<<"   "<<temp2.real()<<"   "<<temp2.imag()   <<"    "<<temp3.real()<<"   "<<temp3.imag()<<   endl;
                file_out2<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<wi<<"   "<<Skw_Mat[nx][ny][wi].real()<<"   "<<Skw_Mat[nx][ny][wi].imag()<<"    "<<squant_kw_Mat[nx][ny][wi].real()<<"   "<<squant_kw_Mat[nx][ny][wi].imag()<<
                           "     "<<Tkw_Mat[nx][ny][wi].real()<<"   "<<Tkw_Mat[nx][ny][wi].imag()<<endl;            }
            file_out2<<endl;
            k_ind +=1;
        }

    }
    else if(use_only_allowed_k != true){
        k_ind=0;
        //From [(pi/2,pi/2)----to----(pi,pi) )

        kx=PI/2 - PI*dk_;
        while(kx<PI){

            ky=kx;

            // kx = (2*nx*PI)/(1.0*Basis.Lx);
            // ky = (2*ny*PI)/(1.0*Basis.Ly);
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(pos_i,pos_j,temp,temp2,temp3)
#endif
            for(int wi=0;wi<n_wpoints;wi++){
                temp3=zero;
                temp2=zero;
                temp=zero;
                for(int x_i=0;x_i<Basis.Lx;x_i++){
                    for(int y_i=0;y_i<Basis.Ly;y_i++){
                        pos_i = y_i*Basis.Lx + x_i;
                        for(int x_j=0;x_j<Basis.Lx;x_j++){
                            for(int y_j=0;y_j<Basis.Ly;y_j++){
                                pos_j = y_j*Basis.Lx + x_j;

                                temp += one*S_rw[pos_i][pos_j][wi]*cos(( (x_j - x_i)*kx +  (y_j - y_i)*ky ) );
                                //temp += S_rw[pos_i][pos_j][wi]*exp(iota*( (x_j - x_i)*kx +  (y_j - y_i)*ky ) );
                                temp2 += one*s_quantum_rw[pos_i][pos_j][wi]*cos(( (x_j - x_i)*kx +  (y_j - y_i)*ky ) );
                                temp3 += one*T_rw[pos_i][pos_j][wi]*cos(( (x_j - x_i)*kx +  (y_j - y_i)*ky ) ) ;

                            }
                        }
                    }
                }
                //file_out2<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<temp.real()<<"   "<<temp.imag()<<"    "<<temp2.real()<<"   "<<temp2.imag()<<"    "<<temp3.real()<<"   "<<temp3.imag()<<endl;
                file_out2<<kx<<"   "<<ky<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<wi<<"   "<<temp.real()<<"   "<<temp.imag()<<"    "<<temp2.real()<<"   "<<temp2.imag()<<"    "<<temp3.real()<<"   "<<temp3.imag()<<endl;
            }
            file_out2<<endl;
            k_ind +=1;
            kx=kx + PI*dk_;
        }




        //From [(pi,pi)----to----(pi,0) )
        kx = PI;
        ky = PI;

        while(ky>0){

            //kx = (2*nx*PI)/(1.0*Basis.Lx);
            //ky = (2*ny*PI)/(1.0*Basis.Ly);
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(pos_i,pos_j,temp,temp2,temp3)
#endif
            for(int wi=0;wi<n_wpoints;wi++){
                temp3=zero;
                temp2=zero;
                temp=zero;
                for(int x_i=0;x_i<Basis.Lx;x_i++){
                    for(int y_i=0;y_i<Basis.Ly;y_i++){
                        pos_i = y_i*Basis.Lx + x_i;
                        for(int x_j=0;x_j<Basis.Lx;x_j++){
                            for(int y_j=0;y_j<Basis.Ly;y_j++){
                                pos_j = y_j*Basis.Lx + x_j;

                                temp += one*S_rw[pos_i][pos_j][wi]*cos(( (x_j - x_i)*kx +  (y_j - y_i)*ky ) );
                                //temp += S_rw[pos_i][pos_j][wi]*exp(iota*( (x_j - x_i)*kx +  (y_j - y_i)*ky ) );
                                temp2 += one*s_quantum_rw[pos_i][pos_j][wi]*cos(( (x_j - x_i)*kx +  (y_j - y_i)*ky ) );
                                temp3 += one*T_rw[pos_i][pos_j][wi]*cos(( (x_j - x_i)*kx +  (y_j - y_i)*ky ) ) ;
                            }
                        }
                    }
                }
                // file_out2<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<temp.real()<<"   "<<temp.imag()<<"   "<<temp2.real()<<"   "<<temp2.imag()   <<"    "<<temp3.real()<<"   "<<temp3.imag()<<   endl;
                file_out2<<kx<<"   "<<ky<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<wi<<"   "<<temp.real()<<"   "<<temp.imag()<<"    "<<temp2.real()<<"   "<<temp2.imag()<<"    "<<temp3.real()<<"   "<<temp3.imag()<<endl;
            }
            file_out2<<endl;
            k_ind +=1;
            ky = ky - PI*dk_;
        }



        //From [(pi,0)----to----(pi/2,pi/2) )
        kx = PI;
        ky = 0;
        while(ky<PI/2.0){
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(pos_i,pos_j,temp,temp2,temp3)
#endif

            for(int wi=0;wi<n_wpoints;wi++){
                temp3=zero;
                temp2=zero;
                temp=zero;
                for(int x_i=0;x_i<Basis.Lx;x_i++){
                    for(int y_i=0;y_i<Basis.Ly;y_i++){
                        pos_i = y_i*Basis.Lx + x_i;
                        for(int x_j=0;x_j<Basis.Lx;x_j++){
                            for(int y_j=0;y_j<Basis.Ly;y_j++){
                                pos_j = y_j*Basis.Lx + x_j;

                                temp += one*S_rw[pos_i][pos_j][wi]*cos(( (x_j - x_i)*kx +  (y_j - y_i)*ky ) );
                                //temp += S_rw[pos_i][pos_j][wi]*exp(iota*( (x_j - x_i)*kx +  (y_j - y_i)*ky ) );
                                temp2 += one*s_quantum_rw[pos_i][pos_j][wi]*cos(( (x_j - x_i)*kx +  (y_j - y_i)*ky ) );
                                temp3 += one*T_rw[pos_i][pos_j][wi]*cos(( (x_j - x_i)*kx +  (y_j - y_i)*ky ) ) ;

                            }
                        }
                    }
                }
                // file_out2<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<temp.real()<<"   "<<temp.imag()<<"   "<<temp2.real()<<"   "<<temp2.imag()   <<"    "<<temp3.real()<<"   "<<temp3.imag()<<endl;
                file_out2<<kx<<"   "<<ky<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<wi<<"   "<<temp.real()<<"   "<<temp.imag()<<"    "<<temp2.real()<<"   "<<temp2.imag()<<"    "<<temp3.real()<<"   "<<temp3.imag()<<endl;

            }
            file_out2<<endl;
            k_ind +=1;
            ky = ky + PI*dk_;
            kx = kx - PI*dk_;
        }

        //From [(pi/2,pi/2)----to----(0,0) )


        kx = PI/2.0;
        ky = PI/2.0;
        while(ky>=0){


#ifdef _OPENMP
#pragma omp parallel for default(shared) private(pos_i,pos_j,temp,temp2,temp3)
#endif
            for(int wi=0;wi<n_wpoints;wi++){
                temp3=zero;
                temp2=zero;
                temp=zero;
                for(int x_i=0;x_i<Basis.Lx;x_i++){
                    for(int y_i=0;y_i<Basis.Ly;y_i++){
                        pos_i = y_i*Basis.Lx + x_i;
                        for(int x_j=0;x_j<Basis.Lx;x_j++){
                            for(int y_j=0;y_j<Basis.Ly;y_j++){
                                pos_j = y_j*Basis.Lx + x_j;

                                temp += one*S_rw[pos_i][pos_j][wi]*cos(( (x_j - x_i)*kx +  (y_j - y_i)*ky ) );
                                //temp += S_rw[pos_i][pos_j][wi]*exp(iota*( (x_j - x_i)*kx +  (y_j - y_i)*ky ) );
                                temp2 += one*s_quantum_rw[pos_i][pos_j][wi]*cos(( (x_j - x_i)*kx +  (y_j - y_i)*ky ) );
                                temp3 += one*T_rw[pos_i][pos_j][wi]*cos(( (x_j - x_i)*kx +  (y_j - y_i)*ky ) ) ;
                            }
                        }
                    }
                }
                // file_out2<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<temp.real()<<"   "<<temp.imag()<<"   "<<temp2.real()<<"   "<<temp2.imag()   <<"    "<<temp3.real()<<"   "<<temp3.imag()<<   endl;
                file_out2<<kx<<"   "<<ky<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<wi<<"   "<<temp.real()<<"   "<<temp.imag()<<"    "<<temp2.real()<<"   "<<temp2.imag()<<"    "<<temp3.real()<<"   "<<temp3.imag()<<endl;
            }
            file_out2<<endl;
            k_ind +=1;
            kx = kx -PI*dk_;
            ky = ky -PI*dk_;
        }




    }
    /*
    for(int nx=0;nx<Basis.Lx;nx++){
        for(int ny=0;ny<Basis.Ly;ny++){
               kx = (2*nx*PI)/(1.0*Basis.Lx);
               ky = (2*ny*PI)/(1.0*Basis.Ly);

            for(int wi=0;wi<n_wpoints;wi++){
                temp=zero;
                for(int x_i=0;x_i<Basis.Lx;x_i++){
                    for(int y_i=0;y_i<Basis.Ly;y_i++){
                        pos_i = y_i*Basis.Lx + x_i;
                        for(int x_j=0;x_j<Basis.Lx;x_j++){
                            for(int y_j=0;y_j<Basis.Ly;y_j++){
                                pos_j = y_j*Basis.Lx + x_j;

                                temp += S_rw[pos_i][pos_j][wi]*exp(iota*( (x_j - x_i)*kx +  (y_j - y_i)*ky ) );


                            }
                        }
                    }
                }
                file_out2<<nx<<"   "<<ny<<"   "<<nx*Basis.Lx + ny<<"   "<<wi*dw<<"   "<<temp.real()<<"   "<<temp.imag()<<endl;
            }
            file_out2<<endl;
        }
    }
*/




}




void SC_SW_ENGINE_ED::Evolve_classical_spins(int ts, MODEL_1_orb_SF &Model,BASIS_1_orb_SF & Basis){


    if(Runge_Kutta_order==1){
        double pi,ti; //phi_i,theta_i,phi_j,theta_j for timestep=ts
        double sy,sx,sz; //Quatum spins for timestep=ts and position=pos

        double derivative_theta, derivative_phi;
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(derivative_theta,derivative_phi, pi, ti, sx, sy, sz)
#endif
        for(int pos=0;pos<=Basis.Max_pos;pos++){
            derivative_theta=0;
            derivative_phi=0;

            pi=Model.Phi[ts][pos];
            ti=Model.Theta[ts][pos];

            sy=Model.quant_s_y[ts][pos];
            sx=Model.quant_s_x[ts][pos];
            sz=Model.quant_s_z[ts][pos];


            Model.Phi[ts+1][pos] = Model.Phi[ts][pos] - (dt_*( Model.Jval*((((sin(pi)*cos(ti))/(sin(ti)))*sy) + (((cos(pi)*cos(ti))/(sin(ti)))*sx) - sz)));

            Model.Theta[ts+1][pos] = Model.Theta[ts][pos] + (dt_*( Model.Jval*((cos(pi)*sy) - (sin(pi)*sx))));

            //Model.Phi[ts+1][pos] = Model.Phi[ts][pos] ;
            //Model.Theta[ts+1][pos] = Model.Theta[ts][pos] ;

            derivative_phi += (-1.0)*( Model.Jval*((((sin(pi)*cos(ti))/(sin(ti)))*sy) + (((cos(pi)*cos(ti))/(sin(ti)))*sx) - sz));
            derivative_theta += ( Model.Jval*((cos(pi)*sy) - (sin(pi)*sx)));

            int pos_ng;
            double pj,tj;
            for(int ng=0;ng<Basis.Neighbour[pos].size();ng++){

                pos_ng = Basis.Neighbour[pos][ng];



                pj=Model.Phi[ts][pos_ng];
                tj=Model.Theta[ts][pos_ng];


                Model.Phi[ts+1][pos] += (-1.0)*(dt_*( (Model.S_mag)*(Model.Jval_p)* (((sin(tj)*cos(ti)*cos(pi - pj))/(sin(ti))) - cos(tj) )));
                derivative_phi += (-1.0)*(( (Model.S_mag)*(Model.Jval_p)*(((sin(tj)*cos(ti)*cos(pi - pj))/(sin(ti))) - cos(tj) ) ));

                Model.Theta[ts+1][pos] += (dt_*(Model.S_mag)*(Model.Jval_p)*(sin(tj)*sin(pj - pi)));
                derivative_theta += (Model.S_mag)*(Model.Jval_p)*(sin(tj)*sin(pj - pi));








            }

            if(Model.Phi[ts+1][pos] > 2*PI){
                Model.Phi[ts+1][pos] += -2*PI;

            }
            if(Model.Phi[ts+1][pos] < 0){
                Model.Phi[ts+1][pos] +=  2*PI;

            }


            if(Model.Theta[ts+1][pos] > PI){
                Model.Theta[ts+1][pos] = Model.Theta[ts+1][pos] -2*PI;
                Model.Phi[ts+1][pos] = fmod( Model.Phi[ts+1][pos] + PI, 2.0*PI );
            }
            if(Model.Theta[ts+1][pos] < 0){
                Model.Theta[ts+1][pos] = - Model.Theta[ts+1][pos];
                Model.Phi[ts+1][pos] = fmod( Model.Phi[ts+1][pos] + PI, 2.0*PI );

            }


            // cout<<pos<<"  "<<derivative_theta<<"  "<<derivative_phi<<endl;

        }
    }

    if(Runge_Kutta_order==4){
        /*
         we let the quantum spins change as well in the intermediate steps of dt interval in RK-4 method,
         here we have to diagonalize the matrix 3 extra times to calculate quantum spins for the changed classical spins.
         */

        double pi,ti; //phi_i,theta_i,phi_j,theta_j for timestep=ts
        double sy,sx,sz; //Quatum spins for timestep=ts and position=pos

        double derivative_theta, derivative_phi;
        Mat_1_doub delta1_theta, delta2_theta, delta3_theta, delta4_theta;
        Mat_1_doub delta1_phi, delta2_phi, delta3_phi, delta4_phi;
        Mat_1_doub quant_s_y_new, quant_s_x_new, quant_s_z_new;

        Mat_2_Complex_doub Hamil_Temp;
        Mat_1_doub Evals_Temp;
        Mat_2_Complex_doub Eigvecs_Temp;

        //Standard RK-4 method Convention is used
        delta1_theta.resize(Basis.Max_pos + 1);delta2_theta.resize(Basis.Max_pos + 1);
        delta3_theta.resize(Basis.Max_pos + 1);delta4_theta.resize(Basis.Max_pos + 1);

        delta1_phi.resize(Basis.Max_pos + 1);delta2_phi.resize(Basis.Max_pos + 1);
        delta3_phi.resize(Basis.Max_pos + 1);delta4_phi.resize(Basis.Max_pos + 1);


        for(int pos=0;pos<=Basis.Max_pos;pos++){
            delta1_theta[pos]=0.0;delta2_theta[pos]=0.0;
            delta3_theta[pos]=0.0;delta4_theta[pos]=0.0;

            delta1_phi[pos]=0.0;delta2_phi[pos]=0.0;
            delta3_phi[pos]=0.0;delta4_phi[pos]=0.0;

        }



        //calculating delta1(2,3,4)_thetha(phi)
        bool intermediate_update = true;
        double mu;
        for(int step_no=0;step_no<4;step_no++){

            //UPDATE QUANTUM SPINS: HERE

            if(step_no==0){
                quant_s_y_new=Model.quant_s_y[ts];
                quant_s_x_new=Model.quant_s_x[ts];
                quant_s_z_new=Model.quant_s_z[ts];
            }

            if(step_no==1){
                Model.Create_Hamil(Basis, Model.Phi[ts],
                                   delta1_phi, Model.Theta[ts],
                                   delta1_theta, Hamil_Temp, 0.5);
                Model.Diagonalize_Hamil(Hamil_Temp, Evals_Temp, Eigvecs_Temp);

                mu=Model.Calculate_mu(Evals_Temp);
                Model.Get_quantum_Spins(quant_s_y_new, quant_s_x_new, quant_s_z_new,
                                        Hamil_Temp, Eigvecs_Temp, Evals_Temp,mu);


            }
            if(step_no==2){
                Model.Create_Hamil(Basis, Model.Phi[ts],
                                   delta2_phi, Model.Theta[ts],
                                   delta2_theta, Hamil_Temp, 0.5);
                Model.Diagonalize_Hamil(Hamil_Temp, Evals_Temp, Eigvecs_Temp);

                mu=Model.Calculate_mu(Evals_Temp);
                Model.Get_quantum_Spins(quant_s_y_new, quant_s_x_new, quant_s_z_new,
                                        Hamil_Temp, Eigvecs_Temp, Evals_Temp,mu);

            }
            if(step_no==3){
                Model.Create_Hamil(Basis, Model.Phi[ts],
                                   delta3_phi, Model.Theta[ts],
                                   delta3_theta, Hamil_Temp, 1.0);
                Model.Diagonalize_Hamil(Hamil_Temp, Evals_Temp, Eigvecs_Temp);

                mu=Model.Calculate_mu(Evals_Temp);
                Model.Get_quantum_Spins(quant_s_y_new, quant_s_x_new, quant_s_z_new,
                                        Hamil_Temp, Eigvecs_Temp, Evals_Temp,mu);


            }
            Hamil_Temp.clear();Eigvecs_Temp.clear();Evals_Temp.clear();


#ifdef _OPENMP
#pragma omp parallel for default(shared) private(derivative_theta,derivative_phi, pi, ti, sx, sy, sz)
#endif
            for(int pos=0;pos<=Basis.Max_pos;pos++){
                derivative_theta=0;
                derivative_phi=0;

                if(step_no==0){
                    pi=Model.Phi[ts][pos];
                    ti=Model.Theta[ts][pos];
                }
                if(step_no==1){
                    pi=Model.Phi[ts][pos] + 0.5*(delta1_phi[pos]);
                    ti=Model.Theta[ts][pos] + 0.5*(delta1_theta[pos]);

                    if(intermediate_update==true){
                        if(pi > 2*PI){
                            pi += -2*PI;

                        }
                        if(pi < 0){
                            pi +=  2*PI;

                        }
                        if(ti > PI){
                            ti = ti -2*PI;
                            pi = fmod(pi + PI, 2.0*PI );
                        }
                        if(ti < 0){
                            ti = - ti;
                            pi = fmod(pi + PI, 2.0*PI );

                        }
                    }
                }

                if(step_no==2){
                    pi=Model.Phi[ts][pos] + 0.5*(delta2_phi[pos]);
                    ti=Model.Theta[ts][pos] + 0.5*(delta2_theta[pos]);

                    if(intermediate_update==true){
                        if(pi > 2*PI){
                            pi += -2*PI;

                        }
                        if(pi < 0){
                            pi +=  2*PI;

                        }
                        if(ti > PI){
                            ti = ti -2*PI;
                            pi = fmod(pi + PI, 2.0*PI );
                        }
                        if(ti < 0){
                            ti = - ti;
                            pi = fmod(pi + PI, 2.0*PI );

                        }
                    }
                }

                if(step_no==3){
                    pi=Model.Phi[ts][pos] + (delta3_phi[pos]);
                    ti=Model.Theta[ts][pos] + (delta3_theta[pos]);

                    if(intermediate_update==true){
                        if(pi > 2*PI){
                            pi += -2*PI;

                        }
                        if(pi < 0){
                            pi +=  2*PI;

                        }
                        if(ti > PI){
                            ti = ti -2*PI;
                            pi = fmod(pi + PI, 2.0*PI );
                        }
                        if(ti < 0){
                            ti = - ti;
                            pi = fmod(pi + PI, 2.0*PI );

                        }
                    }
                }




                sy=quant_s_y_new[pos];
                sx=quant_s_x_new[pos];
                sz=quant_s_z_new[pos];




                derivative_phi = (-1.0)*( Model.Jval*((((sin(pi)*cos(ti))/(sin(ti)))*sy) + (((cos(pi)*cos(ti))/(sin(ti)))*sx) - sz));
                derivative_theta = ( Model.Jval*((cos(pi)*sy) - (sin(pi)*sx)));

                if(step_no==0){
                    delta1_phi[pos] = (dt_*(derivative_phi ));
                    delta1_theta[pos] = (dt_*( derivative_theta));

                }
                if(step_no==1){
                    delta2_phi[pos] = (dt_*(derivative_phi ));
                    delta2_theta[pos] = (dt_*( derivative_theta));

                }
                if(step_no==2){
                    delta3_phi[pos] = (dt_*(derivative_phi ));
                    delta3_theta[pos] = (dt_*( derivative_theta));

                }
                if(step_no==3){
                    delta4_phi[pos] = (dt_*(derivative_phi ));
                    delta4_theta[pos] = (dt_*( derivative_theta));

                }

                //Model.Phi[ts+1][pos] = Model.Phi[ts][pos] ;
                //Model.Theta[ts+1][pos] = Model.Theta[ts][pos] ;




                int pos_ng;
                double pj,tj;
                for(int ng=0;ng<Basis.Neighbour[pos].size();ng++){

                    pos_ng = Basis.Neighbour[pos][ng];


                    if(step_no==0){
                        pj=Model.Phi[ts][pos_ng];
                        tj=Model.Theta[ts][pos_ng];
                    }

                    if(step_no==1){
                        pj=Model.Phi[ts][pos_ng] + 0.5*(delta1_phi[pos_ng]);
                        tj=Model.Theta[ts][pos_ng] + 0.5*(delta1_theta[pos_ng]);

                        if(intermediate_update==true){
                            if(pj > 2*PI){
                                pj += -2*PI;

                            }
                            if(pj < 0){
                                pj +=  2*PI;

                            }
                            if(tj > PI){
                                tj = tj -2*PI;
                                pj = fmod(pj + PI, 2.0*PI );
                            }
                            if(tj < 0){
                                tj = - tj;
                                pj = fmod(pj + PI, 2.0*PI );

                            }
                        }
                    }

                    if(step_no==2){
                        pj=Model.Phi[ts][pos_ng] + 0.5*(delta2_phi[pos_ng]);
                        tj=Model.Theta[ts][pos_ng] + 0.5*(delta2_theta[pos_ng]);

                        if(intermediate_update==true){
                            if(pj > 2*PI){
                                pj += -2*PI;

                            }
                            if(pj < 0){
                                pj +=  2*PI;

                            }
                            if(tj > PI){
                                tj = tj -2*PI;
                                pj = fmod(pj + PI, 2.0*PI );
                            }
                            if(tj < 0){
                                tj = - tj;
                                pj = fmod(pj + PI, 2.0*PI );

                            }
                        }
                    }

                    if(step_no==3){
                        pj=Model.Phi[ts][pos_ng] + (delta3_phi[pos_ng]);
                        tj=Model.Theta[ts][pos_ng] + (delta3_theta[pos_ng]);

                        if(intermediate_update==true){
                            if(pj > 2*PI){
                                pj += -2*PI;

                            }
                            if(pj < 0){
                                pj +=  2*PI;

                            }
                            if(tj > PI){
                                tj = tj -2*PI;
                                pj = fmod(pj + PI, 2.0*PI );
                            }
                            if(tj < 0){
                                tj = - tj;
                                pj = fmod(pj + PI, 2.0*PI );

                            }
                        }
                    }

                    derivative_phi = (-1.0)*(( (Model.S_mag)*(Model.Jval_p)*(((sin(tj)*cos(ti)*cos(pi - pj))/(sin(ti))) - cos(tj) ) ));
                    derivative_theta = (Model.S_mag)*(Model.Jval_p)*(sin(tj)*sin(pj - pi));


                    if(step_no==0){
                        delta1_phi[pos] += (dt_*derivative_phi );
                        delta1_theta[pos] += (dt_*derivative_theta);
                    }

                    if(step_no==1){
                        delta2_phi[pos] += (dt_*derivative_phi );
                        delta2_theta[pos] += (dt_*derivative_theta);
                    }

                    if(step_no==2){
                        delta3_phi[pos] += (dt_*derivative_phi );
                        delta3_theta[pos] += (dt_*derivative_theta);
                    }

                    if(step_no==3){
                        delta4_phi[pos] += (dt_*derivative_phi );
                        delta4_theta[pos] += (dt_*derivative_theta);
                    }


                }

                // cout<<pos<<"  "<<derivative_theta<<"  "<<derivative_phi<<endl;




            }




        }

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
        for(int pos=0;pos<=Basis.Max_pos;pos++){


            Model.Phi[ts+1][pos] = Model.Phi[ts][pos] + (1.0/6.0)*(delta1_phi[pos] + 2.0*delta2_phi[pos] + 2.0*delta3_phi[pos] + delta4_phi[pos]);
            Model.Theta[ts+1][pos] = Model.Theta[ts][pos] + (1.0/6.0)*(delta1_theta[pos] + 2.0*delta2_theta[pos] + 2.0*delta3_theta[pos] + delta4_theta[pos]);


            /*
            Model.Phi[ts+1][pos] = Model.Phi[ts][pos] + (1.0/1.0)*(0.0*delta1_phi[pos] + 1.0*delta2_phi[pos] + 0.0*delta3_phi[pos] + 0.0*delta4_phi[pos]);
            Model.Theta[ts+1][pos] = Model.Theta[ts][pos] + (1.0/1.0)*(0.0*delta1_theta[pos] + 1.0*delta2_theta[pos] + 0.0*delta3_theta[pos] + 0.0*delta4_theta[pos]);
            */

            if(Model.Phi[ts+1][pos] > 2*PI){
                Model.Phi[ts+1][pos] += -2*PI;

            }
            if(Model.Phi[ts+1][pos] < 0){
                Model.Phi[ts+1][pos] +=  2*PI;

            }


            if(Model.Theta[ts+1][pos] > PI){
                Model.Theta[ts+1][pos] = Model.Theta[ts+1][pos] -2*PI;
                Model.Phi[ts+1][pos] = fmod( Model.Phi[ts+1][pos] + PI, 2.0*PI );
            }
            if(Model.Theta[ts+1][pos] < 0){
                Model.Theta[ts+1][pos] = - Model.Theta[ts+1][pos];
                Model.Phi[ts+1][pos] = fmod( Model.Phi[ts+1][pos] + PI, 2.0*PI );

            }

            //cout<<delta1_phi[pos]<<"   "<<delta2_phi[pos]<<"   "<<delta3_phi[pos]<<"   "<<delta4_phi[pos]<<"   "<<(1.0/6.0)*(delta1_phi[pos] + 2.0*delta2_phi[pos] + 2.0*delta3_phi[pos] + delta4_phi[pos])<<endl;
        }





    }


}



void SC_SW_ENGINE_ED::Read_parameters(string filename, MODEL_1_orb_SF &Model,BASIS_1_orb_SF & Basis){




    string beta_, Beta_ = "Beta = ";
    string s_mag_, S_mag_ = "S_mag = ";
    string t_hop_, T_hop_ = "t_hop = ";
    string jval_, Jval_ = "Jval (J_HUND) = ";
    string jval_p_, Jval_p_ = "Jval_p (J_AFM) = ";
    string n_total_, N_total_ = "N_total = ";


    string Geometry_ = "Geometry = ";
    string boundary_conditions_, Boundary_Conditions_ = "Boundary conditions = ";

    string lx_, Lx_ = "Lx = ";
    string ly_, Ly_ = "Ly = ";




    string w_min_, W_Min_ = "w_min = ";
    string w_max_, W_Max_ = "w_max = ";
    string dw_, dW_ = "dw = ";
    string w_conv_, W_conv_ = "w_convolution = ";
    string time_max_, Time_Max_ = "time_max = ";
    string dt__, dT_ = "dt = ";
    string runge_kutta_order_, Runge_Kutta_Order_ = "Runge_Kutta_order = ";
    string Conf_Input_ = "conf_input = ";
    string Spins_R_T_Out_ = "spins_r_t_out = ";
    string SKW_Out_ = "Skw_out = ";
    string SKW_Out_Full_ = "Skw_out_Full = ";

    string no_of_processors_, No_Of_Processors_ = "no_of_threads = ";


    int offset;
    string line;
    ifstream inputfile(filename.c_str());


    if(inputfile.is_open())
    {
        while(!inputfile.eof())
        {
            getline(inputfile,line);


            if ((offset = line.find(N_total_, 0)) != string::npos) {
                n_total_ = line.substr (offset + N_total_.length());		}

            if ((offset = line.find(Jval_p_, 0)) != string::npos) {
                jval_p_ = line.substr (offset + Jval_p_.length());		}

            if ((offset = line.find(W_conv_, 0)) != string::npos) {
                w_conv_ = line.substr (offset + W_conv_.length());		}

            if ((offset = line.find(Jval_, 0)) != string::npos) {
                jval_ = line.substr (offset + Jval_.length());		}

            if ((offset = line.find(T_hop_, 0)) != string::npos) {
                t_hop_ = line.substr (offset + T_hop_.length());		}

            if ((offset = line.find(S_mag_, 0)) != string::npos) {
                s_mag_ = line.substr (offset + S_mag_.length());		}

            if ((offset = line.find(Beta_, 0)) != string::npos) {
                beta_  = line.substr (offset + Beta_.length());		}

            if ((offset = line.find(Boundary_Conditions_, 0)) != string::npos) {
                boundary_conditions_ = line.substr (offset + Boundary_Conditions_.length());		}

            if ((offset = line.find(Geometry_, 0)) != string::npos) {
                Basis.Geometry = line.substr (offset + Geometry_.length());		}

            if ((offset = line.find(Lx_, 0)) != string::npos) {
                lx_ = line.substr (offset + Lx_.length());		}

            if ((offset = line.find(Ly_, 0)) != string::npos) {
                ly_ = line.substr (offset + Ly_.length());		}

            if ((offset = line.find(W_Min_, 0)) != string::npos) {
                w_min_ = line.substr (offset + W_Min_.length());		}

            if ((offset = line.find(W_Max_, 0)) != string::npos) {
                w_max_ = line.substr (offset + W_Max_.length());		}

            if ((offset = line.find(dW_, 0)) != string::npos) {
                dw_= line.substr (offset + dW_.length());		}

            if ((offset = line.find(Time_Max_, 0)) != string::npos) {
                time_max_ = line.substr (offset + Time_Max_.length());		}

            if ((offset = line.find(dT_, 0)) != string::npos) {
                dt__= line.substr (offset + dT_.length());		}

            if ((offset = line.find(Runge_Kutta_Order_, 0)) != string::npos) {
                runge_kutta_order_ = line.substr (offset + Runge_Kutta_Order_.length());		}

            if ((offset = line.find(Conf_Input_, 0)) != string::npos) {
                conf_input = line.substr (offset + Conf_Input_.length());		}

            if ((offset = line.find(Spins_R_T_Out_, 0)) != string::npos) {
                spins_r_t_out = line.substr (offset + Spins_R_T_Out_.length());		}

            if ((offset = line.find(SKW_Out_, 0)) != string::npos) {
                Skw_out = line.substr (offset + SKW_Out_.length());		}

            if ((offset = line.find(SKW_Out_Full_, 0)) != string::npos) {
                Skw_out_full = line.substr (offset + SKW_Out_Full_.length());		}

            if ((offset = line.find(No_Of_Processors_, 0)) != string::npos) {
                no_of_processors_ = line.substr (offset + No_Of_Processors_.length());		}

        }
        inputfile.close();
    }
    else
    {cout<<"Unable to open input file while in the Model class."<<endl;}


    w_min=atof(w_min_.c_str());
    w_max=atof(w_max_.c_str());
    dw=atof(dw_.c_str());
    w_conv=atof(w_conv_.c_str());
    time_max=atof(time_max_.c_str());
    dt_ =atof(dt__.c_str());
    Runge_Kutta_order=atoi(runge_kutta_order_.c_str());


    no_of_processors=atoi(no_of_processors_.c_str());


    cout<<"PARAMETERS::::::::"<<endl;
    cout<<"w_min = "<<w_min<<endl;
    cout<<"w_max = "<<w_max<<endl;
    cout<<"dw = "<<dw<<endl;
    cout<<"w_convolution = "<<w_conv<<endl;
    cout<<"time_max = "<<time_max<<endl;
    cout<<"dt = "<<dt_<<endl;
    cout<<"Runge_Kutta_order = "<<Runge_Kutta_order<<endl;

    if (boundary_conditions_ == "PBC"){
        Basis.PBC =true;
    }

    Basis.Lx=atoi(lx_.c_str());
    Basis.Ly=atoi(ly_.c_str());

    cout<<"Boundary condition = "<<boundary_conditions_<<endl;
    cout<<"Lx = "<<Basis.Lx<<endl;
    cout<<"Ly = "<<Basis.Ly<<endl;



    Model.Jval_p=atof(jval_p_.c_str());
    Model.Jval=atof(jval_.c_str());
    Model.t_hop=atof(t_hop_.c_str());
    Model.S_mag=atof(s_mag_.c_str());
    Model.Beta=atof(beta_.c_str());
    Model.N_total=atoi(n_total_.c_str());

    cout<<"Jval_p = "<<Model.Jval_p<<endl;
    cout<<"Jval = "<<Model.Jval<<endl;
    cout<<"t_hop = "<<Model.t_hop<<endl;
    cout<<"S_mag = "<<Model.S_mag<<endl;
    cout<<"Beta = "<<Model.Beta<<endl;
    cout<<"N_total = "<<Model.N_total<<endl;






}



































