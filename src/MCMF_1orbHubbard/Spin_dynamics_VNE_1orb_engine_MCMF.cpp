#include <stdlib.h>
#include "Spin_dynamics_VNE_1orb_engine_MCMF.h"
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;

void SC_SW_ENGINE_VNE_1orb_MCMF::Initialize_engine(){


    SelfConsistentEvolution = false;

    RK4_type="type_1";   //Not used
    final_conf_out="final_config.txt"; //Not used

    //all above read from input file

    n_wpoints=(int) ((w_max-w_min)/dw + 0.5);

    if(!RESTART){
        Restart_Time=0.0;
    }


    time_steps=(int) ((time_max - Restart_Time)/dt_ + 0.5);



    Theta.resize(4);
    Phi.resize(4);
    Moment_Size.resize(4);

    Theta_time.resize(4);
    Phi_time.resize(4);
    Red_Den_mat_time.resize(4);
    Aux_S_x.resize(4);
    Aux_S_y.resize(4);
    Aux_S_z.resize(4);



    Theta_eq.resize(Parameters_.lx);
    Phi_eq.resize(Parameters_.lx);
    for(int i=0;i<Parameters_.lx;i++){
        Theta_eq[i].resize(Parameters_.ly);
        Phi_eq[i].resize(Parameters_.ly);
    }

    for(int t=0;t<Theta.size();t++){
        Theta[t].resize(Parameters_.lx);
        Phi[t].resize(Parameters_.lx);
        Moment_Size[t].resize(Parameters_.lx);
        Theta_time[t].resize(Parameters_.lx);
        Phi_time[t].resize(Parameters_.lx);
        for(int i=0;i<Parameters_.lx;i++){
            Theta[t][i].resize(Parameters_.ly);
            Phi[t][i].resize(Parameters_.ly);
            Moment_Size[t][i].resize(Parameters_.ly);
            Theta_time[t][i].resize(Parameters_.ly);
            Phi_time[t][i].resize(Parameters_.ly);

        }
    }

    S_kw.resize(Parameters_.ns);
    for(int i=0;i<Parameters_.ns;i++){
        S_kw[i].resize(n_wpoints);

    }




    Red_Den_mat.resize(Parameters_.ns);


    quant_s_x.resize(1);
    quant_s_y.resize(1);
    quant_s_z.resize(1);
    quant_s_x_eq.resize(Parameters_.ns);
    quant_s_y_eq.resize(Parameters_.ns);
    quant_s_z_eq.resize(Parameters_.ns);
    quant_s_x[0].resize(Parameters_.ns);
    quant_s_y[0].resize(Parameters_.ns);
    quant_s_z[0].resize(Parameters_.ns);


    Aux_S_x_eq.resize(Parameters_.ns);
    Aux_S_y_eq.resize(Parameters_.ns);
    Aux_S_z_eq.resize(Parameters_.ns);
    for(int t=0;t<Aux_S_x.size();t++){
        Aux_S_x[t].resize(Parameters_.ns);
        Aux_S_y[t].resize(Parameters_.ns);
        Aux_S_z[t].resize(Parameters_.ns);
    }


    for(int i=0;i<Parameters_.ns;i++){
        Red_Den_mat[i].resize(2);
        for(int s=0;s<2;s++){
            Red_Den_mat[i][s].resize(Parameters_.ns);
            for(int j=0;j<Parameters_.ns;j++){
                Red_Den_mat[i][s][j].resize(2);
            }
        }
    }


    for (int ts=0;ts<Red_Den_mat_time.size();ts++){
        Red_Den_mat_time[ts].resize(Parameters_.ns);
        for(int i=0;i<Parameters_.ns;i++){
            Red_Den_mat_time[ts][i].resize(2);
            for(int s=0;s<2;s++){
                Red_Den_mat_time[ts][i][s].resize(Parameters_.ns);
                for(int j=0;j<Parameters_.ns;j++){
                    Red_Den_mat_time[ts][i][s][j].resize(2);
                }
            }
        }
    }


    center = int((Parameters_.ns)*0.5 + 0.5) - 1;

    Jval_array.resize(Parameters_.ns);
    Bval_array.resize(Parameters_.ns);

    for(int pos=0;pos<Parameters_.ns;pos++){
        Jval_array[pos]=Parameters_.U_COUL*(-2.0);
        Bval_array[pos]=0.0;
    }

    J_Classical = Parameters_.J_Classical;


}



void SC_SW_ENGINE_VNE_1orb_MCMF::Read_equilibrium_configuration(){


    ifstream file_in;
    file_in.open(conf_input.c_str());

    double temp_theta,temp_phi, temp_moment_size, temp_local_den;
    int temp_lx, temp_ly;
    string temp_string;
    getline(file_in, temp_string);


    for(int x=0;x<Parameters_.lx;x++){
        for(int y=0;y<Parameters_.lx;y++){

            file_in>>temp_lx;
            file_in>>temp_ly;
            file_in>>temp_theta;
            file_in>>temp_phi;
            file_in>>temp_moment_size;
            file_in>>temp_local_den;
            Theta[0][temp_lx][temp_ly]=temp_theta;
            Phi[0][temp_lx][temp_ly]=temp_phi;
            Moment_Size[0][temp_lx][temp_ly]=temp_moment_size;

            //cout<<temp_pos<<"  "<<temp_theta<<"   "<<temp_phi<<endl;

            Theta_eq[temp_lx][temp_ly]=temp_theta;
            Phi_eq[temp_lx][temp_ly]=temp_phi;

        }

    }


    cout<<"Inital classical spin configuration is read from given input file "<<conf_input<<endl;





}



void SC_SW_ENGINE_VNE_1orb_MCMF::Start_Engine(){
#ifdef _OPENMP
    omp_set_num_threads(no_of_processors);
    int N_p = omp_get_max_threads();
    cout<<"Max threads which can be used parallely = "<<N_p<<endl;
    cout<<"No. of threads used parallely = "<<no_of_processors<<endl;
#endif

    bool use_only_allowed_k = true;
    string mu_output = "mu.txt";

    complex<double> zero(0,0);
    complex<double> one(1,0);
    complex<double> iota(0,1);
    complex<double> temp;

    int px_,py_;
    double kx, ky;
    double dk_ = 1.0/64.0;

    S_rw.resize(Parameters_.ns);
    for(int pos_i=0;pos_i<Parameters_.ns;pos_i++){
        S_rw[pos_i].resize(Parameters_.ns);
        for(int pos_j=0;pos_j<Parameters_.ns;pos_j++){
            S_rw[pos_i][pos_j].resize(n_wpoints);
            for(int wi=0;wi<n_wpoints;wi++){
                S_rw[pos_i][pos_j][wi]=zero;
            }
        }
    }


    s_quantum_rw.resize(Parameters_.ns);
    for(int pos_i=0;pos_i<Parameters_.ns;pos_i++){
        s_quantum_rw[pos_i].resize(Parameters_.ns);
        for(int pos_j=0;pos_j<Parameters_.ns;pos_j++){
            s_quantum_rw[pos_i][pos_j].resize(n_wpoints);
            for(int wi=0;wi<n_wpoints;wi++){
                s_quantum_rw[pos_i][pos_j][wi]=zero;
            }
        }
    }


    ofstream file_out(spins_r_t_out.c_str());
    file_out<<"# Sz------, Sx-----,Sy-----"<<endl;
    file_out<<scientific<<setprecision(15);

    ofstream file_mu_out(mu_output.c_str());

    for(int ts=0;ts<=time_steps;ts++){


        if(SAVE && (ts==time_steps)){
            Write_final_time_result();
        }


        //create Hamiltonian for time_step=ts
        if(ts==0){

            if(!RESTART){

                MFParams_.Read_classical_DOFs(conf_input);//HERE
                Hamiltonian_.InteractionsCreate();
                //char flag='V';
                Hamiltonian_.Diagonalize('V');
                double mu_ = Hamiltonian_.chemicalpotential(0.5,Parameters_.Fill);
                file_mu_out<<ts<<"     "<<mu_<<endl;
                //Get_quantum_Spins(0);
                Observables_.Get_red_den_mat(Red_Den_mat, mu_);
                Red_Den_mat_time[0]=Red_Den_mat;
                cout<<"Hamiltonian is diagonalized at t=0 "<<endl;
            }
            else{
                cout<<"Restart is done"<<endl;
                Read_Restart_Data();
            }


            for(int pos=0;pos<Parameters_.ns;pos++){
                quant_s_x_eq[pos] = 0.5*real( Red_Den_mat[pos][1][pos][0] +
                        Red_Den_mat[pos][0][pos][1] );
                quant_s_y_eq[pos] = 0.5*imag( Red_Den_mat[pos][1][pos][0] -
                        Red_Den_mat[pos][0][pos][1] );
                quant_s_z_eq[pos] = 0.5*real( Red_Den_mat[pos][0][pos][0] -
                        Red_Den_mat[pos][1][pos][1] );

                quant_s_x[0][pos] = 0.5*real( Red_Den_mat[pos][1][pos][0] +
                        Red_Den_mat[pos][0][pos][1] );
                quant_s_y[0][pos] = 0.5*imag( Red_Den_mat[pos][1][pos][0] -
                        Red_Den_mat[pos][0][pos][1] );
                quant_s_z[0][pos] = 0.5*real( Red_Den_mat[pos][0][pos][0] -
                        Red_Den_mat[pos][1][pos][1] );

                px_ = Coordinates_.indx(pos);
                py_ = Coordinates_.indy(pos);


                if(SelfConsistentEvolution == false){
                    Aux_S_x_eq[pos] = Moment_Size[0][px_][py_]*sin(Theta[0][px_][py_])*cos(Phi[0][px_][py_]);
                    Aux_S_y_eq[pos] = Moment_Size[0][px_][py_]*sin(Theta[0][px_][py_])*sin(Phi[0][px_][py_]);
                    Aux_S_z_eq[pos] = Moment_Size[0][px_][py_]*cos(Theta[0][px_][py_]);

                    Aux_S_x[0][pos] = Moment_Size[0][px_][py_]*sin(Theta[0][px_][py_])*cos(Phi[0][px_][py_]);
                    Aux_S_y[0][pos] = Moment_Size[0][px_][py_]*sin(Theta[0][px_][py_])*sin(Phi[0][px_][py_]);
                    Aux_S_z[0][pos] = Moment_Size[0][px_][py_]*cos(Theta[0][px_][py_]);
                }
                else{
                    Aux_S_x_eq[pos] = quant_s_x_eq[pos];
                    Aux_S_y_eq[pos] = quant_s_y_eq[pos];
                    Aux_S_z_eq[pos] = quant_s_z_eq[pos];

                    Aux_S_x[0][pos] = quant_s_x[0][pos];
                    Aux_S_y[0][pos] = quant_s_y[0][pos];
                    Aux_S_z[0][pos] = quant_s_z[0][pos];

                }

            }



            /*file_out<<(ts*dt_)<<"  ";
            for(int pos=0;pos<=basis.Max_pos;pos++){
                file_out<<S_mag*cos(Theta[0][pos])<<"  ";
            }
            for(int pos=0;pos<=basis.Max_pos;pos++){
                file_out<<S_mag*sin(Theta[0][pos])*cos(Phi[0][pos])<<"  ";
            }
            for(int pos=0;pos<=basis.Max_pos;pos++){
                file_out<<S_mag*sin(Theta[0][pos])*sin(Phi[0][pos])<<"  ";
            }

            for(int pos=0;pos<=basis.Max_pos;pos++){
                file_out<< quant_s_z_eq[pos]<<"  ";
            }
            for(int pos=0;pos<=basis.Max_pos;pos++){
                file_out<< quant_s_x_eq[pos]<<"  ";
            }
            for(int pos=0;pos<=basis.Max_pos;pos++){
                file_out<< quant_s_y_eq[pos]<<"  ";
            }

            file_out<<endl;*/




        }
        else{
        }

        file_out<<(ts*dt_) + Restart_Time<<"  ";

        for(int pos=0;pos<Parameters_.ns;pos++){
            file_out<< Aux_S_z[0][pos]<<"  ";
        }
        for(int pos=0;pos<Parameters_.ns;pos++){
            file_out<< Aux_S_x[0][pos]<<"  ";
        }
        for(int pos=0;pos<Parameters_.ns;pos++){
            file_out<< Aux_S_y[0][pos]<<"  ";
        }

        for(int pos=0;pos<Parameters_.ns;pos++){
            file_out<< quant_s_z[0][pos]<<"  ";
        }
        for(int pos=0;pos<Parameters_.ns;pos++){
            file_out<< quant_s_x[0][pos]<<"  ";
        }
        for(int pos=0;pos<Parameters_.ns;pos++){
            file_out<< quant_s_y[0][pos]<<"  ";
        }

        file_out<<endl;



        if(Insitu_SpaceTimeFourier){
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
            for(int wi=0;wi<n_wpoints;wi++){

                for(int pos_i=0;pos_i<Parameters_.ns;pos_i++){
                    for(int pos_j=0;pos_j<Parameters_.ns;pos_j++){


                        S_rw[pos_i][pos_j][wi] += exp(iota * (-wi * dw) * (ts* dt_))*exp(-0.5*(ts*dt_*w_conv*ts*dt_*w_conv))*dt_*(
                                    ( ( Aux_S_z_eq[pos_i] *  (  1.0*Aux_S_z[0][pos_j] -  0.0*Aux_S_z_eq[pos_j]   )  )  )
                                +
                                (  ( Aux_S_x_eq[pos_i] *  (  1.0*Aux_S_x[0][pos_j] -  0.0*Aux_S_x_eq[pos_j]   )  )  )
                                +
                                (  ( Aux_S_y_eq[pos_i] *  (  1.0*Aux_S_y[0][pos_j] -  0.0*Aux_S_y_eq[pos_j]   )  ) )
                                )
                                ;


                        s_quantum_rw[pos_i][pos_j][wi] += exp(iota * (-wi * dw) * (ts* dt_))*exp(-0.5*(ts*dt_*w_conv*ts*dt_*w_conv))*dt_*

                                ( ((quant_s_x_eq[pos_i])*
                                   ((quant_s_x[0][pos_j]) - 0.0*(quant_s_x_eq[pos_j])    ) )
                                +
                                ((quant_s_y_eq[pos_i])*
                                 ((quant_s_y[0][pos_j]) - 0.0*(quant_s_y_eq[pos_j] )    ) )
                                +
                                ((quant_s_z_eq[pos_i])*
                                 ((quant_s_z[0][pos_j]) - 0.0*(quant_s_z_eq[pos_j] )    ) )
                                )
                                ;

                    }
                }
            }
        }



        if(Predictor_Corrector==false){
            Evolve_classical_spins_Runge_Kutta(0);

            for(int pos=0;pos<Parameters_.ns;pos++){
                quant_s_x[0][pos] = 0.5*real( Red_Den_mat[pos][1][pos][0] + Red_Den_mat[pos][0][pos][1] );
                quant_s_y[0][pos] = 0.5*imag( Red_Den_mat[pos][1][pos][0] - Red_Den_mat[pos][0][pos][1] );
                quant_s_z[0][pos] = 0.5*real( Red_Den_mat[pos][0][pos][0] - Red_Den_mat[pos][1][pos][1] );

                Aux_S_x[0][pos] = Aux_S_x[1][pos];
                Aux_S_y[0][pos] = Aux_S_y[1][pos];
                Aux_S_z[0][pos] = Aux_S_z[1][pos];

            }

        }
        else{
            if(ts<=2){
                Theta_time[3-ts]=Theta[0];
                Phi_time[3-ts]=Phi[0];

                Evolve_classical_spins_Runge_Kutta(0);

                for(int pos=0;pos<Parameters_.ns;pos++){
                    for(int orb =0;orb<Parameters_.orbs;orb++){
                        quant_s_x[0][pos] = 0.5*real( Red_Den_mat[pos][1][pos][0] + Red_Den_mat[pos][0][pos][1] );
                        quant_s_y[0][pos] = 0.5*imag( Red_Den_mat[pos][1][pos][0] - Red_Den_mat[pos][0][pos][1] );
                        quant_s_z[0][pos] = 0.5*real( Red_Den_mat[pos][0][pos][0] - Red_Den_mat[pos][1][pos][1] );
                    }
                }

                Theta[0]=Theta[1];
                Phi[0]=Phi[1];




            }
            else{
                Theta_time[0]=Theta[0];
                Phi_time[0]=Phi[0];

                Evolve_classical_spins_Predictor_Corrector();

                Theta[0]=Theta[1];
                Phi[0]=Phi[1];

                for(int step=0;step<3;step++){
                    Theta_time[step+1]=Theta_time[step];
                    Phi_time[step+1]=Phi_time[step];
                }
            }
        }
    }




    if(Insitu_SpaceTimeFourier){


        ofstream file_out_full(Skw_out_full.c_str());


        ofstream file_out2(Skw_out.c_str());
        file_out2<<"# kx (kx=2*PI*nx/Lx), ky, n(ny*Lx + nx),  w_value, w_no, Skw.real, Skw.imag , skw.real, skw.imag, Tkw.real, Tkw.imag"<<endl;
        file_out_full<<"# kx (kx=2*PI*nx/Lx), ky, n(ny*Lx + nx),  w_value, w_no, Skw.real, Skw.imag , skw.real, skw.imag, Tkw.real, Tkw.imag"<<endl;

        complex<double> temp2, temp3;
        int pos_i, pos_j;
        int nx, ny;
        int k_ind;
        int Lxby2 = (int) ( ((1.0*Parameters_.lx)/2.0) + 0.5);
        int Lyby2 = (int) ( ((1.0*Parameters_.ly)/2.0) + 0.5);
        int Lxby4 = (int) ( ((1.0*Parameters_.lx)/4.0) + 0.5);
        int Lyby4 = (int) ( ((1.0*Parameters_.ly)/4.0) + 0.5);

        k_ind=0;

        Mat_3_Complex_doub Skw_Mat;
        Skw_Mat.resize(Parameters_.lx);

        Mat_3_Complex_doub squant_kw_Mat;
        squant_kw_Mat.resize(Parameters_.lx);

        for(int i=0;i<Parameters_.lx;i++){
            Skw_Mat[i].resize(Parameters_.ly);
            squant_kw_Mat[i].resize(Parameters_.ly);
            for(int j=0;j<Parameters_.ly;j++){
                Skw_Mat[i][j].resize(n_wpoints);
                squant_kw_Mat[i][j].resize(n_wpoints);
            }
        }

        if(use_only_allowed_k != true){

            kx=-(PI*dk_);

            while(kx <  (2*PI + (PI*dk_))  ){
                ky=-(PI*dk_);
                while(ky < (2*PI + (PI*dk_)) ){

                    // kx = (2*nx*PI)/(1.0*Parameters_.lx);
                    // ky = (2*ny*PI)/(1.0*Parameters_.ly);


#ifdef _OPENMP
#pragma omp parallel for default(shared)  private(pos_i,pos_j,temp,temp2,temp3)
#endif
                    for(int wi=0;wi<n_wpoints;wi++){

                        temp3=zero;
                        temp2=zero;
                        temp=zero;
                        for(int x_i=0;x_i<Parameters_.lx;x_i++){
                            for(int y_i=0;y_i<Parameters_.ly;y_i++){
                                pos_i = y_i*Parameters_.lx + x_i;
                                for(int x_j=0;x_j<Parameters_.lx;x_j++){
                                    for(int y_j=0;y_j<Parameters_.ly;y_j++){
                                        pos_j = y_j*Parameters_.lx + x_j;

                                        temp += one*S_rw[pos_i][pos_j][wi]*exp(iota*(( (x_j - x_i)*kx +  (y_j - y_i)*ky ) ));

                                        //                                   temp += S_rw[pos_i][pos_j][wi]*exp(iota*( (x_j - x_i)*kx +  (y_j - y_i)*ky ) );
                                        temp2 += one*s_quantum_rw[pos_i][pos_j][wi]*exp(iota*(( (x_j - x_i)*kx +  (y_j - y_i)*ky ) ));

                                    }
                                }
                            }
                        }
                        file_out_full<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<wi<<"    "<<temp.real()<<"   "<<temp.imag()<<"    "<<temp2.real()<<"   "<<temp2.imag()<<"    "<<endl;
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
            for(nx=0;nx<Parameters_.lx;nx++){
                for(ny=0;ny<Parameters_.ly;ny++){

                    kx = (2*nx*PI)/(1.0*Parameters_.lx);
                    ky = (2*ny*PI)/(1.0*Parameters_.ly);

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(pos_i,pos_j,temp,temp2,temp3)
#endif

                    for(int wi=0;wi<n_wpoints;wi++){
                        temp3=zero;
                        temp2=zero;
                        temp=zero;
                        for(int x_i=0;x_i<Parameters_.lx;x_i++){
                            for(int y_i=0;y_i<Parameters_.ly;y_i++){
                                pos_i = y_i*Parameters_.lx + x_i;
                                for(int x_j=0;x_j<Parameters_.lx;x_j++){
                                    for(int y_j=0;y_j<Parameters_.ly;y_j++){
                                        pos_j = y_j*Parameters_.lx + x_j;

                                        temp += one*S_rw[pos_i][pos_j][wi]*exp(iota*( (x_j - x_i)*kx +  (y_j - y_i)*ky ) );
                                        //temp += S_rw[pos_i][pos_j][wi]*exp(iota*( (x_j - x_i)*kx +  (y_j - y_i)*ky ) );
                                        temp2 += one*s_quantum_rw[pos_i][pos_j][wi]*exp(iota*( (x_j - x_i)*kx +  (y_j - y_i)*ky ) );


                                    }
                                }
                            }
                        }
                        //file_out2<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<temp.real()<<"   "<<temp.imag()<<"    "<<temp2.real()<<"   "<<temp2.imag()<<"    "<<temp3.real()<<"   "<<temp3.imag()<<endl;
                        //file_out_full<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<wi<<"   "<<temp.real()<<"   "<<temp.imag()<<"    "<<"none"<<"   "<<"none"<<"    "<<"none"<<"   "<<"none"<<endl;
                        Skw_Mat[nx][ny][wi]=temp;
                        squant_kw_Mat[nx][ny][wi]=temp2;
                    }
                }
            }


            k_ind=0;
            for(nx=0;nx<Parameters_.lx;nx++){
                for(ny=0;ny<Parameters_.ly;ny++){

                    for(int wi=0;wi<n_wpoints;wi++){

                        //file_out2<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<temp.real()<<"   "<<temp.imag()<<"    "<<temp2.real()<<"   "<<temp2.imag()<<"    "<<temp3.real()<<"   "<<temp3.imag()<<endl;
                        file_out_full<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<wi<<"   "<<Skw_Mat[nx][ny][wi].real()<<"   "<<Skw_Mat[nx][ny][wi].imag()<<"    "<<squant_kw_Mat[nx][ny][wi].real()<<"   "<<squant_kw_Mat[nx][ny][wi].imag()<<
                                       "     "<<endl;
                        //Skw_Mat[nx][ny][wi]=temp;
                    }

                    k_ind +=1;
                    file_out_full<<endl;

                }

            }



        }




        //Only in case of square Lattice
        if(Parameters_.lx==Parameters_.ly){


            if(use_only_allowed_k == true){

                k_ind=0;
                //From [(pi/2,pi/2)----to----(pi,pi) )
                for(nx=Lxby4;nx<Lxby2;nx++){
                    ny=nx;
                    for(int wi=0;wi<n_wpoints;wi++){

                        //file_out2<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<temp.real()<<"   "<<temp.imag()<<"    "<<temp2.real()<<"   "<<temp2.imag()<<"    "<<temp3.real()<<"   "<<temp3.imag()<<endl;
                        file_out2<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<wi<<"   "<<Skw_Mat[nx][ny][wi].real()<<"   "<<Skw_Mat[nx][ny][wi].imag()<<"    "<<squant_kw_Mat[nx][ny][wi].real()<<"   "<<squant_kw_Mat[nx][ny][wi].imag()<<
                                   "     "<<endl;

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
                                   "     "<<endl;
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
                                   "     "<<endl;
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
                                   "     "<<endl;            }
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

                    // kx = (2*nx*PI)/(1.0*Parameters_.lx);
                    // ky = (2*ny*PI)/(1.0*Parameters_.ly);
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(pos_i,pos_j,temp,temp2,temp3)
#endif
                    for(int wi=0;wi<n_wpoints;wi++){
                        temp3=zero;
                        temp2=zero;
                        temp=zero;
                        for(int x_i=0;x_i<Parameters_.lx;x_i++){
                            for(int y_i=0;y_i<Parameters_.ly;y_i++){
                                pos_i = y_i*Parameters_.lx + x_i;
                                for(int x_j=0;x_j<Parameters_.lx;x_j++){
                                    for(int y_j=0;y_j<Parameters_.ly;y_j++){
                                        pos_j = y_j*Parameters_.lx + x_j;

                                        temp += one*S_rw[pos_i][pos_j][wi]*exp(iota*( (x_j - x_i)*kx +  (y_j - y_i)*ky ) );
                                        //temp += S_rw[pos_i][pos_j][wi]*exp(iota*( (x_j - x_i)*kx +  (y_j - y_i)*ky ) );
                                        temp2 += one*s_quantum_rw[pos_i][pos_j][wi]*exp(iota*( (x_j - x_i)*kx +  (y_j - y_i)*ky ) );

                                    }
                                }
                            }
                        }
                        //file_out2<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<temp.real()<<"   "<<temp.imag()<<"    "<<temp2.real()<<"   "<<temp2.imag()<<"    "<<temp3.real()<<"   "<<temp3.imag()<<endl;
                        file_out2<<kx<<"   "<<ky<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<wi<<"   "<<temp.real()<<"   "<<temp.imag()<<"    "<<temp2.real()<<"   "<<temp2.imag()<<"    "<<endl;
                    }
                    file_out2<<endl;
                    k_ind +=1;
                    kx=kx + PI*dk_;
                }




                //From [(pi,pi)----to----(pi,0) )
                kx = PI;
                ky = PI;

                while(ky>0){

                    //kx = (2*nx*PI)/(1.0*Parameters_.lx);
                    //ky = (2*ny*PI)/(1.0*Parameters_.ly);
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(pos_i,pos_j,temp,temp2,temp3)
#endif
                    for(int wi=0;wi<n_wpoints;wi++){
                        temp3=zero;
                        temp2=zero;
                        temp=zero;
                        for(int x_i=0;x_i<Parameters_.lx;x_i++){
                            for(int y_i=0;y_i<Parameters_.ly;y_i++){
                                pos_i = y_i*Parameters_.lx + x_i;
                                for(int x_j=0;x_j<Parameters_.lx;x_j++){
                                    for(int y_j=0;y_j<Parameters_.ly;y_j++){
                                        pos_j = y_j*Parameters_.lx + x_j;

                                        temp += one*S_rw[pos_i][pos_j][wi]*exp(iota*( (x_j - x_i)*kx +  (y_j - y_i)*ky ) );
                                        //temp += S_rw[pos_i][pos_j][wi]*exp(iota*( (x_j - x_i)*kx +  (y_j - y_i)*ky ) );
                                        temp2 += one*s_quantum_rw[pos_i][pos_j][wi]*exp(iota*( (x_j - x_i)*kx +  (y_j - y_i)*ky ) );
                                    }
                                }
                            }
                        }
                        // file_out2<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<temp.real()<<"   "<<temp.imag()<<"   "<<temp2.real()<<"   "<<temp2.imag()   <<"    "<<temp3.real()<<"   "<<temp3.imag()<<   endl;
                        file_out2<<kx<<"   "<<ky<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<wi<<"   "<<temp.real()<<"   "<<temp.imag()<<"    "<<temp2.real()<<"   "<<temp2.imag()<<"    "<<endl;
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
                        for(int x_i=0;x_i<Parameters_.lx;x_i++){
                            for(int y_i=0;y_i<Parameters_.ly;y_i++){
                                pos_i = y_i*Parameters_.lx + x_i;
                                for(int x_j=0;x_j<Parameters_.lx;x_j++){
                                    for(int y_j=0;y_j<Parameters_.ly;y_j++){
                                        pos_j = y_j*Parameters_.lx + x_j;

                                        temp += one*S_rw[pos_i][pos_j][wi]*exp(iota*( (x_j - x_i)*kx +  (y_j - y_i)*ky ) );
                                        //temp += S_rw[pos_i][pos_j][wi]*exp(iota*( (x_j - x_i)*kx +  (y_j - y_i)*ky ) );
                                        temp2 += one*s_quantum_rw[pos_i][pos_j][wi]*exp(iota*( (x_j - x_i)*kx +  (y_j - y_i)*ky ) );

                                    }
                                }
                            }
                        }
                        // file_out2<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<temp.real()<<"   "<<temp.imag()<<"   "<<temp2.real()<<"   "<<temp2.imag()   <<"    "<<temp3.real()<<"   "<<temp3.imag()<<endl;
                        file_out2<<kx<<"   "<<ky<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<wi<<"   "<<temp.real()<<"   "<<temp.imag()<<"    "<<temp2.real()<<"   "<<temp2.imag()<<"    "<<endl;

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
                        for(int x_i=0;x_i<Parameters_.lx;x_i++){
                            for(int y_i=0;y_i<Parameters_.ly;y_i++){
                                pos_i = y_i*Parameters_.lx + x_i;
                                for(int x_j=0;x_j<Parameters_.lx;x_j++){
                                    for(int y_j=0;y_j<Parameters_.ly;y_j++){
                                        pos_j = y_j*Parameters_.lx + x_j;

                                        temp += one*S_rw[pos_i][pos_j][wi]*exp(iota*( (x_j - x_i)*kx +  (y_j - y_i)*ky ) );
                                        //temp += S_rw[pos_i][pos_j][wi]*exp(iota*( (x_j - x_i)*kx +  (y_j - y_i)*ky ) );
                                        temp2 += one*s_quantum_rw[pos_i][pos_j][wi]*exp(iota*( (x_j - x_i)*kx +  (y_j - y_i)*ky ) );
                                    }
                                }
                            }
                        }
                        // file_out2<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<temp.real()<<"   "<<temp.imag()<<"   "<<temp2.real()<<"   "<<temp2.imag()   <<"    "<<temp3.real()<<"   "<<temp3.imag()<<   endl;
                        file_out2<<kx<<"   "<<ky<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<wi<<"   "<<temp.real()<<"   "<<temp.imag()<<"    "<<temp2.real()<<"   "<<temp2.imag()<<"    "<<endl;
                    }
                    file_out2<<endl;
                    k_ind +=1;
                    kx = kx -PI*dk_;
                    ky = ky -PI*dk_;
                }




            }



        }
        /*
    for(int nx=0;nx<Parameters_.lx;nx++){
        for(int ny=0;ny<Parameters_.ly;ny++){
               kx = (2*nx*PI)/(1.0*Parameters_.lx);
               ky = (2*ny*PI)/(1.0*Parameters_.ly);

            for(int wi=0;wi<n_wpoints;wi++){
                temp=zero;
                for(int x_i=0;x_i<Parameters_.lx;x_i++){
                    for(int y_i=0;y_i<Parameters_.ly;y_i++){
                        pos_i = y_i*Parameters_.lx + x_i;
                        for(int x_j=0;x_j<Parameters_.lx;x_j++){
                            for(int y_j=0;y_j<Parameters_.ly;y_j++){
                                pos_j = y_j*Parameters_.lx + x_j;

                                temp += S_rw[pos_i][pos_j][wi]*exp(iota*( (x_j - x_i)*kx +  (y_j - y_i)*ky ) );


                            }
                        }
                    }
                }
                file_out2<<nx<<"   "<<ny<<"   "<<nx*Parameters_.lx + ny<<"   "<<wi*dw<<"   "<<temp.real()<<"   "<<temp.imag()<<endl;
            }
            file_out2<<endl;
        }
    }
*/


    }

}





void SC_SW_ENGINE_VNE_1orb_MCMF::Evolve_classical_spins_Predictor_Corrector(){

    /*

    complex<double> zero(0,0);
    complex<double> one(1,0);
    complex<double> iota(0,1);
    bool EVOLVE_RED_DEN_MAT = false;
    Mat_2_doub Theta_temp, Phi_temp;
    Theta_temp.resize(Parameters_.lx);
    Phi_temp.resize(Parameters_.lx);
    for(int i=0;i<Parameters_.lx;i++){
        Theta_temp[i].resize(Parameters_.ly);
        Phi_temp[i].resize(Parameters_.ly);
    }



    int ts=0;

    int site;

    Mat_6_Complex_doub Red_Den_mat_temp;



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




    //-------Adams-Bashforth Predictor---------//

    double pi,ti; //phi_i,theta_i,phi_j,theta_j for timestep=ts
    vector<double> sy,sx,sz; //Quantum spins for timestep=ts and position=pos
    sx.resize(Parameters_.orbs);sy.resize(Parameters_.orbs);sz.resize(Parameters_.orbs);

    complex<double> derivative_val, temp_val;
    vector<double> coefficients_ABM;
    coefficients_ABM.resize(4);
    coefficients_ABM[0]= 55.0; //19.0 + ((9.0*55.0)*(dt_/24.0));
    coefficients_ABM[1]= -59.0; //-5.0 - ((9.0*59.0)*(dt_/24.0));
    coefficients_ABM[2]= 37.0; //1.0 + ((9.0*37.0)*(dt_/24.0));
    coefficients_ABM[3]= -9.0; //-1.0*((9.0*9.0)*(dt_/24.0));


    Phi[1]= Phi[0];
    Theta[1]= Theta[0];


    for(int step=0;step<4;step++){

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(pi, ti, sx, sy, sz)
#endif
        for(int pos=0;pos<Parameters_.ns;pos++){
            int pos_x= Coordinates_.indx(pos);
            int pos_y= Coordinates_.indy(pos);

            pi=Phi_time[step][pos_x][pos_y];
            ti=Theta_time[step][pos_x][pos_y];

            for(int orb=0;orb<Parameters_.orbs;orb++){
                sx[orb]=0.5*real( Red_Den_mat_time[step][pos][orb][1][pos][orb][0] +
                        Red_Den_mat_time[step][pos][orb][0][pos][orb][1] );
                sy[orb]=0.5*imag( Red_Den_mat_time[step][pos][orb][1][pos][orb][0] -
                        Red_Den_mat_time[step][pos][orb][0][pos][orb][1] );
                sz[orb]=0.5*real( Red_Den_mat_time[step][pos][orb][0][pos][orb][0] -
                        Red_Den_mat_time[step][pos][orb][1][pos][orb][1] );
            }



            for(int orb=0;orb<Parameters_.orbs;orb++){
                Phi[1][pos_x][pos_y] +=  coefficients_ABM[step]*((dt_/24.0)*(  -1.0*(Jval_array[pos]*((((sin(pi)*cos(ti))/(sin(ti)))*sy[orb]) +
                                                                                                      (((cos(pi)*cos(ti))/(sin(ti)))*sx[orb]) - sz[orb]))  -  Bval_array[pos]  )    );

                Theta[1][pos_x][pos_y] +=  coefficients_ABM[step]*((dt_/24.0)*( Jval_array[pos]*((cos(pi)*sy[orb]) -
                                                                                                 (sin(pi)*sx[orb]))));
            }


            int pos_ng;
            double pj,tj;
            double value_;
            for (int ng=0;ng<8;ng++){
                if(ng<2){
                    value_=1.14*Parameters_.J_NN;  // +x, -x
                }
                else if(ng>=2 && ng <4){
                    value_=Parameters_.J_NN;   // +y, -y
                }
                else{
                    value_=Parameters_.J_NNN;  //+x+y, +x-y, -x+y, -x-y
                }
                pos_ng= Coordinates_.neigh(pos,ng);
                int pos_ng_x = Coordinates_.indx(pos_ng);
                int pos_ng_y = Coordinates_.indy(pos_ng);
                pj=Phi_time[step][pos_ng_x][pos_ng_y];
                tj=Theta_time[step][pos_ng_x][pos_ng_y];
                Phi[ts+1][pos_x][pos_y] += coefficients_ABM[step]*(-1.0)*((dt_/24.0)*( (S_mag)*(value_)* (((sin(tj)*cos(ti)*cos(pi - pj))
                                                                                                           /(sin(ti))) - cos(tj) )));
                Theta[ts+1][pos_x][pos_y] += coefficients_ABM[step]*((dt_/24.0)*(S_mag)*(value_)*(sin(tj)*sin(pj - pi)));

            }






            // cout<<pos<<"  "<<derivative_theta<<"  "<<derivative_phi<<endl;

        }
    }

    //-----------------------Corrector------

    Theta_temp=Theta[0];
    Phi_temp=Phi[0];

    coefficients_ABM[0]= 19.0;
    coefficients_ABM[1]= -5.0;
    coefficients_ABM[2]= 1.0;
    coefficients_ABM[3]= 9.0;

    for(int step=0;step<4;step++){

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(pi, ti, sx, sy, sz)
#endif
        for(int pos=0;pos<Parameters_.ns;pos++){
            int pos_x= Coordinates_.indx(pos);
            int pos_y= Coordinates_.indy(pos);

            if(step!=3){
                pi=Phi_time[step][pos_x][pos_y];
                ti=Theta_time[step][pos_x][pos_y];
            }
            else{
                pi=Phi[1][pos_x][pos_y];
                ti=Theta[1][pos_x][pos_y];
            }

            for(int orb=0;orb<Parameters_.orbs;orb++){
                sx[orb]=0.5*real( Red_Den_mat_time[step][pos][orb][1][pos][orb][0] +
                        Red_Den_mat_time[step][pos][orb][0][pos][orb][1] );
                sy[orb]=0.5*imag( Red_Den_mat_time[step][pos][orb][1][pos][orb][0] -
                        Red_Den_mat_time[step][pos][orb][0][pos][orb][1] );
                sz[orb]=0.5*real( Red_Den_mat_time[step][pos][orb][0][pos][orb][0] -
                        Red_Den_mat_time[step][pos][orb][1][pos][orb][1] );
            }



            for(int orb=0;orb<Parameters_.orbs;orb++){
                Phi_temp[pos_x][pos_y] +=  coefficients_ABM[step]*((dt_/24.0)*(  -1.0*(Jval_array[pos]*((((sin(pi)*cos(ti))/(sin(ti)))*sy[orb]) +
                                                                                                        (((cos(pi)*cos(ti))/(sin(ti)))*sx[orb]) - sz[orb]))  -  Bval_array[pos]  )    );

                Theta_temp[pos_x][pos_y] +=  coefficients_ABM[step]*((dt_/24.0)*( Jval_array[pos]*((cos(pi)*sy[orb]) -
                                                                                                   (sin(pi)*sx[orb]))));
            }


            int pos_ng;
            double pj,tj;
            double value_;
            for (int ng=0;ng<8;ng++){
                if(ng<2){
                    value_=1.14*Parameters_.J_NN;  // +x, -x
                }
                else if(ng>=2 && ng <4){
                    value_=Parameters_.J_NN;   // +y, -y
                }
                else{
                    value_=Parameters_.J_NNN;  //+x+y, +x-y, -x+y, -x-y
                }
                pos_ng= Coordinates_.neigh(pos,ng);
                int pos_ng_x = Coordinates_.indx(pos_ng);
                int pos_ng_y = Coordinates_.indy(pos_ng);

                if(step!=3){
                    pj=Phi_time[step][pos_ng_x][pos_ng_y];
                    tj=Theta_time[step][pos_ng_x][pos_ng_y];
                }
                else{
                    pj=Phi[1][pos_ng_x][pos_ng_y];
                    tj=Theta[1][pos_ng_x][pos_ng_y];
                }

                Phi_temp[pos_x][pos_y] += coefficients_ABM[step]*(-1.0)*((dt_/24.0)*( (S_mag)*(value_)* (((sin(tj)*cos(ti)*cos(pi - pj))
                                                                                                          /(sin(ti))) - cos(tj) )));
                Theta_temp[pos_x][pos_y] += coefficients_ABM[step]*((dt_/24.0)*(S_mag)*(value_)*(sin(tj)*sin(pj - pi)));

            }






            // cout<<pos<<"  "<<derivative_theta<<"  "<<derivative_phi<<endl;

        }
    }


    //--------------------


    Phi[1]=Phi_temp;
    Theta[1]=Theta_temp;
    Phi_temp.clear();
    Theta_temp.clear();

    for(int pos=0;pos<Parameters_.ns;pos++){
        int pos_x= Coordinates_.indx(pos);
        int pos_y= Coordinates_.indy(pos);

        if(Phi[ts+1][pos_x][pos_y] > 2*PI){
            Phi[ts+1][pos_x][pos_y] += -2*PI;

        }
        if(Phi[ts+1][pos_x][pos_y] < 0){
            Phi[ts+1][pos_x][pos_y] +=  2*PI;

        }


        if(Theta[ts+1][pos_x][pos_y] > PI){
            Theta[ts+1][pos_x][pos_y] = Theta[ts+1][pos_x][pos_y] -2*PI;
            Phi[ts+1][pos_x][pos_y] = fmod( Phi[ts+1][pos_x][pos_y] + PI, 2.0*PI );
        }
        if(Theta[ts+1][pos_x][pos_y] < 0){
            Theta[ts+1][pos_x][pos_y] = - Theta[ts+1][pos_x][pos_y];
            Phi[ts+1][pos_x][pos_y] = fmod( Phi[ts+1][pos_x][pos_y] + PI, 2.0*PI );

        }

    }



    if(EVOLVE_RED_DEN_MAT==true){
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(derivative_val,temp_val)
#endif
        for(int pos_i=0;pos_i<Parameters_.ns;pos_i++){
            for(int orb_i=0;orb_i<Parameters_.orbs;orb_i++){
                for(int si=0;si<2;si++){
                    for(int pos_j=0;pos_j<Parameters_.ns;pos_j++){
                        for(int orb_j=0;orb_j<Parameters_.orbs;orb_j++){
                            for(int sj=0;sj<2;sj++){




                                int pos_i_x= Coordinates_.indx(pos_i);
                                int pos_i_y= Coordinates_.indy(pos_i);
                                int pos_j_x= Coordinates_.indx(pos_j);
                                int pos_j_y= Coordinates_.indy(pos_j);

                                complex<double> phasex, phasey;
                                int l, spin_l, orb_l;
                                double l_i;
                                double value_temp1;
                                phasex=one_complex;
                                phasey=one_complex;
                                derivative_val=zero_complex;



                                //------------------i connection to l :start---------------------//
                                // Phase from As positions
                                l_i=pow (-1.00, Coordinates_.indx(pos_i) + Coordinates_.indy(pos_i) );

                                spin_l = si;

                                // +x,-x direction Neighbor
                                for (int ng=0;ng<8;ng++){

                                    l = Coordinates_.neigh(pos_i,ng);
                                    for(orb_l=0;orb_l<Parameters_.orbs;orb_l++) {

                                        if(ng==0 || ng==1){
                                            value_temp1=Hamiltonian_.Tx(orb_l,orb_i);
                                        }
                                        if(ng==2 || ng==3){
                                            value_temp1=Hamiltonian_.Ty(orb_l,orb_i);
                                        }
                                        if(ng==4 || ng==6){
                                            value_temp1=Hamiltonian_.Tpxpy(orb_l,orb_i);
                                        }
                                        if(ng==5 || ng==7){
                                            value_temp1=Hamiltonian_.Tpxmy(orb_l,orb_i);
                                        }

                                        if ( (orb_i==2) ^ (orb_l==2) ) {
                                            derivative_val +=iota*complex<double>(1.0*l_i*value_temp1,0.0)*
                                                    Red_Den_mat_temp[l][orb_l][spin_l][pos_j][orb_j][sj];

                                        }
                                        else {
                                            derivative_val +=iota*complex<double>(1.0*value_temp1, 0.0)*
                                                    Red_Den_mat_temp[l][orb_l][spin_l][pos_j][orb_j][sj];

                                        }
                                    }
                                }



                                //------------------i connection to l :done---------------------//


                                //------------------j connection to l :start---------------------//
                                // Phase from As positions
                                l_i=pow (-1.00, Coordinates_.indx(pos_j) + Coordinates_.indy(pos_j) );

                                spin_l = sj;

                                // +x,-x direction Neighbor
                                for (int ng=0;ng<8;ng++){

                                    l = Coordinates_.neigh(pos_j,ng);
                                    for(orb_l=0;orb_l<Parameters_.orbs;orb_l++) {

                                        if(ng==0 || ng==1){
                                            value_temp1=Hamiltonian_.Tx(orb_l,orb_j);

                                        }
                                        if(ng==2 || ng==3){
                                            value_temp1=Hamiltonian_.Ty(orb_l,orb_j);

                                        }
                                        if(ng==4 || ng==6){
                                            value_temp1=Hamiltonian_.Tpxpy(orb_l,orb_j);

                                        }
                                        if(ng==5 || ng==7){
                                            value_temp1=Hamiltonian_.Tpxmy(orb_l,orb_j);

                                        }

                                        if ( (orb_j==2) ^ (orb_l==2) ) {
                                            derivative_val +=iota*complex<double>(1.0*l_i*value_temp1,0.0)*
                                                    Red_Den_mat_temp[pos_i][orb_i][si][l][orb_l][spin_l];

                                        }
                                        else {
                                            derivative_val +=iota*complex<double>(1.0*value_temp1, 0.0)*
                                                    Red_Den_mat_temp[pos_i][orb_i][si][l][orb_l][spin_l];

                                        }
                                    }
                                }

                                // On-Site potential for orb = 2 (xy)
                                derivative_val +=complex<double>(1.0*(Hamiltonian_.potential_local[orb_i] - Hamiltonian_.potential_local[orb_j] ),0.0)*
                                        Red_Den_mat_temp[pos_i][orb_i][si][pos_j][orb_j][sj]  ;



                                //------------------j connection to l :done---------------------//



                                for(int s3=0;s3<2;s3++){
                                    double SX_i,SY_i,SZ_i,SX_j,SY_j,SZ_j;
                                    SX_i = sin(Theta[ts][pos_i_x][pos_i_y])*cos(Phi[ts][pos_i_x][pos_i_y]);
                                    SY_i = sin(Theta[ts][pos_i_x][pos_i_y])*sin(Phi[ts][pos_i_x][pos_i_y]);
                                    SZ_i = cos(Theta[ts][pos_i_x][pos_i_y]);
                                    SX_j = sin(Theta[ts][pos_j_x][pos_j_y])*cos(Phi[ts][pos_j_x][pos_j_y]);
                                    SY_j = sin(Theta[ts][pos_j_x][pos_j_y])*sin(Phi[ts][pos_j_x][pos_j_y]);
                                    SZ_j = cos(Theta[ts][pos_j_x][pos_j_y]);
                                    temp_val=Red_Den_mat_temp[pos_i][orb_i][s3][pos_j][orb_j][sj];
                                    derivative_val+=iota*(-0.5*Jval_array[pos_i])*(S_mag)*(
                                                (SX_i)*Pauli_x[si][s3]
                                                +
                                                (SY_i)*Pauli_y[si][s3]
                                                +
                                                (SZ_i)*Pauli_z[si][s3]
                                                )*temp_val;
                                    temp_val=Red_Den_mat_temp[pos_i][orb_i][si][pos_j][orb_j][s3];
                                    derivative_val +=iota*(0.5*Jval_array[pos_j])*(S_mag)*(
                                                (SX_j)*Pauli_x[s3][sj]
                                                +
                                                (SY_j)*Pauli_y[s3][sj]
                                                +
                                                (SZ_j)*Pauli_z[s3][sj]
                                                )*temp_val;



                                }

                                Red_Den_mat[pos_i][orb_i][si][pos_j][orb_j][sj] += dt_*derivative_val;

                            }
                        }
                    }
                }
            }
        }

    }


    //cout<<"Only RK4 is working for VNE"<<endl;
    //assert (Runge_Kutta_order==4);


    //------------------------------





*/

}

void SC_SW_ENGINE_VNE_1orb_MCMF::Evolve_classical_spins_Runge_Kutta(int ts){


    complex<double> zero(0,0);
    complex<double> one(1,0);
    complex<double> iota(0,1);
    bool EVOLVE_RED_DEN_MAT = true;


    int site;

    Mat_4_Complex_doub Red_Den_mat_temp;
    Red_Den_mat_temp=Red_Den_mat;


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



    if(Runge_Kutta_order==1){

        double Aux_Sx_i,Aux_Sy_i,Aux_Sz_i;
        double Aux_Sx_j, Aux_Sy_j, Aux_Sz_j;
        int pos_neigh;
        double sy,sx,sz; //Quantum spins for timestep=ts and position=pos

        complex<double> derivative_val;

        if(SelfConsistentEvolution==false){
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(Aux_Sx_i, Aux_Sy_i, Aux_Sz_i, sx, sy, sz)
#endif
            for(int pos=0;pos<Parameters_.ns;pos++){

                Aux_Sx_i=Aux_S_x[ts][pos];
                Aux_Sy_i=Aux_S_y[ts][pos];
                Aux_Sz_i=Aux_S_z[ts][pos];

                sx=0.5*real( Red_Den_mat_temp[pos][1][pos][0] + Red_Den_mat_temp[pos][0][pos][1] );
                sy=0.5*imag( Red_Den_mat_temp[pos][1][pos][0] - Red_Den_mat_temp[pos][0][pos][1] );
                sz=0.5*real( Red_Den_mat_temp[pos][0][pos][0] - Red_Den_mat_temp[pos][1][pos][1] );


                Aux_S_x[ts+1][pos] =Aux_S_x[ts][pos];
                Aux_S_y[ts+1][pos] =Aux_S_y[ts][pos];
                Aux_S_z[ts+1][pos] =Aux_S_z[ts][pos];


                Aux_S_x[ts+1][pos] += dt_*(Jval_array[pos])*(sy*Aux_Sz_i - sz*Aux_Sy_i);
                Aux_S_y[ts+1][pos] += dt_*(Jval_array[pos])*(sz*Aux_Sx_i - sx*Aux_Sz_i);
                Aux_S_z[ts+1][pos] += dt_*(Jval_array[pos])*(sx*Aux_Sy_i - sy*Aux_Sx_i);

                // cout<<pos<<"  "<<derivative_theta<<"  "<<derivative_phi<<endl;


                for (int ng=0;ng<4;ng++){
                    pos_neigh = Coordinates_.neigh(pos,ng);

                    Aux_Sx_j=Aux_S_x[ts][pos_neigh];
                    Aux_Sy_j=Aux_S_y[ts][pos_neigh];
                    Aux_Sz_j=Aux_S_z[ts][pos_neigh];

                    Aux_S_x[ts+1][pos] += dt_*J_Classical*(Aux_Sy_j*Aux_Sz_i - Aux_Sz_j*Aux_Sy_i);
                    Aux_S_y[ts+1][pos] += dt_*J_Classical*(Aux_Sz_j*Aux_Sx_i - Aux_Sx_j*Aux_Sz_i);
                    Aux_S_z[ts+1][pos] += dt_*J_Classical*(Aux_Sx_j*Aux_Sy_i - Aux_Sy_j*Aux_Sx_i);

                }

            }

        }


        if(EVOLVE_RED_DEN_MAT==true){
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(derivative_val)
#endif
            for(int pos_i=0;pos_i<Parameters_.ns;pos_i++){
                for(int si=0;si<2;si++){
                    for(int pos_j=0;pos_j<Parameters_.ns;pos_j++){
                        for(int sj=0;sj<2;sj++){

                            complex<double> phasex, phasey;
                            int l, spin_l;
                            phasex=one_complex;
                            phasey=one_complex;
                            derivative_val=zero_complex;



                            //------------------i connection to l :start---------------------//
                            spin_l = si;

                            // +x,-x, +/-y direction Neighbor
                            for (int ng=0;ng<4;ng++){
                                l = Coordinates_.neigh(pos_i,ng);
                                derivative_val +=iota*complex<double>(-1.0*Parameters_.t_hopping, 0.0)*
                                        Red_Den_mat_temp[l][spin_l][pos_j][sj];

                            }
                            //Sign with Parameters_.t_hopping, depends on that Hamiltonian has -t or t, my Hamil is with t


                            //------------------i connection to l :done---------------------//


                            //------------------j connection to l :start---------------------//
                            spin_l = sj;

                            // +x,-x , +/-y direction Neighbor
                            for (int ng=0;ng<4;ng++){
                                l = Coordinates_.neigh(pos_j,ng);
                                derivative_val +=iota*complex<double>(1.0*Parameters_.t_hopping, 0.0)*
                                        Red_Den_mat_temp[pos_i][si][l][spin_l];

                            }

                            //------------------j connection to l :done---------------------//


                            //-------------------coupling with Auxlliary spins--------------//
                            for(int s3=0;s3<2;s3++){
                                double SX_i,SY_i,SZ_i,SX_j,SY_j,SZ_j;
                                SX_i = Aux_S_x[ts][pos_i];
                                SY_i = Aux_S_y[ts][pos_i];
                                SZ_i = Aux_S_z[ts][pos_i];
                                SX_j = Aux_S_x[ts][pos_j];
                                SY_j = Aux_S_y[ts][pos_j];
                                SZ_j = Aux_S_z[ts][pos_j];


                                derivative_val+=iota*(-0.5*Jval_array[pos_i])*(
                                            (SX_i)*Pauli_x[si][s3]
                                            +
                                            (SY_i)*Pauli_y[si][s3]
                                            +
                                            (SZ_i)*Pauli_z[si][s3]
                                            )*Red_Den_mat_temp[pos_i][s3][pos_j][sj];


                                derivative_val +=iota*(0.5*Jval_array[pos_j])*(
                                            (SX_j)*Pauli_x[s3][sj]
                                            +
                                            (SY_j)*Pauli_y[s3][sj]
                                            +
                                            (SZ_j)*Pauli_z[s3][sj]
                                            )*Red_Den_mat_temp[pos_i][si][pos_j][s3];

                            }

                            //---------------------------done----------------------------------//


                            Red_Den_mat[pos_i][si][pos_j][sj] += dt_*derivative_val;

                        }
                    }
                }
            }

        }


        if(SelfConsistentEvolution==true){
            for(int pos=0;pos<Parameters_.ns;pos++){
                sx=0.5*real( Red_Den_mat_temp[pos][1][pos][0] + Red_Den_mat_temp[pos][0][pos][1] );
                sy=0.5*imag( Red_Den_mat_temp[pos][1][pos][0] - Red_Den_mat_temp[pos][0][pos][1] );
                sz=0.5*real( Red_Den_mat_temp[pos][0][pos][0] - Red_Den_mat_temp[pos][1][pos][1] );
                Aux_S_x[ts+1][pos] = sx;
                Aux_S_y[ts+1][pos] = sy;
                Aux_S_z[ts+1][pos] = sz;
            }
        }
        //cout<<"Only RK4 is working for VNE"<<endl;
        //assert (Runge_Kutta_order==4);
    }

    else if(Runge_Kutta_order==4){


        double sy,sx,sz; //Quantum spins for timestep=ts and position=pos
        double Aux_Sx_i, Aux_Sy_i, Aux_Sz_i;
        double Aux_Sx_j, Aux_Sy_j, Aux_Sz_j;


        double derivative_Aux_Sx, derivative_Aux_Sy, derivative_Aux_Sz;
        Mat_1_doub delta1_Aux_Sx, delta2_Aux_Sx, delta3_Aux_Sx, delta4_Aux_Sx;
        Mat_1_doub delta1_Aux_Sy, delta2_Aux_Sy, delta3_Aux_Sy, delta4_Aux_Sy;
        Mat_1_doub delta1_Aux_Sz, delta2_Aux_Sz, delta3_Aux_Sz, delta4_Aux_Sz;

        Mat_4_Complex_doub delta1_Red_Den_mat, delta2_Red_Den_mat, delta3_Red_Den_mat, delta4_Red_Den_mat;

        complex<double> phasex, phasey;
        int l, spin_l, pos_neigh;


        //Standard RK-4 method Convention is used
        delta1_Aux_Sx.resize(Parameters_.ns);delta2_Aux_Sx.resize(Parameters_.ns);
        delta3_Aux_Sx.resize(Parameters_.ns);delta4_Aux_Sx.resize(Parameters_.ns);
        delta1_Aux_Sy.resize(Parameters_.ns);delta2_Aux_Sy.resize(Parameters_.ns);
        delta3_Aux_Sy.resize(Parameters_.ns);delta4_Aux_Sy.resize(Parameters_.ns);
        delta1_Aux_Sz.resize(Parameters_.ns);delta2_Aux_Sz.resize(Parameters_.ns);
        delta3_Aux_Sz.resize(Parameters_.ns);delta4_Aux_Sz.resize(Parameters_.ns);



        delta1_Red_Den_mat.resize(Parameters_.ns);delta2_Red_Den_mat.resize(Parameters_.ns);
        delta3_Red_Den_mat.resize(Parameters_.ns); delta4_Red_Den_mat.resize(Parameters_.ns);

        for(int i=0;i<Parameters_.ns;i++){
            delta1_Red_Den_mat[i].resize(2);delta2_Red_Den_mat[i].resize(2);
            delta3_Red_Den_mat[i].resize(2);delta4_Red_Den_mat[i].resize(2);
            for(int si=0;si<2;si++){
                delta1_Red_Den_mat[i][si].resize(Parameters_.ns);delta2_Red_Den_mat[i][si].resize(Parameters_.ns);
                delta3_Red_Den_mat[i][si].resize(Parameters_.ns);delta4_Red_Den_mat[i][si].resize(Parameters_.ns);
                for(int j=0;j<Parameters_.ns;j++){
                    delta1_Red_Den_mat[i][si][j].resize(2);delta2_Red_Den_mat[i][si][j].resize(2);
                    delta3_Red_Den_mat[i][si][j].resize(2);delta4_Red_Den_mat[i][si][j].resize(2);

                }
            }
        }


        for(int i=0;i<Parameters_.ns;i++){
            delta1_Aux_Sx[i]=0.0;delta2_Aux_Sx[i]=0.0;
            delta3_Aux_Sx[i]=0.0;delta4_Aux_Sx[i]=0.0;
            delta1_Aux_Sy[i]=0.0;delta2_Aux_Sy[i]=0.0;
            delta3_Aux_Sy[i]=0.0;delta4_Aux_Sy[i]=0.0;
            delta1_Aux_Sz[i]=0.0;delta2_Aux_Sz[i]=0.0;
            delta3_Aux_Sz[i]=0.0;delta4_Aux_Sz[i]=0.0;
            for(int si=0;si<2;si++){
                for(int j=0;j<Parameters_.ns;j++){
                    for(int sj=0;sj<2;sj++){
                        delta1_Red_Den_mat[i][si][j][sj]=zero;delta2_Red_Den_mat[i][si][j][sj]=zero;
                        delta3_Red_Den_mat[i][si][j][sj]=zero;delta4_Red_Den_mat[i][si][j][sj]=zero;
                    }
                }
            }
        }



        //calculating delta1(2,3,4)_thetha(phi)
        for(int step_no=0;step_no<4;step_no++){

            //Classical spins time evolution, due to coupling with Quantum spins and Magnetic field
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(derivative_Aux_Sx,derivative_Aux_Sy, derivative_Aux_Sz, Aux_Sx_i, Aux_Sy_i, Aux_Sz_i, sx, sy, sz)
#endif
            for(int pos=0;pos<Parameters_.ns;pos++){

                derivative_Aux_Sx=0;derivative_Aux_Sy=0;derivative_Aux_Sz=0;

                if(step_no==0){
                    Aux_Sx_i=Aux_S_x[ts][pos];
                    Aux_Sy_i=Aux_S_y[ts][pos];
                    Aux_Sz_i=Aux_S_z[ts][pos];
                    sx =0.5*real( Red_Den_mat_temp[pos][1][pos][0] + Red_Den_mat_temp[pos][0][pos][1] );
                    sy =0.5*imag( Red_Den_mat_temp[pos][1][pos][0] - Red_Den_mat_temp[pos][0][pos][1] );
                    sz =0.5*real( Red_Den_mat_temp[pos][0][pos][0] - Red_Den_mat_temp[pos][1][pos][1] );

                }
                if(step_no==1){
                    Aux_Sx_i=Aux_S_x[ts][pos] + 0.5*delta1_Aux_Sx[pos];
                    Aux_Sy_i=Aux_S_y[ts][pos] + 0.5*delta1_Aux_Sy[pos];
                    Aux_Sz_i=Aux_S_z[ts][pos] + 0.5*delta1_Aux_Sz[pos];

                    sx =0.5*real( Red_Den_mat_temp[pos][1][pos][0] + Red_Den_mat_temp[pos][0][pos][1]
                            + 0.5*(delta1_Red_Den_mat[pos][1][pos][0] + delta1_Red_Den_mat[pos][0][pos][1]));
                    sy =0.5*imag( Red_Den_mat_temp[pos][1][pos][0] - Red_Den_mat_temp[pos][0][pos][1]
                            +0.5*(delta1_Red_Den_mat[pos][1][pos][0] - delta1_Red_Den_mat[pos][0][pos][1])
                            );
                    sz =0.5*real( Red_Den_mat_temp[pos][0][pos][0] - Red_Den_mat_temp[pos][1][pos][1]
                            +0.5*(delta1_Red_Den_mat[pos][0][pos][0] - delta1_Red_Den_mat[pos][1][pos][1])
                            );
                }

                if(step_no==2){

                    Aux_Sx_i=Aux_S_x[ts][pos] + 0.5*delta2_Aux_Sx[pos];
                    Aux_Sy_i=Aux_S_y[ts][pos] + 0.5*delta2_Aux_Sy[pos];
                    Aux_Sz_i=Aux_S_z[ts][pos] + 0.5*delta2_Aux_Sz[pos];

                    sx =0.5*real( Red_Den_mat_temp[pos][1][pos][0] + Red_Den_mat_temp[pos][0][pos][1]
                            + 0.5*(delta2_Red_Den_mat[pos][1][pos][0] + delta2_Red_Den_mat[pos][0][pos][1]));
                    sy =0.5*imag( Red_Den_mat_temp[pos][1][pos][0] - Red_Den_mat_temp[pos][0][pos][1]
                            +0.5*(delta2_Red_Den_mat[pos][1][pos][0] - delta2_Red_Den_mat[pos][0][pos][1])
                            );
                    sz =0.5*real( Red_Den_mat_temp[pos][0][pos][0] - Red_Den_mat_temp[pos][1][pos][1]
                            +0.5*(delta2_Red_Den_mat[pos][0][pos][0] - delta2_Red_Den_mat[pos][1][pos][1])
                            );

                }

                if(step_no==3){

                    Aux_Sx_i=Aux_S_x[ts][pos] + delta3_Aux_Sx[pos];
                    Aux_Sy_i=Aux_S_y[ts][pos] + delta3_Aux_Sy[pos];
                    Aux_Sz_i=Aux_S_z[ts][pos] + delta3_Aux_Sz[pos];

                    sx =0.5*real( Red_Den_mat_temp[pos][1][pos][0] + Red_Den_mat_temp[pos][0][pos][1]
                            + (delta3_Red_Den_mat[pos][1][pos][0] + delta3_Red_Den_mat[pos][0][pos][1]));
                    sy =0.5*imag( Red_Den_mat_temp[pos][1][pos][0] - Red_Den_mat_temp[pos][0][pos][1]
                            +(delta3_Red_Den_mat[pos][1][pos][0] - delta3_Red_Den_mat[pos][0][pos][1])
                            );
                    sz =0.5*real( Red_Den_mat_temp[pos][0][pos][0] - Red_Den_mat_temp[pos][1][pos][1]
                            +(delta3_Red_Den_mat[pos][0][pos][0] - delta3_Red_Den_mat[pos][1][pos][1])
                            );

                }

                derivative_Aux_Sx = (Jval_array[pos])*(sy*Aux_Sz_i - sz*Aux_Sy_i);
                derivative_Aux_Sy = (Jval_array[pos])*(sz*Aux_Sx_i - sx*Aux_Sz_i);
                derivative_Aux_Sz = (Jval_array[pos])*(sx*Aux_Sy_i - sy*Aux_Sx_i);


                for (int ng=0;ng<4;ng++){
                    pos_neigh = Coordinates_.neigh(pos,ng);
                    if(step_no==0){
                        Aux_Sx_j=Aux_S_x[ts][pos_neigh];
                        Aux_Sy_j=Aux_S_y[ts][pos_neigh];
                        Aux_Sz_j=Aux_S_z[ts][pos_neigh];
                    }
                    if(step_no==1){
                        Aux_Sx_j=Aux_S_x[ts][pos_neigh] + 0.5*delta1_Aux_Sx[pos_neigh];
                        Aux_Sy_j=Aux_S_y[ts][pos_neigh] + 0.5*delta1_Aux_Sy[pos_neigh];
                        Aux_Sz_j=Aux_S_z[ts][pos_neigh] + 0.5*delta1_Aux_Sz[pos_neigh];
                    }
                    if(step_no==2){
                        Aux_Sx_j=Aux_S_x[ts][pos_neigh] + 0.5*delta2_Aux_Sx[pos_neigh];
                        Aux_Sy_j=Aux_S_y[ts][pos_neigh] + 0.5*delta2_Aux_Sy[pos_neigh];
                        Aux_Sz_j=Aux_S_z[ts][pos_neigh] + 0.5*delta2_Aux_Sz[pos_neigh];
                    }
                    if(step_no==3){
                        Aux_Sx_j=Aux_S_x[ts][pos_neigh] + delta3_Aux_Sx[pos_neigh];
                        Aux_Sy_j=Aux_S_y[ts][pos_neigh] + delta3_Aux_Sy[pos_neigh];
                        Aux_Sz_j=Aux_S_z[ts][pos_neigh] + delta3_Aux_Sz[pos_neigh];
                    }

                    derivative_Aux_Sx +=J_Classical*(Aux_Sy_j*Aux_Sz_i - Aux_Sz_j*Aux_Sy_i);
                    derivative_Aux_Sy +=J_Classical*(Aux_Sz_j*Aux_Sx_i - Aux_Sx_j*Aux_Sz_i);
                    derivative_Aux_Sz +=J_Classical*(Aux_Sx_j*Aux_Sy_i - Aux_Sy_j*Aux_Sx_i);
                }



                if(step_no==0){
                    delta1_Aux_Sx[pos] = dt_*(derivative_Aux_Sx);
                    delta1_Aux_Sy[pos] = dt_*(derivative_Aux_Sy);
                    delta1_Aux_Sz[pos] = dt_*(derivative_Aux_Sz);
                }
                if(step_no==1){
                    delta2_Aux_Sx[pos] = dt_*(derivative_Aux_Sx);
                    delta2_Aux_Sy[pos] = dt_*(derivative_Aux_Sy);
                    delta2_Aux_Sz[pos] = dt_*(derivative_Aux_Sz);

                }
                if(step_no==2){
                    delta3_Aux_Sx[pos] = dt_*(derivative_Aux_Sx);
                    delta3_Aux_Sy[pos] = dt_*(derivative_Aux_Sy);
                    delta3_Aux_Sz[pos] = dt_*(derivative_Aux_Sz);

                }
                if(step_no==3){
                    delta4_Aux_Sx[pos] = dt_*(derivative_Aux_Sx);
                    delta4_Aux_Sy[pos] = dt_*(derivative_Aux_Sy);
                    delta4_Aux_Sz[pos] = dt_*(derivative_Aux_Sz);
                }

            }



            complex<double> temp_val;
            complex<double> derivative_val;

            if(EVOLVE_RED_DEN_MAT==true){
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(derivative_val, temp_val)
#endif

                for(int pos_i=0;pos_i<Parameters_.ns;pos_i++){
                    for(int si=0;si<2;si++){
                        for(int pos_j=0;pos_j<Parameters_.ns;pos_j++){
                            for(int sj=0;sj<2;sj++){


                                phasex=one_complex;
                                phasey=one_complex;
                                derivative_val=zero_complex;



                                //------------------j connection to l :start---------------------//
                                // Phase from As positions

                                spin_l = si;
                                //All neighbours
                                for (int ng=0;ng<4;ng++){
                                    l = Coordinates_.neigh(pos_i,ng);

                                    if(step_no==0){
                                        derivative_val =iota*complex<double>(-1.0*Parameters_.t_hopping, 0.0)*
                                                Red_Den_mat_temp[l][spin_l][pos_j][sj];
                                        delta1_Red_Den_mat[pos_i][si][pos_j][sj] += dt_*derivative_val;
                                    }
                                    if(step_no==1){
                                        temp_val=Red_Den_mat_temp[l][spin_l][pos_j][sj]
                                                + 0.5*delta1_Red_Den_mat[l][spin_l][pos_j][sj];
                                        derivative_val=iota*complex<double>(-1.0*Parameters_.t_hopping,0.0)*temp_val;
                                        delta2_Red_Den_mat[pos_i][si][pos_j][sj] += dt_*derivative_val;
                                    }
                                    if(step_no==2){
                                        temp_val=Red_Den_mat_temp[l][spin_l][pos_j][sj]
                                                + 0.5*delta2_Red_Den_mat[l][spin_l][pos_j][sj];
                                        derivative_val=iota*complex<double>(-1.0*Parameters_.t_hopping,0.0)*temp_val;
                                        delta3_Red_Den_mat[pos_i][si][pos_j][sj] += dt_*derivative_val;
                                    }
                                    if(step_no==3){
                                        temp_val=Red_Den_mat_temp[l][spin_l][pos_j][sj]
                                                + delta3_Red_Den_mat[l][spin_l][pos_j][sj];
                                        derivative_val=iota*complex<double>(-1.0*Parameters_.t_hopping,0.0)*temp_val;
                                        delta4_Red_Den_mat[pos_i][si][pos_j][sj] += dt_*derivative_val;
                                    }

                                }



                                //------------------j connection to l :done---------------------//



                                //------------------l connection to i :start---------------------//
                                spin_l = sj;
                                for (int ng=0;ng<4;ng++){
                                    l = Coordinates_.neigh(pos_j,ng);
                                    if(step_no==0){
                                        derivative_val =iota*complex<double>(1.0*Parameters_.t_hopping, 0.0)*
                                                Red_Den_mat_temp[pos_i][si][l][spin_l];
                                        delta1_Red_Den_mat[pos_i][si][pos_j][sj] += dt_*derivative_val;
                                    }
                                    if(step_no==1){
                                        temp_val=Red_Den_mat_temp[pos_i][si][l][spin_l]
                                                + 0.5*delta1_Red_Den_mat[pos_i][si][l][spin_l];
                                        derivative_val=iota*complex<double>(1.0*Parameters_.t_hopping,0.0)*temp_val;
                                        delta2_Red_Den_mat[pos_i][si][pos_j][sj] += dt_*derivative_val;
                                    }
                                    if(step_no==2){
                                        temp_val=Red_Den_mat_temp[pos_i][si][l][spin_l]
                                                + 0.5*delta2_Red_Den_mat[pos_i][si][l][spin_l];
                                        derivative_val=iota*complex<double>(1.0*Parameters_.t_hopping,0.0)*temp_val;
                                        delta3_Red_Den_mat[pos_i][si][pos_j][sj] += dt_*derivative_val;
                                    }
                                    if(step_no==3){
                                        temp_val=Red_Den_mat_temp[pos_i][si][l][spin_l]
                                                + delta3_Red_Den_mat[pos_i][si][l][spin_l];
                                        derivative_val=iota*complex<double>(1.0*Parameters_.t_hopping,0.0)*temp_val;
                                        delta4_Red_Den_mat[pos_i][si][pos_j][sj] += dt_*derivative_val;
                                    }

                                }

                                //------------------l connection to l :done---------------------//


                                //-------------------coupling with Auxlliary spins--------------//

                                for(int s3=0;s3<2;s3++){
                                    double SX_i,SY_i,SZ_i,SX_j,SY_j,SZ_j;


                                    if(step_no==0){
                                        SX_i = Aux_S_x[ts][pos_i];
                                        SY_i = Aux_S_y[ts][pos_i];
                                        SZ_i = Aux_S_z[ts][pos_i];
                                        SX_j = Aux_S_x[ts][pos_j];
                                        SY_j = Aux_S_y[ts][pos_j];
                                        SZ_j = Aux_S_z[ts][pos_j];
                                        temp_val=Red_Den_mat_temp[pos_i][s3][pos_j][sj];
                                        derivative_val=iota*(-0.5*Jval_array[pos_i])*(
                                                    (SX_i)*Pauli_x[si][s3]
                                                    +
                                                    (SY_i)*Pauli_y[si][s3]
                                                    +
                                                    (SZ_i)*Pauli_z[si][s3]
                                                    )*temp_val;
                                        temp_val=Red_Den_mat_temp[pos_i][si][pos_j][s3];
                                        derivative_val +=iota*(0.5*Jval_array[pos_j])*(
                                                    (SX_j)*Pauli_x[s3][sj]
                                                    +
                                                    (SY_j)*Pauli_y[s3][sj]
                                                    +
                                                    (SZ_j)*Pauli_z[s3][sj]
                                                    )*temp_val;

                                        delta1_Red_Den_mat[pos_i][si][pos_j][sj] += dt_*derivative_val;

                                    }
                                    if(step_no==1){
                                        SX_i = Aux_S_x[ts][pos_i] + 0.5*delta1_Aux_Sx[pos_i];
                                        SY_i = Aux_S_y[ts][pos_i] + 0.5*delta1_Aux_Sy[pos_i];
                                        SZ_i = Aux_S_z[ts][pos_i] + 0.5*delta1_Aux_Sz[pos_i];
                                        SX_j = Aux_S_x[ts][pos_j] + 0.5*delta1_Aux_Sx[pos_j];
                                        SY_j = Aux_S_y[ts][pos_j] + 0.5*delta1_Aux_Sy[pos_j];
                                        SZ_j = Aux_S_z[ts][pos_j] + 0.5*delta1_Aux_Sz[pos_j];

                                        temp_val=Red_Den_mat_temp[pos_i][s3][pos_j][sj] +
                                                0.5*delta1_Red_Den_mat[pos_i][s3][pos_j][sj];

                                        derivative_val=iota*(-0.5*Jval_array[pos_i])*(
                                                    (SX_i)*Pauli_x[si][s3]
                                                    +
                                                    (SY_i)*Pauli_y[si][s3]
                                                    +
                                                    (SZ_i)*Pauli_z[si][s3]
                                                    )*temp_val;

                                        temp_val=Red_Den_mat_temp[pos_i][si][pos_j][s3] +
                                                0.5*delta1_Red_Den_mat[pos_i][si][pos_j][s3];

                                        derivative_val +=iota*(0.5*Jval_array[pos_j])*(
                                                    (SX_j)*Pauli_x[s3][sj]
                                                    +
                                                    (SY_j)*Pauli_y[s3][sj]
                                                    +
                                                    (SZ_j)*Pauli_z[s3][sj]
                                                    )*temp_val;

                                        delta2_Red_Den_mat[pos_i][si][pos_j][sj] += dt_*derivative_val;

                                    }
                                    if(step_no==2){
                                        SX_i = Aux_S_x[ts][pos_i] + 0.5*delta2_Aux_Sx[pos_i];
                                        SY_i = Aux_S_y[ts][pos_i] + 0.5*delta2_Aux_Sy[pos_i];
                                        SZ_i = Aux_S_z[ts][pos_i] + 0.5*delta2_Aux_Sz[pos_i];
                                        SX_j = Aux_S_x[ts][pos_j] + 0.5*delta2_Aux_Sx[pos_j];
                                        SY_j = Aux_S_y[ts][pos_j] + 0.5*delta2_Aux_Sy[pos_j];
                                        SZ_j = Aux_S_z[ts][pos_j] + 0.5*delta2_Aux_Sz[pos_j];

                                        temp_val=Red_Den_mat_temp[pos_i][s3][pos_j][sj] +
                                                0.5*delta2_Red_Den_mat[pos_i][s3][pos_j][sj];

                                        derivative_val=iota*(-0.5*Jval_array[pos_i])*(
                                                    (SX_i)*Pauli_x[si][s3]
                                                    +
                                                    (SY_i)*Pauli_y[si][s3]
                                                    +
                                                    (SZ_i)*Pauli_z[si][s3]
                                                    )*temp_val;

                                        temp_val=Red_Den_mat_temp[pos_i][si][pos_j][s3] +
                                                0.5*delta2_Red_Den_mat[pos_i][si][pos_j][s3];

                                        derivative_val +=iota*(0.5*Jval_array[pos_j])*(
                                                    (SX_j)*Pauli_x[s3][sj]
                                                    +
                                                    (SY_j)*Pauli_y[s3][sj]
                                                    +
                                                    (SZ_j)*Pauli_z[s3][sj]
                                                    )*temp_val;

                                        delta3_Red_Den_mat[pos_i][si][pos_j][sj] += dt_*derivative_val;

                                    }
                                    if(step_no==3){
                                        SX_i = Aux_S_x[ts][pos_i] + delta3_Aux_Sx[pos_i];
                                        SY_i = Aux_S_y[ts][pos_i] + delta3_Aux_Sy[pos_i];
                                        SZ_i = Aux_S_z[ts][pos_i] + delta3_Aux_Sz[pos_i];
                                        SX_j = Aux_S_x[ts][pos_j] + delta3_Aux_Sx[pos_j];
                                        SY_j = Aux_S_y[ts][pos_j] + delta3_Aux_Sy[pos_j];
                                        SZ_j = Aux_S_z[ts][pos_j] + delta3_Aux_Sz[pos_j];

                                        temp_val=Red_Den_mat_temp[pos_i][s3][pos_j][sj] +
                                                delta3_Red_Den_mat[pos_i][s3][pos_j][sj];

                                        derivative_val=iota*(-0.5*Jval_array[pos_i])*(
                                                    (SX_i)*Pauli_x[si][s3]
                                                    +
                                                    (SY_i)*Pauli_y[si][s3]
                                                    +
                                                    (SZ_i)*Pauli_z[si][s3]
                                                    )*temp_val;

                                        temp_val=Red_Den_mat_temp[pos_i][si][pos_j][s3] +
                                                delta3_Red_Den_mat[pos_i][si][pos_j][s3];

                                        derivative_val +=iota*(0.5*Jval_array[pos_j])*(
                                                    (SX_j)*Pauli_x[s3][sj]
                                                    +
                                                    (SY_j)*Pauli_y[s3][sj]
                                                    +
                                                    (SZ_j)*Pauli_z[s3][sj]
                                                    )*temp_val;

                                        delta4_Red_Den_mat[pos_i][si][pos_j][sj] += dt_*derivative_val;

                                    }
                                }
                            }
                        }
                    }
                }
            } //if(EVOLVE_RED_DEN_MAT==true)

        }


#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
        for(int pos=0;pos<Parameters_.ns;pos++){
            Aux_S_x[ts+1][pos] = Aux_S_x[ts][pos] + (1.0/6.0)*(delta1_Aux_Sx[pos] + 2.0*delta2_Aux_Sx[pos] + 2.0*delta3_Aux_Sx[pos] + delta4_Aux_Sx[pos]);
            Aux_S_y[ts+1][pos] = Aux_S_y[ts][pos] + (1.0/6.0)*(delta1_Aux_Sy[pos] + 2.0*delta2_Aux_Sy[pos] + 2.0*delta3_Aux_Sy[pos] + delta4_Aux_Sy[pos]);
            Aux_S_z[ts+1][pos] = Aux_S_z[ts][pos] + (1.0/6.0)*(delta1_Aux_Sz[pos] + 2.0*delta2_Aux_Sz[pos] + 2.0*delta3_Aux_Sz[pos] + delta4_Aux_Sz[pos]);
        }

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
        for(int pos_i=0;pos_i<Parameters_.ns;pos_i++){
            for(int si=0;si<2;si++){
                for(int pos_j=0;pos_j<Parameters_.ns;pos_j++){
                    for(int sj=0;sj<2;sj++){
                        Red_Den_mat[pos_i][si][pos_j][sj] = Red_Den_mat_temp[pos_i][si][pos_j][sj] + (1.0/6.0)*(delta1_Red_Den_mat[pos_i][si][pos_j][sj] + 2.0*delta2_Red_Den_mat[pos_i][si][pos_j][sj] +
                                                                                                                2.0*delta3_Red_Den_mat[pos_i][si][pos_j][sj] + delta4_Red_Den_mat[pos_i][si][pos_j][sj]);
                    }
                }
            }
        }

    }


    /*
    else if(Runge_Kutta_order==5){

        //Butcher’s Fifth Order Runge-Kutta Method is used

        //we let the quantum spins change as well in the intermediate steps of dt interval in RK-4 method,
        //here we not diagonalize the matrix again to calculate quantum spins for the changed classical spins.


        double pi,ti; //phi_i,theta_i,phi_j,theta_j for timestep=ts
        double sy,sx,sz; //Quantum spins for timestep=ts and position=pos

        double derivative_theta, derivative_phi;
        Mat_1_doub delta1_theta, delta2_theta, delta3_theta, delta4_theta, delta5_theta, delta6_theta;
        Mat_1_doub delta1_phi, delta2_phi, delta3_phi, delta4_phi, delta5_phi, delta6_phi;
        Mat_6_Complex_doub delta1_Red_Den_mat, delta2_Red_Den_mat, delta3_Red_Den_mat;
        Mat_6_Complex_doub delta4_Red_Den_mat, delta5_Red_Den_mat, delta6_Red_Den_mat;


        delta1_theta.resize(Parameters_.ns);delta2_theta.resize(Parameters_.ns);
        delta3_theta.resize(Parameters_.ns);delta4_theta.resize(Parameters_.ns);
        delta5_theta.resize(Parameters_.ns);delta6_theta.resize(Parameters_.ns);


        delta1_phi.resize(Parameters_.ns);delta2_phi.resize(Parameters_.ns);
        delta3_phi.resize(Parameters_.ns);delta4_phi.resize(Parameters_.ns);
        delta5_phi.resize(Parameters_.ns);delta6_phi.resize(Parameters_.ns);


        delta1_Red_Den_mat.resize(Parameters_.ns);delta2_Red_Den_mat.resize(Parameters_.ns);
        delta3_Red_Den_mat.resize(Parameters_.ns); delta4_Red_Den_mat.resize(Parameters_.ns);
        delta5_Red_Den_mat.resize(Parameters_.ns); delta6_Red_Den_mat.resize(Parameters_.ns);



        for(int i=0;i<Parameters_.ns;i++){
            delta1_Red_Den_mat[i].resize(Parameters_.orbs);delta2_Red_Den_mat[i].resize(Parameters_.orbs);
            delta3_Red_Den_mat[i].resize(Parameters_.orbs);delta4_Red_Den_mat[i].resize(Parameters_.orbs);
            delta5_Red_Den_mat[i].resize(Parameters_.orbs);delta6_Red_Den_mat[i].resize(Parameters_.orbs);
            for(int orb_i=0;orb_i<Parameters_.orbs;orb_i++){
                delta1_Red_Den_mat[i][orb_i].resize(2);delta2_Red_Den_mat[i][orb_i].resize(2);
                delta3_Red_Den_mat[i][orb_i].resize(2);delta4_Red_Den_mat[i][orb_i].resize(2);
                delta5_Red_Den_mat[i][orb_i].resize(2);delta6_Red_Den_mat[i][orb_i].resize(2);
                for(int si=0;si<2;si++){
                    delta1_Red_Den_mat[i][orb_i][si].resize(Parameters_.ns);delta2_Red_Den_mat[i][orb_i][si].resize(Parameters_.ns);
                    delta3_Red_Den_mat[i][orb_i][si].resize(Parameters_.ns);delta4_Red_Den_mat[i][orb_i][si].resize(Parameters_.ns);
                    delta5_Red_Den_mat[i][orb_i][si].resize(Parameters_.ns);delta6_Red_Den_mat[i][orb_i][si].resize(Parameters_.ns);
                    for(int j=0;j<Parameters_.ns;j++){
                        delta1_Red_Den_mat[i][orb_i][si][j].resize(Parameters_.orbs);delta2_Red_Den_mat[i][orb_i][si][j].resize(Parameters_.orbs);
                        delta3_Red_Den_mat[i][orb_i][si][j].resize(Parameters_.orbs);delta4_Red_Den_mat[i][orb_i][si][j].resize(Parameters_.orbs);
                        delta5_Red_Den_mat[i][orb_i][si][j].resize(Parameters_.orbs);delta6_Red_Den_mat[i][orb_i][si][j].resize(Parameters_.orbs);
                        for(int orb_j=0;orb_j<Parameters_.orbs;orb_j++){
                            delta1_Red_Den_mat[i][orb_i][si][j][orb_j].resize(2);delta2_Red_Den_mat[i][orb_i][si][j][orb_j].resize(2);
                            delta3_Red_Den_mat[i][orb_i][si][j][orb_j].resize(2);delta4_Red_Den_mat[i][orb_i][si][j][orb_j].resize(2);
                            delta5_Red_Den_mat[i][orb_i][si][j][orb_j].resize(2);delta6_Red_Den_mat[i][orb_i][si][j][orb_j].resize(2);
                        }
                    }
                }
            }
        }


        for(int i=0;i<Parameters_.ns;i++){
            delta1_theta[i]=0.0;delta2_theta[i]=0.0;
            delta3_theta[i]=0.0;delta4_theta[i]=0.0;
            delta5_theta[i]=0.0;delta6_theta[i]=0.0;
            delta1_phi[i]=0.0;delta2_phi[i]=0.0;
            delta3_phi[i]=0.0;delta4_phi[i]=0.0;
            delta5_phi[i]=0.0;delta6_phi[i]=0.0;
            for(int orb_i=0;orb_i<Parameters_.orbs;orb_i++){
                for(int si=0;si<2;si++){
                    for(int j=0;j<Parameters_.ns;j++){
                        for(int orb_j=0;orb_j<Parameters_.orbs;orb_j++){
                            for(int sj=0;sj<2;sj++){
                                delta1_Red_Den_mat[i][orb_i][si][j][orb_j][sj]=zero;delta2_Red_Den_mat[i][orb_i][si][j][orb_j][sj]=zero;
                                delta3_Red_Den_mat[i][orb_i][si][j][orb_j][sj]=zero;delta4_Red_Den_mat[i][orb_i][si][j][orb_j][sj]=zero;
                                delta5_Red_Den_mat[i][orb_i][si][j][orb_j][sj]=zero;delta6_Red_Den_mat[i][orb_i][si][j][orb_j][sj]=zero;

                            }
                        }
                    }
                }

            }

        }



        //calculating delta1(2,3,4,5,6)_thetha(phi)
        bool intermediate_update = true;
        for(int step_no=0;step_no<6;step_no++){

            //Classical spins time evolution, due to coupling with Quantum spins and Magnetic field
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(derivative_theta,derivative_phi, pi, ti, sx, sy, sz)
#endif
            for(int pos=0;pos<Parameters_.ns;pos++){

                //if(omp_get_thread_num()==1){
                //cout<<pos<<"in thread = "<<omp_get_thread_num()<<endl;}
                int pos_x = Coordinates_.indx(pos);
                int pos_y = Coordinates_.indy(pos);
                derivative_theta=0;
                derivative_phi=0;

                if(step_no==0){
                    pi=Phi[ts][pos_x][pos_y];
                    ti=Theta[ts][pos_x][pos_y];
                    sx=0.0;sy=0.0;sz=0.0;
                    for(int orb=0;orb<Parameters_.orbs;orb++){
                        sx +=0.5*real( Red_Den_mat_temp[pos][orb][1][pos][orb][0] + Red_Den_mat_temp[pos][orb][0][pos][orb][1] );
                        sy +=0.5*imag( Red_Den_mat_temp[pos][orb][1][pos][orb][0] - Red_Den_mat_temp[pos][orb][0][pos][orb][1] );
                        sz +=0.5*real( Red_Den_mat_temp[pos][orb][0][pos][orb][0] - Red_Den_mat_temp[pos][orb][1][pos][orb][1] );
                    }
                }
                if(step_no==1){
                    pi=Phi[ts][pos_x][pos_y] + (1.0/4.0)*(delta1_phi[pos]);
                    ti=Theta[ts][pos_x][pos_y] + (1.0/4.0)*(delta1_theta[pos]);
                    sx=0.0;sy=0.0;sz=0.0;
                    for(int orb=0;orb<Parameters_.orbs;orb++){
                        sx += 0.5*real( Red_Den_mat_temp[pos][orb][1][pos][orb][0] + Red_Den_mat_temp[pos][orb][0][pos][orb][1]
                                + (1.0/4.0)*(delta1_Red_Den_mat[pos][orb][1][pos][orb][0] + delta1_Red_Den_mat[pos][orb][0][pos][orb][1]));
                        sy +=0.5*imag( Red_Den_mat_temp[pos][orb][1][pos][orb][0] - Red_Den_mat_temp[pos][orb][0][pos][orb][1]
                                + (1.0/4.0)*(delta1_Red_Den_mat[pos][orb][1][pos][orb][0] - delta1_Red_Den_mat[pos][orb][0][pos][orb][1])
                                );
                        sz +=0.5*real( Red_Den_mat_temp[pos][orb][0][pos][orb][0] - Red_Den_mat_temp[pos][orb][1][pos][orb][1]
                                + (1.0/4.0)*(delta1_Red_Den_mat[pos][orb][0][pos][orb][0] - delta1_Red_Den_mat[pos][orb][1][pos][orb][1])
                                );
                    }

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
                    pi=Phi[ts][pos_x][pos_y] + (1.0/8.0)*( (delta1_phi[pos]) + (delta2_phi[pos]) );
                    ti=Theta[ts][pos_x][pos_y] + (1.0/8.0)*( (delta1_theta[pos]) + (delta2_theta[pos]) );

                    sx=0.0;sy=0.0;sz=0.0;
                    for(int orb=0;orb<Parameters_.orbs;orb++){
                        sx +=0.5*real( Red_Den_mat_temp[pos][orb][1][pos][orb][0] + Red_Den_mat_temp[pos][orb][0][pos][orb][1]
                                + ((1.0/8.0)*(delta1_Red_Den_mat[pos][orb][1][pos][orb][0] + delta1_Red_Den_mat[pos][orb][0][pos][orb][1]))
                                + ((1.0/8.0)*(delta2_Red_Den_mat[pos][orb][1][pos][orb][0] + delta2_Red_Den_mat[pos][orb][0][pos][orb][1]))
                                );
                        sy +=0.5*imag( Red_Den_mat_temp[pos][orb][1][pos][orb][0] - Red_Den_mat_temp[pos][orb][0][pos][orb][1]
                                + ((1.0/8.0)*(delta1_Red_Den_mat[pos][orb][1][pos][orb][0] - delta1_Red_Den_mat[pos][orb][0][pos][orb][1]))
                                + ((1.0/8.0)*(delta2_Red_Den_mat[pos][orb][1][pos][orb][0] - delta2_Red_Den_mat[pos][orb][0][pos][orb][1]))
                                );
                        sz +=0.5*real( Red_Den_mat_temp[pos][orb][0][pos][orb][0] - Red_Den_mat_temp[pos][orb][1][pos][orb][1]
                                + ((1.0/8.0)*(delta1_Red_Den_mat[pos][orb][0][pos][orb][0] - delta1_Red_Den_mat[pos][orb][1][pos][orb][1]))
                                + ((1.0/8.0)*(delta2_Red_Den_mat[pos][orb][0][pos][orb][0] - delta2_Red_Den_mat[pos][orb][1][pos][orb][1]))
                                );
                    }

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
                    pi=Phi[ts][pos_x][pos_y] + (delta3_phi[pos]) - ((1.0/2.0)*delta2_phi[pos]);
                    ti=Theta[ts][pos_x][pos_y] + (delta3_theta[pos])  - ((1.0/2.0)*delta2_theta[pos]);

                    sx=0.0;sy=0.0;sz=0.0;
                    for(int orb=0;orb<Parameters_.orbs;orb++){
                        sx +=0.5*real( Red_Den_mat_temp[pos][orb][1][pos][orb][0] + Red_Den_mat_temp[pos][orb][0][pos][orb][1]
                                + (delta3_Red_Den_mat[pos][orb][1][pos][orb][0] + delta3_Red_Den_mat[pos][orb][0][pos][orb][1])
                                - ((1.0/2.0)*(delta2_Red_Den_mat[pos][orb][1][pos][orb][0] + delta2_Red_Den_mat[pos][orb][0][pos][orb][1]))

                                );
                        sy +=0.5*imag( Red_Den_mat_temp[pos][orb][1][pos][orb][0] - Red_Den_mat_temp[pos][orb][0][pos][orb][1]
                                + (delta3_Red_Den_mat[pos][orb][1][pos][orb][0] - delta3_Red_Den_mat[pos][orb][0][pos][orb][1])
                                - ((1.0/2.0)*(delta2_Red_Den_mat[pos][orb][1][pos][orb][0] - delta2_Red_Den_mat[pos][orb][0][pos][orb][1]))
                                );
                        sz +=0.5*real( Red_Den_mat_temp[pos][orb][0][pos][orb][0] - Red_Den_mat_temp[pos][orb][1][pos][orb][1]
                                + (delta3_Red_Den_mat[pos][orb][0][pos][orb][0] - delta3_Red_Den_mat[pos][orb][1][pos][orb][1])
                                - ((1.0/2.0)*(delta2_Red_Den_mat[pos][orb][0][pos][orb][0] - delta2_Red_Den_mat[pos][orb][1][pos][orb][1]))
                                );
                    }


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

                if(step_no==4){
                    pi=Phi[ts][pos_x][pos_y] + ((3.0/16.0)*(delta1_phi[pos])) + ((9.0/16.0)*delta4_phi[pos]);
                    ti=Theta[ts][pos_x][pos_y] + ((3.0/16.0)*delta1_theta[pos])  + ((9.0/16.0)*delta4_theta[pos]);

                    sx=0.0;sy=0.0;sz=0.0;
                    for(int orb=0;orb<Parameters_.orbs;orb++){
                        sx +=0.5*real( Red_Den_mat_temp[pos][orb][1][pos][orb][0] + Red_Den_mat_temp[pos][orb][0][pos][orb][1]
                                + ((3.0/16.0)*(delta1_Red_Den_mat[pos][orb][1][pos][orb][0] + delta1_Red_Den_mat[pos][orb][0][pos][orb][1]))
                                + ((9.0/16.0)*(delta4_Red_Den_mat[pos][orb][1][pos][orb][0] + delta4_Red_Den_mat[pos][orb][0][pos][orb][1]))
                                );
                        sy +=0.5*imag( Red_Den_mat_temp[pos][orb][1][pos][orb][0] - Red_Den_mat_temp[pos][orb][0][pos][orb][1]
                                + ((3.0/16.0)*(delta1_Red_Den_mat[pos][orb][1][pos][orb][0] - delta1_Red_Den_mat[pos][orb][0][pos][orb][1]))
                                + ((9.0/16.0)*(delta4_Red_Den_mat[pos][orb][1][pos][orb][0] - delta4_Red_Den_mat[pos][orb][0][pos][orb][1]))
                                );
                        sz +=0.5*real( Red_Den_mat_temp[pos][orb][0][pos][orb][0] - Red_Den_mat_temp[pos][orb][1][pos][orb][1]
                                + ((3.0/16.0)*(delta1_Red_Den_mat[pos][orb][0][pos][orb][0] - delta1_Red_Den_mat[pos][orb][1][pos][orb][1]))
                                + ((9.0/16.0)*(delta4_Red_Den_mat[pos][orb][0][pos][orb][0] - delta4_Red_Den_mat[pos][orb][1][pos][orb][1]))
                                );
                    }


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

                if(step_no==5){
                    pi=Phi[ts][pos_x][pos_y] + ((-3.0/7.0)*(delta1_phi[pos])) + ((2.0/7.0)*(delta2_phi[pos])) +
                            ((12.0/7.0)*(delta3_phi[pos])) + ((-12.0/7.0)*(delta4_phi[pos]))
                            + ((8.0/7.0)*(delta5_phi[pos]));
                    ti=Theta[ts][pos_x][pos_y] + ((-3.0/7.0)*(delta1_theta[pos])) + ((2.0/7.0)*(delta2_theta[pos])) +
                            ((12.0/7.0)*(delta3_theta[pos])) + ((-12.0/7.0)*(delta4_theta[pos]))
                            + ((8.0/7.0)*(delta5_theta[pos]));;

                    sx=0.0;sy=0.0;sz=0.0;
                    for(int orb=0;orb<Parameters_.orbs;orb++){
                        sx +=0.5*real( Red_Den_mat_temp[pos][orb][1][pos][orb][0] + Red_Den_mat_temp[pos][orb][0][pos][orb][1]
                                + ((-3.0/7.0)*(delta1_Red_Den_mat[pos][orb][1][pos][orb][0] + delta1_Red_Den_mat[pos][orb][0][pos][orb][1]))
                                + ((2.0/7.0)*(delta2_Red_Den_mat[pos][orb][1][pos][orb][0] + delta2_Red_Den_mat[pos][orb][0][pos][orb][1]))
                                + ((12.0/7.0)*(delta3_Red_Den_mat[pos][orb][1][pos][orb][0] + delta3_Red_Den_mat[pos][orb][0][pos][orb][1]))
                                + ((-12.0/7.0)*(delta4_Red_Den_mat[pos][orb][1][pos][orb][0] + delta4_Red_Den_mat[pos][orb][0][pos][orb][1]))
                                + ((8.0/7.0)*(delta5_Red_Den_mat[pos][orb][1][pos][orb][0] + delta5_Red_Den_mat[pos][orb][0][pos][orb][1]))
                                );
                        sy +=0.5*imag( Red_Den_mat_temp[pos][orb][1][pos][orb][0] - Red_Den_mat_temp[pos][orb][0][pos][orb][1]
                                + ((-3.0/7.0)*(delta1_Red_Den_mat[pos][orb][1][pos][orb][0] - delta1_Red_Den_mat[pos][orb][0][pos][orb][1]))
                                + ((2.0/7.0)*(delta2_Red_Den_mat[pos][orb][1][pos][orb][0] - delta2_Red_Den_mat[pos][orb][0][pos][orb][1]))
                                + ((12.0/7.0)*(delta3_Red_Den_mat[pos][orb][1][pos][orb][0] - delta3_Red_Den_mat[pos][orb][0][pos][orb][1]))
                                + ((-12.0/7.0)*(delta4_Red_Den_mat[pos][orb][1][pos][orb][0] - delta4_Red_Den_mat[pos][orb][0][pos][orb][1]))
                                + ((8.0/7.0)*(delta5_Red_Den_mat[pos][orb][1][pos][orb][0] - delta5_Red_Den_mat[pos][orb][0][pos][orb][1]))
                                );
                        sz +=0.5*real( Red_Den_mat_temp[pos][orb][0][pos][orb][0] - Red_Den_mat_temp[pos][orb][1][pos][orb][1]
                                + ((-3.0/7.0)*(delta1_Red_Den_mat[pos][orb][0][pos][orb][0] - delta1_Red_Den_mat[pos][orb][1][pos][orb][1]))
                                + ((2.0/7.0)*(delta2_Red_Den_mat[pos][orb][0][pos][orb][0] - delta2_Red_Den_mat[pos][orb][1][pos][orb][1]))
                                + ((12.0/7.0)*(delta3_Red_Den_mat[pos][orb][0][pos][orb][0] - delta3_Red_Den_mat[pos][orb][1][pos][orb][1]))
                                + ((-12.0/7.0)*(delta4_Red_Den_mat[pos][orb][0][pos][orb][0] - delta4_Red_Den_mat[pos][orb][1][pos][orb][1]))
                                + ((8.0/7.0)*(delta5_Red_Den_mat[pos][orb][0][pos][orb][0] - delta5_Red_Den_mat[pos][orb][1][pos][orb][1]))
                                );
                    }


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


                derivative_phi = (-1.0)*( Jval_array[pos]*((((sin(pi)*cos(ti))/(sin(ti)))*sy) +
                                                           (((cos(pi)*cos(ti))/(sin(ti)))*sx) -  sz ))  -  Bval_array[pos]; //-Bval_array[pos]
                derivative_theta = ( Jval_array[pos]*((cos(pi)*sy) - (sin(pi)*sx)));

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
                if(step_no==4){
                    delta5_phi[pos] = (dt_*(derivative_phi ));
                    delta5_theta[pos] = (dt_*( derivative_theta));

                }
                if(step_no==5){
                    delta6_phi[pos] = (dt_*(derivative_phi ));
                    delta6_theta[pos] = (dt_*( derivative_theta));

                }

                //Phi[ts+1][pos] = Phi[ts][pos] ;
                //Theta[ts+1][pos] = Theta[ts][pos] ;




                int pos_ng;
                double pj,tj, value_;
                for(int ng=0;ng<8;ng++){

                    if(ng<2){
                        value_=1.14*Parameters_.J_NN;  // +x, -x
                    }
                    else if(ng>=2 && ng <4){
                        value_=Parameters_.J_NN;   // +y, -y
                    }
                    else{
                        value_=Parameters_.J_NNN;  //+x+y, +x-y, -x+y, -x-y
                    }

                    pos_ng = Coordinates_.neigh(pos,ng);

                    int pos_ng_x = Coordinates_.indx(pos_ng);
                    int pos_ng_y = Coordinates_.indy(pos_ng);


                    if(step_no==0){
                        pj=Phi[ts][pos_ng_x][pos_ng_y];
                        tj=Theta[ts][pos_ng_x][pos_ng_y];
                    }


                    if(step_no==1){
                        pj=Phi[ts][pos_ng_x][pos_ng_y] + (1.0/4.0)*(delta1_phi[pos_ng]);
                        tj=Theta[ts][pos_ng_x][pos_ng_y] + (1.0/4.0)*(delta1_theta[pos_ng]);

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
                        pj=Phi[ts][pos_ng_x][pos_ng_y] + (1.0/8.0)*(delta1_phi[pos_ng] + delta2_phi[pos_ng]);
                        tj=Theta[ts][pos_ng_x][pos_ng_y] + (1.0/8.0)*(delta1_theta[pos_ng] +delta2_theta[pos_ng] );

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
                        pj=Phi[ts][pos_ng_x][pos_ng_y] + (delta3_phi[pos_ng] - ((1.0/2.0)*delta2_phi[pos_ng]));
                        tj=Theta[ts][pos_ng_x][pos_ng_y] + (delta3_theta[pos_ng] - ((1.0/2.0)*delta2_theta[pos_ng]));

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

                    if(step_no==4){
                        pj=Phi[ts][pos_ng_x][pos_ng_y] + (((3.0/16.0)*delta1_phi[pos_ng]) + ((9.0/16.0)*delta4_phi[pos_ng]));
                        tj=Theta[ts][pos_ng_x][pos_ng_y] + (((3.0/16.0)*delta1_theta[pos_ng]) + ((9.0/16.0)*delta4_theta[pos_ng]));

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

                    if(step_no==5){
                        pj=Phi[ts][pos_ng_x][pos_ng_y] + ((-3.0/7.0)*delta1_phi[pos_ng]) + ((2.0/7.0)*delta2_phi[pos_ng]) +
                                ((12.0/7.0)*delta3_phi[pos_ng])
                                + ((-12.0/7.0)*delta4_phi[pos_ng]) + ((8.0/7.0)*delta5_phi[pos_ng]);
                        tj=Theta[ts][pos_ng_x][pos_ng_y] + ((-3.0/7.0)*delta1_theta[pos_ng]) + ((2.0/7.0)*delta2_theta[pos_ng]) +
                                ((12.0/7.0)*delta3_theta[pos_ng])
                                + ((-12.0/7.0)*delta4_theta[pos_ng]) + ((8.0/7.0)*delta5_theta[pos_ng]);

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

                    derivative_phi = (-1.0)*(( (S_mag)*(value_)*(((sin(tj)*cos(ti)*cos(pi - pj))/(sin(ti))) - cos(tj) ) ));
                    derivative_theta = (S_mag)*(value_)*(sin(tj)*sin(pj - pi));


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

                    if(step_no==4){
                        delta5_phi[pos] += (dt_*derivative_phi );
                        delta5_theta[pos] += (dt_*derivative_theta);
                    }

                    if(step_no==5){
                        delta6_phi[pos] += (dt_*derivative_phi );
                        delta6_theta[pos] += (dt_*derivative_theta);
                    }

                }

            }

        }


#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
        for(int pos=0;pos<Parameters_.ns;pos++){

            int pos_x= Coordinates_.indx(pos);
            int pos_y= Coordinates_.indy(pos);


            Phi[ts+1][pos_x][pos_y] = Phi[ts][pos_x][pos_y]  + (1.0/90.0)*( ((7.0)*delta1_phi[pos]) + ((32.0)*delta3_phi[pos]) + ((12.0)*delta4_phi[pos])
                                                                            + ((32.0)*delta5_phi[pos]) + ((7.0)*delta6_phi[pos]) );
            Theta[ts+1][pos_x][pos_y]  = Theta[ts][pos_x][pos_y]  + (1.0/90.0)*( ((7.0)*delta1_theta[pos]) + ((32.0)*delta3_theta[pos]) + ((12.0)*delta4_theta[pos])
                                                                                 + ((32.0)*delta5_theta[pos]) + ((7.0)*delta6_theta[pos]) );



            //Phi[ts+1][pos] = Phi[ts][pos] + (1.0/1.0)*(0.0*delta1_phi[pos] + 1.0*delta2_phi[pos] + 0.0*delta3_phi[pos] + 0.0*delta4_phi[pos]);
            //Theta[ts+1][pos] = Theta[ts][pos] + (1.0/1.0)*(0.0*delta1_theta[pos] + 1.0*delta2_theta[pos] + 0.0*delta3_theta[pos] + 0.0*delta4_theta[pos]);


            if(Phi[ts+1][pos_x][pos_y]  > 2*PI){
                Phi[ts+1][pos_x][pos_y]  += -2*PI;

            }
            if(Phi[ts+1][pos_x][pos_y]  < 0){
                Phi[ts+1][pos_x][pos_y]  +=  2*PI;

            }


            if(Theta[ts+1][pos_x][pos_y]  > PI){
                Theta[ts+1][pos_x][pos_y]  = Theta[ts+1][pos_x][pos_y]  -2*PI;
                Phi[ts+1][pos_x][pos_y]  = fmod( Phi[ts+1][pos_x][pos_y]  + PI, 2.0*PI );
            }
            if(Theta[ts+1][pos_x][pos_y]  < 0){
                Theta[ts+1][pos_x][pos_y]  = - Theta[ts+1][pos_x][pos_y] ;
                Phi[ts+1][pos_x][pos_y]  = fmod( Phi[ts+1][pos_x][pos_y]  + PI, 2.0*PI );

            }

            //cout<<delta1_phi[pos]<<"   "<<delta2_phi[pos]<<"   "<<delta3_phi[pos]<<"   "<<delta4_phi[pos]<<"   "<<(1.0/6.0)*(delta1_phi[pos] + 2.0*delta2_phi[pos] + 2.0*delta3_phi[pos] + delta4_phi[pos])<<endl;
        }

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
        for(int pos_i=0;pos_i<Parameters_.ns;pos_i++){
            for (int orb_i=0;orb_i<Parameters_.orbs;orb_i++){
                for(int si=0;si<2;si++){
                    for(int pos_j=0;pos_j<Parameters_.ns;pos_j++){
                        for (int orb_j=0;orb_j<Parameters_.orbs;orb_j++){
                            for(int sj=0;sj<2;sj++){
                                Red_Den_mat[pos_i][orb_i][si][pos_j][orb_j][sj] = Red_Den_mat_temp[pos_i][orb_i][si][pos_j][orb_j][sj] + (1.0/90.0)*(7.0*delta1_Red_Den_mat[pos_i][orb_i][si][pos_j][orb_j][sj] + 32.0*delta3_Red_Den_mat[pos_i][orb_i][si][pos_j][orb_j][sj] +
                                                                                                                                                     12.0*delta4_Red_Den_mat[pos_i][orb_i][si][pos_j][orb_j][sj] + 32.0*delta5_Red_Den_mat[pos_i][orb_i][si][pos_j][orb_j][sj]
                                                                                                                                                     + 7.0*delta6_Red_Den_mat[pos_i][orb_i][si][pos_j][orb_j][sj]);
                            }
                        }
                    }
                }
            }
        }

    }

*/


    /*
    else if(Runge_Kutta_order==6){

        //Sixth Order Runge-Kutta Method is used
        //From "https://www.ams.org/journals/mcom/1968-22-102/S0025-5718-68-99876-1/S0025-5718-68-99876-1.pdf"
        //or google "An Explicit Sixth-Order Runge-Kutta Formula By H. A. Luther"

        double nu=1.0;
        double sign_=-1.0;

        //we let the quantum spins change as well in the intermediate steps of dt interval in RK-4 method,
        //here we not diagonalize the matrix again to calculate quantum spins for the changed classical spins.


        double pi,ti; //phi_i,theta_i,phi_j,theta_j for timestep=ts
        double sy,sx,sz; //Quantum spins for timestep=ts and position=pos

        double derivative_theta, derivative_phi;
        Mat_1_doub delta1_theta, delta2_theta, delta3_theta, delta4_theta;
        Mat_1_doub delta5_theta, delta6_theta, delta7_theta;
        Mat_1_doub delta1_phi, delta2_phi, delta3_phi, delta4_phi;
        Mat_1_doub delta5_phi, delta6_phi, delta7_phi;
        Mat_6_Complex_doub delta1_Red_Den_mat, delta2_Red_Den_mat, delta3_Red_Den_mat;
        Mat_6_Complex_doub delta4_Red_Den_mat, delta5_Red_Den_mat, delta6_Red_Den_mat, delta7_Red_Den_mat;


        delta1_theta.resize(Parameters_.ns);delta2_theta.resize(Parameters_.ns);
        delta3_theta.resize(Parameters_.ns);delta4_theta.resize(Parameters_.ns);
        delta5_theta.resize(Parameters_.ns);delta6_theta.resize(Parameters_.ns);
        delta7_theta.resize(Parameters_.ns);


        delta1_phi.resize(Parameters_.ns);delta2_phi.resize(Parameters_.ns);
        delta3_phi.resize(Parameters_.ns);delta4_phi.resize(Parameters_.ns);
        delta5_phi.resize(Parameters_.ns);delta6_phi.resize(Parameters_.ns);
        delta7_phi.resize(Parameters_.ns);


        delta1_Red_Den_mat.resize(Parameters_.ns); delta2_Red_Den_mat.resize(Parameters_.ns);
        delta3_Red_Den_mat.resize(Parameters_.ns); delta4_Red_Den_mat.resize(Parameters_.ns);
        delta5_Red_Den_mat.resize(Parameters_.ns); delta6_Red_Den_mat.resize(Parameters_.ns);
        delta7_Red_Den_mat.resize(Parameters_.ns);



        for(int i=0;i<Parameters_.ns;i++){
            delta1_Red_Den_mat[i].resize(Parameters_.orbs);delta2_Red_Den_mat[i].resize(Parameters_.orbs);
            delta3_Red_Den_mat[i].resize(Parameters_.orbs);delta4_Red_Den_mat[i].resize(Parameters_.orbs);
            delta5_Red_Den_mat[i].resize(Parameters_.orbs);delta6_Red_Den_mat[i].resize(Parameters_.orbs);
            delta7_Red_Den_mat[i].resize(Parameters_.orbs);
            for(int orb_i=0;orb_i<Parameters_.orbs;orb_i++){
                delta1_Red_Den_mat[i][orb_i].resize(2);delta2_Red_Den_mat[i][orb_i].resize(2);
                delta3_Red_Den_mat[i][orb_i].resize(2);delta4_Red_Den_mat[i][orb_i].resize(2);
                delta5_Red_Den_mat[i][orb_i].resize(2);delta6_Red_Den_mat[i][orb_i].resize(2);
                delta7_Red_Den_mat[i][orb_i].resize(2);
                for(int si=0;si<2;si++){
                    delta1_Red_Den_mat[i][orb_i][si].resize(Parameters_.ns);delta2_Red_Den_mat[i][orb_i][si].resize(Parameters_.ns);
                    delta3_Red_Den_mat[i][orb_i][si].resize(Parameters_.ns);delta4_Red_Den_mat[i][orb_i][si].resize(Parameters_.ns);
                    delta5_Red_Den_mat[i][orb_i][si].resize(Parameters_.ns);delta6_Red_Den_mat[i][orb_i][si].resize(Parameters_.ns);
                    delta7_Red_Den_mat[i][orb_i][si].resize(Parameters_.ns);
                    for(int j=0;j<Parameters_.ns;j++){
                        delta1_Red_Den_mat[i][orb_i][si][j].resize(Parameters_.orbs);delta2_Red_Den_mat[i][orb_i][si][j].resize(Parameters_.orbs);
                        delta3_Red_Den_mat[i][orb_i][si][j].resize(Parameters_.orbs);delta4_Red_Den_mat[i][orb_i][si][j].resize(Parameters_.orbs);
                        delta5_Red_Den_mat[i][orb_i][si][j].resize(Parameters_.orbs);delta6_Red_Den_mat[i][orb_i][si][j].resize(Parameters_.orbs);
                        delta7_Red_Den_mat[i][orb_i][si][j].resize(Parameters_.orbs);
                        for(int orb_j=0;orb_j<Parameters_.orbs;orb_j++){
                            delta1_Red_Den_mat[i][orb_i][si][j][orb_j].resize(2);delta2_Red_Den_mat[i][orb_i][si][j][orb_j].resize(2);
                            delta3_Red_Den_mat[i][orb_i][si][j][orb_j].resize(2);delta4_Red_Den_mat[i][orb_i][si][j][orb_j].resize(2);
                            delta5_Red_Den_mat[i][orb_i][si][j][orb_j].resize(2);delta6_Red_Den_mat[i][orb_i][si][j][orb_j].resize(2);
                            delta7_Red_Den_mat[i][orb_i][si][j][orb_j].resize(2);
                        }
                    }
                }
            }
        }


        for(int i=0;i<Parameters_.ns;i++){
            delta1_theta[i]=0.0;delta2_theta[i]=0.0;
            delta3_theta[i]=0.0;delta4_theta[i]=0.0;
            delta5_theta[i]=0.0;delta6_theta[i]=0.0;delta7_theta[i]=0.0;
            delta1_phi[i]=0.0;delta2_phi[i]=0.0;
            delta3_phi[i]=0.0;delta4_phi[i]=0.0;
            delta5_phi[i]=0.0;delta6_phi[i]=0.0;delta7_phi[i]=0.0;
            for(int orb_i=0;orb_i<Parameters_.orbs;orb_i++){
                for(int si=0;si<2;si++){
                    for(int j=0;j<Parameters_.ns;j++){
                        for(int orb_j=0;orb_j<Parameters_.orbs;orb_j++){
                            for(int sj=0;sj<2;sj++){
                                delta1_Red_Den_mat[i][orb_i][si][j][orb_j][sj]=zero;delta2_Red_Den_mat[i][orb_i][si][j][orb_j][sj]=zero;
                                delta3_Red_Den_mat[i][orb_i][si][j][orb_j][sj]=zero;delta4_Red_Den_mat[i][orb_i][si][j][orb_j][sj]=zero;
                                delta5_Red_Den_mat[i][orb_i][si][j][orb_j][sj]=zero;delta6_Red_Den_mat[i][orb_i][si][j][orb_j][sj]=zero;
                                delta7_Red_Den_mat[i][orb_i][si][j][orb_j][sj]=zero;

                            }
                        }
                    }
                }

            }

        }



        //calculating delta1(2,3,4,5,6,7)_thetha(phi)
        bool intermediate_update = true;
        for(int step_no=0;step_no<7;step_no++){

            //Classical spins time evolution, due to coupling with Quantum spins and Magnetic field
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(derivative_theta,derivative_phi, pi, ti, sx, sy, sz)
#endif
            for(int pos=0;pos<Parameters_.ns;pos++){

                //if(omp_get_thread_num()==1){
                //cout<<pos<<"in thread = "<<omp_get_thread_num()<<endl;}
                int pos_x = Coordinates_.indx(pos);
                int pos_y = Coordinates_.indy(pos);
                derivative_theta=0;
                derivative_phi=0;

                if(step_no==0){
                    pi=Phi[ts][pos_x][pos_y];
                    ti=Theta[ts][pos_x][pos_y];
                    sx=0.0;sy=0.0;sz=0.0;
                    for(int orb=0;orb<Parameters_.orbs;orb++){
                        sx +=0.5*real( Red_Den_mat_temp[pos][orb][1][pos][orb][0] + Red_Den_mat_temp[pos][orb][0][pos][orb][1] );
                        sy +=0.5*imag( Red_Den_mat_temp[pos][orb][1][pos][orb][0] - Red_Den_mat_temp[pos][orb][0][pos][orb][1] );
                        sz +=0.5*real( Red_Den_mat_temp[pos][orb][0][pos][orb][0] - Red_Den_mat_temp[pos][orb][1][pos][orb][1] );
                    }
                }
                if(step_no==1){
                    pi=Phi[ts][pos_x][pos_y] + (nu)*(delta1_phi[pos]);
                    ti=Theta[ts][pos_x][pos_y] + (nu)*(delta1_theta[pos]);
                    sx=0.0;sy=0.0;sz=0.0;
                    for(int orb=0;orb<Parameters_.orbs;orb++){
                        sx += 0.5*real( Red_Den_mat_temp[pos][orb][1][pos][orb][0] + Red_Den_mat_temp[pos][orb][0][pos][orb][1]
                                + (nu)*(delta1_Red_Den_mat[pos][orb][1][pos][orb][0] + delta1_Red_Den_mat[pos][orb][0][pos][orb][1]));
                        sy +=0.5*imag( Red_Den_mat_temp[pos][orb][1][pos][orb][0] - Red_Den_mat_temp[pos][orb][0][pos][orb][1]
                                + (nu)*(delta1_Red_Den_mat[pos][orb][1][pos][orb][0] - delta1_Red_Den_mat[pos][orb][0][pos][orb][1])
                                );
                        sz +=0.5*real( Red_Den_mat_temp[pos][orb][0][pos][orb][0] - Red_Den_mat_temp[pos][orb][1][pos][orb][1]
                                + (nu)*(delta1_Red_Den_mat[pos][orb][0][pos][orb][0] - delta1_Red_Den_mat[pos][orb][1][pos][orb][1])
                                );
                    }

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
                    pi=Phi[ts][pos_x][pos_y] + ( ((4*nu) - 1.0)*(delta1_phi[pos]) + (1.0)*(delta2_phi[pos]) )*(1.0/(8.0*nu));
                    ti=Theta[ts][pos_x][pos_y] + ( ((4*nu) - 1.0)*(delta1_theta[pos]) + (1.0)*(delta2_theta[pos]) )*(1.0/(8.0*nu));

                    sx=0.0;sy=0.0;sz=0.0;
                    for(int orb=0;orb<Parameters_.orbs;orb++){
                        sx +=0.5*real( Red_Den_mat_temp[pos][orb][1][pos][orb][0] + Red_Den_mat_temp[pos][orb][0][pos][orb][1]
                                +
                                ( (((4*nu) - 1.0)*(delta1_Red_Den_mat[pos][orb][1][pos][orb][0] + delta1_Red_Den_mat[pos][orb][0][pos][orb][1]))
                                + ((1.0)*(delta2_Red_Den_mat[pos][orb][1][pos][orb][0] + delta2_Red_Den_mat[pos][orb][0][pos][orb][1]))
                                )*(1.0/(8.0*nu))
                                );
                        sy +=0.5*imag( Red_Den_mat_temp[pos][orb][1][pos][orb][0] - Red_Den_mat_temp[pos][orb][0][pos][orb][1]
                                +
                                ( (((4*nu) - 1.0)*(delta1_Red_Den_mat[pos][orb][1][pos][orb][0] - delta1_Red_Den_mat[pos][orb][0][pos][orb][1]))
                                + ((1.0)*(delta2_Red_Den_mat[pos][orb][1][pos][orb][0] - delta2_Red_Den_mat[pos][orb][0][pos][orb][1]))
                                )*(1.0/(8.0*nu))
                                );
                        sz +=0.5*real( Red_Den_mat_temp[pos][orb][0][pos][orb][0] - Red_Den_mat_temp[pos][orb][1][pos][orb][1]
                                +
                                ( (((4*nu) - 1.0)*(delta1_Red_Den_mat[pos][orb][0][pos][orb][0] - delta1_Red_Den_mat[pos][orb][1][pos][orb][1]))
                                + ((1.0)*(delta2_Red_Den_mat[pos][orb][0][pos][orb][0] - delta2_Red_Den_mat[pos][orb][1][pos][orb][1]))
                                )*(1.0/(8.0*nu))
                                );
                    }

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
                    pi=Phi[ts][pos_x][pos_y] + (1.0/(27*nu))*( (( (10.0*nu) - 2.0)*delta1_phi[pos]) +
                                                               ((2.0)*delta2_phi[pos]) + ((8*nu)*delta3_phi[pos]) );
                    ti=Theta[ts][pos_x][pos_y] + (1.0/(27*nu))*( (( (10.0*nu) - 2.0)*delta1_theta[pos]) +
                                                                 ((2.0)*delta2_theta[pos]) + ((8*nu)*delta3_theta[pos]) );

                    sx=0.0;sy=0.0;sz=0.0;
                    for(int orb=0;orb<Parameters_.orbs;orb++){
                        sx +=0.5*real( Red_Den_mat_temp[pos][orb][1][pos][orb][0] + Red_Den_mat_temp[pos][orb][0][pos][orb][1]
                                +
                                (1.0/(27*nu))*(
                                    ((10.0*nu) - 2.0)*(delta1_Red_Den_mat[pos][orb][1][pos][orb][0] + delta1_Red_Den_mat[pos][orb][0][pos][orb][1])
                                + (2.0)*(delta2_Red_Den_mat[pos][orb][1][pos][orb][0] + delta2_Red_Den_mat[pos][orb][0][pos][orb][1])
                                + (8.0*nu)*(delta3_Red_Den_mat[pos][orb][1][pos][orb][0] + delta3_Red_Den_mat[pos][orb][0][pos][orb][1])
                                )
                                );
                        sy +=0.5*imag( Red_Den_mat_temp[pos][orb][1][pos][orb][0] - Red_Den_mat_temp[pos][orb][0][pos][orb][1]
                                +
                                (1.0/(27*nu))*(
                                    ((10.0*nu) - 2.0)*(delta1_Red_Den_mat[pos][orb][1][pos][orb][0] - delta1_Red_Den_mat[pos][orb][0][pos][orb][1])
                                + (2.0)*(delta2_Red_Den_mat[pos][orb][1][pos][orb][0] - delta2_Red_Den_mat[pos][orb][0][pos][orb][1])
                                + (8.0*nu)*(delta3_Red_Den_mat[pos][orb][1][pos][orb][0] - delta3_Red_Den_mat[pos][orb][0][pos][orb][1])
                                )
                                );
                        sz +=0.5*real( Red_Den_mat_temp[pos][orb][0][pos][orb][0] - Red_Den_mat_temp[pos][orb][1][pos][orb][1]
                                +
                                (1.0/(27*nu))*(
                                    ((10.0*nu) - 2.0)*(delta1_Red_Den_mat[pos][orb][0][pos][orb][0] - delta1_Red_Den_mat[pos][orb][1][pos][orb][1])
                                + (2.0)*(delta2_Red_Den_mat[pos][orb][0][pos][orb][0] - delta2_Red_Den_mat[pos][orb][1][pos][orb][1])
                                + (8.0*nu)*(delta3_Red_Den_mat[pos][orb][0][pos][orb][0] - delta3_Red_Den_mat[pos][orb][1][pos][orb][1])
                                )
                                );
                    }


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


                if(step_no==4){
                    pi=Phi[ts][pos_x][pos_y] + (1.0/(392.0*nu))*( (-1.0*( (77.0*nu - 56.0) + ((17.0*nu - 8.0)*sqrt(21.0)*sign_))*(delta1_phi[pos]))
                                                                  + (-8.0*(7.0 + sqrt(21.0)*sign_)*delta2_phi[pos] )
                                                                  + (48.0*nu*(7.0 + sqrt(21.0)*sign_)*delta3_phi[pos] )
                                                                  + (-3.0*nu*(21.0 + sqrt(21.0)*sign_)*delta4_phi[pos] )
                                                                  );
                    ti=Theta[ts][pos_x][pos_y] + (1.0/(392.0*nu))*( (-1.0*( (77.0*nu - 56.0) + ((17.0*nu - 8.0)*sqrt(21.0)*sign_))*(delta1_theta[pos]))
                                                                    + (-8.0*(7.0 + sqrt(21.0)*sign_)*delta2_theta[pos] )
                                                                    + (48.0*nu*(7.0 + sqrt(21.0)*sign_)*delta3_theta[pos] )
                                                                    + (-3.0*nu*(21.0 + sqrt(21.0)*sign_)*delta4_theta[pos] )
                                                                    );

                    sx=0.0;sy=0.0;sz=0.0;
                    for(int orb=0;orb<Parameters_.orbs;orb++){
                        sx +=0.5*real( Red_Den_mat_temp[pos][orb][1][pos][orb][0] + Red_Den_mat_temp[pos][orb][0][pos][orb][1]
                                +
                                (1.0/(392.0*nu))*(
                                    ((-1.0*( (77.0*nu - 56.0) + ((17.0*nu - 8.0)*sqrt(21.0)*sign_)))*(delta1_Red_Den_mat[pos][orb][1][pos][orb][0] + delta1_Red_Den_mat[pos][orb][0][pos][orb][1]))
                                + ((-8.0*(7.0 + sqrt(21.0)*sign_))*(delta2_Red_Den_mat[pos][orb][1][pos][orb][0] + delta2_Red_Den_mat[pos][orb][0][pos][orb][1]))
                                + ((48.0*nu*(7.0 + sqrt(21.0)*sign_))*(delta3_Red_Den_mat[pos][orb][1][pos][orb][0] + delta3_Red_Den_mat[pos][orb][0][pos][orb][1]))
                                + ((-3.0*nu*(21.0 + sqrt(21.0)*sign_))*(delta4_Red_Den_mat[pos][orb][1][pos][orb][0] + delta4_Red_Den_mat[pos][orb][0][pos][orb][1]))
                                )
                                );
                        sy +=0.5*imag( Red_Den_mat_temp[pos][orb][1][pos][orb][0] - Red_Den_mat_temp[pos][orb][0][pos][orb][1]
                                +
                                (1.0/(392.0*nu))*(
                                    ((-1.0*( (77.0*nu - 56.0) + ((17.0*nu - 8.0)*sqrt(21.0)*sign_)))*(delta1_Red_Den_mat[pos][orb][1][pos][orb][0] - delta1_Red_Den_mat[pos][orb][0][pos][orb][1]))
                                + ((-8.0*(7.0 + sqrt(21.0)*sign_))*(delta2_Red_Den_mat[pos][orb][1][pos][orb][0] - delta2_Red_Den_mat[pos][orb][0][pos][orb][1]))
                                + ((48.0*nu*(7.0 + sqrt(21.0)*sign_))*(delta3_Red_Den_mat[pos][orb][1][pos][orb][0] - delta3_Red_Den_mat[pos][orb][0][pos][orb][1]))
                                + ((-3.0*nu*(21.0 + sqrt(21.0)*sign_))*(delta4_Red_Den_mat[pos][orb][1][pos][orb][0] - delta4_Red_Den_mat[pos][orb][0][pos][orb][1]))
                                )
                                );
                        sz +=0.5*real( Red_Den_mat_temp[pos][orb][0][pos][orb][0] - Red_Den_mat_temp[pos][orb][1][pos][orb][1]
                                +
                                (1.0/(392.0*nu))*(
                                    ((-1.0*( (77.0*nu - 56.0) + ((17.0*nu - 8.0)*sqrt(21.0)*sign_)))*(delta1_Red_Den_mat[pos][orb][0][pos][orb][0] - delta1_Red_Den_mat[pos][orb][1][pos][orb][1]))
                                + ((-8.0*(7.0 + sqrt(21.0)*sign_))*(delta2_Red_Den_mat[pos][orb][0][pos][orb][0] - delta2_Red_Den_mat[pos][orb][1][pos][orb][1]))
                                + ((48.0*nu*(7.0 + sqrt(21.0)*sign_))*(delta3_Red_Den_mat[pos][orb][0][pos][orb][0] - delta3_Red_Den_mat[pos][orb][1][pos][orb][1]))
                                + ((-3.0*nu*(21.0 + sqrt(21.0)*sign_))*(delta4_Red_Den_mat[pos][orb][0][pos][orb][0] - delta4_Red_Den_mat[pos][orb][1][pos][orb][1]))
                                )
                                );
                    }


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

                if(step_no==5){
                    pi=Phi[ts][pos_x][pos_y] + (1.0/(1960.0*nu))*( (-5.0*( (287.0*nu - 56.0) - ((59.0*nu - 8.0)*sqrt(21.0)*sign_))*(delta1_phi[pos]))
                                                                   + (-40.0*(7.0 - sqrt(21.0)*sign_)*delta2_phi[pos] )
                                                                   + (320.0*nu*(sqrt(21.0)*sign_)*delta3_phi[pos] )
                                                                   + (3.0*nu*(21.0 - 121.0*sqrt(21.0)*sign_)*delta4_phi[pos] )
                                                                   + (392.0*nu*(6.0 - sqrt(21.0)*sign_)*delta5_phi[pos] )
                                                                   );
                    ti=Theta[ts][pos_x][pos_y] + (1.0/(1960.0*nu))*( (-5.0*( (287.0*nu - 56.0) - ((59.0*nu - 8.0)*sqrt(21.0)*sign_))*(delta1_theta[pos]))
                                                                     + (-40.0*(7.0 - sqrt(21.0)*sign_)*delta2_theta[pos] )
                                                                     + (320.0*nu*(sqrt(21.0)*sign_)*delta3_theta[pos] )
                                                                     + (3.0*nu*(21.0 - 121.0*sqrt(21.0)*sign_)*delta4_theta[pos] )
                                                                     + (392.0*nu*(6.0 - sqrt(21.0)*sign_)*delta5_theta[pos] )
                                                                     );


                    sx=0.0;sy=0.0;sz=0.0;
                    for(int orb=0;orb<Parameters_.orbs;orb++){
                        sx +=0.5*real( Red_Den_mat_temp[pos][orb][1][pos][orb][0] + Red_Den_mat_temp[pos][orb][0][pos][orb][1]
                                +
                                (1.0/(1960.0*nu))*(
                                    ((-5.0*( (287.0*nu - 56.0) - ((59.0*nu - 8.0)*sqrt(21.0)*sign_)))*(delta1_Red_Den_mat[pos][orb][1][pos][orb][0] + delta1_Red_Den_mat[pos][orb][0][pos][orb][1]))
                                + ((-40.0*(7.0 - sqrt(21.0)*sign_))*(delta2_Red_Den_mat[pos][orb][1][pos][orb][0] + delta2_Red_Den_mat[pos][orb][0][pos][orb][1]))
                                + ((320.0*nu*(sqrt(21.0)*sign_))*(delta3_Red_Den_mat[pos][orb][1][pos][orb][0] + delta3_Red_Den_mat[pos][orb][0][pos][orb][1]))
                                + ((3.0*nu*(21.0 - 121.0*sqrt(21.0)*sign_))*(delta4_Red_Den_mat[pos][orb][1][pos][orb][0] + delta4_Red_Den_mat[pos][orb][0][pos][orb][1]))
                                + ((392.0*nu*(6.0 - sqrt(21.0)*sign_))*(delta5_Red_Den_mat[pos][orb][1][pos][orb][0] + delta5_Red_Den_mat[pos][orb][0][pos][orb][1]))
                                )
                                );
                        sy +=0.5*imag( Red_Den_mat_temp[pos][orb][1][pos][orb][0] - Red_Den_mat_temp[pos][orb][0][pos][orb][1]
                                +
                                (1.0/(1960.0*nu))*(
                                    ((-5.0*( (287.0*nu - 56.0) - ((59.0*nu - 8.0)*sqrt(21.0)*sign_)))*(delta1_Red_Den_mat[pos][orb][1][pos][orb][0] - delta1_Red_Den_mat[pos][orb][0][pos][orb][1]))
                                + ((-40.0*(7.0 - sqrt(21.0)*sign_))*(delta2_Red_Den_mat[pos][orb][1][pos][orb][0] - delta2_Red_Den_mat[pos][orb][0][pos][orb][1]))
                                + ((320.0*nu*(sqrt(21.0)*sign_))*(delta3_Red_Den_mat[pos][orb][1][pos][orb][0] - delta3_Red_Den_mat[pos][orb][0][pos][orb][1]))
                                + ((3.0*nu*(21.0 - 121.0*sqrt(21.0)*sign_))*(delta4_Red_Den_mat[pos][orb][1][pos][orb][0] - delta4_Red_Den_mat[pos][orb][0][pos][orb][1]))
                                + ((392.0*nu*(6.0 - sqrt(21.0)*sign_))*(delta5_Red_Den_mat[pos][orb][1][pos][orb][0] - delta5_Red_Den_mat[pos][orb][0][pos][orb][1]))
                                )
                                );
                        sz +=0.5*real( Red_Den_mat_temp[pos][orb][0][pos][orb][0] - Red_Den_mat_temp[pos][orb][1][pos][orb][1]
                                +
                                (1.0/(1960.0*nu))*(
                                    ((-5.0*( (287.0*nu - 56.0) - ((59.0*nu - 8.0)*sqrt(21.0)*sign_)))*(delta1_Red_Den_mat[pos][orb][0][pos][orb][0] - delta1_Red_Den_mat[pos][orb][1][pos][orb][1]))
                                + ((-40.0*(7.0 - sqrt(21.0)*sign_))*(delta2_Red_Den_mat[pos][orb][0][pos][orb][0] - delta2_Red_Den_mat[pos][orb][1][pos][orb][1]))
                                + ((320.0*nu*(sqrt(21.0)*sign_))*(delta3_Red_Den_mat[pos][orb][0][pos][orb][0] - delta3_Red_Den_mat[pos][orb][1][pos][orb][1]))
                                + ((3.0*nu*(21.0 - 121.0*sqrt(21.0)*sign_))*(delta4_Red_Den_mat[pos][orb][0][pos][orb][0] - delta4_Red_Den_mat[pos][orb][1][pos][orb][1]))
                                + ((392.0*nu*(6.0 - sqrt(21.0)*sign_))*(delta5_Red_Den_mat[pos][orb][0][pos][orb][0] - delta5_Red_Den_mat[pos][orb][1][pos][orb][1]))
                                )
                                );
                    }


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

                if(step_no==6){
                    pi=Phi[ts][pos_x][pos_y] + (1.0/(180.0*nu))*( (15.0*( (30.0*nu - 8.0) - ((7.0*nu)*sqrt(21.0)*sign_))*(delta1_phi[pos]))
                                                                  + ((120.0)*delta2_phi[pos] )
                                                                  + (-40.0*nu*(5.0 + (7.0*(sqrt(21.0)*sign_)) )*delta3_phi[pos] )
                                                                  + (63.0*nu*(2.0 + 3.0*sqrt(21.0)*sign_)*delta4_phi[pos] )
                                                                  + (-14.0*nu*(49.0 - (9.0*sqrt(21.0)*sign_) )*delta5_phi[pos] )
                                                                  + (70.0*nu*(7.0 + (sqrt(21.0)*sign_) )*delta6_phi[pos] )
                                                                  );
                    ti=Theta[ts][pos_x][pos_y] + (1.0/(180.0*nu))*( (15.0*( (30.0*nu - 8.0) - ((7.0*nu)*sqrt(21.0)*sign_))*(delta1_theta[pos]))
                                                                    + ((120.0)*delta2_theta[pos] )
                                                                    + (-40.0*nu*(5.0 + (7.0*(sqrt(21.0)*sign_)) )*delta3_theta[pos] )
                                                                    + (63.0*nu*(2.0 + 3.0*sqrt(21.0)*sign_)*delta4_theta[pos] )
                                                                    + (-14.0*nu*(49.0 - (9.0*sqrt(21.0)*sign_) )*delta5_theta[pos] )
                                                                    + (70.0*nu*(7.0 + (sqrt(21.0)*sign_) )*delta6_theta[pos] )
                                                                    );


                    sx=0.0;sy=0.0;sz=0.0;
                    for(int orb=0;orb<Parameters_.orbs;orb++){
                        sx +=0.5*real( Red_Den_mat_temp[pos][orb][1][pos][orb][0] + Red_Den_mat_temp[pos][orb][0][pos][orb][1]
                                +
                                (1.0/(180.0*nu))*(
                                    (15.0*( (30.0*nu - 8.0) - ((7.0*nu)*sqrt(21.0)*sign_))*(delta1_Red_Den_mat[pos][orb][1][pos][orb][0] + delta1_Red_Den_mat[pos][orb][0][pos][orb][1]))
                                + ((120.0)*(delta2_Red_Den_mat[pos][orb][1][pos][orb][0] + delta2_Red_Den_mat[pos][orb][0][pos][orb][1]))
                                + (-40.0*nu*(5.0 + (7.0*(sqrt(21.0)*sign_)) )*(delta3_Red_Den_mat[pos][orb][1][pos][orb][0] + delta3_Red_Den_mat[pos][orb][0][pos][orb][1]))
                                + (63.0*nu*(2.0 + 3.0*sqrt(21.0)*sign_)*(delta4_Red_Den_mat[pos][orb][1][pos][orb][0] + delta4_Red_Den_mat[pos][orb][0][pos][orb][1]))
                                + (-14.0*nu*(49.0 - (9.0*sqrt(21.0)*sign_) )*(delta5_Red_Den_mat[pos][orb][1][pos][orb][0] + delta5_Red_Den_mat[pos][orb][0][pos][orb][1]))
                                + (70.0*nu*(7.0 + (sqrt(21.0)*sign_) )*(delta6_Red_Den_mat[pos][orb][1][pos][orb][0] + delta6_Red_Den_mat[pos][orb][0][pos][orb][1]))
                                )
                                );
                        sy +=0.5*imag( Red_Den_mat_temp[pos][orb][1][pos][orb][0] - Red_Den_mat_temp[pos][orb][0][pos][orb][1]
                                +
                                (1.0/(180.0*nu))*(
                                    (15.0*( (30.0*nu - 8.0) - ((7.0*nu)*sqrt(21.0)*sign_))*(delta1_Red_Den_mat[pos][orb][1][pos][orb][0] - delta1_Red_Den_mat[pos][orb][0][pos][orb][1]))
                                + ((120.0)*(delta2_Red_Den_mat[pos][orb][1][pos][orb][0] - delta2_Red_Den_mat[pos][orb][0][pos][orb][1]))
                                + (-40.0*nu*(5.0 + (7.0*(sqrt(21.0)*sign_)) )*(delta3_Red_Den_mat[pos][orb][1][pos][orb][0] - delta3_Red_Den_mat[pos][orb][0][pos][orb][1]))
                                + (63.0*nu*(2.0 + 3.0*sqrt(21.0)*sign_)*(delta4_Red_Den_mat[pos][orb][1][pos][orb][0] - delta4_Red_Den_mat[pos][orb][0][pos][orb][1]))
                                + (-14.0*nu*(49.0 - (9.0*sqrt(21.0)*sign_) )*(delta5_Red_Den_mat[pos][orb][1][pos][orb][0] - delta5_Red_Den_mat[pos][orb][0][pos][orb][1]))
                                + (70.0*nu*(7.0 + (sqrt(21.0)*sign_) )*(delta6_Red_Den_mat[pos][orb][1][pos][orb][0] - delta6_Red_Den_mat[pos][orb][0][pos][orb][1]))
                                )
                                );
                        sz +=0.5*real( Red_Den_mat_temp[pos][orb][0][pos][orb][0] - Red_Den_mat_temp[pos][orb][1][pos][orb][1]
                                +
                                (1.0/(180.0*nu))*(
                                    (15.0*( (30.0*nu - 8.0) - ((7.0*nu)*sqrt(21.0)*sign_))*(delta1_Red_Den_mat[pos][orb][0][pos][orb][0] - delta1_Red_Den_mat[pos][orb][1][pos][orb][1]))
                                + ((120.0)*(delta2_Red_Den_mat[pos][orb][0][pos][orb][0] - delta2_Red_Den_mat[pos][orb][1][pos][orb][1]))
                                + (-40.0*nu*(5.0 + (7.0*(sqrt(21.0)*sign_)) )*(delta3_Red_Den_mat[pos][orb][0][pos][orb][0] - delta3_Red_Den_mat[pos][orb][1][pos][orb][1]))
                                + (63.0*nu*(2.0 + 3.0*sqrt(21.0)*sign_)*(delta4_Red_Den_mat[pos][orb][0][pos][orb][0] - delta4_Red_Den_mat[pos][orb][1][pos][orb][1]))
                                + (-14.0*nu*(49.0 - (9.0*sqrt(21.0)*sign_) )*(delta5_Red_Den_mat[pos][orb][0][pos][orb][0] - delta5_Red_Den_mat[pos][orb][1][pos][orb][1]))
                                + (70.0*nu*(7.0 + (sqrt(21.0)*sign_) )*(delta6_Red_Den_mat[pos][orb][0][pos][orb][0] - delta6_Red_Den_mat[pos][orb][1][pos][orb][1]))
                                )
                                );
                    }


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


                derivative_phi = (-1.0)*( Jval_array[pos]*((((sin(pi)*cos(ti))/(sin(ti)))*sy) +
                                                           (((cos(pi)*cos(ti))/(sin(ti)))*sx) -  sz ))  -  Bval_array[pos]; //-Bval_array[pos]
                derivative_theta = ( Jval_array[pos]*((cos(pi)*sy) - (sin(pi)*sx)));

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
                if(step_no==4){
                    delta5_phi[pos] = (dt_*(derivative_phi ));
                    delta5_theta[pos] = (dt_*( derivative_theta));

                }
                if(step_no==5){
                    delta6_phi[pos] = (dt_*(derivative_phi ));
                    delta6_theta[pos] = (dt_*( derivative_theta));

                }
                if(step_no==6){
                    delta7_phi[pos] = (dt_*(derivative_phi ));
                    delta7_theta[pos] = (dt_*( derivative_theta));

                }

                //Phi[ts+1][pos] = Phi[ts][pos] ;
                //Theta[ts+1][pos] = Theta[ts][pos] ;




                int pos_ng;
                double pj,tj, value_;
                for(int ng=0;ng<8;ng++){

                    if(ng<2){
                        value_=1.14*Parameters_.J_NN;  // +x, -x
                    }
                    else if(ng>=2 && ng <4){
                        value_=Parameters_.J_NN;   // +y, -y
                    }
                    else{
                        value_=Parameters_.J_NNN;  //+x+y, +x-y, -x+y, -x-y
                    }

                    pos_ng = Coordinates_.neigh(pos,ng);

                    int pos_ng_x = Coordinates_.indx(pos_ng);
                    int pos_ng_y = Coordinates_.indy(pos_ng);


                    if(step_no==0){
                        pj=Phi[ts][pos_ng_x][pos_ng_y];
                        tj=Theta[ts][pos_ng_x][pos_ng_y];
                    }


                    if(step_no==1){
                        pj=Phi[ts][pos_ng_x][pos_ng_y] + (nu)*(delta1_phi[pos_ng]);
                        tj=Theta[ts][pos_ng_x][pos_ng_y] + (nu)*(delta1_theta[pos_ng]);

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
                        pj=Phi[ts][pos_ng_x][pos_ng_y] + ( ((4*nu) - 1.0)*(delta1_phi[pos_ng]) + (1.0)*(delta2_phi[pos_ng]) )*(1.0/(8.0*nu));
                        tj=Theta[ts][pos_ng_x][pos_ng_y] + ( ((4*nu) - 1.0)*(delta1_theta[pos_ng]) + (1.0)*(delta2_theta[pos_ng]) )*(1.0/(8.0*nu));

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
                        pj=Phi[ts][pos_ng_x][pos_ng_y] + (1.0/(27*nu))*( (( (10.0*nu) - 2.0)*delta1_phi[pos_ng]) +
                                                                         ((2.0)*delta2_phi[pos_ng]) + ((8*nu)*delta3_phi[pos_ng]) );
                        tj=Theta[ts][pos_ng_x][pos_ng_y] + (1.0/(27*nu))*( (( (10.0*nu) - 2.0)*delta1_theta[pos_ng]) +
                                                                           ((2.0)*delta2_theta[pos_ng]) + ((8*nu)*delta3_theta[pos_ng]) );

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

                    if(step_no==4){
                        pj=Phi[ts][pos_ng_x][pos_ng_y] + (1.0/(392.0*nu))*( (-1.0*( (77.0*nu - 56.0) + ((17.0*nu - 8.0)*sqrt(21.0)*sign_))*(delta1_phi[pos_ng]))
                                                                            + (-8.0*(7.0 + sqrt(21.0)*sign_)*delta2_phi[pos_ng] )
                                                                            + (48.0*nu*(7.0 + sqrt(21.0)*sign_)*delta3_phi[pos_ng] )
                                                                            + (-3.0*nu*(21.0 + sqrt(21.0)*sign_)*delta4_phi[pos_ng] )
                                                                            );
                        tj=Theta[ts][pos_ng_x][pos_ng_y] + (1.0/(392.0*nu))*( (-1.0*( (77.0*nu - 56.0) + ((17.0*nu - 8.0)*sqrt(21.0)*sign_))*(delta1_theta[pos_ng]))
                                                                              + (-8.0*(7.0 + sqrt(21.0)*sign_)*delta2_theta[pos_ng] )
                                                                              + (48.0*nu*(7.0 + sqrt(21.0)*sign_)*delta3_theta[pos_ng] )
                                                                              + (-3.0*nu*(21.0 + sqrt(21.0)*sign_)*delta4_theta[pos_ng] )
                                                                              );

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

                    if(step_no==5){
                        pj=Phi[ts][pos_ng_x][pos_ng_y] + (1.0/(1960.0*nu))*( (-5.0*( (287.0*nu - 56.0) - ((59.0*nu - 8.0)*sqrt(21.0)*sign_))*(delta1_phi[pos_ng]))
                                                                             + (-40.0*(7.0 - sqrt(21.0)*sign_)*delta2_phi[pos_ng] )
                                                                             + (320.0*nu*(sqrt(21.0)*sign_)*delta3_phi[pos_ng] )
                                                                             + (3.0*nu*(21.0 - 121.0*sqrt(21.0)*sign_)*delta4_phi[pos_ng] )
                                                                             + (392.0*nu*(6.0 - sqrt(21.0)*sign_)*delta5_phi[pos_ng] )
                                                                             );
                        tj=Theta[ts][pos_ng_x][pos_ng_y] + (1.0/(1960.0*nu))*( (-5.0*( (287.0*nu - 56.0) - ((59.0*nu - 8.0)*sqrt(21.0)*sign_))*(delta1_theta[pos_ng]))
                                                                               + (-40.0*(7.0 - sqrt(21.0)*sign_)*delta2_theta[pos_ng] )
                                                                               + (320.0*nu*(sqrt(21.0)*sign_)*delta3_theta[pos_ng] )
                                                                               + (3.0*nu*(21.0 - 121.0*sqrt(21.0)*sign_)*delta4_theta[pos_ng] )
                                                                               + (392.0*nu*(6.0 - sqrt(21.0)*sign_)*delta5_theta[pos_ng] )
                                                                               );

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

                    if(step_no==6){
                        pj=Phi[ts][pos_ng_x][pos_ng_y] + (1.0/(180.0*nu))*( (15.0*( (30.0*nu - 8.0) - ((7.0*nu)*sqrt(21.0)*sign_))*(delta1_phi[pos_ng]))
                                                                            + ((120.0)*delta2_phi[pos_ng] )
                                                                            + (-40.0*nu*(5.0 + (7.0*(sqrt(21.0)*sign_)) )*delta3_phi[pos_ng] )
                                                                            + (63.0*nu*(2.0 + 3.0*sqrt(21.0)*sign_)*delta4_phi[pos_ng] )
                                                                            + (-14.0*nu*(49.0 - (9.0*sqrt(21.0)*sign_) )*delta5_phi[pos_ng] )
                                                                            + (70.0*nu*(7.0 + (sqrt(21.0)*sign_) )*delta6_phi[pos_ng] )
                                                                            );
                        tj=Theta[ts][pos_ng_x][pos_ng_y] + (1.0/(180.0*nu))*( (15.0*( (30.0*nu - 8.0) - ((7.0*nu)*sqrt(21.0)*sign_))*(delta1_theta[pos_ng]))
                                                                              + ((120.0)*delta2_theta[pos_ng] )
                                                                              + (-40.0*nu*(5.0 + (7.0*(sqrt(21.0)*sign_)) )*delta3_theta[pos_ng] )
                                                                              + (63.0*nu*(2.0 + 3.0*sqrt(21.0)*sign_)*delta4_theta[pos_ng] )
                                                                              + (-14.0*nu*(49.0 - (9.0*sqrt(21.0)*sign_) )*delta5_theta[pos_ng] )
                                                                              + (70.0*nu*(7.0 + (sqrt(21.0)*sign_) )*delta6_theta[pos_ng] )
                                                                              );

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

                    derivative_phi = (-1.0)*(( (S_mag)*(value_)*(((sin(tj)*cos(ti)*cos(pi - pj))/(sin(ti))) - cos(tj) ) ));
                    derivative_theta = (S_mag)*(value_)*(sin(tj)*sin(pj - pi));


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

                    if(step_no==4){
                        delta5_phi[pos] += (dt_*derivative_phi );
                        delta5_theta[pos] += (dt_*derivative_theta);
                    }

                    if(step_no==5){
                        delta6_phi[pos] += (dt_*derivative_phi );
                        delta6_theta[pos] += (dt_*derivative_theta);
                    }

                    if(step_no==6){
                        delta7_phi[pos] += (dt_*derivative_phi );
                        delta7_theta[pos] += (dt_*derivative_theta);
                    }

                }

            }

        }


#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
        for(int pos=0;pos<Parameters_.ns;pos++){

            int pos_x= Coordinates_.indx(pos);
            int pos_y= Coordinates_.indy(pos);


            Phi[ts+1][pos_x][pos_y] = Phi[ts][pos_x][pos_y]  + (1.0/180.0)*( ((9.0)*delta1_phi[pos]) + ((64.0)*delta3_phi[pos]) + ((49.0)*delta5_phi[pos])
                                                                             + ((49.0)*delta6_phi[pos]) + ((9.0)*delta7_phi[pos]) );
            Theta[ts+1][pos_x][pos_y]  = Theta[ts][pos_x][pos_y]  + (1.0/180.0)*( ((9.0)*delta1_theta[pos]) + ((64.0)*delta3_theta[pos]) + ((49.0)*delta5_theta[pos])
                                                                                  + ((49.0)*delta6_theta[pos]) + ((9.0)*delta7_theta[pos]) );



            //Phi[ts+1][pos] = Phi[ts][pos] + (1.0/1.0)*(0.0*delta1_phi[pos] + 1.0*delta2_phi[pos] + 0.0*delta3_phi[pos] + 0.0*delta4_phi[pos]);
            //Theta[ts+1][pos] = Theta[ts][pos] + (1.0/1.0)*(0.0*delta1_theta[pos] + 1.0*delta2_theta[pos] + 0.0*delta3_theta[pos] + 0.0*delta4_theta[pos]);


            if(Phi[ts+1][pos_x][pos_y]  > 2*PI){
                Phi[ts+1][pos_x][pos_y]  += -2*PI;

            }
            if(Phi[ts+1][pos_x][pos_y]  < 0){
                Phi[ts+1][pos_x][pos_y]  +=  2*PI;

            }


            if(Theta[ts+1][pos_x][pos_y]  > PI){
                Theta[ts+1][pos_x][pos_y]  = Theta[ts+1][pos_x][pos_y]  -2*PI;
                Phi[ts+1][pos_x][pos_y]  = fmod( Phi[ts+1][pos_x][pos_y]  + PI, 2.0*PI );
            }
            if(Theta[ts+1][pos_x][pos_y]  < 0){
                Theta[ts+1][pos_x][pos_y]  = - Theta[ts+1][pos_x][pos_y] ;
                Phi[ts+1][pos_x][pos_y]  = fmod( Phi[ts+1][pos_x][pos_y]  + PI, 2.0*PI );

            }

            //cout<<delta1_phi[pos]<<"   "<<delta2_phi[pos]<<"   "<<delta3_phi[pos]<<"   "<<delta4_phi[pos]<<"   "<<(1.0/6.0)*(delta1_phi[pos] + 2.0*delta2_phi[pos] + 2.0*delta3_phi[pos] + delta4_phi[pos])<<endl;
        }

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
        for(int pos_i=0;pos_i<Parameters_.ns;pos_i++){
            for (int orb_i=0;orb_i<Parameters_.orbs;orb_i++){
                for(int si=0;si<2;si++){
                    for(int pos_j=0;pos_j<Parameters_.ns;pos_j++){
                        for (int orb_j=0;orb_j<Parameters_.orbs;orb_j++){
                            for(int sj=0;sj<2;sj++){
                                Red_Den_mat[pos_i][orb_i][si][pos_j][orb_j][sj] = Red_Den_mat_temp[pos_i][orb_i][si][pos_j][orb_j][sj] + (1.0/180.0)*(9.0*delta1_Red_Den_mat[pos_i][orb_i][si][pos_j][orb_j][sj] + 64.0*delta3_Red_Den_mat[pos_i][orb_i][si][pos_j][orb_j][sj] +
                                                                                                                                                      49.0*delta5_Red_Den_mat[pos_i][orb_i][si][pos_j][orb_j][sj] + 49.0*delta6_Red_Den_mat[pos_i][orb_i][si][pos_j][orb_j][sj]
                                                                                                                                                      + 9.0*delta7_Red_Den_mat[pos_i][orb_i][si][pos_j][orb_j][sj]);
                            }
                        }
                    }
                }
            }
        }

    }

    */

    Red_Den_mat_temp.clear();
}





void SC_SW_ENGINE_VNE_1orb_MCMF::Write_final_time_result(){



    ofstream file_out_DM(Save_DM_file.c_str());

    for(int pi=0;pi<Red_Den_mat.size();pi++){
        for(int si=0;si<Red_Den_mat[pi].size();si++){

            for(int pj=0;pj<Red_Den_mat[pi][si].size();pj++){
                for(int sj=0;sj<Red_Den_mat[pi][si][pj].size();sj++){
                    file_out_DM<<pi<<"   "<<si<<"    "<<pj<<"   "<<sj<<"    "<<
                                 Red_Den_mat[pi][si][pj][sj].real()<<"  "<<
                                 Red_Den_mat[pi][si][pj][sj].imag()<<endl;

                }
            }
        }

    }

    ofstream file_out_Classical(Save_Classical_file.c_str());

    for(int px=0;px<Parameters_.lx;px++){
        for(int py=0;py<Parameters_.ly;py++){

            file_out_Classical<<px<<"    "<<py<<"   "<<Theta[0][px][py]<<"    "<<Phi[0][px][py]<<endl;
        }
    }


}



void SC_SW_ENGINE_VNE_1orb_MCMF::Read_Restart_Data(){

    complex<double> one(1,0);
    complex<double> iota(0,1);

    ifstream file_in_DM(Restart_DM_file.c_str());

    int temp;
    double temp_real, temp_imag;
    for(int pi=0;pi<Red_Den_mat.size();pi++){
        for(int si=0;si<2;si++){

            for(int pj=0;pj<Red_Den_mat.size();pj++){
                for(int sj=0;sj<2;sj++){
                    file_in_DM>>temp>>temp>>temp>>temp>>temp_real>>temp_imag;
                    Red_Den_mat[pi][si][pj][sj]=(one*temp_real) + (iota*temp_imag);


                }
            }
        }



    }

    ifstream file_in_Classical(Restart_Classical_file.c_str());

    int temp_lx, temp_ly;
    double temp_theta, temp_phi;
    for(int x=0;x<Parameters_.lx;x++){
        for(int y=0;y<Parameters_.lx;y++){

            file_in_Classical>>temp_lx>>temp_ly>>temp_theta>>temp_phi;
            Theta[0][temp_lx][temp_ly]=temp_theta;
            Phi[0][temp_lx][temp_ly]=temp_phi;

        }

    }





}



void SC_SW_ENGINE_VNE_1orb_MCMF::Read_parameters(string filename){


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


    string restart_, Restart_ = "Restart = ";
    string restart_time_, Restart_Time_ = "Restart_time = ";
    string restart_dm_file_, Restart_Dm_File_ = "Restart_DM_Configuration = ";
    string restart_classical_file_, Restart_Classical_File_ = "Restart_Classical_angles = ";


    string predictor_corrector_ ,Predictor_Corrector_ = "Predictor_Corrector = ";
    string use_fft_, Use_FFT_ = "Use_FFT = ";

    string save_final_time_, Save_Final_Time_ = "Save_Final_Time_Results = ";
    string Save_Dm_File_ = "Save_DM_Configuration = ";
    string Save_Classical_File_ = "Save_Classical_angles = ";


    string insitu_spacetimefourier_, Insitu_SpaceTimeFourier_ = "Insitu_SpaceTimeFourier = ";


    /*
        bool RESTART;
        double Restart_Time;
        string Restart_DM_file;
        string Restart_Classical_file;
             */

    int offset;
    string line;
    ifstream inputfile(filename.c_str());


    if(inputfile.is_open())
    {
        while(!inputfile.eof())
        {
            getline(inputfile,line);


            if ((offset = line.find(Restart_, 0)) != string::npos) {
                restart_ = line.substr (offset + Restart_.length());		}

            if ((offset = line.find(Insitu_SpaceTimeFourier_, 0)) != string::npos) {
                insitu_spacetimefourier_ = line.substr (offset + Insitu_SpaceTimeFourier_.length());		}

            if ((offset = line.find(Restart_Time_, 0)) != string::npos) {
                restart_time_ = line.substr (offset + Restart_Time_.length());		}

            if ((offset = line.find(Restart_Dm_File_, 0)) != string::npos) {
                restart_dm_file_ = line.substr (offset + Restart_Dm_File_.length());		}

            if ((offset = line.find(Restart_Classical_File_, 0)) != string::npos) {
                restart_classical_file_ = line.substr (offset + Restart_Classical_File_.length());		}

            if ((offset = line.find(Save_Final_Time_, 0)) != string::npos) {
                save_final_time_= line.substr (offset + Save_Final_Time_.length());		}


            if ((offset = line.find(Predictor_Corrector_, 0)) != string::npos) {
                predictor_corrector_= line.substr (offset + Predictor_Corrector_.length());		}


            if ((offset = line.find(Use_FFT_, 0)) != string::npos) {
                use_fft_= line.substr (offset + Use_FFT_.length());		}



            if ((offset = line.find(Save_Dm_File_, 0)) != string::npos) {
                Save_DM_file = line.substr (offset + Save_Dm_File_.length());		}

            if ((offset = line.find(Save_Classical_File_, 0)) != string::npos) {
                Save_Classical_file = line.substr (offset + Save_Classical_File_.length());		}



            if ((offset = line.find(W_conv_, 0)) != string::npos) {
                w_conv_ = line.substr (offset + W_conv_.length());		}



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
    {cout<<"Unable to open input file while in the model class."<<endl;}


    w_min=atof(w_min_.c_str());
    w_max=atof(w_max_.c_str());
    Restart_Time=atof(restart_time_.c_str());
    dw=atof(dw_.c_str());
    w_conv=atof(w_conv_.c_str());
    time_max=atof(time_max_.c_str());
    dt_ =atof(dt__.c_str());
    Runge_Kutta_order=atoi(runge_kutta_order_.c_str());


    no_of_processors=atoi(no_of_processors_.c_str());




    if (restart_ == "true"){
        RESTART =true;
    }
    else{
        RESTART=false;

    }

    if (insitu_spacetimefourier_ == "true" ){
        Insitu_SpaceTimeFourier=true;
    }
    else{
        Insitu_SpaceTimeFourier=false;
    }

    Restart_DM_file=restart_dm_file_;
    Restart_Classical_file=restart_classical_file_;

    if(save_final_time_=="true"){
        SAVE=true;
    }
    else{
        SAVE=false;
    }

    if(predictor_corrector_=="true"){
        Predictor_Corrector=true;
    }
    else{
        Predictor_Corrector=false;
    }



    if(use_fft_=="true"){
        Use_FFT=true;
    }
    else{
        Use_FFT=false;
    }



    RANDOM_NO_SEED=Parameters_.RandomSeed;

    if(RESTART){
        cout<<"This is a Restart from time = "<<Restart_Time<<endl;
        cout<<"Files "<< Restart_DM_file<<", "<<Restart_Classical_file<<" are used"<<endl;
    }
    {
        cout <<"This starts from equibrium configuration 0 time using file "<<conf_input<<endl;
    }

    cout<<"PARAMETERS::::::::"<<endl;

    cout<<"Lx = "<<Parameters_.lx<<endl;
    cout<<"Ly = "<<Parameters_.ly<<endl;
    cout<<"w_min = "<<w_min<<endl;
    cout<<"w_max = "<<w_max<<endl;
    cout<<"dw = "<<dw<<endl;
    cout<<"w_convolution = "<<w_conv<<endl;
    cout<<"time_max = "<<time_max<<endl;
    cout<<"dt = "<<dt_<<endl;
    cout<<"Runge_Kutta_order = "<<Runge_Kutta_order<<endl;
    cout<<"U_COUL = "<<Parameters_.U_COUL<<endl;
    cout<<"Temperature = "<<Parameters_.temp<<endl;
    cout<<"N_total = "<<Parameters_.Fill*Parameters_.ns*Parameters_.orbs*2<<endl;
    cout<<"RANDOM_NO_SEED = "<<Parameters_.RandomSeed<<endl;

    if(SAVE){
        cout<<"Final Results for time = "<<time_max<<endl;
        cout<<"will be written in "<< Save_DM_file<<", "<<Save_Classical_file<<endl;
    }








}
















































