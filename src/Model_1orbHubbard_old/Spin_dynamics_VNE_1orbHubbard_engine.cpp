#include "Spin_dynamics_VNE_1orbHubbard_engine.h"
#include <stdlib.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;



void SC_SW_ENGINE_VNE_1orbHubbard::Initialize_engine(MODEL_1_orb_SF &model, BASIS_1_orb_SF &basis){


    //******************To reproduce New J. Phys. 17 (2015) 113058********//
    model.Single_classical_S_in_B=false;
    //-------------------------------------------------------//


    //******************To do 1D with given initial condition for classical spins********//
    model.Decide_initial_cond_1D=true;
    //-------------------------------------------------------//

    RK4_type="type_1";   //Not used
    final_conf_out="final_config.txt"; //Not used

    //all above read from input file

    n_wpoints=(int) ((w_max-w_min)/dw + 0.5);

    if(!RESTART){
        Restart_Time=0.0;
    }


    time_steps=(int) ((time_max - Restart_Time)/dt_ + 0.5);



    model.Theta.resize(Runge_Kutta_order+1);
    model.Phi.resize(Runge_Kutta_order+1);


    model.Theta_eq.resize(basis.Max_pos + 1);
    model.Phi_eq.resize(basis.Max_pos + 1);

    for(int i=0;i<model.Theta.size();i++){

        model.Theta[i].resize(basis.Max_pos + 1);
        model.Phi[i].resize(basis.Max_pos + 1);
    }

    S_kw.resize(basis.Max_pos + 1 );
    for(int i=0;i<=basis.Max_pos;i++){
        S_kw[i].resize(n_wpoints);

    }



    model.quant_s_x.resize(1);
    model.quant_s_y.resize(1);
    model.quant_s_z.resize(1);

    model.Red_Den_mat.resize(basis.Max_pos + 1);
    model.quant_s_x_eq.resize(basis.Max_pos + 1);
    model.quant_s_y_eq.resize(basis.Max_pos + 1);
    model.quant_s_z_eq.resize(basis.Max_pos + 1);
    model.quant_s_x[0].resize(basis.Max_pos + 1);
    model.quant_s_y[0].resize(basis.Max_pos + 1);
    model.quant_s_z[0].resize(basis.Max_pos + 1);
    for(int i=0;i<=basis.Max_pos;i++){
        model.Red_Den_mat[i].resize(2);
        for(int si=0;si<2;si++){
            model.Red_Den_mat[i][si].resize(basis.Max_pos + 1);
            for(int j=0;j<=basis.Max_pos;j++){
                model.Red_Den_mat[i][si][j].resize(2);

            }

        }

    }





    model.center = int((basis.Max_pos + 1)*0.5 + 0.5) - 1;



    model.Jval_array.resize(basis.Max_pos + 1);
    model.Bval_array.resize(basis.Max_pos + 1);


    if(model.Single_classical_S_in_B==true){
        cout<<"only one classical spin is coupled with quantum bath"<<endl;
        for(int pos=0;pos<basis.Max_pos + 1;pos++){

            if(pos== model.center){
                model.Jval_array[pos]=model.Jval;
                model.Bval_array[pos]=model.B_mag;
            }
            else{
                model.Jval_array[pos]=0;
                model.Bval_array[pos]=0;
            }
        }
        cout<<"Magnetic Field at site = "<<model.center<<endl;
    }
    else{
        for(int pos=0;pos<basis.Max_pos + 1;pos++){
            model.Jval_array[pos]=model.Jval;
            model.Bval_array[pos]=0.0;
        }

    }






}



void SC_SW_ENGINE_VNE_1orbHubbard::Read_equilibrium_configuration(MODEL_1_orb_SF & model, BASIS_1_orb_SF &basis){


    ifstream file_in(conf_input.c_str());

    double temp_theta,temp_phi;
    int temp_pos;
    //Later on read from file
    for(int x=0;x<basis.Lx;x++){
        for(int y=0;y<basis.Ly;y++){

            //int pos= y*basis.Lx + x;

            file_in>>temp_pos>>temp_theta>>temp_phi;
            model.Theta[0][temp_pos-1]=temp_theta*PI;
            model.Phi[0][temp_pos-1]=temp_phi*PI;

            //cout<<temp_pos<<"  "<<temp_theta<<"   "<<temp_phi<<endl;

            model.Theta_eq[temp_pos-1]=temp_theta*PI;
            model.Phi_eq[temp_pos-1]=temp_phi*PI;

        }

    }


    double temp, noise;
    noise=0.5;
    if(model.Decide_initial_cond_1D==true){
        cout<<"Initial classical spin configuration is decided internally by code"<<endl;

        if(model.Single_classical_S_in_B==true){
            //Little tilted away from +z axis.
            model.Phi[0][model.center]=0;
            model.Theta[0][model.center]=PI*(1.0 - (1.0/50.0));
        }
        else{//pi order in classical spins with some noise
            srand(RANDOM_NO_SEED);

            for(int pos=0;pos<basis.Max_pos + 1;pos++){
                temp=(rand()%RAND_MAX);
                temp=temp/RAND_MAX;
                model.Phi[0][pos]=0;
                model.Phi_eq[pos]=0;
                if(pos%2==0){
                    model.Theta[0][pos]=noise*temp;
                    model.Theta_eq[pos]=noise*temp;
                }
                else{
                    model.Theta[0][pos]=(PI) - noise*temp;
                    model.Theta_eq[pos]=(PI) - noise*temp;
                }
            }
        }
    }
    else{
        cout<<"Inital classical spin configuration is read from given input file"<<endl;
    }





    string file_out_Sk_name = "Sk_static.txt";
    ofstream file_out_Sk(file_out_Sk_name.c_str());

    file_out_Sk<<"#nx    ny      S_classical(nx,ny)"<<endl;
    double kx,ky;
    int pos_i,pos_j;
    for(int nx=0;nx<basis.Lx;nx++){
        for(int ny=0;ny<basis.Ly;ny++){

            kx = (2*nx*PI)/(1.0*basis.Lx);
            ky = (2*ny*PI)/(1.0*basis.Ly);
            temp=0;
            for(int x_i=0;x_i<basis.Lx;x_i++){
                for(int y_i=0;y_i<basis.Ly;y_i++){
                    pos_i = y_i*basis.Lx + x_i;
                    for(int x_j=0;x_j<basis.Lx;x_j++){
                        for(int y_j=0;y_j<basis.Ly;y_j++){
                            pos_j = y_j*basis.Lx + x_j;

                            temp += model.S_mag*model.S_mag*(
                                        ( (cos(model.Theta_eq[pos_i])*  ( cos(model.Theta_eq[pos_j])   )  )  )
                                        +
                                        (  (sin(model.Theta_eq[pos_i])*cos(model.Phi_eq[pos_i])*   ( sin(model.Theta_eq[pos_j])*cos(model.Phi_eq[pos_j])     )    )   )
                                        +
                                        (  (sin(model.Theta_eq[pos_i])*sin(model.Phi_eq[pos_i])*    (sin(model.Theta_eq[pos_j])*sin(model.Phi_eq[pos_j])      )      )  )
                                        )*cos(( (x_j - x_i)*kx +  (y_j - y_i)*ky ) );


                        }
                    }
                }
            }
            file_out_Sk<<nx<<"    "<<ny<<"      "<<temp<<endl;


        }
    }





}



void SC_SW_ENGINE_VNE_1orbHubbard::Start_Engine(MODEL_1_orb_SF & model, BASIS_1_orb_SF &basis){
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

    S_rw.resize(basis.Max_pos +1);
    for(int pos_i=0;pos_i<basis.Max_pos +1;pos_i++){
        S_rw[pos_i].resize(basis.Max_pos +1);
        for(int pos_j=0;pos_j<basis.Max_pos +1;pos_j++){
            S_rw[pos_i][pos_j].resize(n_wpoints);
            for(int wi=0;wi<n_wpoints;wi++){
                S_rw[pos_i][pos_j][wi]=zero;
            }
        }
    }


    s_quantum_rw.resize(basis.Max_pos +1);
    for(int pos_i=0;pos_i<basis.Max_pos +1;pos_i++){
        s_quantum_rw[pos_i].resize(basis.Max_pos +1);
        for(int pos_j=0;pos_j<basis.Max_pos +1;pos_j++){
            s_quantum_rw[pos_i][pos_j].resize(n_wpoints);
            for(int wi=0;wi<n_wpoints;wi++){
                s_quantum_rw[pos_i][pos_j][wi]=zero;
            }
        }
    }


    T_rw.resize(basis.Max_pos +1);
    for(int pos_i=0;pos_i<basis.Max_pos +1;pos_i++){
        T_rw[pos_i].resize(basis.Max_pos +1);
        for(int pos_j=0;pos_j<basis.Max_pos +1;pos_j++){
            T_rw[pos_i][pos_j].resize(n_wpoints);
            for(int wi=0;wi<n_wpoints;wi++){
                T_rw[pos_i][pos_j][wi]=zero;
            }
        }
    }


    ofstream file_out(spins_r_t_out.c_str());
    file_out<<"# Sz------, Sx-----,Sy-----"<<endl;

    ofstream file_mu_out(mu_output.c_str());

    for(int ts=0;ts<=time_steps;ts++){


        if(SAVE && (ts==time_steps)){
            Write_final_time_result(model, basis);
        }


        //create Hamiltonian for time_step=ts
        if(ts==0){

            if(!RESTART){
                model.Create_Hamil(basis, 0);
                model.Diagonalize_Hamil();
                model.Mu_=model.Calculate_mu(model.Evals);
                file_mu_out<<ts<<"     "<<model.Mu_<<endl;
                //model.Get_quantum_Spins(0);
                model.Get_red_den_mat(0);
                cout<<"Hamiltonian is diagonalized at t=0 "<<endl;
            }
            else{
                cout<<"Restart is done"<<endl;
                Read_Restart_Data(model, basis);
            }

            for(int pos=0;pos<basis.Max_pos +1;pos++){
                model.quant_s_x_eq[pos] = 0.5*real( model.Red_Den_mat[pos][1][pos][0] + model.Red_Den_mat[pos][0][pos][1] );
                model.quant_s_y_eq[pos] = 0.5*imag( model.Red_Den_mat[pos][1][pos][0] - model.Red_Den_mat[pos][0][pos][1] );
                model.quant_s_z_eq[pos] = 0.5*real( model.Red_Den_mat[pos][0][pos][0] - model.Red_Den_mat[pos][1][pos][1] );
            }

            for(int pos=0;pos<basis.Max_pos +1;pos++){
                model.quant_s_x[0][pos] = 0.5*real( model.Red_Den_mat[pos][1][pos][0] + model.Red_Den_mat[pos][0][pos][1] );
                model.quant_s_y[0][pos] = 0.5*imag( model.Red_Den_mat[pos][1][pos][0] - model.Red_Den_mat[pos][0][pos][1] );
                model.quant_s_z[0][pos] = 0.5*real( model.Red_Den_mat[pos][0][pos][0] - model.Red_Den_mat[pos][1][pos][1] );
            }



            /*file_out<<(ts*dt_)<<"  ";
            for(int pos=0;pos<=basis.Max_pos;pos++){
                file_out<<model.S_mag*cos(model.Theta[0][pos])<<"  ";
            }
            for(int pos=0;pos<=basis.Max_pos;pos++){
                file_out<<model.S_mag*sin(model.Theta[0][pos])*cos(model.Phi[0][pos])<<"  ";
            }
            for(int pos=0;pos<=basis.Max_pos;pos++){
                file_out<<model.S_mag*sin(model.Theta[0][pos])*sin(model.Phi[0][pos])<<"  ";
            }

            for(int pos=0;pos<=basis.Max_pos;pos++){
                file_out<< model.quant_s_z_eq[pos]<<"  ";
            }
            for(int pos=0;pos<=basis.Max_pos;pos++){
                file_out<< model.quant_s_x_eq[pos]<<"  ";
            }
            for(int pos=0;pos<=basis.Max_pos;pos++){
                file_out<< model.quant_s_y_eq[pos]<<"  ";
            }

            file_out<<endl;*/




        }
        else{
        }

        file_out<<(ts*dt_) + Restart_Time<<"  ";
        for(int pos=0;pos<=basis.Max_pos;pos++){
            file_out<<model.S_mag*cos(model.Theta[0][pos])<<"  ";
        }
        for(int pos=0;pos<=basis.Max_pos;pos++){
            file_out<<model.S_mag*sin(model.Theta[0][pos])*cos(model.Phi[0][pos])<<"  ";
        }
        for(int pos=0;pos<=basis.Max_pos;pos++){
            file_out<<model.S_mag*sin(model.Theta[0][pos])*sin(model.Phi[0][pos])<<"  ";
        }

        for(int pos=0;pos<=basis.Max_pos;pos++){
            file_out<< model.quant_s_z[0][pos]<<"  ";
        }
        for(int pos=0;pos<=basis.Max_pos;pos++){
            file_out<< model.quant_s_x[0][pos]<<"  ";
        }
        for(int pos=0;pos<=basis.Max_pos;pos++){
            file_out<< model.quant_s_y[0][pos]<<"  ";
        }

        file_out<<endl;


        if(Insitu_SpaceTimeFourier){
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
            for(int wi=0;wi<n_wpoints;wi++){

                for(int pos_i=0;pos_i<basis.Max_pos +1;pos_i++){
                    for(int pos_j=0;pos_j<basis.Max_pos +1;pos_j++){



                        /*               S_rw[pos_i][pos_j][wi] += one*cos((wi * dw) * (ts* dt_))*exp(-0.5*(ts*dt_*w_conv*ts*dt_*w_conv))*dt_*model.S_mag*model.S_mag*(  (cos(model.Theta_eq[pos_i])*cos(model.Theta[0][pos_j]))
                            +
                            (sin(model.Theta_eq[pos_i])*cos(model.Phi_eq[pos_i])*sin(model.Theta[0][pos_j])*cos(model.Phi[0][pos_j]) )
                            +
                            (sin(model.Theta_eq[pos_i])*sin(model.Phi_eq[pos_i])*sin(model.Theta[0][pos_j])*sin(model.Phi[0][pos_j]) )
                            );
*/

                        S_rw[pos_i][pos_j][wi] += exp(iota * (-wi * dw) * (ts* dt_))*exp(-0.5*(ts*dt_*w_conv*ts*dt_*w_conv))*dt_*model.S_mag*model.S_mag*(
                                    ( (cos(model.Theta_eq[pos_i])*  (  0.0*cos(model.Theta[0][pos_j]) +  1.0*cos(model.Theta_eq[pos_j])   )  )  )
                                +
                                (  (sin(model.Theta_eq[pos_i])*cos(model.Phi_eq[pos_i])*   ( 0.0*sin(model.Theta[0][pos_j])*cos(model.Phi[0][pos_j])  +   1.0*sin(model.Theta_eq[pos_j])*cos(model.Phi_eq[pos_j])     )    )   )
                                +
                                (  (sin(model.Theta_eq[pos_i])*sin(model.Phi_eq[pos_i])*    (0.0*sin(model.Theta[0][pos_j])*sin(model.Phi[0][pos_j])   +   1.0*sin(model.Theta_eq[pos_j])*sin(model.Phi_eq[pos_j])      )      )  )
                                )





                                /*-
                                                                                                    0.5*exp(iota * (-wi * dw) * (ts* dt_))*exp(-0.5*(ts*dt_*w_conv*ts*dt_*w_conv))*dt_*model.S_mag*model.S_mag*(
                                                                                                        ( (cos(model.Theta_eq[pos_i])*cos(model.Theta_eq[pos_j]))  )
                                                                                                        +
                                                                                                        (  (sin(model.Theta_eq[pos_i])*cos(model.Phi_eq[pos_i])*sin(model.Theta_eq[pos_j])*cos(model.Phi_eq[pos_j]) )   )
                                                                                                        +
                                                                                                        (  (sin(model.Theta_eq[pos_i])*sin(model.Phi_eq[pos_i])*sin(model.Theta_eq[pos_j])*sin(model.Phi_eq[pos_j]) )  )
                                                                                                        )*/


                                ;


                        s_quantum_rw[pos_i][pos_j][wi] += exp(iota * (-wi * dw) * (ts* dt_))*exp(-0.5*(ts*dt_*w_conv*ts*dt_*w_conv))*dt_*

                                ( (model.quant_s_x_eq[pos_i]*     (model.quant_s_x[0][pos_j]    -  1.0*model.quant_s_x_eq[pos_j]) )
                                +
                                (model.quant_s_y_eq[pos_i]*    ( model.quant_s_y[0][pos_j]  -  1.0*model.quant_s_y_eq[pos_j]) )
                                +
                                (model.quant_s_z_eq[pos_i]*    (model.quant_s_z[0][pos_j] -  1.0*model.quant_s_z_eq[pos_j]  ) )


                                )

                                /*-

                                                                                                    0.5*exp(iota * (-wi * dw) * (ts* dt_))*exp(-0.5*(ts*dt_*w_conv*ts*dt_*w_conv))*dt_*

                                                                                                    ( (model.quant_s_x_eq[pos_i]*model.quant_s_x_eq[pos_j]  )
                                                                                                      +
                                                                                                      (model.quant_s_y_eq[pos_i]*model.quant_s_y_eq[pos_j]  )
                                                                                                      +
                                                                                                      (model.quant_s_z_eq[pos_i]*model.quant_s_z_eq[pos_j] )


                                                                                                      )*/

                                ;




                        T_rw[pos_i][pos_j][wi] +=exp(iota * (-wi * dw) * (ts* dt_))*exp(-0.5*(ts*dt_*w_conv*ts*dt_*w_conv))*dt_*  ( model.S_mag*(
                                                                                                                                        (
                                                                                                                                            (cos(model.Theta_eq[pos_i])*model.quant_s_z[0][pos_j])
                                                                                                                                        +
                                                                                                                                        (sin(model.Theta_eq[pos_i])*cos(model.Phi_eq[pos_i])*model.quant_s_x[0][pos_j])
                                                                                                                                    +
                                                                                                                                    (sin(model.Theta_eq[pos_i])*sin(model.Phi_eq[pos_i])*model.quant_s_y[0][pos_j])

                                )
                                +
                                (
                                    (model.quant_s_z_eq[pos_i]*cos(model.Theta[0][pos_j]))
                                +
                                (model.quant_s_x_eq[pos_i]*sin(model.Theta[0][pos_j])*cos(model.Phi[0][pos_j]) )
                                +
                                (model.quant_s_y_eq[pos_i]*sin(model.Theta[0][pos_j])*sin(model.Phi[0][pos_j]) )

                                )


                                )

                                +


                                (
                                    (model.quant_s_x_eq[pos_i]*model.quant_s_x[0][pos_j])
                                +
                                (model.quant_s_y_eq[pos_i]*model.quant_s_y[0][pos_j])
                                +
                                (model.quant_s_z_eq[pos_i]*model.quant_s_z[0][pos_j])


                                )

                                +

                                (
                                    model.S_mag*model.S_mag*(  (cos(model.Theta_eq[pos_i])*cos(model.Theta[0][pos_j]))
                                    +
                                    (sin(model.Theta_eq[pos_i])*cos(model.Phi_eq[pos_i])*sin(model.Theta[0][pos_j])*cos(model.Phi[0][pos_j]) )
                                +
                                (sin(model.Theta_eq[pos_i])*sin(model.Phi_eq[pos_i])*sin(model.Theta[0][pos_j])*sin(model.Phi[0][pos_j]) )
                                )

                                )



                                );



                        /*  S_rw[pos_i][pos_j][wi] += exp(iota * (wi * dw) * (ts* dt_))*dt_*(  cos(model.Theta[0][pos_i])*cos(model.Theta[0][pos_j])
                                +
                                sin(model.Theta[0][pos_i])*cos(model.Phi[0][pos_i])*sin(model.Theta[0][pos_j])*cos(model.Phi[0][pos_j])
                                +
                                sin(model.Theta[0][pos_i])*sin(model.Phi[0][pos_i])*sin(model.Theta[0][pos_j])*sin(model.Phi[0][pos_j])
                                );

                        */




                    }
                }
            }


        }


        Evolve_classical_spins(0, model, basis);

        for(int pos=0;pos<basis.Max_pos +1;pos++){
            model.quant_s_x[0][pos] = 0.5*real( model.Red_Den_mat[pos][1][pos][0] + model.Red_Den_mat[pos][0][pos][1] );
            model.quant_s_y[0][pos] = 0.5*imag( model.Red_Den_mat[pos][1][pos][0] - model.Red_Den_mat[pos][0][pos][1] );
            model.quant_s_z[0][pos] = 0.5*real( model.Red_Den_mat[pos][0][pos][0] - model.Red_Den_mat[pos][1][pos][1] );
        }

        model.Theta[0]=model.Theta[1];
        model.Phi[0]=model.Phi[1];


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
        int Lxby2 = (int) ( ((1.0*basis.Lx)/2.0) + 0.5);
        int Lyby2 = (int) ( ((1.0*basis.Ly)/2.0) + 0.5);
        int Lxby4 = (int) ( ((1.0*basis.Lx)/4.0) + 0.5);
        int Lyby4 = (int) ( ((1.0*basis.Ly)/4.0) + 0.5);

        k_ind=0;

        Mat_3_Complex_doub Skw_Mat;
        Skw_Mat.resize(basis.Lx);

        Mat_3_Complex_doub squant_kw_Mat;
        squant_kw_Mat.resize(basis.Lx);

        Mat_3_Complex_doub Tkw_Mat;
        Tkw_Mat.resize(basis.Lx);

        for(int i=0;i<basis.Lx;i++){
            Skw_Mat[i].resize(basis.Ly);
            squant_kw_Mat[i].resize(basis.Ly);
            Tkw_Mat[i].resize(basis.Ly);
            for(int j=0;j<basis.Ly;j++){
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

                    // kx = (2*nx*PI)/(1.0*basis.Lx);
                    // ky = (2*ny*PI)/(1.0*basis.Ly);


#ifdef _OPENMP
#pragma omp parallel for default(shared)  private(pos_i,pos_j,temp,temp2,temp3)
#endif

                    for(int wi=0;wi<n_wpoints;wi++){
                        temp3=zero;
                        temp2=zero;
                        temp=zero;
                        for(int x_i=0;x_i<basis.Lx;x_i++){
                            for(int y_i=0;y_i<basis.Ly;y_i++){
                                pos_i = y_i*basis.Lx + x_i;
                                for(int x_j=0;x_j<basis.Lx;x_j++){
                                    for(int y_j=0;y_j<basis.Ly;y_j++){
                                        pos_j = y_j*basis.Lx + x_j;

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
            for(nx=0;nx<basis.Lx;nx++){
                for(ny=0;ny<basis.Ly;ny++){

                    kx = (2*nx*PI)/(1.0*basis.Lx);
                    ky = (2*ny*PI)/(1.0*basis.Ly);

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(pos_i,pos_j,temp,temp2,temp3)
#endif

                    for(int wi=0;wi<n_wpoints;wi++){
                        temp3=zero;
                        temp2=zero;
                        temp=zero;
                        for(int x_i=0;x_i<basis.Lx;x_i++){
                            for(int y_i=0;y_i<basis.Ly;y_i++){
                                pos_i = y_i*basis.Lx + x_i;
                                for(int x_j=0;x_j<basis.Lx;x_j++){
                                    for(int y_j=0;y_j<basis.Ly;y_j++){
                                        pos_j = y_j*basis.Lx + x_j;

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
            for(nx=0;nx<basis.Lx;nx++){
                for(ny=0;ny<basis.Ly;ny++){

                    for(int wi=0;wi<n_wpoints;wi++){

                        //file_out2<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<temp.real()<<"   "<<temp.imag()<<"    "<<temp2.real()<<"   "<<temp2.imag()<<"    "<<temp3.real()<<"   "<<temp3.imag()<<endl;
                        file_out_full<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<wi<<"   "<<Skw_Mat[nx][ny][wi].real()<<"   "<<Skw_Mat[nx][ny][wi].imag()<<"    "<<squant_kw_Mat[nx][ny][wi].real()<<"   "<<squant_kw_Mat[nx][ny][wi].imag()<<
                                       "     "<<Tkw_Mat[nx][ny][wi].real()<<"   "<<Tkw_Mat[nx][ny][wi].imag()<<endl;
                        //Skw_Mat[nx][ny][wi]=temp;
                    }

                    k_ind +=1;
                    file_out_full<<endl;

                }

            }



        }




        //Only in case of square Lattice
        if(basis.Lx==basis.Ly){


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

                    // kx = (2*nx*PI)/(1.0*basis.Lx);
                    // ky = (2*ny*PI)/(1.0*basis.Ly);
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(pos_i,pos_j,temp,temp2,temp3)
#endif
                    for(int wi=0;wi<n_wpoints;wi++){
                        temp3=zero;
                        temp2=zero;
                        temp=zero;
                        for(int x_i=0;x_i<basis.Lx;x_i++){
                            for(int y_i=0;y_i<basis.Ly;y_i++){
                                pos_i = y_i*basis.Lx + x_i;
                                for(int x_j=0;x_j<basis.Lx;x_j++){
                                    for(int y_j=0;y_j<basis.Ly;y_j++){
                                        pos_j = y_j*basis.Lx + x_j;

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

                    //kx = (2*nx*PI)/(1.0*basis.Lx);
                    //ky = (2*ny*PI)/(1.0*basis.Ly);
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(pos_i,pos_j,temp,temp2,temp3)
#endif
                    for(int wi=0;wi<n_wpoints;wi++){
                        temp3=zero;
                        temp2=zero;
                        temp=zero;
                        for(int x_i=0;x_i<basis.Lx;x_i++){
                            for(int y_i=0;y_i<basis.Ly;y_i++){
                                pos_i = y_i*basis.Lx + x_i;
                                for(int x_j=0;x_j<basis.Lx;x_j++){
                                    for(int y_j=0;y_j<basis.Ly;y_j++){
                                        pos_j = y_j*basis.Lx + x_j;

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
                        for(int x_i=0;x_i<basis.Lx;x_i++){
                            for(int y_i=0;y_i<basis.Ly;y_i++){
                                pos_i = y_i*basis.Lx + x_i;
                                for(int x_j=0;x_j<basis.Lx;x_j++){
                                    for(int y_j=0;y_j<basis.Ly;y_j++){
                                        pos_j = y_j*basis.Lx + x_j;

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
                        for(int x_i=0;x_i<basis.Lx;x_i++){
                            for(int y_i=0;y_i<basis.Ly;y_i++){
                                pos_i = y_i*basis.Lx + x_i;
                                for(int x_j=0;x_j<basis.Lx;x_j++){
                                    for(int y_j=0;y_j<basis.Ly;y_j++){
                                        pos_j = y_j*basis.Lx + x_j;

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



        }
        /*
    for(int nx=0;nx<basis.Lx;nx++){
        for(int ny=0;ny<basis.Ly;ny++){
               kx = (2*nx*PI)/(1.0*basis.Lx);
               ky = (2*ny*PI)/(1.0*basis.Ly);

            for(int wi=0;wi<n_wpoints;wi++){
                temp=zero;
                for(int x_i=0;x_i<basis.Lx;x_i++){
                    for(int y_i=0;y_i<basis.Ly;y_i++){
                        pos_i = y_i*basis.Lx + x_i;
                        for(int x_j=0;x_j<basis.Lx;x_j++){
                            for(int y_j=0;y_j<basis.Ly;y_j++){
                                pos_j = y_j*basis.Lx + x_j;

                                temp += S_rw[pos_i][pos_j][wi]*exp(iota*( (x_j - x_i)*kx +  (y_j - y_i)*ky ) );


                            }
                        }
                    }
                }
                file_out2<<nx<<"   "<<ny<<"   "<<nx*basis.Lx + ny<<"   "<<wi*dw<<"   "<<temp.real()<<"   "<<temp.imag()<<endl;
            }
            file_out2<<endl;
        }
    }
*/


    }

}




void SC_SW_ENGINE_VNE_1orbHubbard::Evolve_classical_spins(int ts, MODEL_1_orb_SF & model, BASIS_1_orb_SF &basis){

    complex<double> zero(0,0);
    complex<double> one(1,0);
    complex<double> iota(0,1);


    Mat_4_Complex_doub Red_Den_mat_temp;
    Red_Den_mat_temp=model.Red_Den_mat;


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

        double pi,ti; //phi_i,theta_i,phi_j,theta_j for timestep=ts
        double sy,sx,sz; //Quantum spins for timestep=ts and position=pos


        complex<double> derivative_val, temp_val;

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(pi, ti, sx, sy, sz)
#endif
        for(int pos=0;pos<=basis.Max_pos;pos++){


            pi=model.Phi[ts][pos];
            ti=model.Theta[ts][pos];

            sx=0.5*real( Red_Den_mat_temp[pos][1][pos][0] + Red_Den_mat_temp[pos][0][pos][1] );
            sy=0.5*imag( Red_Den_mat_temp[pos][1][pos][0] - Red_Den_mat_temp[pos][0][pos][1] );
            sz=0.5*real( Red_Den_mat_temp[pos][0][pos][0] - Red_Den_mat_temp[pos][1][pos][1] );


            model.Phi[ts+1][pos] = model.Phi[ts][pos] + (dt_*(  -1.0*(model.Jval_array[pos]*((((sin(pi)*cos(ti))/(sin(ti)))*sy) + (((cos(pi)*cos(ti))/(sin(ti)))*sx) - sz))  -  model.Bval_array[pos]  )    );

            model.Theta[ts+1][pos] = model.Theta[ts][pos] + (dt_*( model.Jval_array[pos]*((cos(pi)*sy) - (sin(pi)*sx))));

            //model.Phi[ts+1][pos] = model.Phi[ts][pos] ;
            //model.Theta[ts+1][pos] = model.Theta[ts][pos] ;

            //derivative_phi += (-1.0)*( model.Jval*((((sin(pi)*cos(ti))/(sin(ti)))*sy) + (((cos(pi)*cos(ti))/(sin(ti)))*sx) - sz));
            //derivative_theta += ( model.Jval*((cos(pi)*sy) - (sin(pi)*sx)));

            int pos_ng;
            double pj,tj;
            for(int ng=0;ng<basis.Neighbour[pos].size();ng++){

                pos_ng = basis.Neighbour[pos][ng];

                pj=model.Phi[ts][pos_ng];
                tj=model.Theta[ts][pos_ng];


                model.Phi[ts+1][pos] += (-1.0)*(dt_*( (model.S_mag)*(model.Jval_p)* (((sin(tj)*cos(ti)*cos(pi - pj))/(sin(ti))) - cos(tj) )));
                //derivative_phi += (-1.0)*(( (model.S_mag)*(model.Jval_p)*(((sin(tj)*cos(ti)*cos(pi - pj))/(sin(ti))) - cos(tj) ) ));

                model.Theta[ts+1][pos] += (dt_*(model.S_mag)*(model.Jval_p)*(sin(tj)*sin(pj - pi)));
                //derivative_theta += (model.S_mag)*(model.Jval_p)*(sin(tj)*sin(pj - pi));


            }

            if(model.Phi[ts+1][pos] > 2*PI){
                model.Phi[ts+1][pos] += -2*PI;

            }
            if(model.Phi[ts+1][pos] < 0){
                model.Phi[ts+1][pos] +=  2*PI;

            }


            if(model.Theta[ts+1][pos] > PI){
                model.Theta[ts+1][pos] = model.Theta[ts+1][pos] -2*PI;
                model.Phi[ts+1][pos] = fmod( model.Phi[ts+1][pos] + PI, 2.0*PI );
            }
            if(model.Theta[ts+1][pos] < 0){
                model.Theta[ts+1][pos] = - model.Theta[ts+1][pos];
                model.Phi[ts+1][pos] = fmod( model.Phi[ts+1][pos] + PI, 2.0*PI );

            }


            // cout<<pos<<"  "<<derivative_theta<<"  "<<derivative_phi<<endl;

        }


#ifdef _OPENMP
#pragma omp parallel for default(shared) private(derivative_val,temp_val)
#endif
        for(int pos_i=0;pos_i<=basis.Max_pos;pos_i++){
            for(int si=0;si<2;si++){
                for(int pos_j=0;pos_j<=basis.Max_pos;pos_j++){
                    for(int sj=0;sj<2;sj++){


                        int pos_ng;
                        for(int ng=0;ng<basis.Neighbour[pos_i].size();ng++){
                            pos_ng=basis.Neighbour[pos_i][ng];

                            derivative_val=iota*model.t_hop*Red_Den_mat_temp[pos_ng][si][pos_j][sj];
                            model.Red_Den_mat[pos_i][si][pos_j][sj] += dt_*derivative_val;

                        }

                        for(int ng=0;ng<basis.Neighbour[pos_j].size();ng++){
                            pos_ng=basis.Neighbour[pos_j][ng];

                            derivative_val=iota*(-1.0)*model.t_hop*Red_Den_mat_temp[pos_i][si][pos_ng][sj];
                            model.Red_Den_mat[pos_i][si][pos_j][sj] += dt_*derivative_val;

                        }

                        for(int s3=0;s3<2;s3++){
                            double SX_i,SY_i,SZ_i,SX_j,SY_j,SZ_j;
                            SX_i = sin(model.Theta[ts][pos_i])*cos(model.Phi[ts][pos_i]);
                            SY_i = sin(model.Theta[ts][pos_i])*sin(model.Phi[ts][pos_i]);
                            SZ_i = cos(model.Theta[ts][pos_i]);
                            SX_j = sin(model.Theta[ts][pos_j])*cos(model.Phi[ts][pos_j]);
                            SY_j = sin(model.Theta[ts][pos_j])*sin(model.Phi[ts][pos_j]);
                            SZ_j = cos(model.Theta[ts][pos_j]);
                            temp_val=Red_Den_mat_temp[pos_i][s3][pos_j][sj];
                            derivative_val=iota*(-0.5*model.Jval_array[pos_i])*(model.S_mag)*(
                                        (SX_i)*Pauli_x[si][s3]
                                        +
                                        (SY_i)*Pauli_y[si][s3]
                                        +
                                        (SZ_i)*Pauli_z[si][s3]
                                        )*temp_val;
                            temp_val=Red_Den_mat_temp[pos_i][si][pos_j][s3];
                            derivative_val +=iota*(0.5*model.Jval_array[pos_j])*(model.S_mag)*(
                                        (SX_j)*Pauli_x[s3][sj]
                                        +
                                        (SY_j)*Pauli_y[s3][sj]
                                        +
                                        (SZ_j)*Pauli_z[s3][sj]
                                        )*temp_val;

                            model.Red_Den_mat[pos_i][si][pos_j][sj] += dt_*derivative_val;

                        }

                    }
                }
            }
        }


        //cout<<"Only RK4 is working for VNE"<<endl;
        //assert (Runge_Kutta_order==4);
    }

    if(Runge_Kutta_order==4){
        /*
         we let the quantum spins change as well in the intermediate steps of dt interval in RK-4 method,
         here we not diagonalize the matrix again calculate quantum spins for the changed classical spins.
         */

        double pi,ti; //phi_i,theta_i,phi_j,theta_j for timestep=ts
        double sy,sx,sz; //Quantum spins for timestep=ts and position=pos

        double derivative_theta, derivative_phi;
        Mat_1_doub delta1_theta, delta2_theta, delta3_theta, delta4_theta;
        Mat_1_doub delta1_phi, delta2_phi, delta3_phi, delta4_phi;
        Mat_4_Complex_doub delta1_Red_Den_mat, delta2_Red_Den_mat, delta3_Red_Den_mat, delta4_Red_Den_mat;


        //Standard RK-4 method Convention is used
        delta1_theta.resize(basis.Max_pos + 1);delta2_theta.resize(basis.Max_pos + 1);
        delta3_theta.resize(basis.Max_pos + 1);delta4_theta.resize(basis.Max_pos + 1);

        delta1_phi.resize(basis.Max_pos + 1);delta2_phi.resize(basis.Max_pos + 1);
        delta3_phi.resize(basis.Max_pos + 1);delta4_phi.resize(basis.Max_pos + 1);


        delta1_Red_Den_mat.resize(basis.Max_pos + 1);delta2_Red_Den_mat.resize(basis.Max_pos + 1);
        delta3_Red_Den_mat.resize(basis.Max_pos + 1); delta4_Red_Den_mat.resize(basis.Max_pos + 1);

        for(int i=0;i<=basis.Max_pos;i++){
            delta1_Red_Den_mat[i].resize(2);delta2_Red_Den_mat[i].resize(2);
            delta3_Red_Den_mat[i].resize(2);delta4_Red_Den_mat[i].resize(2);
            for(int si=0;si<2;si++){
                delta1_Red_Den_mat[i][si].resize(basis.Max_pos + 1);delta2_Red_Den_mat[i][si].resize(basis.Max_pos + 1);
                delta3_Red_Den_mat[i][si].resize(basis.Max_pos + 1);delta4_Red_Den_mat[i][si].resize(basis.Max_pos + 1);
                for(int j=0;j<=basis.Max_pos;j++){
                    delta1_Red_Den_mat[i][si][j].resize(2);delta2_Red_Den_mat[i][si][j].resize(2);
                    delta3_Red_Den_mat[i][si][j].resize(2);delta4_Red_Den_mat[i][si][j].resize(2);
                }
            }
        }


        for(int pos_i=0;pos_i<=basis.Max_pos;pos_i++){
            delta1_theta[pos_i]=0.0;delta2_theta[pos_i]=0.0;
            delta3_theta[pos_i]=0.0;delta4_theta[pos_i]=0.0;
            delta1_phi[pos_i]=0.0;delta2_phi[pos_i]=0.0;
            delta3_phi[pos_i]=0.0;delta4_phi[pos_i]=0.0;
            for(int si=0;si<2;si++){
                for(int pos_j=0;pos_j<=basis.Max_pos;pos_j++){
                    for(int sj=0;sj<2;sj++){
                        delta1_Red_Den_mat[pos_i][si][pos_j][sj]=zero;delta2_Red_Den_mat[pos_i][si][pos_j][sj]=zero;
                        delta3_Red_Den_mat[pos_i][si][pos_j][sj]=zero;delta4_Red_Den_mat[pos_i][si][pos_j][sj]=zero;
                    }
                }

            }

        }



        //calculating delta1(2,3,4)_thetha(phi)
        bool intermediate_update = true;
        double mu;
        for(int step_no=0;step_no<4;step_no++){


#ifdef _OPENMP
#pragma omp parallel for default(shared) private(derivative_theta,derivative_phi, pi, ti, sx, sy, sz)
#endif

            for(int pos=0;pos<=basis.Max_pos;pos++){

                derivative_theta=0;
                derivative_phi=0;

                if(step_no==0){
                    pi=model.Phi[ts][pos];
                    ti=model.Theta[ts][pos];
                    sx=0.5*real( Red_Den_mat_temp[pos][1][pos][0] + Red_Den_mat_temp[pos][0][pos][1] );
                    sy=0.5*imag( Red_Den_mat_temp[pos][1][pos][0] - Red_Den_mat_temp[pos][0][pos][1] );
                    sz=0.5*real( Red_Den_mat_temp[pos][0][pos][0] - Red_Den_mat_temp[pos][1][pos][1] );
                }
                if(step_no==1){
                    pi=model.Phi[ts][pos] + 0.5*(delta1_phi[pos]);
                    ti=model.Theta[ts][pos] + 0.5*(delta1_theta[pos]);
                    sx=0.5*real( Red_Den_mat_temp[pos][1][pos][0] + Red_Den_mat_temp[pos][0][pos][1]
                            + 0.5*(delta1_Red_Den_mat[pos][1][pos][0] + delta1_Red_Den_mat[pos][0][pos][1]));
                    sy=0.5*imag( Red_Den_mat_temp[pos][1][pos][0] - Red_Den_mat_temp[pos][0][pos][1]
                            +0.5*(delta1_Red_Den_mat[pos][1][pos][0] - delta1_Red_Den_mat[pos][0][pos][1])
                            );
                    sz=0.5*real( Red_Den_mat_temp[pos][0][pos][0] - Red_Den_mat_temp[pos][1][pos][1]
                            +0.5*(delta1_Red_Den_mat[pos][0][pos][0] - delta1_Red_Den_mat[pos][1][pos][1])
                            );

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
                    pi=model.Phi[ts][pos] + 0.5*(delta2_phi[pos]);
                    ti=model.Theta[ts][pos] + 0.5*(delta2_theta[pos]);

                    sx=0.5*real( Red_Den_mat_temp[pos][1][pos][0] + Red_Den_mat_temp[pos][0][pos][1]
                            + 0.5*(delta2_Red_Den_mat[pos][1][pos][0] + delta2_Red_Den_mat[pos][0][pos][1]));
                    sy=0.5*imag( Red_Den_mat_temp[pos][1][pos][0] - Red_Den_mat_temp[pos][0][pos][1]
                            +0.5*(delta2_Red_Den_mat[pos][1][pos][0] - delta2_Red_Den_mat[pos][0][pos][1])
                            );
                    sz=0.5*real( Red_Den_mat_temp[pos][0][pos][0] - Red_Den_mat_temp[pos][1][pos][1]
                            +0.5*(delta2_Red_Den_mat[pos][0][pos][0] - delta2_Red_Den_mat[pos][1][pos][1])
                            );

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
                    pi=model.Phi[ts][pos] + (delta3_phi[pos]);
                    ti=model.Theta[ts][pos] + (delta3_theta[pos]);

                    sx=0.5*real( Red_Den_mat_temp[pos][1][pos][0] + Red_Den_mat_temp[pos][0][pos][1]
                            + (delta3_Red_Den_mat[pos][1][pos][0] + delta3_Red_Den_mat[pos][0][pos][1]));
                    sy=0.5*imag( Red_Den_mat_temp[pos][1][pos][0] - Red_Den_mat_temp[pos][0][pos][1]
                            + (delta3_Red_Den_mat[pos][1][pos][0] - delta3_Red_Den_mat[pos][0][pos][1])
                            );
                    sz=0.5*real( Red_Den_mat_temp[pos][0][pos][0] - Red_Den_mat_temp[pos][1][pos][1]
                            + (delta3_Red_Den_mat[pos][0][pos][0] - delta3_Red_Den_mat[pos][1][pos][1])
                            );


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


                derivative_phi = (-1.0)*( model.Jval_array[pos]*((((sin(pi)*cos(ti))/(sin(ti)))*sy) + (((cos(pi)*cos(ti))/(sin(ti)))*sx) -  sz ))  -  model.Bval_array[pos]; //-model.Bval_array[pos]
                derivative_theta = ( model.Jval_array[pos]*((cos(pi)*sy) - (sin(pi)*sx)));

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

                //model.Phi[ts+1][pos] = model.Phi[ts][pos] ;
                //model.Theta[ts+1][pos] = model.Theta[ts][pos] ;




                int pos_ng;
                double pj,tj;
                for(int ng=0;ng<basis.Neighbour[pos].size();ng++){

                    pos_ng = basis.Neighbour[pos][ng];


                    if(step_no==0){
                        pj=model.Phi[ts][pos_ng];
                        tj=model.Theta[ts][pos_ng];
                    }

                    if(step_no==1){
                        pj=model.Phi[ts][pos_ng] + 0.5*(delta1_phi[pos_ng]);
                        tj=model.Theta[ts][pos_ng] + 0.5*(delta1_theta[pos_ng]);

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
                        pj=model.Phi[ts][pos_ng] + 0.5*(delta2_phi[pos_ng]);
                        tj=model.Theta[ts][pos_ng] + 0.5*(delta2_theta[pos_ng]);

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
                        pj=model.Phi[ts][pos_ng] + (delta3_phi[pos_ng]);
                        tj=model.Theta[ts][pos_ng] + (delta3_theta[pos_ng]);

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

                    derivative_phi = (-1.0)*(( (model.S_mag)*(model.Jval_p)*(((sin(tj)*cos(ti)*cos(pi - pj))/(sin(ti))) - cos(tj) ) ));
                    derivative_theta = (model.S_mag)*(model.Jval_p)*(sin(tj)*sin(pj - pi));


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

            }



            complex<double> temp_val;
            complex<double> derivative_val;
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(derivative_val, temp_val)
#endif
            for(int pos_i=0;pos_i<=basis.Max_pos;pos_i++){
                for(int si=0;si<2;si++){
                    for(int pos_j=0;pos_j<=basis.Max_pos;pos_j++){
                        for(int sj=0;sj<2;sj++){

                            int pos_ng;
                            for(int ng=0;ng<basis.Neighbour[pos_i].size();ng++){
                                pos_ng=basis.Neighbour[pos_i][ng];
                                if(step_no==0){
                                    temp_val=Red_Den_mat_temp[pos_ng][si][pos_j][sj];
                                    derivative_val=iota*model.t_hop*temp_val;
                                    delta1_Red_Den_mat[pos_i][si][pos_j][sj] += dt_*derivative_val;
                                }

                                if(step_no==1){
                                    temp_val=Red_Den_mat_temp[pos_ng][si][pos_j][sj]
                                            + 0.5*delta1_Red_Den_mat[pos_ng][si][pos_j][sj];
                                    derivative_val=iota*model.t_hop*temp_val;
                                    delta2_Red_Den_mat[pos_i][si][pos_j][sj] += dt_*derivative_val;
                                }
                                if(step_no==2){
                                    temp_val=Red_Den_mat_temp[pos_ng][si][pos_j][sj]
                                            + 0.5*delta2_Red_Den_mat[pos_ng][si][pos_j][sj];
                                    derivative_val=iota*model.t_hop*temp_val;
                                    delta3_Red_Den_mat[pos_i][si][pos_j][sj] += dt_*derivative_val;
                                }
                                if(step_no==3){
                                    temp_val=Red_Den_mat_temp[pos_ng][si][pos_j][sj]
                                            + delta3_Red_Den_mat[pos_ng][si][pos_j][sj];
                                    derivative_val=iota*model.t_hop*temp_val;
                                    delta4_Red_Den_mat[pos_i][si][pos_j][sj] += dt_*derivative_val;
                                }

                            }

                            for(int ng=0;ng<basis.Neighbour[pos_j].size();ng++){
                                pos_ng=basis.Neighbour[pos_j][ng];
                                if(step_no==0){
                                    temp_val=Red_Den_mat_temp[pos_i][si][pos_ng][sj];
                                    derivative_val=iota*(-1.0*model.t_hop)*temp_val;
                                    delta1_Red_Den_mat[pos_i][si][pos_j][sj] += dt_*derivative_val;
                                }
                                if(step_no==1){
                                    temp_val=Red_Den_mat_temp[pos_i][si][pos_ng][sj]
                                            +0.5*delta1_Red_Den_mat[pos_i][si][pos_ng][sj];
                                    derivative_val=iota*(-1.0*model.t_hop)*temp_val;
                                    delta2_Red_Den_mat[pos_i][si][pos_j][sj] += dt_*derivative_val;
                                }
                                if(step_no==2){
                                    temp_val=Red_Den_mat_temp[pos_i][si][pos_ng][sj]
                                            +0.5*delta2_Red_Den_mat[pos_i][si][pos_ng][sj];
                                    derivative_val=iota*(-1.0*model.t_hop)*temp_val;
                                    delta3_Red_Den_mat[pos_i][si][pos_j][sj] += dt_*derivative_val;
                                }
                                if(step_no==3){
                                    temp_val=Red_Den_mat_temp[pos_i][si][pos_ng][sj]
                                            + delta3_Red_Den_mat[pos_i][si][pos_ng][sj];
                                    derivative_val=iota*(-1.0*model.t_hop)*temp_val;
                                    delta4_Red_Den_mat[pos_i][si][pos_j][sj] += dt_*derivative_val;
                                }

                            }


                            for(int s3=0;s3<2;s3++){
                                double SX_i,SY_i,SZ_i,SX_j,SY_j,SZ_j;
                                if(step_no==0){
                                    SX_i = sin(model.Theta[ts][pos_i])*cos(model.Phi[ts][pos_i]);
                                    SY_i = sin(model.Theta[ts][pos_i])*sin(model.Phi[ts][pos_i]);
                                    SZ_i = cos(model.Theta[ts][pos_i]);
                                    SX_j = sin(model.Theta[ts][pos_j])*cos(model.Phi[ts][pos_j]);
                                    SY_j = sin(model.Theta[ts][pos_j])*sin(model.Phi[ts][pos_j]);
                                    SZ_j = cos(model.Theta[ts][pos_j]);
                                    temp_val=Red_Den_mat_temp[pos_i][s3][pos_j][sj];
                                    derivative_val=iota*(-0.5*model.Jval_array[pos_i])*(model.S_mag)*(
                                                (SX_i)*Pauli_x[si][s3]
                                                +
                                                (SY_i)*Pauli_y[si][s3]
                                                +
                                                (SZ_i)*Pauli_z[si][s3]
                                                )*temp_val;
                                    temp_val=Red_Den_mat_temp[pos_i][si][pos_j][s3];
                                    derivative_val +=iota*(0.5*model.Jval_array[pos_j])*(model.S_mag)*(
                                                (SX_j)*Pauli_x[s3][sj]
                                                +
                                                (SY_j)*Pauli_y[s3][sj]
                                                +
                                                (SZ_j)*Pauli_z[s3][sj]
                                                )*temp_val;

                                    delta1_Red_Den_mat[pos_i][si][pos_j][sj] += dt_*derivative_val;

                                }
                                if(step_no==1){
                                    SX_i = sin(model.Theta[ts][pos_i] + 0.5*delta1_theta[pos_i])*
                                            cos(model.Phi[ts][pos_i] + 0.5*delta1_phi[pos_i]);
                                    SY_i = sin(model.Theta[ts][pos_i] + 0.5*delta1_theta[pos_i])*
                                            sin(model.Phi[ts][pos_i] + 0.5*delta1_phi[pos_i]);
                                    SZ_i = cos(model.Theta[ts][pos_i] + 0.5*delta1_theta[pos_i]);

                                    SX_j = sin(model.Theta[ts][pos_j] + 0.5*delta1_theta[pos_j])*
                                            cos(model.Phi[ts][pos_j] + 0.5*delta1_phi[pos_j]);
                                    SY_j = sin(model.Theta[ts][pos_j] + 0.5*delta1_theta[pos_j])*
                                            sin(model.Phi[ts][pos_j] + 0.5*delta1_phi[pos_j]);
                                    SZ_j = cos(model.Theta[ts][pos_j] + 0.5*delta1_theta[pos_j]);

                                    temp_val=Red_Den_mat_temp[pos_i][s3][pos_j][sj] +
                                            0.5*delta1_Red_Den_mat[pos_i][s3][pos_j][sj];

                                    derivative_val=iota*(-0.5*model.Jval_array[pos_i])*(model.S_mag)*(
                                                (SX_i)*Pauli_x[si][s3]
                                                +
                                                (SY_i)*Pauli_y[si][s3]
                                                +
                                                (SZ_i)*Pauli_z[si][s3]
                                                )*temp_val;

                                    temp_val=Red_Den_mat_temp[pos_i][si][pos_j][s3] +
                                            0.5*delta1_Red_Den_mat[pos_i][si][pos_j][s3];

                                    derivative_val +=iota*(0.5*model.Jval_array[pos_j])*(model.S_mag)*(
                                                (SX_j)*Pauli_x[s3][sj]
                                                +
                                                (SY_j)*Pauli_y[s3][sj]
                                                +
                                                (SZ_j)*Pauli_z[s3][sj]
                                                )*temp_val;

                                    delta2_Red_Den_mat[pos_i][si][pos_j][sj] += dt_*derivative_val;

                                }
                                if(step_no==2){
                                    SX_i = sin(model.Theta[ts][pos_i] + 0.5*delta2_theta[pos_i])*
                                            cos(model.Phi[ts][pos_i] + 0.5*delta2_phi[pos_i]);
                                    SY_i = sin(model.Theta[ts][pos_i] + 0.5*delta2_theta[pos_i])*
                                            sin(model.Phi[ts][pos_i] + 0.5*delta2_phi[pos_i]);
                                    SZ_i = cos(model.Theta[ts][pos_i] + 0.5*delta2_theta[pos_i]);

                                    SX_j = sin(model.Theta[ts][pos_j] + 0.5*delta2_theta[pos_j])*
                                            cos(model.Phi[ts][pos_j] + 0.5*delta2_phi[pos_j]);
                                    SY_j = sin(model.Theta[ts][pos_j] + 0.5*delta2_theta[pos_j])*
                                            sin(model.Phi[ts][pos_j] + 0.5*delta2_phi[pos_j]);
                                    SZ_j = cos(model.Theta[ts][pos_j] + 0.5*delta2_theta[pos_j]);

                                    temp_val=Red_Den_mat_temp[pos_i][s3][pos_j][sj] +
                                            0.5*delta2_Red_Den_mat[pos_i][s3][pos_j][sj];

                                    derivative_val=iota*(-0.5*model.Jval_array[pos_i])*(model.S_mag)*(
                                                (SX_i)*Pauli_x[si][s3]
                                                +
                                                (SY_i)*Pauli_y[si][s3]
                                                +
                                                (SZ_i)*Pauli_z[si][s3]
                                                )*temp_val;

                                    temp_val=Red_Den_mat_temp[pos_i][si][pos_j][s3] +
                                            0.5*delta2_Red_Den_mat[pos_i][si][pos_j][s3];

                                    derivative_val +=iota*(0.5*model.Jval_array[pos_j])*(model.S_mag)*(
                                                (SX_j)*Pauli_x[s3][sj]
                                                +
                                                (SY_j)*Pauli_y[s3][sj]
                                                +
                                                (SZ_j)*Pauli_z[s3][sj]
                                                )*temp_val;

                                    delta3_Red_Den_mat[pos_i][si][pos_j][sj] += dt_*derivative_val;

                                }
                                if(step_no==3){
                                    SX_i = sin(model.Theta[ts][pos_i] + 0.5*delta3_theta[pos_i])*
                                            cos(model.Phi[ts][pos_i] + 0.5*delta3_phi[pos_i]);
                                    SY_i = sin(model.Theta[ts][pos_i] + 0.5*delta3_theta[pos_i])*
                                            sin(model.Phi[ts][pos_i] + 0.5*delta3_phi[pos_i]);
                                    SZ_i = cos(model.Theta[ts][pos_i] + 0.5*delta3_theta[pos_i]);

                                    SX_j = sin(model.Theta[ts][pos_j] + 0.5*delta3_theta[pos_j])*
                                            cos(model.Phi[ts][pos_j] + 0.5*delta3_phi[pos_j]);
                                    SY_j = sin(model.Theta[ts][pos_j] + 0.5*delta3_theta[pos_j])*
                                            sin(model.Phi[ts][pos_j] + 0.5*delta3_phi[pos_j]);
                                    SZ_j = cos(model.Theta[ts][pos_j] + 0.5*delta3_theta[pos_j]);

                                    temp_val=Red_Den_mat_temp[pos_i][s3][pos_j][sj] +
                                            0.5*delta3_Red_Den_mat[pos_i][s3][pos_j][sj];

                                    derivative_val=iota*(-0.5*model.Jval_array[pos_i])*(model.S_mag)*(
                                                (SX_i)*Pauli_x[si][s3]
                                                +
                                                (SY_i)*Pauli_y[si][s3]
                                                +
                                                (SZ_i)*Pauli_z[si][s3]
                                                )*temp_val;

                                    temp_val=Red_Den_mat_temp[pos_i][si][pos_j][s3] +
                                            0.5*delta3_Red_Den_mat[pos_i][si][pos_j][s3];

                                    derivative_val +=iota*(0.5*model.Jval_array[pos_j])*(model.S_mag)*(
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


        }


#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
        for(int pos=0;pos<=basis.Max_pos;pos++){


            model.Phi[ts+1][pos] = model.Phi[ts][pos] + (1.0/6.0)*(delta1_phi[pos] + 2.0*delta2_phi[pos] + 2.0*delta3_phi[pos] + delta4_phi[pos]);
            model.Theta[ts+1][pos] = model.Theta[ts][pos] + (1.0/6.0)*(delta1_theta[pos] + 2.0*delta2_theta[pos] + 2.0*delta3_theta[pos] + delta4_theta[pos]);


            /*
            model.Phi[ts+1][pos] = model.Phi[ts][pos] + (1.0/1.0)*(0.0*delta1_phi[pos] + 1.0*delta2_phi[pos] + 0.0*delta3_phi[pos] + 0.0*delta4_phi[pos]);
            model.Theta[ts+1][pos] = model.Theta[ts][pos] + (1.0/1.0)*(0.0*delta1_theta[pos] + 1.0*delta2_theta[pos] + 0.0*delta3_theta[pos] + 0.0*delta4_theta[pos]);
            */

            if(model.Phi[ts+1][pos] > 2*PI){
                model.Phi[ts+1][pos] += -2*PI;

            }
            if(model.Phi[ts+1][pos] < 0){
                model.Phi[ts+1][pos] +=  2*PI;

            }


            if(model.Theta[ts+1][pos] > PI){
                model.Theta[ts+1][pos] = model.Theta[ts+1][pos] -2*PI;
                model.Phi[ts+1][pos] = fmod( model.Phi[ts+1][pos] + PI, 2.0*PI );
            }
            if(model.Theta[ts+1][pos] < 0){
                model.Theta[ts+1][pos] = - model.Theta[ts+1][pos];
                model.Phi[ts+1][pos] = fmod( model.Phi[ts+1][pos] + PI, 2.0*PI );

            }

            //cout<<delta1_phi[pos]<<"   "<<delta2_phi[pos]<<"   "<<delta3_phi[pos]<<"   "<<delta4_phi[pos]<<"   "<<(1.0/6.0)*(delta1_phi[pos] + 2.0*delta2_phi[pos] + 2.0*delta3_phi[pos] + delta4_phi[pos])<<endl;
        }

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
        for(int pos_i=0;pos_i<=basis.Max_pos;pos_i++){
            for(int si=0;si<2;si++){
                for(int pos_j=0;pos_j<=basis.Max_pos;pos_j++){
                    for(int sj=0;sj<2;sj++){
                        model.Red_Den_mat[pos_i][si][pos_j][sj] = Red_Den_mat_temp[pos_i][si][pos_j][sj] + (1.0/6.0)*(delta1_Red_Den_mat[pos_i][si][pos_j][sj] + 2.0*delta2_Red_Den_mat[pos_i][si][pos_j][sj] +
                                                                                                                      2.0*delta3_Red_Den_mat[pos_i][si][pos_j][sj] + delta4_Red_Den_mat[pos_i][si][pos_j][sj]);
                    }
                }
            }
        }



    }

    Red_Den_mat_temp.clear();
}





void SC_SW_ENGINE_VNE_1orbHubbard::Write_final_time_result(MODEL_1_orb_SF & model, BASIS_1_orb_SF &basis){



    ofstream file_out_DM(Save_DM_file.c_str());

    for(int pi=0;pi<model.Red_Den_mat.size();pi++){
        for(int si=0;si<model.Red_Den_mat[pi].size();si++){

            for(int pj=0;pj<model.Red_Den_mat[pi][si].size();pj++){
                for(int sj=0;sj<model.Red_Den_mat[pi][si][pj].size();sj++){
                    file_out_DM<<pi<<"    "<<si<<"    "<<pj<<"    "<<sj<<"    "<<
                                 model.Red_Den_mat[pi][si][pj][sj].real()<<"  "<<model.Red_Den_mat[pi][si][pj][sj].imag()<<endl;

                }
            }

        }
    }

    ofstream file_out_Classical(Save_Classical_file.c_str());

    for(int pi=0;pi<model.Theta[0].size();pi++){

        file_out_Classical<<pi<<"    "<<model.Theta[0][pi]<<"    "<<model.Phi[0][pi]<<endl;
    }


}



void SC_SW_ENGINE_VNE_1orbHubbard::Read_Restart_Data(MODEL_1_orb_SF & model, BASIS_1_orb_SF &basis){

    complex<double> one(1,0);
    complex<double> iota(0,1);

    ifstream file_in_DM(Restart_DM_file.c_str());

    int temp;
    double temp_real, temp_imag;
    for(int pi=0;pi<model.Red_Den_mat.size();pi++){
        for(int si=0;si<model.Red_Den_mat[pi].size();si++){

            for(int pj=0;pj<model.Red_Den_mat[pi][si].size();pj++){
                for(int sj=0;sj<model.Red_Den_mat[pi][si][pj].size();sj++){
                    file_in_DM>>temp>>temp>>temp>>temp>>temp_real>>temp_imag;
                    model.Red_Den_mat[pi][si][pj][sj]=(one*temp_real) + (iota*temp_imag);


                }
            }

        }
    }

    ifstream file_in_Classical(Restart_Classical_file.c_str());

    for(int pi=0;pi<model.Theta[0].size();pi++){
        file_in_Classical>> temp >>model.Theta[0][pi]>>model.Phi[0][pi];
    }





}



void SC_SW_ENGINE_VNE_1orbHubbard::Read_parameters(string filename, MODEL_1_orb_SF & model, BASIS_1_orb_SF &basis){




    string temperature_, Temperature_ = "Temperature = ";
    string s_mag_, S_mag_ = "S_mag = ";
    string t_hop_, T_hop_ = "t_hop = ";
    string jvalh_, JvalH_ = "J_HUND = ";
    string jval_nn_, Jval_NN_ = "J_NN = ";
    string jval_nnn_, Jval_NNN_ = "J_NNN = ";
    string n_total_, N_total_ = "N_total = ";


    string Geometry_ = "Geometry = ";
    string boundary_conditions_, Boundary_Conditions_ = "Boundary conditions = ";

    string lx_, Lx_ = "Lx = ";
    string ly_, Ly_ = "Ly = ";

    string random_no_seed_, RANDOM_NO_SEED_  = "RANDOM_NO_SEED = ";



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



            if ((offset = line.find(Save_Dm_File_, 0)) != string::npos) {
                Save_DM_file = line.substr (offset + Save_Dm_File_.length());		}

            if ((offset = line.find(Save_Classical_File_, 0)) != string::npos) {
                Save_Classical_file = line.substr (offset + Save_Classical_File_.length());		}


            if ((offset = line.find(N_total_, 0)) != string::npos) {
                n_total_ = line.substr (offset + N_total_.length());		}



            if ((offset = line.find(RANDOM_NO_SEED_, 0)) != string::npos) {
                random_no_seed_ = line.substr (offset + RANDOM_NO_SEED_.length());		}

            if ((offset = line.find(Jval_NN_, 0)) != string::npos) {
                jval_nn_ = line.substr (offset + Jval_NN_.length());		}

            if ((offset = line.find(Jval_NNN_, 0)) != string::npos) {
                jval_nnn_ = line.substr (offset + Jval_NNN_.length());		}

            if ((offset = line.find(W_conv_, 0)) != string::npos) {
                w_conv_ = line.substr (offset + W_conv_.length());		}

            if ((offset = line.find(JvalH_, 0)) != string::npos) {
                jvalh_ = line.substr (offset + JvalH_.length());		}

            if ((offset = line.find(T_hop_, 0)) != string::npos) {
                t_hop_ = line.substr (offset + T_hop_.length());		}

            if ((offset = line.find(S_mag_, 0)) != string::npos) {
                s_mag_ = line.substr (offset + S_mag_.length());		}

            if ((offset = line.find(Temperature_, 0)) != string::npos) {
                temperature_  = line.substr (offset + Temperature_.length());		}

            if ((offset = line.find(Boundary_Conditions_, 0)) != string::npos) {
                boundary_conditions_ = line.substr (offset + Boundary_Conditions_.length());		}

            if ((offset = line.find(Geometry_, 0)) != string::npos) {
                basis.Geometry = line.substr (offset + Geometry_.length());		}

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




    if (boundary_conditions_ == "PBC"){
        basis.PBC =true;
    }
    else{
       basis.PBC = false;
    }



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

    basis.Lx=atoi(lx_.c_str());
    basis.Ly=atoi(ly_.c_str());


    RANDOM_NO_SEED=atoi(random_no_seed_.c_str());


    model.Jval_p=atof(jval_nn_.c_str());
    model.Jval=atof(jvalh_.c_str());
    model.t_hop=atof(t_hop_.c_str());
    model.S_mag=atof(s_mag_.c_str());
    model.Beta=double(11604.0/ (atof(temperature_.c_str())));
    model.N_total=atoi(n_total_.c_str());

    model.B_mag=1.0;




    if(RESTART){
        cout<<"This is a Restart from time = "<<Restart_Time<<endl;
        cout<<"Files "<< Restart_DM_file<<", "<<Restart_Classical_file<<" are used"<<endl;
    }

    cout<<"PARAMETERS::::::::"<<endl;
    cout<<"Boundary condition = "<<boundary_conditions_<<endl;
    cout<<"Lx = "<<basis.Lx<<endl;
    cout<<"Ly = "<<basis.Ly<<endl;
    cout<<"w_min = "<<w_min<<endl;
    cout<<"w_max = "<<w_max<<endl;
    cout<<"dw = "<<dw<<endl;
    cout<<"w_convolution = "<<w_conv<<endl;
    cout<<"time_max = "<<time_max<<endl;
    cout<<"dt = "<<dt_<<endl;
    cout<<"Runge_Kutta_order = "<<Runge_Kutta_order<<endl;
    cout<<"Jval_p = "<<model.Jval_p<<endl;
    cout<<"Jval = "<<model.Jval<<endl;
    cout<<"t_hop = "<<model.t_hop<<endl;
    cout<<"S_mag = "<<model.S_mag<<endl;
    cout<<"Beta = "<<model.Beta<<endl;
    cout<<"N_total = "<<model.N_total<<endl;
    cout<<"RANDOM_NO_SEED = "<<RANDOM_NO_SEED<<endl;

    if(SAVE){
        cout<<"Final Results for time = "<<time_max<<endl;
        cout<<"will be written in "<< Save_DM_file<<", "<<Save_Classical_file<<endl;
    }








}

































