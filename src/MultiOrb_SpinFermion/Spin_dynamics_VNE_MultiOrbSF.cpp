#include <stdlib.h>
#include "Spin_dynamics_VNE_MultiOrbSF.h"
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;

void SC_SW_ENGINE_VNE_MultiOrbSF::Initialize_engine(){


    SelfConsistentEvolution = false;

    RK4_type="type_1";   //Not used
    final_conf_out="final_config.txt"; //Not used

    //all above read from input file

    n_wpoints=(int) ((w_max-w_min)/dw + 0.5);

    if(!RESTART){
        Restart_Time=0.0;
    }


    time_steps=(int) ((time_max - Restart_Time)/(fabs(dt_)) + 0.5);



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
        Red_Den_mat[i].resize(Parameters_.n_orbs*2);
        for(int s=0;s<Parameters_.n_orbs*2;s++){
            Red_Den_mat[i][s].resize(Parameters_.ns);
            for(int j=0;j<Parameters_.ns;j++){
                Red_Den_mat[i][s][j].resize(Parameters_.n_orbs*2);
            }
        }
    }


    for (int ts=0;ts<Red_Den_mat_time.size();ts++){
        Red_Den_mat_time[ts].resize(Parameters_.ns);
        for(int i=0;i<Parameters_.ns;i++){
            Red_Den_mat_time[ts][i].resize(Parameters_.n_orbs*2);
            for(int s=0;s<Parameters_.n_orbs*2;s++){
                Red_Den_mat_time[ts][i][s].resize(Parameters_.ns);
                for(int j=0;j<Parameters_.ns;j++){
                    Red_Den_mat_time[ts][i][s][j].resize(Parameters_.n_orbs*2);
                }
            }
        }
    }


    center = int((Parameters_.ns)*0.5 + 0.5) - 1;

    Jval_array.resize(Parameters_.ns*Parameters_.n_orbs);
    Bval_array.resize(Parameters_.ns);

    int pos_x,pos_y;
    for(int pos=0;pos<Parameters_.ns;pos++){
        pos_x = Coordinates_.indx_cellwise(pos);
        pos_y = Coordinates_.indy_cellwise(pos);
        Bval_array[pos]=0.0;
        for(int orb=0;orb<Parameters_.n_orbs;orb++){
            Jval_array[Coordinates_.Nbasis(pos_x,pos_y,orb)] = Parameters_.J_Hund[orb];
        }
    }

    J_Classical = Parameters_.K1x;





    Pauli_x.resize(2);Pauli_y.resize(2);Pauli_z.resize(2);

    for(int i=0;i<2;i++){
        Pauli_x[i].resize(2); Pauli_y[i].resize(2);Pauli_z[i].resize(2);
        for(int j=0;j<2;j++){
            Pauli_x[i][j]=0;Pauli_y[i][j]=0;Pauli_z[i][j]=0;
        }
    }

    Pauli_x[0][1]=1.0;Pauli_x[1][0]=1.0;
    Pauli_y[0][1]=-1.0*iota_complex;Pauli_y[1][0]=1.0*iota_complex;
    Pauli_z[0][0]=1.0;Pauli_z[1][1]=-1.0;


}



void SC_SW_ENGINE_VNE_MultiOrbSF::Read_equilibrium_configuration(){


    ifstream file_in;
    file_in.open(conf_input.c_str());

    double temp_theta,temp_phi, temp_moment_size;
    int temp_lx, temp_ly;
    string temp_string;
    getline(file_in, temp_string);


    for(int x=0;x<Parameters_.lx;x++){
        for(int y=0;y<Parameters_.ly;y++){

            file_in>>temp_lx;
            file_in>>temp_ly;
            file_in>>temp_theta;
            file_in>>temp_phi;
            file_in>>temp_moment_size;
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



void SC_SW_ENGINE_VNE_MultiOrbSF::Start_Engine(){
#ifdef _OPENMP
    omp_set_num_threads(no_of_processors);
    int N_p = omp_get_max_threads();
    cout<<"Max threads which can be used parallely = "<<N_p<<endl;
    cout<<"No. of threads used parallely = "<<no_of_processors<<endl;
#endif

    string mu_output = "mu.txt";
    double initial_mu_guess;
    int n_states_occupied_zeroT;
    int px_,py_;


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

                MFParams_.Read_classical_DOFs(conf_input);
                Hamiltonian_.InteractionsCreate();
                //char flag='V';
                Hamiltonian_.Diagonalize('V');

                n_states_occupied_zeroT = Coordinates_.nbasis_*(Parameters_.Fill/(Parameters_.n_orbs));
                if(!Parameters_.fix_mu){
                    initial_mu_guess = 0.5 * (Hamiltonian_.eigs_[n_states_occupied_zeroT - 1] + Hamiltonian_.eigs_[n_states_occupied_zeroT]);
                }
                else{
                    initial_mu_guess=Parameters_.fixed_mu_value;
                }
                double mu_ = Hamiltonian_.chemicalpotential(initial_mu_guess,Parameters_.Fill/(Parameters_.n_orbs*2.0));
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

                px_ = Coordinates_.indx_cellwise(pos);
                py_ = Coordinates_.indy_cellwise(pos);

                Aux_S_x[0][pos] = Moment_Size[0][px_][py_]*sin(Theta[0][px_][py_])*cos(Phi[0][px_][py_]);
                Aux_S_y[0][pos] = Moment_Size[0][px_][py_]*sin(Theta[0][px_][py_])*sin(Phi[0][px_][py_]);
                Aux_S_z[0][pos] = Moment_Size[0][px_][py_]*cos(Theta[0][px_][py_]);

            }

            Map_Variables_to_Y(Red_Den_mat, Aux_S_x[0], Aux_S_y[0], Aux_S_z[0], YVec0);
        }


        file_out<<(ts*dt_) + Restart_Time<<"  ";

        for(int pos=0;pos<Parameters_.ns;pos++){
            file_out<< YVec0[AuxSz_to_index[pos]].real()<<"  ";
        }
        for(int pos=0;pos<Parameters_.ns;pos++){
            file_out<< YVec0[AuxSx_to_index[pos]].real()<<"  ";
        }
        for(int pos=0;pos<Parameters_.ns;pos++){
            file_out<< YVec0[AuxSy_to_index[pos]].real()<<"  ";
        }

        file_out<<endl;


        Evolve_classical_spins_Runge_Kutta(0);
        YVec0=YVec1;
    }


}



void SC_SW_ENGINE_VNE_MultiOrbSF::Evolve_classical_spins_Predictor_Corrector(){

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


void SC_SW_ENGINE_VNE_MultiOrbSF::Evolve_classical_spins_Runge_Kutta(int ts){

    Mat_4_Complex_doub Red_Den_mat_temp;
    Red_Den_mat_temp=Red_Den_mat;


    if(Runge_Kutta_order==1){
        RungeKuttaOne(YVec0, YVec1);
    }

    else if(Runge_Kutta_order==4){
        RungeKuttaFour(YVec0, YVec1);
    }


    /*
    else if(Runge_Kutta_order==5){

        //Butcher’s Fifth Order Runge-Kutta Method is used
            }
*/



    else if(Runge_Kutta_order==6){
        RungeKuttaSix(YVec0, YVec1);
    }



    Red_Den_mat_temp.clear();
}




void SC_SW_ENGINE_VNE_MultiOrbSF::Map_Variables_to_Y(Mat_4_Complex_doub & Red_Den_mat_temp, Mat_1_doub & AuxSx,
                                                     Mat_1_doub & AuxSy, Mat_1_doub & AuxSz,  Mat_1_Complex_doub & Y_)
{

    int Y_size = 3*Parameters_.ns + Parameters_.n_orbs*2*Parameters_.n_orbs*2*Parameters_.ns*Parameters_.ns;
    Y_.resize(Y_size);
    //Convention Aux_Sx--Aux_Sy---Aux_Sz-- Red_Den_mat[][][][]

    int index=0;
    for(int i=0;i<Parameters_.ns;i++){
        Y_[index]=complex<double>(AuxSx[i],0.0);
        index++;
    }
    for(int i=0;i<Parameters_.ns;i++){
        Y_[index]=complex<double>(AuxSy[i],0.0);
        index++;
    }
    for(int i=0;i<Parameters_.ns;i++){
        Y_[index]=complex<double>(AuxSz[i],0.0);
        index++;
    }

    for(int pos_i=0;pos_i<Parameters_.ns;pos_i++){
        for(int si=0;si<Parameters_.n_orbs*2;si++){
            for(int pos_j=0;pos_j<Parameters_.ns;pos_j++){
                for(int sj=0;sj<Parameters_.n_orbs*2;sj++){
                    Y_[index]=Red_Den_mat_temp[pos_i][si][pos_j][sj];
                    index++;
                }
            }
        }
    }

    assert(index==Y_size);


}



void SC_SW_ENGINE_VNE_MultiOrbSF::Map_Y_to_Variables(Mat_1_Complex_doub & Y_, Mat_4_Complex_doub & Red_Den_mat_temp, Mat_1_doub & AuxSx,
                                                     Mat_1_doub & AuxSy, Mat_1_doub & AuxSz)
{

    int Y_size = 3*Parameters_.ns + Parameters_.n_orbs*2*Parameters_.n_orbs*2*Parameters_.ns*Parameters_.ns;

    AuxSx.resize(Parameters_.ns);
    AuxSy.resize(Parameters_.ns);
    AuxSz.resize(Parameters_.ns);

    Red_Den_mat_temp.resize(Parameters_.ns);
    for(int i=0;i<Parameters_.ns;i++){
        Red_Den_mat_temp[i].resize(Parameters_.n_orbs*2);
        for(int si=0;si<Parameters_.n_orbs*2;si++){
            Red_Den_mat_temp[i][si].resize(Parameters_.ns);
            for(int j=0;j<Parameters_.ns;j++){
                Red_Den_mat_temp[i][si][j].resize(Parameters_.n_orbs*2);
            }
        }
    }


    //Convention Aux_Sx--Aux_Sy---Aux_Sz-- Red_Den_mat[][][][]


    int index=0;
    for(int i=0;i<Parameters_.ns;i++){
        AuxSx[i]=Y_[index].real();
        index++;
    }
    for(int i=0;i<Parameters_.ns;i++){
        AuxSy[i]=Y_[index].real();
        index++;
    }
    for(int i=0;i<Parameters_.ns;i++){
        AuxSz[i]=Y_[index].real();
        index++;
    }

    for(int pos_i=0;pos_i<Parameters_.ns;pos_i++){
        for(int si=0;si<Parameters_.n_orbs*2;si++){
            for(int pos_j=0;pos_j<Parameters_.ns;pos_j++){
                for(int sj=0;sj<Parameters_.n_orbs*2;sj++){
                    Red_Den_mat_temp[pos_i][si][pos_j][sj]=Y_[index];
                    index++;
                }
            }
        }
    }

    assert(index==Y_size);


}


void SC_SW_ENGINE_VNE_MultiOrbSF::IndexMapping_bw_Y_and_Variables()
{

    int Y_size = 3*Parameters_.ns + 4*Parameters_.n_orbs*Parameters_.n_orbs*Parameters_.ns*Parameters_.ns;
    index_to_AuxSx.resize(Y_size);
    index_to_AuxSy.resize(Y_size);
    index_to_AuxSz.resize(Y_size);
    index_to_RedDen.resize(Y_size);



    AuxSx_to_index.resize(Parameters_.ns);
    AuxSy_to_index.resize(Parameters_.ns);
    AuxSz_to_index.resize(Parameters_.ns);
    RedDen_to_index.resize(Parameters_.ns);
    for(int i=0;i<Parameters_.ns;i++){
        RedDen_to_index[i].resize(Parameters_.n_orbs*2);
        for(int si=0;si<Parameters_.n_orbs*2;si++){
            RedDen_to_index[i][si].resize(Parameters_.ns);
            for(int j=0;j<Parameters_.ns;j++){
                RedDen_to_index[i][si][j].resize(Parameters_.n_orbs*2);

            }
        }
    }




    //Convention Aux_Sx--Aux_Sy---Aux_Sz-- Red_Den_mat[][][][]

    int index=0;
    for(int i=0;i<Parameters_.ns;i++){
        AuxSx_to_index[i] = index;
        index_to_AuxSx[index] =i;
        index++;
    }
    for(int i=0;i<Parameters_.ns;i++){
        AuxSy_to_index[i] = index;
        index_to_AuxSy[index] =i;
        index++;
    }
    for(int i=0;i<Parameters_.ns;i++){
        AuxSz_to_index[i] = index;
        index_to_AuxSz[index] =i;
        index++;
    }


    for(int pos_i=0;pos_i<Parameters_.ns;pos_i++){
        for(int si=0;si<Parameters_.n_orbs*2;si++){
            for(int pos_j=0;pos_j<Parameters_.ns;pos_j++){
                for(int sj=0;sj<Parameters_.n_orbs*2;sj++){
                    RedDen_to_index[pos_i][si][pos_j][sj]=index;
                    index_to_RedDen[index].first=pos_i;
                    index_to_RedDen[index].second=si;
                    index_to_RedDen[index].third=pos_j;
                    index_to_RedDen[index].fourth=sj;
                    index++;
                }
            }
        }
    }

    assert(index==Y_size);


}



void SC_SW_ENGINE_VNE_MultiOrbSF::Derivative(Mat_1_Complex_doub & Y_, Mat_1_Complex_doub & dYbydt){



    double sign_check=1.0; // In Batista's paper sign_check=-1 [may be that is typo]
    dYbydt.resize(Y_.size());
    for(int i=0;i<Y_.size();i++){
        dYbydt[i]=zero_complex;
    }

    double sx, sy, sz;
    int pos_neigh;




#ifdef _OPENMP
#pragma omp parallel for default(shared) private(sx, sy, sz, pos_neigh)
#endif
    for(int pos=0;pos<Parameters_.ns;pos++){


        for(int alpha=0;alpha<Parameters_.n_orbs;alpha++){
            sx =0.5*real( Y_[RedDen_to_index[pos][alpha+Parameters_.n_orbs*1][pos][alpha+Parameters_.n_orbs*0]] + Y_[RedDen_to_index[pos][alpha+Parameters_.n_orbs*0][pos][alpha+Parameters_.n_orbs*1]] );
            sy =0.5*imag( Y_[RedDen_to_index[pos][alpha+Parameters_.n_orbs*1][pos][alpha+Parameters_.n_orbs*0]] - Y_[RedDen_to_index[pos][alpha+Parameters_.n_orbs*0][pos][alpha+Parameters_.n_orbs*1]] );
            sz =0.5*real( Y_[RedDen_to_index[pos][alpha+Parameters_.n_orbs*0][pos][alpha+Parameters_.n_orbs*0]] - Y_[RedDen_to_index[pos][alpha+Parameters_.n_orbs*1][pos][alpha+Parameters_.n_orbs*1]] );


            dYbydt[AuxSx_to_index[pos]] += (Parameters_.J_Hund[alpha])*(sy*Y_[AuxSz_to_index[pos]] - sz*Y_[AuxSy_to_index[pos]]);
            dYbydt[AuxSy_to_index[pos]] += (Parameters_.J_Hund[alpha])*(sz*Y_[AuxSx_to_index[pos]] - sx*Y_[AuxSz_to_index[pos]]);
            dYbydt[AuxSz_to_index[pos]] += (Parameters_.J_Hund[alpha])*(sx*Y_[AuxSy_to_index[pos]] - sy*Y_[AuxSx_to_index[pos]]);
        }

        for (int ng=0;ng<4;ng++){
            pos_neigh = Coordinates_.neigh(pos,ng);
            dYbydt[AuxSx_to_index[pos]] +=J_Classical*(Y_[AuxSy_to_index[pos_neigh]]*Y_[AuxSz_to_index[pos]] - Y_[AuxSz_to_index[pos_neigh]]*Y_[AuxSy_to_index[pos]]);
            dYbydt[AuxSy_to_index[pos]] +=J_Classical*(Y_[AuxSz_to_index[pos_neigh]]*Y_[AuxSx_to_index[pos]] - Y_[AuxSx_to_index[pos_neigh]]*Y_[AuxSz_to_index[pos]]);
            dYbydt[AuxSz_to_index[pos]] +=J_Classical*(Y_[AuxSx_to_index[pos_neigh]]*Y_[AuxSy_to_index[pos]] - Y_[AuxSy_to_index[pos_neigh]]*Y_[AuxSx_to_index[pos]]);
        }


    }


    int k;
    double hopping_val;


#ifdef _OPENMP
#pragma omp parallel for default(shared) private(k, hopping_val)
#endif
    for(int pos_i=0;pos_i<Parameters_.ns;pos_i++){
        for(int si=0;si<2;si++){
            for(int alpha=0;alpha<Parameters_.n_orbs;alpha++){
                for(int pos_j=0;pos_j<Parameters_.ns;pos_j++){
                    for(int sj=0;sj<2;sj++){
                        for(int beta=0;beta<Parameters_.n_orbs;beta++){


                            //------------------i-k connection :start---------------------//
                            //All neighbours
                            for (int ng=0;ng<4;ng++){
                                k = Coordinates_.neigh(pos_i,ng);

                                for(int eta=0;eta<Parameters_.n_orbs;eta++){
                                if(ng==Coordinates_.PX || ng==Coordinates_.MX){
                                    hopping_val = Parameters_.hopping_NN_X(eta,alpha);
                                }
                                if(ng==Coordinates_.PY || ng==Coordinates_.MY){
                                   hopping_val = Parameters_.hopping_NN_Y(eta,alpha);
                                }

                                dYbydt[RedDen_to_index[pos_i][alpha+Parameters_.n_orbs*si][pos_j][beta+Parameters_.n_orbs*sj]] +=
                                        iota_complex*sign_check*complex<double>(1.0*hopping_val, 0.0)*
                                        Y_[RedDen_to_index[k][eta+Parameters_.n_orbs*si][pos_j][beta+Parameters_.n_orbs*sj]];
                            }
                            }

                            //------------------i-k connection :done---------------------//


                            //------------------j-k connection:start---------------------//
                            for (int ng=0;ng<4;ng++){
                                k = Coordinates_.neigh(pos_j,ng);

                                for(int eta=0;eta<Parameters_.n_orbs;eta++){
                                if(ng==Coordinates_.PX || ng==Coordinates_.MX){
                                    hopping_val = Parameters_.hopping_NN_X(beta,eta);
                                }
                                if(ng==Coordinates_.PY || ng==Coordinates_.MY){
                                   hopping_val = Parameters_.hopping_NN_Y(beta,eta);
                                }

                                dYbydt[RedDen_to_index[pos_i][alpha+Parameters_.n_orbs*si][pos_j][beta+Parameters_.n_orbs*sj]] +=
                                        iota_complex*sign_check*complex<double>(-1.0*hopping_val, 0.0)*
                                        Y_[RedDen_to_index[pos_i][alpha+Parameters_.n_orbs*si][k][eta+Parameters_.n_orbs*sj]];
                            }
                            }

                            //------------------j-k connection :done---------------------//




                            //-------------------coupling with Auxlliary spins--------------//

                            for(int tau=0;tau<2;tau++){

                                dYbydt[RedDen_to_index[pos_i][alpha+Parameters_.n_orbs*si][pos_j][beta+Parameters_.n_orbs*sj]] +=
                                        iota_complex*sign_check*(0.5*Jval_array[pos_i])*(
                                            (Y_[AuxSx_to_index[pos_i]])*Pauli_x[si][tau]
                                        +
                                        (Y_[AuxSy_to_index[pos_i]])*Pauli_y[si][tau]
                                        +
                                        (Y_[AuxSz_to_index[pos_i]])*Pauli_z[si][tau]
                                        )*Y_[RedDen_to_index[pos_i][alpha+Parameters_.n_orbs*tau][pos_j][beta+Parameters_.n_orbs*sj]];

                                dYbydt[RedDen_to_index[pos_i][alpha+Parameters_.n_orbs*si][pos_j][beta+Parameters_.n_orbs*sj]] +=
                                        iota_complex*sign_check*(-0.5*Jval_array[pos_j])*(
                                            (Y_[AuxSx_to_index[pos_j]])*Pauli_x[tau][sj]
                                        +
                                        (Y_[AuxSy_to_index[pos_j]])*Pauli_y[tau][sj]
                                        +
                                        (Y_[AuxSz_to_index[pos_j]])*Pauli_z[tau][sj]
                                        )*Y_[RedDen_to_index[pos_i][alpha+Parameters_.n_orbs*si][pos_j][beta+Parameters_.n_orbs*tau]];


                            }

                            //-------------------coupling with Auxlliary spins: done--------------//


                        }
                    }
                }
            }
        }
    }





}



void SC_SW_ENGINE_VNE_MultiOrbSF::RungeKuttaOne(Mat_1_Complex_doub & Yn, Mat_1_Complex_doub & Ynp1){


    Ynp1.resize(Yn.size());
    Mat_1_Complex_doub K1;

    //step_no==0
    Derivative(Yn,K1);

    for(int i=0;i<Yn.size();i++){
        Ynp1[i] = Yn[i] + (dt_)*(K1[i]);

    }


}
void SC_SW_ENGINE_VNE_MultiOrbSF::RungeKuttaFour(Mat_1_Complex_doub & Yn, Mat_1_Complex_doub & Ynp1){



    Mat_1_Complex_doub Y_temp;
    Y_temp.resize(Yn.size());

    Ynp1.resize(Yn.size());
    Mat_1_Complex_doub K1, K2, K3, K4;


    //step_no==0
    Derivative(Yn,K1);
    for(int i=0;i<Yn.size();i++){
        Y_temp[i] = Yn[i]  + (K1[i]*dt_*0.5);
    }

    //step_no==1
    Derivative(Y_temp,K2);
    for(int i=0;i<Yn.size();i++){
        Y_temp[i] = Yn[i]  + (K2[i]*dt_*0.5);
    }

    //step_no==2
    Derivative(Y_temp,K3);
    for(int i=0;i<Yn.size();i++){
        Y_temp[i] = Yn[i]  + (K3[i]*dt_);
    }

    //step_no==3
    Derivative(Y_temp,K4);


    for(int i=0;i<Yn.size();i++){
        Ynp1[i] = Yn[i] + (dt_/6.0)*(K1[i] + 2.0*K2[i] + 2.0*K3[i] + K4[i]);

    }


}


void SC_SW_ENGINE_VNE_MultiOrbSF::RungeKuttaSix(Mat_1_Complex_doub & Yn, Mat_1_Complex_doub & Ynp1){



    //Sixth Order Runge-Kutta Method is used
    //From "https://www.ams.org/journals/mcom/1968-22-102/S0025-5718-68-99876-1/S0025-5718-68-99876-1.pdf"
    //or google "An Explicit Sixth-Order Runge-Kutta Formula By H. A. Luther"



    double nu, surd;

    nu=0.9;
    surd=-1.0*sqrt(21.0);

    Mat_1_Complex_doub Y_temp;
    Y_temp.resize(Yn.size());

    Ynp1.resize(Yn.size());
    Mat_1_Complex_doub K1, K2, K3, K4, K5, K6, K7;


    //step_no==0
    Derivative(Yn,K1);


    //step_no==1
    for(int i=0;i<Yn.size();i++){
        Y_temp[i] = Yn[i]  + (K1[i]*dt_*nu);
    }
    Derivative(Y_temp,K2);


    //step_no=2
    for(int i=0;i<Yn.size();i++){
        Y_temp[i] = Yn[i]  + ( ((4.0*nu-1.0)*K1[i]) +  K2[i])*dt_*(1.0/(8.0*nu));
    }
    Derivative(Y_temp,K3);


    //step_no=3
    for(int i=0;i<Yn.size();i++){
        Y_temp[i] = Yn[i]  + ( ((10.0*nu-2.0)*K1[i]) +  (2.0*K2[i])   +   (8.0*nu*K3[i]) )*dt_*(1.0/(27.0*nu));
    }
    Derivative(Y_temp,K4);


    //step_no=4
    for(int i=0;i<Yn.size();i++){
        Y_temp[i] = Yn[i]  + ( (-1.0*(  (77.0*nu - 56.0)  + ((17.0*nu - 8.0)*surd)   )*K1[i]) +
                               (-8.0*(7.0 + surd)*K2[i])   +   (48.0*(7.0 + surd)*nu*K3[i])
                               + (-3.0*(21.0 + surd)*nu*K4[i]) )*dt_*(1.0/(392.0*nu));
    }
    Derivative(Y_temp,K5);


    //step_no=5
    for(int i=0;i<Yn.size();i++){
        Y_temp[i] = Yn[i]  +  (  (-5.0*( (287.0*nu - 56.0) - ((59.0*nu -8.0)*surd) )*K1[i])
                                 +  (-40.0*(7.0 - surd)*K2[i])
                                 +  (320.0*surd*nu*K3[i])
                                 +  (3.0*(21.0 - 121.0*surd)*nu*K4[i])
                                 +  (392.0*(6.0 - surd)*nu*K5[i])  )*dt_*(1.0/(1960.0*nu));
    }
    Derivative(Y_temp,K6);


    //step_no=6
    for(int i=0;i<Yn.size();i++){
        Y_temp[i] = Yn[i]  +  (  (15.0*( (30.0*nu - 8.0)  -  (7.0*nu*surd) )*K1[i])
                                 +  (120.0*K2[i])
                                 +  (-40.0*(5.0 + 7.0*surd)*nu*K3[i])
                                 +  (63.0*(2.0 + 3.0*surd)*nu*K4[i])
                                 +  (-14.0*(49.0 - 9.0*surd)*nu*K5[i])
                                 + (70.0*(7.0 + surd)*nu*K6[i]) )*dt_*(1.0/(180.0*nu));
    }
    Derivative(Y_temp,K7);




    for(int i=0;i<Yn.size();i++){
        Ynp1[i] = Yn[i] + (dt_/180.0)*(9.0*K1[i] + 64.0*K3[i] + 49.0*K5[i] + 49.0*K6[i]  +  9.0*K7[i]);

    }


}


void SC_SW_ENGINE_VNE_MultiOrbSF::Write_final_time_result(){



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



void SC_SW_ENGINE_VNE_MultiOrbSF::Read_Restart_Data(){

    complex<double> one(1,0);
    complex<double> iota(0,1);

    ifstream file_in_DM(Restart_DM_file.c_str());

    int temp;
    double temp_real, temp_imag;
    for(int pi=0;pi<Red_Den_mat.size();pi++){
        for(int si=0;si<Red_Den_mat[pi].size();si++){

            for(int pj=0;pj<Red_Den_mat.size();pj++){
                for(int sj=0;sj<Red_Den_mat[pj].size();sj++){
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



void SC_SW_ENGINE_VNE_MultiOrbSF::Read_parameters(string filename){


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
    cout<<"Temperature = "<<Parameters_.temp<<endl;
    cout<<"N_total = "<<Parameters_.Fill*Parameters_.ns<<endl;
    cout<<"RANDOM_NO_SEED = "<<Parameters_.RandomSeed<<endl;

    if(SAVE){
        cout<<"Final Results for time = "<<time_max<<endl;
        cout<<"will be written in "<< Save_DM_file<<", "<<Save_Classical_file<<endl;
    }








}
















































