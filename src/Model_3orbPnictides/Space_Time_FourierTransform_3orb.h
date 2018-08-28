#include <iostream>
#include <math.h>  //fabs(double x) =|x|
#include <algorithm>
#include <stdlib.h>  //for div(q,n).rem(quot),rand
#include <assert.h>
#include <time.h>
#include <fstream>
#include <sstream>
#include <limits>
#include <iomanip>
#include <stdio.h>
#include "Coordinates.h"
#include "Hamiltonian.h"
#include "Matrix.h"
#include "MFParams.h"
#include "Observables.h"
#include "ParametersEngine.h"
#include "Spin_dynamics_VNE_3orbPnictides_engine.h"
#include "../../tensor_type.h"
#define PI acos(-1.0)

using namespace std;

#ifndef ST_Fourier_3orb_3orb
#define ST_Fourier_3orb_3orb

class ST_Fourier_3orb{

public:
    ST_Fourier_3orb(Parameters& Parameters__, Coordinates& Coordinates__,
                                       MFParams& MFParams__, Hamiltonian& Hamiltonian__,
                    Observables& Observables__, SC_SW_ENGINE_VNE_3orbPnictides& SC_SW_ENGINE_VNE_3orbPnictides__)
            : Parameters_(Parameters__), Coordinates_(Coordinates__), MFParams_(MFParams__),
              Hamiltonian_(Hamiltonian__), Observables_(Observables__), SC_SW_ENGINE_VNE_3orbPnictides_(SC_SW_ENGINE_VNE_3orbPnictides__)
    {

    }

double w_min,w_max,dw;
double w_conv;
int n_wpoints;

double time_max;
double dt_;
int time_steps;

int no_of_processors;

Mat_1_string conf_inputs;
int No_Of_Inputs;


string Skw_out;
string Skw_out_full;


/*
Mat_3_doub SS_rt;
Mat_2_doub S_rt;


Mat_3_doub sqsq_rt;
Mat_2_doub sq_rt;

Mat_3_doub TT_rt;
Mat_2_doub T_rt;
*/

Mat_2_doub S_tr;

Mat_3_doub S_rw; //
Mat_3_doub S_kw; //

Mat_3_doub s_quantum_rw; //
Mat_3_doub s_quantum_kw; //

Mat_3_doub T_rw; //
Mat_3_doub T_kw; //


Mat_1_doub Sz_eq,Sx_eq,Sy_eq;
Mat_1_doub sz_eq,sx_eq,sy_eq;
Mat_1_doub Tz_eq,Tx_eq,Ty_eq;


Mat_1_doub Sz_t,Sx_t,Sy_t;
Mat_1_doub sz_t,sx_t,sy_t;
Mat_1_doub Tz_t,Tx_t,Ty_t;

Parameters& Parameters_;
Coordinates& Coordinates_;
MFParams& MFParams_;
Hamiltonian& Hamiltonian_;
Observables& Observables_;
SC_SW_ENGINE_VNE_3orbPnictides& SC_SW_ENGINE_VNE_3orbPnictides_;




void Initialize_engine();

void Read_parameters();

void Perform_Averaging_on_one_point();
void Perform_Smarter_Averaging_on_one_point();


void Calculate_Skw_from_Srt_file(string filename, string fileout);

};




void ST_Fourier_3orb::Initialize_engine(){

    n_wpoints=(int) ((w_max-w_min)/dw + 0.5);
    cout<<"n_wpoints = "<<n_wpoints<<endl;



    ifstream check_file;
    check_file.open(conf_inputs[0].c_str());
    string line_temp, line_last;



    double t0, t1;
    for(int curLine = 0; getline(check_file, line_temp); curLine++) {

        if(curLine==1){
            stringstream line_temp_ss(line_temp, stringstream::in);
            line_temp_ss>>t0;
        }
        if(curLine==2){
            stringstream line_temp_ss(line_temp, stringstream::in);
            line_temp_ss>>t1;
        }
        line_last=line_temp;
    }


    dt_ = t1 - t0;


    stringstream line_temp_ss(line_last, stringstream::in);
    line_temp_ss>>time_max;

    //cout<<line_temp<<endl;
    time_steps=(int) ((time_max - t0)/dt_ + 0.5);

    if(t0 != 0){
        cout<<endl;
        cout<<"ERROR:"<<endl;
        cout<<"Initial time [1st row] have to be 0.0 in the file \""<<conf_inputs[0]<<"\""<<endl<<endl;
        assert (false);
    }

    cout <<"dt_ = "<<dt_<<endl;
    cout<< "time_max = "<<time_max<<endl;
    cout<< "time_steps = "<<time_steps<<endl;


    check_file.clear();
    check_file.seekg(0, ios::beg);


}
void ST_Fourier_3orb::Perform_Averaging_on_one_point(){

    string Temp_file_Srt="Temp_file_Srt.txt";
    ofstream Temp_file_Srt_out(Temp_file_Srt.c_str());


    Temp_file_Srt_out<<"#time   Sz   Sx  Sy  sz  sx  sy"<<endl;

    string line_temp;
    double row_time, double_temp;
    int row_ts;

    Mat_1_doub Sr;
    Sr.resize(6*Parameters_.ns);


    for(int ts=0;ts<=time_steps;ts++){

        for(int r=0;r<Sr.size();r++){
            Sr[r]=0.0;
        }

        Temp_file_Srt_out<<ts*dt_<<"  ";

        for(int conf=0;conf<No_Of_Inputs;conf++){
            ifstream Specific_conf_in;
            Specific_conf_in.open(conf_inputs[conf].c_str());

            int stop_here = ts+2; //First line is #Sz Sx ..........

            for(int curLine = 0;curLine < stop_here ; curLine++) {

                getline(Specific_conf_in, line_temp);
                //cout<<line_temp<<endl;

            }

            // int check;
            //cin>>check;
            stringstream line_temp_ss(line_temp, stringstream::in);
            line_temp_ss>>row_time;



            row_ts = (int) ((row_time)/dt_ + 0.5);

            if(row_ts != ts){
                cout<<endl;
                cout<<"Error;"<<endl;
                cout<< row_ts <<"  "<<ts<<endl;
                cout<< "Problem while reading "<<conf_inputs[conf]<<endl<<endl;
                assert(false);
            }

            for(int r=0;r<6*Parameters_.ns;r++){
                line_temp_ss>>double_temp;

                Sr[r] += double_temp;
            }

            Specific_conf_in.clear();
            Specific_conf_in.seekg(0, ios::beg);


        }

        for(int r=0;r<6*Parameters_.ns;r++){
            Temp_file_Srt_out<<Sr[r]/(No_Of_Inputs)<<"  ";
        }
        Temp_file_Srt_out<<endl;

    }


}


void ST_Fourier_3orb::Perform_Smarter_Averaging_on_one_point(){

    string Temp_file_Srt="Temp_file_Srt.txt";
    ofstream Temp_file_Srt_out(Temp_file_Srt.c_str());


    Temp_file_Srt_out<<"#time   Sz   Sx  Sy  sz  sx  sy"<<endl;

    string line_temp;
    double row_time, double_temp;
    int row_ts;


    S_tr.resize(time_steps+1);

    for(int ts=0;ts<=time_steps;ts++){
        S_tr[ts].resize(6*Parameters_.ns);
        for(int r=0;r<6*Parameters_.ns;r++){
            S_tr[ts][r]=0;
        }

    }




    for(int conf=0;conf<No_Of_Inputs;conf++){


        ifstream Specific_conf_in;
        Specific_conf_in.open(conf_inputs[conf].c_str());
        getline(Specific_conf_in, line_temp);

        for(int ts=0;ts<=time_steps;ts++){

            getline(Specific_conf_in, line_temp);
            stringstream line_temp_ss(line_temp, stringstream::in);

            line_temp_ss>>row_time;
            row_ts = (int) ((row_time)/dt_ + 0.5);
            if(row_ts != ts){
                cout<<endl;
                cout<<"Error;"<<endl;
                cout<< row_time<<"   "<<row_ts <<"  "<<ts<<endl;
                cout<< "Problem while reading "<<conf_inputs[conf]<<endl<<endl;
                assert(false);
            }

            for(int r=0;r<6*Parameters_.ns;r++){
                line_temp_ss>>double_temp;
                S_tr[ts][r]+= double_temp;
            }


        }

    }



    for(int ts=0;ts<=time_steps;ts++){
        Temp_file_Srt_out<<ts*dt_<<"  ";
        for(int r=0;r<6*Parameters_.ns;r++){
            Temp_file_Srt_out<<S_tr[ts][r]/(No_Of_Inputs)<<"  ";
        }
        Temp_file_Srt_out<<endl;
    }


}

void ST_Fourier_3orb::Calculate_Skw_from_Srt_file( string filename, string fileout){


    S_rw.resize(Parameters_.ns);
    s_quantum_rw.resize(Parameters_.ns);
    T_rw.resize(Parameters_.ns);

    for(int pos_i=0;pos_i<Parameters_.ns;pos_i++){
        S_rw[pos_i].resize(Parameters_.ns);
        s_quantum_rw[pos_i].resize(Parameters_.ns);
        T_rw[pos_i].resize(Parameters_.ns);

        for(int pos_j=0;pos_j<Parameters_.ns;pos_j++){
            S_rw[pos_i][pos_j].resize(n_wpoints);
            s_quantum_rw[pos_i][pos_j].resize(n_wpoints);
            T_rw[pos_i][pos_j].resize(n_wpoints);

            for(int wi=0;wi<n_wpoints;wi++){
                S_rw[pos_i][pos_j][wi]=0;
                s_quantum_rw[pos_i][pos_j][wi]=0;
                T_rw[pos_i][pos_j][wi]=0;

            }
        }
    }

    Sz_eq.resize(Parameters_.ns);Sx_eq.resize(Parameters_.ns);Sy_eq.resize(Parameters_.ns);
    sz_eq.resize(Parameters_.ns);sx_eq.resize(Parameters_.ns);sy_eq.resize(Parameters_.ns);

    Sz_t.resize(Parameters_.ns);Sx_t.resize(Parameters_.ns);Sy_t.resize(Parameters_.ns);
    sz_t.resize(Parameters_.ns);sx_t.resize(Parameters_.ns);sy_t.resize(Parameters_.ns);


    string line_temp;
    double temp_waste;
    ifstream file_Srt_in;
    file_Srt_in.open(filename.c_str());

    getline(file_Srt_in, line_temp);//First line is a waste
    getline(file_Srt_in, line_temp);
    stringstream line_temp_ss(line_temp, stringstream::in);

    line_temp_ss>>temp_waste;

    for(int pos_i=0;pos_i<Parameters_.ns;pos_i++){
        line_temp_ss>>Sz_eq[pos_i];
    }
    for(int pos_i=0;pos_i<Parameters_.ns;pos_i++){
        line_temp_ss>>Sx_eq[pos_i];
    }
    for(int pos_i=0;pos_i<Parameters_.ns;pos_i++){
        line_temp_ss>>Sy_eq[pos_i];
    }
    for(int pos_i=0;pos_i<Parameters_.ns;pos_i++){
        line_temp_ss>>sz_eq[pos_i];
    }
    for(int pos_i=0;pos_i<Parameters_.ns;pos_i++){
        line_temp_ss>>sx_eq[pos_i];
    }
    for(int pos_i=0;pos_i<Parameters_.ns;pos_i++){
        line_temp_ss>>sy_eq[pos_i];
    }

    file_Srt_in.clear();
    file_Srt_in.seekg(0, ios::beg);





    //open again
    ifstream file_Srt_in2;
    file_Srt_in2.open(filename.c_str());

    getline(file_Srt_in2, line_temp);//First line is a waste

//cout<<line_temp<<endl;
    for(int ts=0;ts<=time_steps;ts++){
        getline(file_Srt_in2, line_temp);
  //      cout<<line_temp<<endl;
        stringstream line_temp2_ss(line_temp, stringstream::in);


        line_temp2_ss>>temp_waste;
        for(int pos_i=0;pos_i<Parameters_.ns;pos_i++){
            line_temp2_ss>>Sz_t[pos_i];
        }
        for(int pos_i=0;pos_i<Parameters_.ns;pos_i++){
            line_temp2_ss>>Sx_t[pos_i];
        }
        for(int pos_i=0;pos_i<Parameters_.ns;pos_i++){
            line_temp2_ss>>Sy_t[pos_i];
        }
        for(int pos_i=0;pos_i<Parameters_.ns;pos_i++){
            line_temp2_ss>>sz_t[pos_i];
        }
        for(int pos_i=0;pos_i<Parameters_.ns;pos_i++){
            line_temp2_ss>>sx_t[pos_i];
        }
        for(int pos_i=0;pos_i<Parameters_.ns;pos_i++){
            line_temp2_ss>>sy_t[pos_i];
        }






#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
        for(int wi=0;wi<n_wpoints;wi++){

            for(int pos_i=0;pos_i<Parameters_.ns;pos_i++){
                for(int pos_j=0;pos_j<Parameters_.ns;pos_j++){



                    S_rw[pos_i][pos_j][wi] += cos((-wi * dw) * (ts* dt_))*exp(-0.5*(ts*dt_*w_conv*ts*dt_*w_conv))*dt_*(
                                ( (Sz_eq[pos_i]*Sz_t[pos_j]  )  )
                                +
                                ( (Sx_eq[pos_i]*Sx_t[pos_j]  )  )
                                +
                                ( (Sy_eq[pos_i]*Sy_t[pos_j]  )  )
                                );


                    s_quantum_rw[pos_i][pos_j][wi] += cos( (-wi * dw) * (ts* dt_))*exp(-0.5*(ts*dt_*w_conv*ts*dt_*w_conv))*dt_*(
                                ( (sz_eq[pos_i]*( sz_t[pos_j] - 0.0*sz_eq[pos_j])  )  )
                                +
                                ( (sx_eq[pos_i]*( sx_t[pos_j] - 0.0*sx_eq[pos_j])  )  )
                                +
                                ( (sy_eq[pos_i]*( sy_t[pos_j] - 0.0*sy_eq[pos_j])  )  )
                                );




                    /*    T_rw[pos_i][pos_j][wi] +=exp(iota * (-wi * dw) * (ts* dt_))*exp(-0.5*(ts*dt_*w_conv*ts*dt_*w_conv))*dt_* (


                            );*/






                }
            }
        }

    }







    //Now "rw - space" to "kw - space"

    S_kw.resize(Parameters_.lx);
    s_quantum_kw.resize(Parameters_.lx);
    T_kw.resize(Parameters_.lx);


    for(int x_i=0;x_i<Parameters_.lx;x_i++){

        S_kw[x_i].resize(Parameters_.ly);
        s_quantum_kw[x_i].resize(Parameters_.ly);
        T_kw[x_i].resize(Parameters_.ly);

        for(int y_j=0;y_j<Parameters_.ly;y_j++){

            S_kw[x_i][y_j].resize(n_wpoints);
            s_quantum_kw[x_i][y_j].resize(n_wpoints);
            T_kw[x_i][y_j].resize(n_wpoints);

        }
    }





    double kx, ky;
    int pos_i, pos_j;
    for(int nx=0;nx<Parameters_.lx;nx++){
        for(int ny=0;ny<Parameters_.ly;ny++){

            kx = (2*nx*PI)/(1.0*Parameters_.lx);
            ky = (2*ny*PI)/(1.0*Parameters_.ly);


#ifdef _OPENMP
#pragma omp parallel for default(shared) private(pos_i,pos_j)
#endif

            for(int wi=0;wi<n_wpoints;wi++){
                double temp3=0;
                double temp2=0;
                double temp=0;
                for(int x_i=0;x_i<Parameters_.lx;x_i++){
                    for(int y_i=0;y_i<Parameters_.ly;y_i++){

                        pos_i = y_i*Parameters_.lx + x_i;

                        for(int x_j=0;x_j<Parameters_.lx;x_j++){
                            for(int y_j=0;y_j<Parameters_.ly;y_j++){

                                pos_j = y_j*Parameters_.lx + x_j;

                                temp += S_rw[pos_i][pos_j][wi]*cos(( (x_j - x_i)*kx +  (y_j - y_i)*ky ) );
                                //temp += S_rw[pos_i][pos_j][wi]*exp(iota*( (x_j - x_i)*kx +  (y_j - y_i)*ky ) );
                                temp2 += s_quantum_rw[pos_i][pos_j][wi]*cos(( (x_j - x_i)*kx +  (y_j - y_i)*ky ) );
                                temp3 += T_rw[pos_i][pos_j][wi]*cos(( (x_j - x_i)*kx +  (y_j - y_i)*ky ) ) ;


                            }
                        }
                    }
                }
                //file_out2<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<temp.real()<<"   "<<temp.imag()<<"    "<<temp2.real()<<"   "<<temp2.imag()<<"    "<<temp3.real()<<"   "<<temp3.imag()<<endl;
                //file_out_full<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<wi<<"   "<<temp.real()<<"   "<<temp.imag()<<"    "<<"none"<<"   "<<"none"<<"    "<<"none"<<"   "<<"none"<<endl;
                S_kw[nx][ny][wi]=temp;
                s_quantum_kw[nx][ny][wi]=temp2;
                T_kw[nx][ny][wi]=temp3;
            }



        }

    }





    ofstream file_out_full(fileout.c_str());

    file_out_full<<"#nx   ny   k_ind   wi*dw   wi   S_kw[nx][ny][wi]    s_quantum_kw[nx][ny][wi]          T_kw[nx][ny][wi]"<<endl;
    int k_ind=0;
    for(int nx=0;nx<Parameters_.lx;nx++){
        for(int ny=0;ny<Parameters_.ly;ny++){

            for(int wi=0;wi<n_wpoints;wi++){

                //file_out2<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<temp.real()<<"   "<<temp.imag()<<"    "<<temp2.real()<<"   "<<temp2.imag()<<"    "<<temp3.real()<<"   "<<temp3.imag()<<endl;
                file_out_full<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<wi<<"   "<<S_kw[nx][ny][wi]<<"    "<<s_quantum_kw[nx][ny][wi]<<
                               "     "<<T_kw[nx][ny][wi]<<endl;
                //Skw_Mat[nx][ny][wi]=temp;
            }

            k_ind +=1;
            file_out_full<<endl;

        }

    }








}

void ST_Fourier_3orb::Read_parameters(){


    w_min=SC_SW_ENGINE_VNE_3orbPnictides_.w_min;
    w_max=SC_SW_ENGINE_VNE_3orbPnictides_.w_max;

    dw=SC_SW_ENGINE_VNE_3orbPnictides_.dw;
    w_conv=SC_SW_ENGINE_VNE_3orbPnictides_.w_conv;


    no_of_processors=SC_SW_ENGINE_VNE_3orbPnictides_.no_of_processors;


}


























#endif


