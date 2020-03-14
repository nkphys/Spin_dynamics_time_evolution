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
#include "Coordinates_MCMF.h"
#include "Hamiltonian_MCMF.h"
#include "../../Matrix.h"
#include "MFParams_MCMF.h"
#include "Observables_MCMF.h"
#include "ParametersEngine_MCMF.h"
#include "Spin_dynamics_VNE_1orb_engine_MCMF.h"
#include "../../tensor_type.h"
#include "../../FftComplex.h"
#define PI acos(-1.0)

using namespace std;

#ifndef ST_Fourier_1orbMCMF
#define ST_Fourier_1orbMCMF

class ST_Fourier_1orb_MCMF{

public:
    ST_Fourier_1orb_MCMF(Parameters_MCMF& Parameters__, Coordinates_MCMF& Coordinates__,
                         MFParams_MCMF& MFParams__, Hamiltonian_MCMF& Hamiltonian__,
                         Observables_MCMF& Observables__, SC_SW_ENGINE_VNE_1orb_MCMF& SC_SW_ENGINE_VNE_1orb_MCMF__)
        : Parameters_(Parameters__), Coordinates_(Coordinates__), MFParams_(MFParams__),
          Hamiltonian_(Hamiltonian__), Observables_(Observables__), SC_SW_ENGINE_VNE_1orb_MCMF_(SC_SW_ENGINE_VNE_1orb_MCMF__)
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


    bool Use_FFT;


    /*
Mat_3_doub SS_rt;
Mat_2_doub S_rt;


Mat_3_doub sqsq_rt;
Mat_2_doub sq_rt;

Mat_3_doub TT_rt;
Mat_2_doub T_rt;
*/

    Mat_2_doub S_tr;
    Mat_3_doub C_tr; //Space-time dispaced correlation, avg over Ensemble

    Mat_3_Complex_doub S_rw; //
    Mat_3_Complex_doub S_kw; //

    Mat_3_Complex_doub s_quantum_rw; //
    Mat_3_Complex_doub s_quantum_kw; //

    Mat_1_doub Sz_eq,Sx_eq,Sy_eq;
    Mat_1_doub sz_eq,sx_eq,sy_eq;

    Mat_1_doub Sz_t,Sx_t,Sy_t;
    Mat_1_doub sz_t,sx_t,sy_t;

    Parameters_MCMF& Parameters_;
    Coordinates_MCMF& Coordinates_;
    MFParams_MCMF& MFParams_;
    Hamiltonian_MCMF& Hamiltonian_;
    Observables_MCMF& Observables_;
    SC_SW_ENGINE_VNE_1orb_MCMF& SC_SW_ENGINE_VNE_1orb_MCMF_;




    void Initialize_engine();

    void Read_parameters();

    void Perform_Averaging_on_one_point();
    void Perform_Smarter_Averaging_on_one_point();


    void Calculate_Skw_from_Srt_file(string filename, string fileout);
    void Calculate_Skw_from_Crt(string fileout);


    void Calculate_SpaceTimeDisplacedCorrelations(string STdisplaced_Crt_fileout);

};




void ST_Fourier_1orb_MCMF::Initialize_engine(){



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

    if(Use_FFT==true){
        cout<<"FFT is used, so dw, w_min are fixed, and w_max is used from input file."<<endl;
        dw = (2.0*PI)/ time_max;
        w_min = 0.0;
        cout<<"dw = 2*pi/T_max = "<< dw<<endl;
        cout<<"w_min = "<< w_min<<endl;
        cout<<"w_max = "<< w_max<<endl;
        n_wpoints=(int) ((w_max-w_min)/dw + 0.5);
        cout<<"n_wpoints = "<<n_wpoints<<endl;
    }
    else{
        cout<<"dw, w_min, w_max are used from input file."<<endl;
        cout<<"dw = "<< dw<<endl;
        cout<<"w_min = "<< w_min<<endl;
        cout<<"w_max = "<< w_max<<endl;
        n_wpoints=(int) ((w_max-w_min)/dw + 0.5);
        cout<<"n_wpoints = "<<n_wpoints<<endl;

    }


    check_file.clear();
    check_file.seekg(0, ios::beg);


}
void ST_Fourier_1orb_MCMF::Perform_Averaging_on_one_point(){

    string Temp_file_Srt="Average_Srt.txt";
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


void ST_Fourier_1orb_MCMF::Perform_Smarter_Averaging_on_one_point(){


    int N_p;
    int no_threads_used;
    double begin_time, end_time;
    clock_t oprt_SB_time = clock();

    cout<<"<S(r,t)> is being calculated"<<endl;
#ifdef _OPENMP
    N_p = omp_get_max_threads();
    cout<<"Maximum threads can be used parallely = "<<N_p<<endl;
    cout<<"No. of threads you asked for = "<<no_of_processors<<endl;

    begin_time = omp_get_wtime();
#endif




    ostringstream ostr_w_conv;
    ostr_w_conv << w_conv;
    string string_w_conv = ostr_w_conv.str();

    string Temp_file_Srt="Average_Srt_w_conv"+string_w_conv+".txt";
    ofstream Temp_file_Srt_out(Temp_file_Srt.c_str());


    Temp_file_Srt_out<<"#time   Sz   Sx  Sy  sz  sx  sy"<<endl;


    S_tr.resize(time_steps+1);


#ifdef _OPENMP
#pragma omp parallel for default(shared) //private()
#endif
    for(int ts=0;ts<=time_steps;ts++){
        S_tr[ts].resize(6*Parameters_.ns);
        for(int r=0;r<6*Parameters_.ns;r++){
            S_tr[ts][r]=0;
        }
    }


#ifdef _OPENMP
    no_threads_used = min(no_of_processors, No_Of_Inputs);
    omp_set_num_threads(no_threads_used);
    N_p = omp_get_max_threads();
    cout<<"threads being used parallely = "<<N_p<<endl;
    cout<<"No. of threads you asked for = "<<no_of_processors<<endl;
#pragma omp parallel for default(shared) //private()
#endif
    for(int conf=0;conf<No_Of_Inputs;conf++){

#ifdef _OPENMP
        //cout<< "using thread no "<<omp_get_thread_num()<<" out of total "<<omp_get_num_threads()<<"threads"<<endl;
#endif

        string line_temp;
        double row_time, double_temp;
        int row_ts;
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
            S_tr[ts][r] = S_tr[ts][r]/(No_Of_Inputs);
            Temp_file_Srt_out<<S_tr[ts][r]<<"  ";
        }
        Temp_file_Srt_out<<endl;
    }


#ifdef _OPENMP
    end_time = omp_get_wtime();
    cout<<"Time to Calculate <S(r,t)>, <s(r,t)> [actual time, using OMP]: "<<double(end_time - begin_time)<<endl;//cout<<"here"<<endl;

#endif

    cout<<"Time to Calculate <S(r,t)>, <s(r,t)> : "<<double( clock() - oprt_SB_time ) / (double)CLOCKS_PER_SEC<<endl;//cout<<"here"<<endl;





}


void ST_Fourier_1orb_MCMF::Calculate_SpaceTimeDisplacedCorrelations(string STdisplaced_Crt_fileout){


    int N_p;
    int no_threads_used;
    double begin_time, end_time;
    clock_t oprt_SB_time = clock();


#ifdef _OPENMP
    begin_time = omp_get_wtime();
#endif


    Mat_2_doub S_r_t0; //[Conf_no][r \in [0,6*N-1] ]
    Mat_2_doub S_r_t; //[Conf_no][r \in [0,6*N-1] ]

    S_r_t0.resize(conf_inputs.size());
    S_r_t.resize(conf_inputs.size());

#ifdef _OPENMP
#pragma omp parallel for default(shared) //private()
#endif
    for(int i =0;i<S_r_t0.size();i++){
        S_r_t0[i].resize(6*Parameters_.ns);
        S_r_t[i].resize(6*Parameters_.ns);
    }


    C_tr.resize(time_steps+1);

#ifdef _OPENMP
#pragma omp parallel for default(shared) //private()
#endif
    for(int ts=0;ts<=time_steps;ts++){
        C_tr[ts].resize(6*Parameters_.ns);
        for(int r=0;r<6*Parameters_.ns;r++){
            C_tr[ts][r].resize(6*Parameters_.ns);
            for(int rp=0;rp<6*Parameters_.ns;rp++){
                C_tr[ts][r][rp]=0;
            }
        }
    }


    string line_temp;
    double row_time, double_temp;
    int row_ts;

    ostringstream ostr_w_conv;
    ostr_w_conv << w_conv;
    string string_w_conv = ostr_w_conv.str();
    string SpaceTimeDisplaced_Crt_r0_file = "SpaceTimeDisplaced_Crt_r0_w_conv" + string_w_conv + ".txt";

    ofstream STdisplaced_Crt_r0_Fileout_(SpaceTimeDisplaced_Crt_r0_file.c_str());
    STdisplaced_Crt_r0_Fileout_<<"#time   C(r=0,t)_Sz   C(r=0,t)_Sx  C(r=0,t)_Sy  C(r=0,t)_sz   C(r=0,t)_sx  C(r=0,t)_sy"<<endl;


    ofstream STdisplaced_Crt_Fileout_(STdisplaced_Crt_fileout.c_str());
    STdisplaced_Crt_Fileout_<<"#time   C(rt)_Sz   C(rt)_Sx  C(rt)_Sy  C(rt)_sz   C(rt)_sx  C(rt)_sy"<<endl;




#ifdef _OPENMP
    no_threads_used = min(no_of_processors, No_Of_Inputs);
    omp_set_num_threads(no_threads_used);
    N_p = omp_get_max_threads();
    cout<<"threads being used parallely = "<<N_p<<endl;
    cout<<"No. of threads you asked for = "<<no_of_processors<<endl;
#pragma omp parallel for default(shared) private(line_temp, row_time, row_ts, double_temp)
#endif
    for(int conf=0;conf<No_Of_Inputs;conf++){

        ifstream Specific_conf_in;
        Specific_conf_in.open(conf_inputs[conf].c_str());
        getline(Specific_conf_in, line_temp);

        for(int ts=0;ts<time_steps;ts++){

            getline(Specific_conf_in, line_temp);
            stringstream line_temp_ss(line_temp, stringstream::in);

            line_temp_ss>>row_time;
            row_ts = (int) ((row_time)/dt_ + 0.5);
            if(row_ts != ts){
                cout<<endl;
                cout<<"Error;"<<endl;
                cout<< row_time<<"   "<<row_ts <<"  "<<ts<<endl;
                cout<< "Problem while reading from "<<conf_inputs[conf]<<endl<<endl;
                assert(false);
            }

            if(ts==0){
                for(int r=0;r<6*Parameters_.ns;r++){
                    line_temp_ss>>double_temp;
                    S_r_t0[conf][r] = double_temp;
                }
                for(int r=0;r<Parameters_.ns;r++){
                    for(int rp=0;rp<Parameters_.ns;rp++){

                        //Classical
                        C_tr[ts][r][rp] += (S_r_t0[conf][r]*S_r_t0[conf][rp]);
                        C_tr[ts][r+(Parameters_.ns)][rp+(Parameters_.ns)] +=
                                (S_r_t0[conf][r+Parameters_.ns]*S_r_t0[conf][rp+Parameters_.ns]);
                        C_tr[ts][r+(2*Parameters_.ns)][rp+(2*Parameters_.ns)] +=
                                (S_r_t0[conf][r+(2*Parameters_.ns)]*S_r_t0[conf][rp+(2*Parameters_.ns)]);


                        //Quantum
                        C_tr[ts][r+(3*Parameters_.ns)][rp+(3*Parameters_.ns)] +=
                                (S_r_t0[conf][r+(3*Parameters_.ns)]*S_r_t0[conf][rp+(3*Parameters_.ns)]);
                        C_tr[ts][r+(4*Parameters_.ns)][rp+(4*Parameters_.ns)] +=
                                (S_r_t0[conf][r+(4*Parameters_.ns)]*S_r_t0[conf][rp+(4*Parameters_.ns)]);
                        C_tr[ts][r+(5*Parameters_.ns)][rp+(5*Parameters_.ns)] +=
                                (S_r_t0[conf][r+(5*Parameters_.ns)]*S_r_t0[conf][rp+(5*Parameters_.ns)]);

                    }
                }
            }
            else{
                for(int r=0;r<6*Parameters_.ns;r++){
                    line_temp_ss>>double_temp;
                    S_r_t[conf][r] = double_temp;
                }
                for(int r=0;r<Parameters_.ns;r++){
                    for(int rp=0;rp<Parameters_.ns;rp++){

                        //Classical
                        C_tr[ts][r][rp] += (S_r_t[conf][r]*S_r_t0[conf][rp]);
                        C_tr[ts][r+(Parameters_.ns)][rp+(Parameters_.ns)] +=
                                (S_r_t[conf][r+Parameters_.ns]*S_r_t0[conf][rp+Parameters_.ns]);
                        C_tr[ts][r+(2*Parameters_.ns)][rp+(2*Parameters_.ns)] +=
                                (S_r_t[conf][r+(2*Parameters_.ns)]*S_r_t0[conf][rp+(2*Parameters_.ns)]);


                        //Quantum
                        C_tr[ts][r+(3*Parameters_.ns)][rp+(3*Parameters_.ns)] +=
                                (S_r_t[conf][r+(3*Parameters_.ns)]*S_r_t0[conf][rp+(3*Parameters_.ns)]);
                        C_tr[ts][r+(4*Parameters_.ns)][rp+(4*Parameters_.ns)] +=
                                (S_r_t[conf][r+(4*Parameters_.ns)]*S_r_t0[conf][rp+(4*Parameters_.ns)]);
                        C_tr[ts][r+(5*Parameters_.ns)][rp+(5*Parameters_.ns)] +=
                                (S_r_t[conf][r+(5*Parameters_.ns)]*S_r_t0[conf][rp+(5*Parameters_.ns)]);
                    }
                }
            }
        }
    }


        cout <<"<S_tr . S_t=0,rp> is done"<<endl;


#ifdef _OPENMP
    no_threads_used = min(no_of_processors, time_steps);
    omp_set_num_threads(no_threads_used);
    N_p = omp_get_max_threads();
    cout<<"threads being used parallely = "<<N_p<<endl;
    cout<<"No. of threads you asked for = "<<no_of_processors<<endl;
#pragma omp parallel for default(shared) //private()
#endif
    for(int ts=0;ts<time_steps;ts++){
        for(int r=0;r<Parameters_.ns;r++){
            for(int rp=0;rp<Parameters_.ns;rp++){
                C_tr[ts][r][rp] = (C_tr[ts][r][rp]/No_Of_Inputs)  - ( S_tr[ts][r]*(S_tr[0][rp]) );
                C_tr[ts][r+(Parameters_.ns)][rp+(Parameters_.ns)] =
                        (C_tr[ts][r+(Parameters_.ns)][rp+(Parameters_.ns)]/No_Of_Inputs) -
                        ( S_tr[ts][r+(Parameters_.ns)]*(S_tr[0][rp+(Parameters_.ns)]) );
                C_tr[ts][r+(2*Parameters_.ns)][rp+(2*Parameters_.ns)] =
                        (C_tr[ts][r+(2*Parameters_.ns)][rp+(2*Parameters_.ns)]/No_Of_Inputs) -
                        ( S_tr[ts][r+(2*Parameters_.ns)]*(S_tr[0][rp+(2*Parameters_.ns)]) );

                C_tr[ts][r+(3*Parameters_.ns)][rp+(3*Parameters_.ns)] =
                        (C_tr[ts][r+(3*Parameters_.ns)][rp+(3*Parameters_.ns)]/No_Of_Inputs) -
                        ( S_tr[ts][r+(3*Parameters_.ns)]*(S_tr[0][rp+(3*Parameters_.ns)]) );
                C_tr[ts][r+(4*Parameters_.ns)][rp+(4*Parameters_.ns)] =
                        (C_tr[ts][r+(4*Parameters_.ns)][rp+(4*Parameters_.ns)]/No_Of_Inputs) -
                        ( S_tr[ts][r+(4*Parameters_.ns)]*(S_tr[0][rp+(4*Parameters_.ns)]) );
                C_tr[ts][r+(5*Parameters_.ns)][rp+(5*Parameters_.ns)] =
                        (C_tr[ts][r+(5*Parameters_.ns)][rp+(5*Parameters_.ns)]/No_Of_Inputs) -
                        ( S_tr[ts][r+(5*Parameters_.ns)]*(S_tr[0][rp+(5*Parameters_.ns)]) );

            }
        }
    }


    for(int ts=0;ts<time_steps;ts++){
        STdisplaced_Crt_Fileout_<<ts*dt_<<"  ";
        STdisplaced_Crt_r0_Fileout_<<ts*dt_<<"  ";
        for(int r=0;r<Parameters_.ns;r++){
            for(int rp=0;rp<Parameters_.ns;rp++){

                STdisplaced_Crt_Fileout_<< C_tr[ts][r][rp]*exp(-0.5*(ts*dt_*w_conv*ts*dt_*w_conv)) <<"  ";
                if(r==rp && r==0){
                    STdisplaced_Crt_r0_Fileout_<< C_tr[ts][r][rp]*exp(-0.5*(ts*dt_*w_conv*ts*dt_*w_conv));
                }

            }
        }
        STdisplaced_Crt_Fileout_<<endl;
        STdisplaced_Crt_r0_Fileout_<<endl;
    }


#ifdef _OPENMP
    end_time = omp_get_wtime();
    cout<<"Time to Calculate <C(r,t)> [actual time, using OMP]: "<<double(end_time - begin_time)<<endl;//cout<<"here"<<endl;

#endif

    cout<<"Time to Calculate <C(r,t)> : "<<double( clock() - oprt_SB_time ) / (double)CLOCKS_PER_SEC<<endl;//cout<<"here"<<endl;


}


void ST_Fourier_1orb_MCMF::Calculate_Skw_from_Crt(string fileout){



    int N_p;
    int no_threads_used;
    Fft DO_Fft;
    DO_Fft.PI_EFF=PI;


    complex<double> iota(0,1);
    S_rw.resize(Parameters_.ns);
    s_quantum_rw.resize(Parameters_.ns);

    for(int pos_i=0;pos_i<Parameters_.ns;pos_i++){
        S_rw[pos_i].resize(Parameters_.ns);
        s_quantum_rw[pos_i].resize(Parameters_.ns);

        for(int pos_j=0;pos_j<Parameters_.ns;pos_j++){
            S_rw[pos_i][pos_j].resize(n_wpoints);
            s_quantum_rw[pos_i][pos_j].resize(n_wpoints);

            for(int wi=0;wi<n_wpoints;wi++){
                S_rw[pos_i][pos_j][wi]=0;
                s_quantum_rw[pos_i][pos_j][wi]=0;

            }
        }
    }

    Sz_eq.resize(Parameters_.ns);Sx_eq.resize(Parameters_.ns);Sy_eq.resize(Parameters_.ns);
    sz_eq.resize(Parameters_.ns);sx_eq.resize(Parameters_.ns);sy_eq.resize(Parameters_.ns);

    Sz_t.resize(Parameters_.ns);Sx_t.resize(Parameters_.ns);Sy_t.resize(Parameters_.ns);
    sz_t.resize(Parameters_.ns);sx_t.resize(Parameters_.ns);sy_t.resize(Parameters_.ns);


    string line_temp;
    double temp_waste;


    //cout<<line_temp<<endl;


    double begin_time, end_time;
    clock_t oprt_SB_time;

#ifdef _OPENMP
    begin_time = omp_get_wtime();
#endif

    oprt_SB_time = clock();


    if(!Use_FFT){
        for(int ts=0;ts<time_steps;ts++){

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
            for(int wi=0;wi<n_wpoints;wi++){

                for(int pos_i=0;pos_i<Parameters_.ns;pos_i++){
                    for(int pos_j=0;pos_j<Parameters_.ns;pos_j++){


                        //exp(iota*(-wi * dw) * (ts* dt_))
                        S_rw[pos_i][pos_j][wi] += exp(iota*(-wi * dw) * (ts* dt_))*exp(-0.5*(ts*dt_*w_conv*ts*dt_*w_conv))*dt_*(
                                    ( ( C_tr[ts][pos_i][pos_j]  )  )
                                    +
                                    ( ( C_tr[ts][pos_i+(Parameters_.ns)][pos_j+(Parameters_.ns)]  )  )
                                +
                                ( ( C_tr[ts][pos_i+(2*Parameters_.ns)][pos_j+(2*Parameters_.ns)] )  )
                                );


                        s_quantum_rw[pos_i][pos_j][wi] += exp(iota*(-wi * dw) * (ts* dt_))*exp(-0.5*(ts*dt_*w_conv*ts*dt_*w_conv))*dt_*(
                                    ( (( C_tr[ts][pos_i+(3*Parameters_.ns)][pos_j+(3*Parameters_.ns)] )  )  )
                                +
                                ( (( C_tr[ts][pos_i+(4*Parameters_.ns)][pos_j+(4*Parameters_.ns)] )  )  )
                                +
                                ( (( C_tr[ts][pos_i+(5*Parameters_.ns)][pos_j+(5*Parameters_.ns)] )  )  )
                                );




                        /*    T_rw[pos_i][pos_j][wi] +=exp(iota * (-wi * dw) * (ts* dt_))*exp(-0.5*(ts*dt_*w_conv*ts*dt_*w_conv))*dt_* (


                            );*/






                    }
                }
            }

        }
    }

    else{
        //USING FFT



#ifdef _OPENMP
        no_threads_used = min(no_of_processors, Parameters_.ns);
        omp_set_num_threads(no_threads_used);
        N_p = omp_get_max_threads();
        cout<<"threads being used parallely = "<<N_p<<endl;
        cout<<"No. of threads you asked for = "<<no_of_processors<<endl;
#pragma omp parallel for default(shared) //private()
#endif
        for(int pos_i=0;pos_i<Parameters_.ns;pos_i++){

            Mat_1_Complex_doub Vec_1, Vec_2;
            Vec_1.resize(time_steps);Vec_2.resize(time_steps);

            for(int pos_j=0;pos_j<Parameters_.ns;pos_j++){
                for(int ts=0;ts<time_steps;ts++){
                    Vec_1[ts] = exp(-0.5*(ts*dt_*w_conv*ts*dt_*w_conv))*dt_*(
                                ( ( C_tr[ts][pos_i][pos_j]  )  )
                                +
                                ( ( C_tr[ts][pos_i+(Parameters_.ns)][pos_j+(Parameters_.ns)]  )  )
                            +
                            ( ( C_tr[ts][pos_i+(2*Parameters_.ns)][pos_j+(2*Parameters_.ns)] )  )
                            );
                    Vec_2[ts] = exp(-0.5*(ts*dt_*w_conv*ts*dt_*w_conv))*dt_*(
                                ( ( C_tr[ts][pos_i+(3*Parameters_.ns)][pos_j+(3*Parameters_.ns)]  )  )
                            +
                            ( ( C_tr[ts][pos_i+(4*Parameters_.ns)][pos_j+(4*Parameters_.ns)]  )  )
                            +
                            ( ( C_tr[ts][pos_i+(5*Parameters_.ns)][pos_j+(5*Parameters_.ns)] )  )
                            );
                }

                DO_Fft.transform(Vec_1);
                DO_Fft.transform(Vec_2);


                for(int wi=0;wi<n_wpoints;wi++){
                    S_rw[pos_i][pos_j][wi] = Vec_1[wi];
                    s_quantum_rw[pos_i][pos_j][wi] =Vec_2[wi];
                }
            }

            vector< complex<double> >().swap( Vec_1 );
            vector< complex<double> >().swap( Vec_1 );
        }

    }



    //Now "rw - space" to "kw - space"
    S_kw.resize(Parameters_.lx);
    s_quantum_kw.resize(Parameters_.lx);

    for(int x_i=0;x_i<Parameters_.lx;x_i++){
        S_kw[x_i].resize(Parameters_.ly);
        s_quantum_kw[x_i].resize(Parameters_.ly);
        for(int y_j=0;y_j<Parameters_.ly;y_j++){
            S_kw[x_i][y_j].resize(n_wpoints);
            s_quantum_kw[x_i][y_j].resize(n_wpoints);
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
                complex<double> temp2(0,0);
                complex<double> temp(0,0);
                for(int x_i=0;x_i<Parameters_.lx;x_i++){
                    for(int y_i=0;y_i<Parameters_.ly;y_i++){

                        pos_i = y_i*Parameters_.lx + x_i;

                        for(int x_j=0;x_j<Parameters_.lx;x_j++){
                            for(int y_j=0;y_j<Parameters_.ly;y_j++){

                                pos_j = y_j*Parameters_.lx + x_j;

                                //expd(iota*( (x_j - x_i)*kx +  (y_j - y_i)*ky ) )
                                temp += S_rw[pos_i][pos_j][wi]*exp(iota*( (x_j - x_i)*kx +  (y_j - y_i)*ky ) );
                                //temp += S_rw[pos_i][pos_j][wi]*cos(( (x_j - x_i)*kx +  (y_j - y_i)*ky ) );
                                temp2 += s_quantum_rw[pos_i][pos_j][wi]*exp(iota*( (x_j - x_i)*kx +  (y_j - y_i)*ky ) );


                            }
                        }
                    }
                }
                //file_out2<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<temp.real()<<"   "<<temp.imag()<<"    "<<temp2.real()<<"   "<<temp2.imag()<<"    "<<temp3.real()<<"   "<<temp3.imag()<<endl;
                //file_out_full<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<wi<<"   "<<temp.real()<<"   "<<temp.imag()<<"    "<<"none"<<"   "<<"none"<<"    "<<"none"<<"   "<<"none"<<endl;
                S_kw[nx][ny][wi]=temp;
                s_quantum_kw[nx][ny][wi]=temp2;
            }



        }

    }



#ifdef _OPENMP
    end_time = omp_get_wtime();
    cout<<"Time to Calculate <SS(kx,ky,w)> [actual time, using OMP]: "<<double(end_time - begin_time)<<endl;//cout<<"here"<<endl;

#endif

    cout<<"Time to Calculate <SS(kx,ky,w)> : "<<double( clock() - oprt_SB_time ) / (double)CLOCKS_PER_SEC<<endl;//cout<<"here"<<endl;




    ofstream file_out_full(fileout.c_str());

    file_out_full<<"#nx   ny   k_ind   wi*dw   wi   S_kw[nx][ny][wi]    s_quantum_kw[nx][ny][wi]          T_kw[nx][ny][wi]"<<endl;
    int k_ind=0;
    for(int nx=0;nx<Parameters_.lx;nx++){
        for(int ny=0;ny<Parameters_.ly;ny++){

            for(int wi=0;wi<n_wpoints;wi++){

                //file_out2<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<temp.real()<<"   "<<temp.imag()<<"    "<<temp2.real()<<"   "<<temp2.imag()<<"    "<<temp3.real()<<"   "<<temp3.imag()<<endl;
                file_out_full<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<wi<<"   "<<S_kw[nx][ny][wi].real()<<"   "<<S_kw[nx][ny][wi].imag()<<"    "<<s_quantum_kw[nx][ny][wi].real()<<
                               "    "<<s_quantum_kw[nx][ny][wi].imag()<<endl;
                //Skw_Mat[nx][ny][wi]=temp;
            }

            k_ind +=1;
            file_out_full<<endl;

        }

    }



}

void ST_Fourier_1orb_MCMF::Calculate_Skw_from_Srt_file( string filename, string fileout){

    complex<double> iota(0,1);
    S_rw.resize(Parameters_.ns);
    s_quantum_rw.resize(Parameters_.ns);

    for(int pos_i=0;pos_i<Parameters_.ns;pos_i++){
        S_rw[pos_i].resize(Parameters_.ns);
        s_quantum_rw[pos_i].resize(Parameters_.ns);

        for(int pos_j=0;pos_j<Parameters_.ns;pos_j++){
            S_rw[pos_i][pos_j].resize(n_wpoints);
            s_quantum_rw[pos_i][pos_j].resize(n_wpoints);

            for(int wi=0;wi<n_wpoints;wi++){
                S_rw[pos_i][pos_j][wi]=0;
                s_quantum_rw[pos_i][pos_j][wi]=0;

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


                    //exp(iota*(-wi * dw) * (ts* dt_))
                    S_rw[pos_i][pos_j][wi] += cos((wi * dw) * (ts* dt_))*exp(-0.5*(ts*dt_*w_conv*ts*dt_*w_conv))*dt_*(
                                ( (Sz_eq[pos_i]*Sz_t[pos_j]  )  )
                                +
                                ( (Sx_eq[pos_i]*Sx_t[pos_j]  )  )
                                +
                                ( (Sy_eq[pos_i]*Sy_t[pos_j]  )  )
                                );


                    s_quantum_rw[pos_i][pos_j][wi] += cos((wi * dw) * (ts* dt_))*exp(-0.5*(ts*dt_*w_conv*ts*dt_*w_conv))*dt_*(
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


    for(int x_i=0;x_i<Parameters_.lx;x_i++){

        S_kw[x_i].resize(Parameters_.ly);
        s_quantum_kw[x_i].resize(Parameters_.ly);

        for(int y_j=0;y_j<Parameters_.ly;y_j++){

            S_kw[x_i][y_j].resize(n_wpoints);
            s_quantum_kw[x_i][y_j].resize(n_wpoints);

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
                complex<double> temp3(0,0);
                complex<double> temp2(0,0);
                complex<double> temp(0,0);
                for(int x_i=0;x_i<Parameters_.lx;x_i++){
                    for(int y_i=0;y_i<Parameters_.ly;y_i++){

                        pos_i = y_i*Parameters_.lx + x_i;

                        for(int x_j=0;x_j<Parameters_.lx;x_j++){
                            for(int y_j=0;y_j<Parameters_.ly;y_j++){

                                pos_j = y_j*Parameters_.lx + x_j;

                                //expd(iota*( (x_j - x_i)*kx +  (y_j - y_i)*ky ) )
                                temp += S_rw[pos_i][pos_j][wi]*cos(( (x_j - x_i)*kx +  (y_j - y_i)*ky ) );
                                //temp += S_rw[pos_i][pos_j][wi]*cos(( (x_j - x_i)*kx +  (y_j - y_i)*ky ) );
                                temp2 += s_quantum_rw[pos_i][pos_j][wi]*cos(( (x_j - x_i)*kx +  (y_j - y_i)*ky ) );


                            }
                        }
                    }
                }
                //file_out2<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<temp.real()<<"   "<<temp.imag()<<"    "<<temp2.real()<<"   "<<temp2.imag()<<"    "<<temp3.real()<<"   "<<temp3.imag()<<endl;
                //file_out_full<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<wi<<"   "<<temp.real()<<"   "<<temp.imag()<<"    "<<"none"<<"   "<<"none"<<"    "<<"none"<<"   "<<"none"<<endl;
                S_kw[nx][ny][wi]=temp;
                s_quantum_kw[nx][ny][wi]=temp2;
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
                file_out_full<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<wi<<"   "<<S_kw[nx][ny][wi].real()<<"   "<<S_kw[nx][ny][wi].imag()<<"    "<<s_quantum_kw[nx][ny][wi].real()<<
                               "    "<<s_quantum_kw[nx][ny][wi].imag()<<endl;
                //Skw_Mat[nx][ny][wi]=temp;
            }

            k_ind +=1;
            file_out_full<<endl;

        }

    }








}

void ST_Fourier_1orb_MCMF::Read_parameters(){


    w_min=SC_SW_ENGINE_VNE_1orb_MCMF_.w_min;
    w_max=SC_SW_ENGINE_VNE_1orb_MCMF_.w_max;

    dw=SC_SW_ENGINE_VNE_1orb_MCMF_.dw;
    w_conv=SC_SW_ENGINE_VNE_1orb_MCMF_.w_conv;
    Use_FFT=SC_SW_ENGINE_VNE_1orb_MCMF_.Use_FFT;

    no_of_processors=SC_SW_ENGINE_VNE_1orb_MCMF_.no_of_processors;


}


























#endif


