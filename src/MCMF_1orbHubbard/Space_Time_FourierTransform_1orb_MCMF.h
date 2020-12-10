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

    Mat_1_string conf_inputs, F_wq_inputs, Aq_inputs;
    int No_Of_Inputs;


    string Skw_out;
    string Skw_out_full;

    bool Space_Fourier_using_single_S;
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
    Mat_2_doub C_tr_; //Space-time dispaced correlation, avg over Ensemble
    Mat_2_doub C_Quantum_tr_;
    Mat_2_doub C_Classical_tr_;
    Mat_3_doub C_tr; //Space-time dispaced correlation, avg over Ensemble
    Mat_3_doub C_Quantum_tr;
    Mat_3_doub C_Classical_tr;

    Mat_2_Complex_doub F_rw, F_qw, S_qw, D_qw, S_qw_conv;
    Mat_3_Complex_doub F_qw_;
    Mat_1_Complex_doub Aq, delAq, Aq_avg, Aq2_avg;


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

    void Calculate_Skw_from_Crt_(string fileout);
    void Calculate_SpaceTimeDisplacedCorrelations_Smarter(string STdisplaced_Crt_fileout);

    void Calculate_Fw_and_Aq(string fileout, string fileout_Aq);

    void Calculate_Sqw_using_Aq_Fwq(string Sqw_file, string Dqw_file);

    void Calculate_Sqw_using_Fwq(string Sqw_file, string Dqw_file);
    void Calculate_Sqw_using_Fwq_with_negativeandpositive_Time(string fileout, string Dqw_file);
    void Sum_Rule_For_Way1(string Sqw_file_in, string Sq_static_out, string Sq_dynamical_out);

    void Convolute_the_spectrum(string Sqw_file_in, string Sqw_file_out);


};




void ST_Fourier_1orb_MCMF::Initialize_engine(){



    ifstream check_file;
    check_file.open(conf_inputs[0].c_str());
    string line_temp, line_last;


    double t0, t1;
    //    for(int curLine = 0; getline(check_file, line_temp); curLine++) {

    //        if(curLine==1){
    //            stringstream line_temp_ss(line_temp, stringstream::in);
    //            line_temp_ss>>t0;
    //        }
    //        if(curLine==2){
    //            stringstream line_temp_ss(line_temp, stringstream::in);
    //            line_temp_ss>>t1;
    //        }
    //        line_last=line_temp;

    //        if(curLine%100==0){
    //            cout<<"test read of "<<conf_inputs[0]<<", line no = "<<curLine<<endl;
    //        }
    //    }


    getline(check_file, line_temp); //line=0

    getline(check_file, line_temp);  //line=1
    stringstream line_temp_ss0(line_temp, stringstream::in);
    line_temp_ss0>>t0;

    getline(check_file, line_temp);  //line=2
    stringstream line_temp_ss1(line_temp, stringstream::in);
    line_temp_ss1>>t1;

    if( (dt_  - (t1 -t0)) > 10e-10  ){
        cout <<"dt from file and input does not match"<<endl;
        assert( (dt_  - (t1 -t0)) < 10e-10 );
    }

    //dt_ = t1 - t0;


    stringstream line_temp_ss(line_last, stringstream::in);
    line_temp_ss>>time_max;

    //cout<<line_temp<<endl;
    time_steps=(int) ((time_max - t0)/(fabs(dt_)) + 0.5);

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
        cout<<"FFT is used, so dw, w_min are fixed."<<endl;
        dw = (2.0*PI)/ time_max;
        w_min = 0.0;
        n_wpoints=(int) ((w_max-w_min)/dw + 0.5);
        //n_wpoints = time_steps/2;
        //w_max = w_min + n_wpoints*(dw);
        cout<<"dw = 2*pi/T_max = "<< dw<<endl;
        cout<<"w_min = "<< w_min<<endl;
        cout<<"w_max = "<< w_max<<endl;
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


    Space_Fourier_using_single_S=false;


    cout<<"ST_Fourier engine initialized"<<endl;


}
void ST_Fourier_1orb_MCMF::Perform_Averaging_on_one_point(){

    string Temp_file_Srt="Average_Srt.txt";
    ofstream Temp_file_Srt_out(Temp_file_Srt.c_str());


    Temp_file_Srt_out<<"#time   Sz   Sx  Sy"<<endl;

    string line_temp;
    double row_time, double_temp;
    int row_ts;

    Mat_1_doub Sr;
    Sr.resize(3*Parameters_.ns);


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

            for(int r=0;r<3*Parameters_.ns;r++){
                line_temp_ss>>double_temp;

                Sr[r] += double_temp;
            }

            Specific_conf_in.clear();
            Specific_conf_in.seekg(0, ios::beg);


        }

        for(int r=0;r<3*Parameters_.ns;r++){
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
    //#ifdef _OPENMP
    //    N_p = omp_get_max_threads();
    //    cout<<"Maximum threads can be used parallely = "<<N_p<<endl;
    //    cout<<"No. of threads you asked for = "<<no_of_processors<<endl;

    //    begin_time = omp_get_wtime();
    //    omp_set_num_threads(1);
    //#endif


    ostringstream ostr_w_conv;
    ostr_w_conv << w_conv;
    string string_w_conv = ostr_w_conv.str();

    string Temp_file_Srt="Average_Srt_w_conv"+string_w_conv+".txt";
    ofstream Temp_file_Srt_out(Temp_file_Srt.c_str());


    Temp_file_Srt_out<<"#time   Sz   Sx  Sy  sz  sx  sy"<<endl;


    S_tr.resize(time_steps+1);


    //#ifdef _OPENMP
    //#pragma omp parallel for default(shared) //private()
    //#endif
    for(int ts=0;ts<=time_steps;ts++){
        S_tr[ts].resize(3*Parameters_.ns);
        for(int r=0;r<3*Parameters_.ns;r++){
            S_tr[ts][r]=0;
        }
    }


    //#ifdef _OPENMP
    //    no_threads_used = min(no_of_processors, No_Of_Inputs);
    //    omp_set_num_threads(no_threads_used);
    //    N_p = omp_get_max_threads();
    //    cout<<"threads being used parallely = "<<N_p<<endl;
    //    cout<<"No. of threads you asked for = "<<no_of_processors<<endl;
    //#pragma omp parallel for default(shared) //private()
    //#endif
    for(int conf=0;conf<No_Of_Inputs;conf++){

        //#ifdef _OPENMP
        //        //cout<< "using thread no "<<omp_get_thread_num()<<" out of total "<<omp_get_num_threads()<<"threads"<<endl;
        //#endif

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

            for(int r=0;r<3*Parameters_.ns;r++){
                line_temp_ss>>double_temp;
                S_tr[ts][r]+= double_temp;
            }
        }
    }



    for(int ts=0;ts<=time_steps;ts++){
        Temp_file_Srt_out<<ts*dt_<<"  ";
        for(int r=0;r<3*Parameters_.ns;r++){
            S_tr[ts][r] = S_tr[ts][r]/(No_Of_Inputs);
            Temp_file_Srt_out<<S_tr[ts][r]<<"  ";
        }
        Temp_file_Srt_out<<endl;
    }


    //#ifdef _OPENMP
    //    end_time = omp_get_wtime();
    //    cout<<"Time to Calculate <S(r,t)>, <s(r,t)> [actual time, using OMP]: "<<double(end_time - begin_time)<<endl;//cout<<"here"<<endl;

    //#endif

    cout<<"Time to Calculate <S(r,t)>, <s(r,t)> : "<<double( clock() - oprt_SB_time ) / (double)CLOCKS_PER_SEC<<endl;//cout<<"here"<<endl;





}


void ST_Fourier_1orb_MCMF::Calculate_SpaceTimeDisplacedCorrelations(string STdisplaced_Crt_fileout){

    //TO DO: Parallelize this routine as well, later on : 16th March-2020

    int N_p;
    int no_threads_used;
    double begin_time, end_time;
    clock_t oprt_SB_time = clock();


    //#ifdef _OPENMP
    //    omp_set_num_threads(1);
    //    begin_time = omp_get_wtime();
    //#endif


    Mat_2_doub S_r_t0; //[Conf_no][r \in [0,6*N-1] ]
    Mat_2_doub S_r_t; //[Conf_no][r \in [0,6*N-1] ]

    S_r_t0.resize(conf_inputs.size());
    S_r_t.resize(conf_inputs.size());

    //#ifdef _OPENMP
    //#pragma omp parallel for default(shared) //private()
    //#endif
    for(int i =0;i<S_r_t0.size();i++){
        S_r_t0[i].resize(3*Parameters_.ns);
        S_r_t[i].resize(3*Parameters_.ns);
    }


    C_Quantum_tr.resize(time_steps+1);
    C_Classical_tr.resize(time_steps+1);

    //#ifdef _OPENMP
    //#pragma omp parallel for default(shared) //private()
    //#endif
    for(int ts=0;ts<=time_steps;ts++){
        C_Quantum_tr[ts].resize(Parameters_.ns);
        C_Classical_tr[ts].resize(Parameters_.ns);
        for(int r=0;r<Parameters_.ns;r++){
            C_Quantum_tr[ts][r].resize(Parameters_.ns);
            C_Classical_tr[ts][r].resize(Parameters_.ns);
            for(int rp=0;rp<Parameters_.ns;rp++){
                C_Quantum_tr[ts][r][rp]=0;
                C_Classical_tr[ts][r][rp]=0;
            }
        }
    }



    ostringstream ostr_w_conv;
    ostr_w_conv << w_conv;
    string string_w_conv = ostr_w_conv.str();
    string SpaceTimeDisplaced_Crt_r0_file = "SpaceTimeDisplaced_Crt_r0_w_conv" + string_w_conv + ".txt";

    ofstream STdisplaced_Crt_r0_Fileout_(SpaceTimeDisplaced_Crt_r0_file.c_str());
    STdisplaced_Crt_r0_Fileout_<<"#time   C(r=0,t)_S  C(r=0,t)_s"<<endl;


    ofstream STdisplaced_Crt_Fileout_(STdisplaced_Crt_fileout.c_str());
    STdisplaced_Crt_Fileout_<<"#time   C(rt)_S   C(rt)_s"<<endl;


    //This parallelization gives noise, why??
    //#ifdef _OPENMP
    //    no_threads_used = min(no_of_processors, No_Of_Inputs);
    //    omp_set_num_threads(no_threads_used);
    //    N_p = omp_get_max_threads();
    //    cout<<"threads being used parallely = "<<N_p<<endl;
    //    cout<<"No. of threads you asked for = "<<no_of_processors<<endl;
    //#pragma omp parallel for default(shared) //private()
    //#endif
    for(int conf=0;conf<No_Of_Inputs;conf++){

        string line_temp2;

        ifstream Specific_conf_in;
        Specific_conf_in.open(conf_inputs[conf].c_str());
        getline(Specific_conf_in, line_temp2);

        int ts=0;
        string line_temp;
        double row_time, double_temp;
        int row_ts;
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

        for(int r=0;r<3*Parameters_.ns;r++){
            line_temp_ss>>double_temp;
            S_r_t0[conf][r] = double_temp;
        }

        for(int r=0;r<Parameters_.ns;r++){
            for(int rp=0;rp<Parameters_.ns;rp++){

                //Classical
                C_Classical_tr[ts][r][rp] +=
                        (S_r_t0[conf][r]*S_r_t0[conf][rp]) +
                        (S_r_t0[conf][r+Parameters_.ns]*S_r_t0[conf][rp+Parameters_.ns]) +
                        (S_r_t0[conf][r+(2*Parameters_.ns)]*S_r_t0[conf][rp+(2*Parameters_.ns)]);


            }
        }
    }


    //-----------------------------REPEAT for ts!=0---------------//
    //#ifdef _OPENMP
    //    no_threads_used = min(no_of_processors, No_Of_Inputs);
    //    omp_set_num_threads(no_threads_used);
    //    N_p = omp_get_max_threads();
    //#pragma omp parallel for default(shared) //private()
    //#endif
    for(int conf=0;conf<No_Of_Inputs;conf++){

        string line_temp2;

        ifstream Specific_conf_in;
        Specific_conf_in.open(conf_inputs[conf].c_str());
        getline(Specific_conf_in, line_temp2);
        getline(Specific_conf_in, line_temp2);

        for(int ts=1;ts<time_steps;ts++){
            string line_temp;
            double row_time, double_temp;
            int row_ts;
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

            for(int r=0;r<3*Parameters_.ns;r++){
                line_temp_ss>>double_temp;
                S_r_t[conf][r] = double_temp;
            }

            for(int r=0;r<Parameters_.ns;r++){
                for(int rp=0;rp<Parameters_.ns;rp++){

                    //Classical
                    C_Classical_tr[ts][r][rp] +=
                            (S_r_t[conf][r]*S_r_t0[conf][rp]) +
                            (S_r_t[conf][r+Parameters_.ns]*S_r_t0[conf][rp+Parameters_.ns]) +
                            (S_r_t[conf][r+(2*Parameters_.ns)]*S_r_t0[conf][rp+(2*Parameters_.ns)]);

                }
            }
        }

    }

    //-------------------------------------------------------------//


    cout <<"<S_tr . S_t=0,rp> is done"<<endl;


    for(int ts=0;ts<time_steps;ts++){
        for(int r=0;r<Parameters_.ns;r++){
            for(int rp=0;rp<Parameters_.ns;rp++){
                C_Classical_tr[ts][r][rp] =
                        (C_Classical_tr[ts][r][rp]/No_Of_Inputs)
                        - ( S_tr[ts][r]*(S_tr[0][rp]) )
                        - ( S_tr[ts][r+(Parameters_.ns)]*(S_tr[0][rp+(Parameters_.ns)]) )
                        - ( S_tr[ts][r+(2*Parameters_.ns)]*(S_tr[0][rp+(2*Parameters_.ns)]) );

            }
        }
    }


    for(int ts=0;ts<time_steps;ts++){
        STdisplaced_Crt_Fileout_<<ts*dt_<<"  ";
        STdisplaced_Crt_r0_Fileout_<<ts*dt_<<"  ";
        for(int r=0;r<Parameters_.ns;r++){
            for(int rp=0;rp<Parameters_.ns;rp++){

                STdisplaced_Crt_Fileout_<< C_Classical_tr[ts][r][rp]*exp(-0.5*(ts*dt_*w_conv*ts*dt_*w_conv)) <<"  ";
                if(r==rp && r==0){
                    STdisplaced_Crt_r0_Fileout_<< C_Classical_tr[ts][r][rp]*exp(-0.5*(ts*dt_*w_conv*ts*dt_*w_conv))<<"   ";
                }

            }
        }
        STdisplaced_Crt_Fileout_<<endl;
        STdisplaced_Crt_r0_Fileout_<<endl;
    }


    //#ifdef _OPENMP
    //    end_time = omp_get_wtime();
    //    cout<<"Time to Calculate <C(r,t)> [actual time, using OMP]: "<<double(end_time - begin_time)<<endl;//cout<<"here"<<endl;

    //#endif

    cout<<"Time to Calculate <C(r,t)> : "<<double( clock() - oprt_SB_time ) / (double)CLOCKS_PER_SEC<<endl;//cout<<"here"<<endl;


}




void ST_Fourier_1orb_MCMF::Calculate_SpaceTimeDisplacedCorrelations_Smarter(string STdisplaced_Crt_fileout){

    clock_t oprt_SB_time = clock();

    int r_new, rx_new, ry_new;
    int rx, ry, rpx, rpy;
    int rp_max;
    rp_max = (Parameters_.ns*int(!Space_Fourier_using_single_S));

    Mat_2_doub S_r_t0; //[Conf_no][r \in [0,6*N-1] ]
    Mat_2_doub S_r_t; //[Conf_no][r \in [0,6*N-1] ]

    S_r_t0.resize(conf_inputs.size());
    S_r_t.resize(conf_inputs.size());

    for(int i =0;i<S_r_t0.size();i++){
        S_r_t0[i].resize(3*Parameters_.ns);
        S_r_t[i].resize(3*Parameters_.ns);
    }

    C_Quantum_tr_.resize(time_steps+1);
    C_Classical_tr_.resize(time_steps+1);

    //#ifdef _OPENMP
    //#pragma omp parallel for default(shared) //private()
    //#endif
    for(int ts=0;ts<=time_steps;ts++){
        C_Quantum_tr_[ts].resize(Parameters_.ns);
        C_Classical_tr_[ts].resize(Parameters_.ns);
    }



    ostringstream ostr_w_conv;
    ostr_w_conv << w_conv;
    string string_w_conv = ostr_w_conv.str();
    string SpaceTimeDisplaced_Crt_r0_file = "SpaceTimeDisplaced_Crt_r0_w_conv" + string_w_conv + ".txt";

    ofstream STdisplaced_Crt_r0_Fileout_(SpaceTimeDisplaced_Crt_r0_file.c_str());
    STdisplaced_Crt_r0_Fileout_<<"#time   C(r=0,t)_S "<<endl;


    ofstream STdisplaced_Crt_Fileout_(STdisplaced_Crt_fileout.c_str());
    STdisplaced_Crt_Fileout_<<"#time   C(rt)_S "<<endl;


    //This parallelization gives noise, why??
    //#ifdef _OPENMP
    //    no_threads_used = min(no_of_processors, No_Of_Inputs);
    //    omp_set_num_threads(no_threads_used);
    //    N_p = omp_get_max_threads();
    //    cout<<"threads being used parallely = "<<N_p<<endl;
    //    cout<<"No. of threads you asked for = "<<no_of_processors<<endl;
    //#pragma omp parallel for default(shared) //private()
    //#endif
    for(int conf=0;conf<No_Of_Inputs;conf++){

        string line_temp2;

        ifstream Specific_conf_in;
        Specific_conf_in.open(conf_inputs[conf].c_str());
        getline(Specific_conf_in, line_temp2);

        int ts=0;
        string line_temp;
        double row_time, double_temp;
        int row_ts;
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

        for(int r=0;r<3*Parameters_.ns;r++){
            line_temp_ss>>double_temp;
            S_r_t0[conf][r] = double_temp;
        }

        for(int r=0;r<Parameters_.ns;r++){
            rx=Coordinates_.indx(r);
            ry=Coordinates_.indy(r);

            for(int rp=0;rp<rp_max;rp++){
                rpx=Coordinates_.indx(rp);
                rpy=Coordinates_.indy(rp);

                rx_new = (rpx+rx)%Parameters_.lx;
                ry_new = (rpy+ry)%Parameters_.ly;
                r_new = Coordinates_.Nc(rx_new, ry_new);

                //Classical
                C_Classical_tr_[ts][r] +=
                        (S_r_t0[conf][rp]*S_r_t0[conf][r_new]) +
                        (S_r_t0[conf][rp+Parameters_.ns]*S_r_t0[conf][r_new+Parameters_.ns]) +
                        (S_r_t0[conf][rp+(2*Parameters_.ns)]*S_r_t0[conf][r_new+(2*Parameters_.ns)]);


            }
        }
    }


    //-----------------------------REPEAT for ts!=0---------------//
    //#ifdef _OPENMP
    //    no_threads_used = min(no_of_processors, No_Of_Inputs);
    //    omp_set_num_threads(no_threads_used);
    //    N_p = omp_get_max_threads();
    //#pragma omp parallel for default(shared) //private()
    //#endif
    for(int conf=0;conf<No_Of_Inputs;conf++){

        string line_temp2;

        ifstream Specific_conf_in;
        Specific_conf_in.open(conf_inputs[conf].c_str());
        getline(Specific_conf_in, line_temp2);
        getline(Specific_conf_in, line_temp2);

        for(int ts=1;ts<time_steps;ts++){
            string line_temp;
            double row_time, double_temp;
            int row_ts;
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

            for(int r=0;r<3*Parameters_.ns;r++){
                line_temp_ss>>double_temp;
                S_r_t[conf][r] = double_temp;
            }

            for(int r=0;r<Parameters_.ns;r++){
                rx=Coordinates_.indx(r);
                ry=Coordinates_.indy(r);

                for(int rp=0;rp<rp_max;rp++){
                    rpx=Coordinates_.indx(rp);
                    rpy=Coordinates_.indy(rp);

                    rx_new = (rpx+rx)%Parameters_.lx;
                    ry_new = (rpy+ry)%Parameters_.ly;
                    r_new = Coordinates_.Nc(rx_new, ry_new);

                    //Classical
                    C_Classical_tr_[ts][r] +=
                            (S_r_t[conf][rp]*S_r_t0[conf][r_new]) +
                            (S_r_t[conf][rp+Parameters_.ns]*S_r_t0[conf][r_new+Parameters_.ns]) +
                            (S_r_t[conf][rp+(2*Parameters_.ns)]*S_r_t0[conf][r_new+(2*Parameters_.ns)]);

                }
            }
        }

    }

    //-------------------------------------------------------------//


    cout <<"<S_tr . S_t=0,rp> is done"<<endl;


    for(int ts=0;ts<time_steps;ts++){

        for(int r=0;r<Parameters_.ns;r++){
            rx=Coordinates_.indx(r);
            ry=Coordinates_.indy(r);

            C_Classical_tr_[ts][r] = (C_Classical_tr_[ts][r]/No_Of_Inputs);
            C_Quantum_tr_[ts][r] = (C_Quantum_tr_[ts][r]/No_Of_Inputs);

            if(No_Of_Inputs>1){
                for(int rp=0;rp<Parameters_.ns;rp++){
                    rpx=Coordinates_.indx(rp);
                    rpy=Coordinates_.indy(rp);

                    rx_new = (rpx+rx)%Parameters_.lx;
                    ry_new = (rpy+ry)%Parameters_.ly;
                    r_new = Coordinates_.Nc(rx_new, ry_new);


                    C_Classical_tr_[ts][r] +=
                            - ( S_tr[ts][rp]*(S_tr[0][r_new]) )
                            - ( S_tr[ts][rp+(Parameters_.ns)]*(S_tr[0][r_new+(Parameters_.ns)]) )
                            - ( S_tr[ts][rp+(2*Parameters_.ns)]*(S_tr[0][r_new+(2*Parameters_.ns)]) );


                }

            }
            C_Classical_tr_[ts][r] = C_Classical_tr_[ts][r]*(1.0/(1.0*Parameters_.ns));
            C_Quantum_tr_[ts][r] = C_Quantum_tr_[ts][r]*(1.0/(1.0*Parameters_.ns));

        }
    }

    if(No_Of_Inputs==1){
        cout<<"Only <S(t,r) . S(t=0,rp)>  is used"<<endl;
    }


    for(int ts=0;ts<time_steps;ts++){
        STdisplaced_Crt_Fileout_<<ts*dt_<<"  ";
        STdisplaced_Crt_r0_Fileout_<<ts*dt_<<"  ";
        for(int r=0;r<Parameters_.ns;r++){

            STdisplaced_Crt_Fileout_<< C_Classical_tr_[ts][r]*exp(-0.5*(ts*dt_*w_conv*ts*dt_*w_conv)) <<"  ";
            if(r==0){
                STdisplaced_Crt_r0_Fileout_<< C_Classical_tr_[ts][r]*exp(-0.5*(ts*dt_*w_conv*ts*dt_*w_conv))<<"   ";
            }
        }
        STdisplaced_Crt_Fileout_<<endl;
        STdisplaced_Crt_r0_Fileout_<<endl;
    }

    cout<<"Time to Calculate <C(r,t)> : "<<double( clock() - oprt_SB_time ) / (double)CLOCKS_PER_SEC<<endl;//cout<<"here"<<endl;


}



void ST_Fourier_1orb_MCMF::Calculate_Skw_from_Crt(string fileout){



    int N_p;
    int no_threads_used;
    no_threads_used = min(no_of_processors, No_Of_Inputs);
#ifdef _OPENMP
    //            no_threads_used = min(no_of_processors, No_Of_Inputs);
    //            omp_set_num_threads(no_threads_used);
    //            N_p = omp_get_max_threads();
    //            cout<<"threads being used parallely = "<<N_p<<endl;
    //            cout<<"No. of threads you asked for = "<<no_of_processors<<endl;
#endif

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
            no_threads_used = min(no_of_processors, No_Of_Inputs);
            omp_set_num_threads(no_threads_used);
            N_p = omp_get_max_threads();
            cout<<"threads being used parallely = "<<N_p<<endl;
            cout<<"No. of threads you asked for = "<<no_of_processors<<endl;
#pragma omp parallel for default(shared)
#endif
            for(int wi=0;wi<n_wpoints;wi++){

                for(int pos_i=0;pos_i<Parameters_.ns;pos_i++){
                    for(int pos_j=0;pos_j<Parameters_.ns;pos_j++){


                        //exp(iota*(-wi * dw) * (ts* dt_))
                        S_rw[pos_i][pos_j][wi] += exp(iota*(-wi * dw) * (ts* dt_))*exp(-0.5*(ts*dt_*w_conv*ts*dt_*w_conv))*dt_*(
                                    ( ( C_Classical_tr[ts][pos_i][pos_j]  )  )
                                    );

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
                                ( ( C_Classical_tr[ts][pos_i][pos_j]  )  )
                                );
                    Vec_2[ts] = exp(-0.5*(ts*dt_*w_conv*ts*dt_*w_conv))*dt_*(
                                ( ( C_Quantum_tr[ts][pos_i][pos_j]  )  )
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


void ST_Fourier_1orb_MCMF::Calculate_Skw_from_Crt_(string fileout){

    int N_p;
    int no_threads_used;
    no_threads_used = min(no_of_processors, No_Of_Inputs);

    Fft DO_Fft;
    DO_Fft.PI_EFF=PI;

    int pos_i, pos_j;

    complex<double> iota(0,1);
    S_rw.resize(Parameters_.lx);
    s_quantum_rw.resize(Parameters_.lx);

    for(int pos_ix=0;pos_ix<Parameters_.lx;pos_ix++){
        S_rw[pos_ix].resize(Parameters_.ly);
        s_quantum_rw[pos_ix].resize(Parameters_.ly);

        for(int pos_iy=0;pos_iy<Parameters_.ly;pos_iy++){
            S_rw[pos_ix][pos_iy].resize(n_wpoints);
            s_quantum_rw[pos_ix][pos_iy].resize(n_wpoints);

            for(int wi=0;wi<n_wpoints;wi++){
                S_rw[pos_ix][pos_iy][wi]=0;
                s_quantum_rw[pos_ix][pos_iy][wi]=0;
            }
        }
    }

    //    Sz_eq.resize(Parameters_.ns);Sx_eq.resize(Parameters_.ns);Sy_eq.resize(Parameters_.ns);
    //    sz_eq.resize(Parameters_.ns);sx_eq.resize(Parameters_.ns);sy_eq.resize(Parameters_.ns);

    //    Sz_t.resize(Parameters_.ns);Sx_t.resize(Parameters_.ns);Sy_t.resize(Parameters_.ns);
    //    sz_t.resize(Parameters_.ns);sx_t.resize(Parameters_.ns);sy_t.resize(Parameters_.ns);


    string line_temp;
    double temp_waste;

    double begin_time, end_time;
    clock_t oprt_SB_time;

#ifdef _OPENMP
    begin_time = omp_get_wtime();
#endif

    oprt_SB_time = clock();

    if(!Use_FFT){
        for(int ts=0;ts<time_steps;ts++){

#ifdef _OPENMP
            no_threads_used = min(no_of_processors, No_Of_Inputs);
            omp_set_num_threads(no_threads_used);
            N_p = omp_get_max_threads();
            cout<<"threads being used parallely = "<<N_p<<endl;
            cout<<"No. of threads you asked for = "<<no_of_processors<<endl;
#pragma omp parallel for default(shared) private(pos_i)
#endif
            for(int wi=0;wi<n_wpoints;wi++){

                for(int pos_ix=0;pos_ix<Parameters_.lx;pos_ix++){
                    for(int pos_iy=0;pos_iy<Parameters_.ly;pos_iy++){

                        pos_i = Coordinates_.Nc(pos_ix,pos_iy);

                        //exp(iota*(-wi * dw) * (ts* dt_))
                        S_rw[pos_ix][pos_iy][wi] += exp(iota*(-wi * dw) * (ts* dt_))*exp(-0.5*(ts*dt_*w_conv*ts*dt_*w_conv))*dt_*(
                                    (( C_Classical_tr_[ts][pos_i] ))
                                    );


                        s_quantum_rw[pos_ix][pos_iy][wi] += exp(iota*(-wi * dw) * (ts* dt_))*exp(-0.5*(ts*dt_*w_conv*ts*dt_*w_conv))*dt_*(
                                    ( ( C_Quantum_tr_[ts][pos_i]  )  )
                                    );

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
#pragma omp parallel for default(shared) private(pos_i)
#endif
        for(int pos_ix=0;pos_ix<Parameters_.lx;pos_ix++){

            Mat_1_Complex_doub Vec_1, Vec_2;
            Vec_1.resize(time_steps);Vec_2.resize(time_steps);

            for(int pos_iy=0;pos_iy<Parameters_.ly;pos_iy++){

                pos_i = Coordinates_.Nc(pos_ix,pos_iy);

                for(int ts=0;ts<time_steps;ts++){

                    Vec_1[ts] = exp(-0.5*(ts*dt_*w_conv*ts*dt_*w_conv))*dt_*(
                                ( ( C_Classical_tr_[ts][pos_i] )  )
                                );
                    Vec_2[ts] = exp(-0.5*(ts*dt_*w_conv*ts*dt_*w_conv))*dt_*(
                                ( ( C_Quantum_tr_[ts][pos_i] )  )
                                );
                }

                DO_Fft.transform(Vec_1);
                DO_Fft.transform(Vec_2);


                for(int wi=0;wi<n_wpoints;wi++){
                    S_rw[pos_ix][pos_iy][wi] = Vec_1[wi];
                    s_quantum_rw[pos_ix][pos_iy][wi] =Vec_2[wi];
                }
            }

            vector< complex<double> >().swap( Vec_1 );
            vector< complex<double> >().swap( Vec_2 );
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

                        pos_i = Coordinates_.Nc(x_i,y_i);

                        temp += S_rw[x_i][y_i][wi]*exp(iota*( (x_i)*kx +  (y_i)*ky ) );
                        //temp += S_rw[pos_i][pos_j][wi]*cos(( (x_j - x_i)*kx +  (y_j - y_i)*ky ) );
                        temp2 += s_quantum_rw[x_i][y_i][wi]*exp(iota*( (x_i)*kx +  (y_i)*ky ) );

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


void ST_Fourier_1orb_MCMF::Calculate_Fw_and_Aq(string fileout, string fileout_Aq){



    if(!Use_FFT){
        cout<<"FFT is NOT used"<<endl;
    }
    else{
        cout<<"FFT is used"<<endl;
    }

    int GaussianCenteredAtTmaxby2=0;
    int GCATm2=GaussianCenteredAtTmaxby2;
    double begin_time, end_time;
    clock_t oprt_SB_time;
    double expnt;

#ifdef _OPENMP
    begin_time = omp_get_wtime();
#endif

    oprt_SB_time = clock();


    int N_p;
    int no_threads_used;
    no_threads_used = min(no_of_processors, No_Of_Inputs);

    Fft DO_Fft;
    complex<double> iota(0,1);


    S_tr.resize(time_steps);
    for(int i=0;i<time_steps;i++){
        S_tr[i].resize(3*Parameters_.ns);
    }

    F_rw.resize(3*Parameters_.ns);
    F_qw.resize(3*Parameters_.ns);
    for(int i=0;i<3*Parameters_.ns;i++){
        F_rw[i].resize(n_wpoints);
        F_qw[i].resize(n_wpoints);
    }

    Aq.resize(3*Parameters_.ns);

    string line_temp2;

    ifstream Specific_conf_in;
    Specific_conf_in.open(conf_inputs[0].c_str());
    getline(Specific_conf_in, line_temp2);

    for(int ts=0;ts<time_steps;ts++)
    {
        string line_temp;
        double row_time, double_temp;
        int row_ts;
        getline(Specific_conf_in, line_temp);
        stringstream line_temp_ss(line_temp, stringstream::in);

        line_temp_ss>>row_time;
        row_ts = (int) ((row_time)/dt_ + 0.5);
        if(row_ts != ts){
            cout<<endl;
            cout<<"Error;"<<endl;
            cout<< row_time<<"   "<<row_ts <<"  "<<ts<<endl;
            cout<< "Problem while reading from "<<conf_inputs[0]<<endl<<endl;
            assert(false);
        }

        for(int r=0;r<3*Parameters_.ns;r++){
            line_temp_ss>>double_temp;
            S_tr[ts][r] = double_temp;
        }

        if(ts%100==0){
            cout<<"Reading Str upto time slice i.e. t/dt="<<ts<<endl;
        }

    }


    cout<<"Reading S[t][r] is completed"<<endl;



    int pos_i, index;

#ifdef _OPENMP
    no_threads_used = min(no_of_processors, Parameters_.ns);
    omp_set_num_threads(no_threads_used);
    N_p = omp_get_max_threads();
    cout<<"threads being used parallely = "<<N_p<<endl;
    cout<<"No. of threads you asked for = "<<no_of_processors<<endl;
#pragma omp parallel for default(shared) private(pos_i, index)
#endif
    for(int pos_ix=0;pos_ix<Parameters_.lx;pos_ix++){

        Mat_1_Complex_doub Vec_1;
        Vec_1.resize(time_steps);

        for(int pos_iy=0;pos_iy<Parameters_.ly;pos_iy++){

            for(int type=0;type<3;type++){

                pos_i = Coordinates_.Nc(pos_ix,pos_iy);
                index  = pos_i + (type*Parameters_.ns);


                if(!Use_FFT){

                    for(int wi=0;wi<n_wpoints;wi++){

                        for(int ts=0;ts<time_steps;ts++){

                            expnt = (ts - (0.5*GCATm2*time_steps))*dt_*w_conv;
                            expnt = expnt*expnt;
                            F_rw[index][wi] += exp(iota*(-wi * dw) * (ts* dt_))*exp(-0.5*expnt)*dt_*(
                                        ((  S_tr[ts][index] ))
                                        );
                        }

                    }
                }

                else{

                    for(int ts=0;ts<time_steps;ts++){

                        expnt = (ts - (0.5*GCATm2*time_steps))*dt_*w_conv;
                        expnt = expnt*expnt;
                        Vec_1[ts] = exp(-0.5*expnt)*fabs(dt_)*(( ( S_tr[ts][index] )  ));
                    }

                    DO_Fft.PI_EFF=PI*(dt_/fabs(dt_));
                    DO_Fft.transform(Vec_1);

                    for(int wi=0;wi<n_wpoints;wi++){
                        //F_rw[index][wi] = Vec_1[wi].real();
                        F_rw[index][wi] = Vec_1[wi];

                    }
                }

            }

        }

        vector< complex<double> >().swap( Vec_1 );

        cout<<"FFT for S[t][rx="   << pos_ix <<  ", ry] is completed for all ry, and all 6 types of S"<<endl;

    }



    double kx, ky;
    int k_index;
    complex<double> temp;

    for(int type=0;type<3;type++){
        for(int nx=0;nx<Parameters_.lx;nx++){
            for(int ny=0;ny<Parameters_.ly;ny++){

                kx = (2*nx*PI)/(1.0*Parameters_.lx);
                ky = (2*ny*PI)/(1.0*Parameters_.ly);
                k_index = Coordinates_.Nc(nx,ny) + (type*Parameters_.ns);

                temp = complex<double>(0,0);
                for(int x_i=0;x_i<Parameters_.lx;x_i++){
                    for(int y_i=0;y_i<Parameters_.ly;y_i++){

                        pos_i = Coordinates_.Nc(x_i,y_i);
                        index  = pos_i + (type*Parameters_.ns);
                        temp += S_tr[(GCATm2*time_steps)/2][index]*exp(iota*(1.0*( (x_i)*kx +  (y_i)*ky ) ) );

                    }
                }
                Aq[k_index] = temp;

            }
        }
    }




    for(int type=0;type<3;type++){
        for(int nx=0;nx<Parameters_.lx;nx++){
            for(int ny=0;ny<Parameters_.ly;ny++){

                kx = (2*nx*PI)/(1.0*Parameters_.lx);
                ky = (2*ny*PI)/(1.0*Parameters_.ly);
                k_index = Coordinates_.Nc(nx,ny) + (type*Parameters_.ns);

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(pos_i, index, temp)
#endif

                for(int wi=0;wi<n_wpoints;wi++){

                    temp = complex<double> (0,0);
                    for(int x_i=0;x_i<Parameters_.lx;x_i++){
                        for(int y_i=0;y_i<Parameters_.ly;y_i++){

                            pos_i = Coordinates_.Nc(x_i,y_i);
                            index  = pos_i + (type*Parameters_.ns);

                            temp += F_rw[index][wi]*exp(iota*(-1.0*( (x_i)*kx +  (y_i)*ky ) ) );

                        }
                    }
                    //file_out2<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<temp.real()<<"   "<<temp.imag()<<"    "<<temp2.real()<<"   "<<temp2.imag()<<"    "<<temp3.real()<<"   "<<temp3.imag()<<endl;
                    //file_out_full<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<wi<<"   "<<temp.real()<<"   "<<temp.imag()<<"    "<<"none"<<"   "<<"none"<<"    "<<"none"<<"   "<<"none"<<endl;
                    F_qw[k_index][wi]=temp*(1.0/(Parameters_.lx*Parameters_.ly));
                }
            }
        }

    }




#ifdef _OPENMP
    end_time = omp_get_wtime();
    cout<<"Time to Calculate F(kx,ky,w) [actual time, using OMP]: "<<double(end_time - begin_time)<<endl;//cout<<"here"<<endl;

#endif

    cout<<"Time to Calculate F(kx,ky,w) : "<<double( clock() - oprt_SB_time ) / (double)CLOCKS_PER_SEC<<endl;//cout<<"here"<<endl;





    cout<< "Printing results in files"<<endl;

    ofstream file_out_full(fileout.c_str());

    file_out_full<<"#w_details:  "<<"   "<<w_max<<"  "<<w_min<<"  "<<dw<<"  "<<n_wpoints<<endl;
    file_out_full<<"#nx   ny   k_ind   wi*dw   wi   F_qw[k_ind][wi].real()   F_qw[k_ind][wi].imag()"<<endl;
    for(int nx=0;nx<Parameters_.lx;nx++){
        for(int ny=0;ny<Parameters_.ly;ny++){
            k_index=Coordinates_.Nc(nx,ny);

            for(int wi=0;wi<n_wpoints;wi++){

                //file_out2<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<temp.real()<<"   "<<temp.imag()<<"    "<<temp2.real()<<"   "<<temp2.imag()<<"    "<<temp3.real()<<"   "<<temp3.imag()<<endl;
                file_out_full<<nx<<"   "<<ny<<"   "<<k_index<<"   "<<wi*dw<<"   "<<wi<<"   ";

                for(int type=0;type<3;type++){

                    file_out_full<<F_qw[k_index + (type*Parameters_.ns)][wi].real()<<"   "<<F_qw[k_index + (type*Parameters_.ns)][wi].imag()<<"   ";
                }

                file_out_full<<endl;

            }

            file_out_full<<endl;

        }
    }



    ofstream file_out_Aq(fileout_Aq.c_str());
    file_out_Aq<<"#nx   ny   k_ind     Aq.real()    Aq.imag()"<<endl;

    for(int nx=0;nx<Parameters_.lx;nx++){
        for(int ny=0;ny<Parameters_.ly;ny++){
            k_index=Coordinates_.Nc(nx,ny);

            file_out_Aq<<nx<<"   "<<ny<<"   "<<k_index<<"   ";

            for(int type=0;type<3;type++){
                file_out_Aq<<Aq[k_index + (type*Parameters_.ns)].real()<<"   "<<Aq[k_index + (type*Parameters_.ns)].imag()<<"    ";

            }
            file_out_Aq<<endl;

        }
        file_out_Aq<<endl;
    }


    cout<<"Calculate_Fw_and_Aq completed"<<endl;

}


void ST_Fourier_1orb_MCMF::Calculate_Sqw_using_Fwq_with_negativeandpositive_Time(string fileout, string Dqw_file){


    F_qw_.resize(No_Of_Inputs);

    string line_temp2;
    ifstream file_Fwq(F_wq_inputs[0].c_str());
    getline(file_Fwq, line_temp2);
    stringstream line_temp2_ss(line_temp2, stringstream::in);

    //"#w_details:  "<<"   "<<w_max<<"  "<<w_min<<"  "<<dw<<"  "<<n_wpoints<<endl;

    line_temp2_ss>>line_temp2;
    line_temp2_ss>>w_max;
    line_temp2_ss>>w_min;
    line_temp2_ss>>dw;
    line_temp2_ss>>n_wpoints;


    S_qw.resize(3*Parameters_.ns);
    D_qw.resize(3*Parameters_.ns);
    for(int i=0;i<3*Parameters_.ns;i++){
        S_qw[i].resize(n_wpoints);
        D_qw[i].resize(n_wpoints);
        for(int n=0;n<n_wpoints;n++){
            S_qw[i][n]=complex<double>(0,0);
            D_qw[i][n]=complex<double>(0,0);
        }
    }



    for(int j=0;j<F_qw_.size();j++){
        F_qw_[j].resize(3*Parameters_.ns);
        for(int i=0;i<3*Parameters_.ns;i++){
            F_qw_[j][i].resize(n_wpoints);
        }
    }



    string line_temp;
    int nx_temp, ny_temp, k_index_temp, k_index, wi_temp;
    double double_temp, temp1, temp2;


    for(int ms=0;ms<(No_Of_Inputs/2);ms++){


        ifstream file_Fwq_in(F_wq_inputs[ms].c_str());
        getline(file_Fwq_in, line_temp);
        getline(file_Fwq_in, line_temp);
        //stringstream line_temp_ss(line_temp, stringstream::in);

        // file_out_full<<"#nx   ny   k_ind   wi*dw   wi   F_qw[k_ind][wi].real()   F_qw[k_ind][wi].imag()"<<endl;
        for(int nx=0;nx<Parameters_.lx;nx++){
            for(int ny=0;ny<Parameters_.ly;ny++){
                k_index=Coordinates_.Nc(nx,ny);

                for(int wi=0;wi<n_wpoints;wi++){

                    // file_out_full<<nx<<"   "<<ny<<"   "<<k_index<<"   "<<wi*dw<<"   "<<wi<<"   ";
                    getline(file_Fwq_in, line_temp);
                    stringstream line_temp_ss(line_temp, stringstream::in);
                    line_temp_ss>>nx_temp>>ny_temp>>k_index_temp;
                    line_temp_ss>>double_temp>>wi_temp;
                    assert((nx_temp==nx) && (ny_temp==ny) && (k_index_temp==k_index) && (wi==wi_temp) );

                    for(int type=0;type<3;type++){

                        //file_out_full<<F_qw[k_index + (type*Parameters_.ns)][wi].real()<<"   "
                        //<<F_qw[k_index + (type*Parameters_.ns)][wi].imag()<<"   ";

                        line_temp_ss>>temp1>>temp2;
                        F_qw_[ms][k_index + (type*Parameters_.ns)][wi]=complex<double>(temp1, temp2);


                        S_qw[k_index + (type*Parameters_.ns)][wi] += (  (F_qw_[ms][k_index + (type*Parameters_.ns)][wi]  +  F_qw_[ms+(No_Of_Inputs/2)][k_index + (type*Parameters_.ns)][wi])
                                *(  conj(F_qw_[ms][k_index + (type*Parameters_.ns)][wi])  + conj(F_qw_[ms+(No_Of_Inputs/2)][k_index + (type*Parameters_.ns)][wi]) )
                                )
                                *(1.0/(0.5*No_Of_Inputs));



                        D_qw[k_index + (type*Parameters_.ns)][wi] += (  (F_qw_[ms][k_index + (type*Parameters_.ns)][wi]  +  F_qw_[ms+(No_Of_Inputs/2)][k_index + (type*Parameters_.ns)][wi])
                                *(  conj(F_qw_[ms][k_index + (type*Parameters_.ns)][wi])  + conj(F_qw_[ms+(No_Of_Inputs/2)][k_index + (type*Parameters_.ns)][wi]) )
                                )
                                *(1.0/(0.25*No_Of_Inputs*No_Of_Inputs));;

                    }

                    //file_out_full<<endl;

                }

                //file_out_full<<endl;
                getline(file_Fwq_in, line_temp);

            }
        }


    }



    cout<<"Now Sqw is being calculated :"<<endl;


    for(int i=0;i<3*Parameters_.ns;i++){
        for(int wi=0;wi<n_wpoints;wi++){

            for(int ms=0;ms<No_Of_Inputs/2;ms++){
                for(int ns=0;ns<No_Of_Inputs/2;ns++){

                    S_qw[i][wi] = S_qw[i][wi] -
                            ((     (  (F_qw_[ms][i][wi]  +  F_qw_[ms+(No_Of_Inputs/2)][i][wi])
                                   *(  conj(F_qw_[ns][i][wi])  + conj(F_qw_[ns+(No_Of_Inputs/2)][i][wi]) )
                             )
                            )
                            *(1.0/(0.25*No_Of_Inputs*No_Of_Inputs))   );

                }
            }
        }
    }



    ofstream file_out_full(fileout.c_str());

    file_out_full<<"#nx   ny   k_ind   wi*dw   wi   S_qw[k_ind][wi].real()   S_qw[k_ind][wi].imag()"<<endl;
    for(int nx=0;nx<Parameters_.lx;nx++){
        for(int ny=0;ny<Parameters_.ly;ny++){
            k_index=Coordinates_.Nc(nx,ny);

            for(int wi=0;wi<n_wpoints;wi++){

                //file_out2<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<temp.real()<<"   "<<temp.imag()<<"    "<<temp2.real()<<"   "<<temp2.imag()<<"    "<<temp3.real()<<"   "<<temp3.imag()<<endl;
                file_out_full<<nx<<"   "<<ny<<"   "<<k_index<<"   "<<wi*dw<<"   "<<wi<<"   ";

                for(int type=0;type<3;type++){

                    file_out_full<<S_qw[k_index + (type*Parameters_.ns)][wi].real()<<"   "<<S_qw[k_index + (type*Parameters_.ns)][wi].imag()<<"   ";
                }

                file_out_full<<endl;

            }

            file_out_full<<endl;

        }
    }





    ofstream file_out(Dqw_file.c_str());

    file_out<<"#nx   ny   k_ind   wi*dw   wi   |D_qw[k_ind][wi]|   (D_qw[k_ind][wi].real())  (D_qw[k_ind][wi].imag()) "<<endl;
    for(int nx=0;nx<Parameters_.lx;nx++){
        for(int ny=0;ny<Parameters_.ly;ny++){
            k_index=Coordinates_.Nc(nx,ny);

            for(int wi=0;wi<n_wpoints;wi++){

                //file_out2<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<temp.real()<<"   "<<temp.imag()<<"    "<<temp2.real()<<"   "<<temp2.imag()<<"    "<<temp3.real()<<"   "<<temp3.imag()<<endl;
                file_out<<nx<<"   "<<ny<<"   "<<k_index<<"   "<<wi*dw<<"   "<<wi<<"   ";

                for(int type=0;type<3;type++){
                    file_out<< abs(D_qw[k_index + (type*Parameters_.ns)][wi])<<"    ";
                    file_out<< D_qw[k_index + (type*Parameters_.ns)][wi].real()<<"    ";
                    file_out<< D_qw[k_index + (type*Parameters_.ns)][wi].imag()<<"    ";
                }

                file_out<<endl;
            }

            file_out<<endl;
        }
    }



}


void ST_Fourier_1orb_MCMF::Calculate_Sqw_using_Fwq(string fileout, string Dqw_file){


   // F_qw_.resize(1);

    string line_temp2;
    ifstream file_Fwq(F_wq_inputs[0].c_str());
    getline(file_Fwq, line_temp2);
    stringstream line_temp2_ss(line_temp2, stringstream::in);

    //"#w_details:  "<<"   "<<w_max<<"  "<<w_min<<"  "<<dw<<"  "<<n_wpoints<<endl;

    line_temp2_ss>>line_temp2;
    line_temp2_ss>>w_max;
    line_temp2_ss>>w_min;
    line_temp2_ss>>dw;
    line_temp2_ss>>n_wpoints;


    S_qw.resize(3*Parameters_.ns);
    //D_qw.resize(3*Parameters_.ns);
    for(int i=0;i<3*Parameters_.ns;i++){
        S_qw[i].resize(n_wpoints);
        //D_qw[i].resize(n_wpoints);
        for(int n=0;n<n_wpoints;n++){
            S_qw[i][n]=complex<double>(0,0);
           // D_qw[i][n]=complex<double>(0,0);
        }
    }




//    for(int j=0;j<F_qw_.size();j++){
//        F_qw_[j].resize(3*Parameters_.ns);
//        for(int i=0;i<3*Parameters_.ns;i++){
//            F_qw_[j][i].resize(n_wpoints);
//        }
//    }

    complex<double> fqw_temp;

     cout<<"F_qw_ initializing"<<endl;
    F_qw.resize(3*Parameters_.ns);
    for(int i=0;i<3*Parameters_.ns;i++){
        F_qw[i].resize(n_wpoints);
    }
    cout<<"done"<<endl;

    string line_temp;
    int nx_temp, ny_temp, k_index_temp, k_index, wi_temp;
    double double_temp, temp1, temp2;


    for(int ms=0;ms<(No_Of_Inputs);ms++){


        ifstream file_Fwq_in(F_wq_inputs[ms].c_str());
        cout<<"'"<<F_wq_inputs[ms]<<"'"<<endl;
        getline(file_Fwq_in, line_temp);
        getline(file_Fwq_in, line_temp);
        //stringstream line_temp_ss(line_temp, stringstream::in);

        // file_out_full<<"#nx   ny   k_ind   wi*dw   wi   F_qw[k_ind][wi].real()   F_qw[k_ind][wi].imag()"<<endl;
        for(int nx=0;nx<Parameters_.lx;nx++){
            for(int ny=0;ny<Parameters_.ly;ny++){
                k_index=Coordinates_.Nc(nx,ny);

                for(int wi=0;wi<n_wpoints;wi++){

                    // file_out_full<<nx<<"   "<<ny<<"   "<<k_index<<"   "<<wi*dw<<"   "<<wi<<"   ";
                    getline(file_Fwq_in, line_temp);
                    stringstream line_temp_ss(line_temp, stringstream::in);
                    line_temp_ss>>nx_temp>>ny_temp>>k_index_temp;
                    line_temp_ss>>double_temp>>wi_temp;
                    assert((nx_temp==nx) && (ny_temp==ny) && (k_index_temp==k_index) && (wi==wi_temp) );

                    for(int type=0;type<3;type++){

                        //file_out_full<<F_qw[k_index + (type*Parameters_.ns)][wi].real()<<"   "
                        //<<F_qw[k_index + (type*Parameters_.ns)][wi].imag()<<"   ";

                        line_temp_ss>>temp1>>temp2;
                        fqw_temp=complex<double>(temp1, temp2);
                        F_qw[k_index + (type*Parameters_.ns)][wi] +=(1.0/(1.0*No_Of_Inputs))*complex<double>(temp1, temp2);

                        S_qw[k_index + (type*Parameters_.ns)][wi] += (fqw_temp*conj(fqw_temp))*(1.0/(1.0*No_Of_Inputs));
                       // D_qw[k_index + (type*Parameters_.ns)][wi] += (fqw_temp*conj(fqw_temp))*(1.0/(1.0*No_Of_Inputs*No_Of_Inputs));

                    }


                    //file_out_full<<endl;

                }

                //file_out_full<<endl;
                getline(file_Fwq_in, line_temp);

            }
        }

        cout<<"Microstate "<<ms<<" done"<<endl;
    }






    ofstream file_out(Dqw_file.c_str());

    file_out<<"#nx   ny   k_ind   wi*dw   wi   |D_qw[k_ind][wi]|   (D_qw[k_ind][wi].real())  (D_qw[k_ind][wi].imag()) "<<endl;
    for(int nx=0;nx<Parameters_.lx;nx++){
        for(int ny=0;ny<Parameters_.ly;ny++){
            k_index=Coordinates_.Nc(nx,ny);

            for(int wi=0;wi<n_wpoints;wi++){

                //file_out2<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<temp.real()<<"   "<<temp.imag()<<"    "<<temp2.real()<<"   "<<temp2.imag()<<"    "<<temp3.real()<<"   "<<temp3.imag()<<endl;
                file_out<<nx<<"   "<<ny<<"   "<<k_index<<"   "<<wi*dw<<"   "<<wi<<"   ";

                for(int type=0;type<3;type++){
                    file_out<< abs(S_qw[k_index + (type*Parameters_.ns)][wi])*(1.0/(1.0*No_Of_Inputs))<<"    ";
                    file_out<< S_qw[k_index + (type*Parameters_.ns)][wi].real()*(1.0/(1.0*No_Of_Inputs))<<"    ";
                    file_out<< S_qw[k_index + (type*Parameters_.ns)][wi].imag()*(1.0/(1.0*No_Of_Inputs))<<"    ";
                }

                file_out<<endl;
            }

            file_out<<endl;
        }
    }





    cout<<"Now Sqw is being calculated :"<<endl;


    for(int i=0;i<3*Parameters_.ns;i++){
        for(int wi=0;wi<n_wpoints;wi++){

            S_qw[i][wi] = S_qw[i][wi] -
                    ((F_qw[i][wi]*conj(F_qw[i][wi]) )   );

        }
    }



    ofstream file_out_full(fileout.c_str());

    file_out_full<<"#nx   ny   k_ind   wi*dw   wi   S_qw[k_ind][wi].real()   S_qw[k_ind][wi].imag()"<<endl;
    for(int nx=0;nx<Parameters_.lx;nx++){
        for(int ny=0;ny<Parameters_.ly;ny++){
            k_index=Coordinates_.Nc(nx,ny);

            for(int wi=0;wi<n_wpoints;wi++){

                //file_out2<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<temp.real()<<"   "<<temp.imag()<<"    "<<temp2.real()<<"   "<<temp2.imag()<<"    "<<temp3.real()<<"   "<<temp3.imag()<<endl;
                file_out_full<<nx<<"   "<<ny<<"   "<<k_index<<"   "<<wi*dw<<"   "<<wi<<"   ";

                for(int type=0;type<3;type++){

                    file_out_full<<S_qw[k_index + (type*Parameters_.ns)][wi].real()<<"   "<<S_qw[k_index + (type*Parameters_.ns)][wi].imag()<<"   ";
                }

                file_out_full<<endl;

            }

            file_out_full<<endl;

        }
    }









}



void ST_Fourier_1orb_MCMF::Calculate_Sqw_using_Aq_Fwq(string fileout, string Dqw_file){


    string line_temp2;
    ifstream file_Fwq(F_wq_inputs[0].c_str());
    getline(file_Fwq, line_temp2);
    stringstream line_temp2_ss(line_temp2, stringstream::in);

    //"#w_details:  "<<"   "<<w_max<<"  "<<w_min<<"  "<<dw<<"  "<<n_wpoints<<endl;

    line_temp2_ss>>line_temp2;
    line_temp2_ss>>w_max;
    line_temp2_ss>>w_min;
    line_temp2_ss>>dw;
    line_temp2_ss>>n_wpoints;


    S_qw.resize(3*Parameters_.ns);
    D_qw.resize(3*Parameters_.ns);
    for(int i=0;i<3*Parameters_.ns;i++){
        S_qw[i].resize(n_wpoints);
        D_qw[i].resize(n_wpoints);
        for(int n=0;n<n_wpoints;n++){
            S_qw[i][n]=complex<double>(0,0);
            D_qw[i][n]=complex<double>(0,0);
        }
    }

    F_qw.resize(3*Parameters_.ns);
    for(int i=0;i<3*Parameters_.ns;i++){
        F_qw[i].resize(n_wpoints);
    }

    Aq.resize(3*Parameters_.ns);
    delAq.resize(3*Parameters_.ns);
    Aq_avg.resize(3*Parameters_.ns);
    for(int i=0;i<3*Parameters_.ns;i++){
        Aq_avg[i]=complex<double>(0,0);
    }




    string line_temp;
    int nx_temp, ny_temp, k_index_temp, k_index, wi_temp;
    double double_temp, temp1, temp2;


    for(int ms=0;ms<No_Of_Inputs;ms++){

        ifstream file_Aq_in(Aq_inputs[ms].c_str());
        getline(file_Aq_in, line_temp);

        //file_out_Aq<<"#nx   ny   k_ind     Aq.real()    Aq.imag()"<<endl;

        for(int nx=0;nx<Parameters_.lx;nx++){
            for(int ny=0;ny<Parameters_.ly;ny++){
                k_index=Coordinates_.Nc(nx,ny);

                getline(file_Aq_in, line_temp);
                stringstream line_temp_ss(line_temp, stringstream::in);
                line_temp_ss>>nx_temp>>ny_temp>>k_index_temp;
                assert((nx_temp==nx) && (ny_temp==ny) && (k_index==k_index_temp));
                // file_out_Aq<<nx<<"   "<<ny<<"   "<<k_index<<"   ";


                for(int type=0;type<3;type++){
                    //file_out_Aq<<Aq[k_index + (type*Parameters_.ns)].real()<<"   "<<Aq[k_index + (type*Parameters_.ns)].imag()<<"    ";

                    line_temp_ss>>temp1>>temp2;
                    //Aq[k_index + (type*Parameters_.ns)]=complex<double> (temp1, temp2);
                    Aq_avg[k_index + (type*Parameters_.ns)] +=complex<double> (temp1, temp2);

                }
                //file_out_Aq<<endl;

            }
            //file_out_Aq<<endl;
            getline(file_Aq_in, line_temp);
        }
    }


    for(int i=0;i<3*Parameters_.ns;i++){
        Aq_avg[i]=Aq_avg[i]*(1.0/(No_Of_Inputs*1.0));
    }



    for(int ms=0;ms<No_Of_Inputs;ms++){



        ifstream file_Aq_in(Aq_inputs[ms].c_str());
        getline(file_Aq_in, line_temp);

        //file_out_Aq<<"#nx   ny   k_ind     Aq.real()    Aq.imag()"<<endl;

        for(int nx=0;nx<Parameters_.lx;nx++){
            for(int ny=0;ny<Parameters_.ly;ny++){
                k_index=Coordinates_.Nc(nx,ny);

                getline(file_Aq_in, line_temp);
                stringstream line_temp_ss(line_temp, stringstream::in);
                line_temp_ss>>nx_temp>>ny_temp>>k_index_temp;
                assert((nx_temp==nx) && (ny_temp==ny) && (k_index==k_index_temp));
                // file_out_Aq<<nx<<"   "<<ny<<"   "<<k_index<<"   ";


                for(int type=0;type<3;type++){
                    //file_out_Aq<<Aq[k_index + (type*Parameters_.ns)].real()<<"   "<<Aq[k_index + (type*Parameters_.ns)].imag()<<"    ";

                    line_temp_ss>>temp1>>temp2;
                    //Aq[k_index + (type*Parameters_.ns)]=complex<double> (temp1, temp2);
                    delAq[k_index + (type*Parameters_.ns)] =complex<double> (temp1, temp2) -  Aq_avg[k_index + (type*Parameters_.ns)];

                }
                //file_out_Aq<<endl;

            }
            //file_out_Aq<<endl;
            getline(file_Aq_in, line_temp);
        }






        ifstream file_Fwq_in(F_wq_inputs[ms].c_str());
        getline(file_Fwq_in, line_temp);
        getline(file_Fwq_in, line_temp);
        //stringstream line_temp_ss(line_temp, stringstream::in);

        // file_out_full<<"#nx   ny   k_ind   wi*dw   wi   F_qw[k_ind][wi].real()   F_qw[k_ind][wi].imag()"<<endl;
        for(int nx=0;nx<Parameters_.lx;nx++){
            for(int ny=0;ny<Parameters_.ly;ny++){
                k_index=Coordinates_.Nc(nx,ny);

                for(int wi=0;wi<n_wpoints;wi++){

                    // file_out_full<<nx<<"   "<<ny<<"   "<<k_index<<"   "<<wi*dw<<"   "<<wi<<"   ";
                    getline(file_Fwq_in, line_temp);
                    stringstream line_temp_ss(line_temp, stringstream::in);
                    line_temp_ss>>nx_temp>>ny_temp>>k_index_temp;
                    line_temp_ss>>double_temp>>wi_temp;
                    assert((nx_temp==nx) && (ny_temp==ny) && (k_index_temp==k_index) && (wi==wi_temp) );

                    for(int type=0;type<3;type++){

                        //file_out_full<<F_qw[k_index + (type*Parameters_.ns)][wi].real()<<"   "
                        //<<F_qw[k_index + (type*Parameters_.ns)][wi].imag()<<"   ";

                        line_temp_ss>>temp1>>temp2;
                        F_qw[k_index + (type*Parameters_.ns)][wi]=complex<double>(temp1, temp2);

                        S_qw[k_index + (type*Parameters_.ns)][wi] += F_qw[k_index + (type*Parameters_.ns)][wi]*delAq[k_index + (type*Parameters_.ns)];

                        D_qw[k_index + (type*Parameters_.ns)][wi] += F_qw[k_index + (type*Parameters_.ns)][wi];
                    }

                    //file_out_full<<endl;

                }

                //file_out_full<<endl;
                getline(file_Fwq_in, line_temp);

            }
        }


    }



    for(int i=0;i<3*Parameters_.ns;i++){
        for(int n=0;n<n_wpoints;n++){
            S_qw[i][n]=S_qw[i][n]*(1.0/(1.0*No_Of_Inputs*Parameters_.ns));
            D_qw[i][n]=D_qw[i][n]*(1.0/(1.0*No_Of_Inputs*Parameters_.ns*Parameters_.ns));
        }
    }





    ofstream file_out_full(fileout.c_str());

    file_out_full<<"#nx   ny   k_ind   wi*dw   wi   S_qw[k_ind][wi].real()   S_qw[k_ind][wi].imag()"<<endl;
    for(int nx=0;nx<Parameters_.lx;nx++){
        for(int ny=0;ny<Parameters_.ly;ny++){
            k_index=Coordinates_.Nc(nx,ny);

            for(int wi=0;wi<n_wpoints;wi++){

                //file_out2<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<temp.real()<<"   "<<temp.imag()<<"    "<<temp2.real()<<"   "<<temp2.imag()<<"    "<<temp3.real()<<"   "<<temp3.imag()<<endl;
                file_out_full<<nx<<"   "<<ny<<"   "<<k_index<<"   "<<wi*dw<<"   "<<wi<<"   ";

                for(int type=0;type<3;type++){

                    file_out_full<<S_qw[k_index + (type*Parameters_.ns)][wi].real()<<"   "<<S_qw[k_index + (type*Parameters_.ns)][wi].imag()<<"   ";
                }

                file_out_full<<endl;

            }

            file_out_full<<endl;

        }
    }





    ofstream file_out(Dqw_file.c_str());

    file_out<<"#nx   ny   k_ind   wi*dw   wi   |D_qw[k_ind][wi]|^2   (D_qw[k_ind][wi].real())^2  "<<endl;
    for(int nx=0;nx<Parameters_.lx;nx++){
        for(int ny=0;ny<Parameters_.ly;ny++){
            k_index=Coordinates_.Nc(nx,ny);

            for(int wi=0;wi<n_wpoints;wi++){

                //file_out2<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<temp.real()<<"   "<<temp.imag()<<"    "<<temp2.real()<<"   "<<temp2.imag()<<"    "<<temp3.real()<<"   "<<temp3.imag()<<endl;
                file_out<<nx<<"   "<<ny<<"   "<<k_index<<"   "<<wi*dw<<"   "<<wi<<"   ";

                for(int type=0;type<3;type++){
                    file_out<< abs(D_qw[k_index + (type*Parameters_.ns)][wi] * conj(D_qw[k_index + (type*Parameters_.ns)][wi]))<<"    ";
                    file_out<< D_qw[k_index + (type*Parameters_.ns)][wi].real() * D_qw[k_index + (type*Parameters_.ns)][wi].real()<<"    ";
                }

                file_out<<endl;
            }

            file_out<<endl;
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

void ST_Fourier_1orb_MCMF::Convolute_the_spectrum(string Sqw_file_in, string Sqw_file_out){



    string line_temp;

    ifstream Sqw_in;
    Sqw_in.open(Sqw_file_in.c_str());

    //stringstream line_temp_ss(line_temp, stringstream::in);

    int nx, ny, k_ind, omega_ind;
    double omega_val;
    double real_temp, imag_temp;

    S_qw.resize(3*Parameters_.ns);
    for(int i=0;i<S_qw.size();i++){
        S_qw[i].clear();
    }






    for(int curLine = 0; getline(Sqw_in, line_temp); curLine++) {

        if( curLine>0 && ( !line_temp.empty() ) ){

            stringstream line_temp_ss(line_temp, stringstream::in);

            line_temp_ss >> nx >> ny >> k_ind >> omega_val >> omega_ind;

            for(int type=0;type<3;type++){
                line_temp_ss >> real_temp >> imag_temp;
                S_qw[k_ind + (type*Parameters_.ns)].push_back( complex<double> (real_temp, imag_temp) );
            }


            w_max = omega_val;
            n_wpoints = omega_ind + 1;

        }

    }

    dw = w_max/(1.0*(n_wpoints - 1));

    cout<<"dw = "<<dw<<endl;
    cout<<"w_max = "<<w_max<<endl;
    cout<<"n_wpoints = "<<n_wpoints<<endl;

    S_qw_conv.resize(3*Parameters_.ns);
    for(int i=0;i<S_qw_conv.size();i++){
        S_qw_conv[i].resize(n_wpoints);
    }


    double wj_val, wi_val;


    for(int i=0;i<3*Parameters_.ns;i++){
        for(int wj=0;wj<n_wpoints;wj++){
            wj_val = 1.0*wj*dw;
            S_qw_conv[i][wj] = zero_complex;

            for(int wi =0;wi<n_wpoints;wi++){
                wi_val = 1.0*wi*dw;

                S_qw_conv[i][wj] += (dw/w_conv)*S_qw[i][wi]*exp( (-1.0*(pow(wj_val - wi_val,2.0) )) / (2.0*w_conv*w_conv)  );

            }

        }
    }


    cout<<"Smoothening done"<<endl;




    ofstream Sqw_out(Sqw_file_out.c_str());

    int k_index;
    Sqw_out<<"#nx   ny   k_ind   wi*dw   wi   S_qw[k_ind][wi].real()   S_qw[k_ind][wi].imag()"<<endl;
    for(int nx=0;nx<Parameters_.lx;nx++){
        for(int ny=0;ny<Parameters_.ly;ny++){
            k_index=Coordinates_.Nc(nx,ny);

            for(int wi=0;wi<n_wpoints;wi++){

                //file_out2<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<temp.real()<<"   "<<temp.imag()<<"    "<<temp2.real()<<"   "<<temp2.imag()<<"    "<<temp3.real()<<"   "<<temp3.imag()<<endl;
                Sqw_out<<nx<<"   "<<ny<<"   "<<k_index<<"   "<<wi*dw<<"   "<<wi<<"   ";

                for(int type=0;type<3;type++){

                    Sqw_out<<S_qw_conv[k_index + (type*Parameters_.ns)][wi].real()<<"   "<<S_qw_conv[k_index + (type*Parameters_.ns)][wi].imag()<<"   ";
                }

                Sqw_out<<endl;

            }

            Sqw_out<<endl;

        }
    }

    Sqw_out<<"#w_conv = "<<w_conv<<endl;


}


void ST_Fourier_1orb_MCMF::Sum_Rule_For_Way1(string Sqw_file_in, string Sq_static_out, string Sq_dynamical_out){


    string line_temp;

    ifstream Sqw_in;
    Sqw_in.open(Sqw_file_in.c_str());

    //stringstream line_temp_ss(line_temp, stringstream::in);

    int nx, ny, k_ind, omega_ind;
    double omega_val;
    double real_temp, imag_temp;

    Mat_2_Complex_doub Sq_dyn, Sq_static;

    Sq_dyn.resize(Parameters_.lx);
    Sq_static.resize(Parameters_.lx);
    for(int ix=0;ix<Parameters_.lx;ix++){
        Sq_dyn[ix].resize(Parameters_.ly);
        Sq_static[ix].resize(Parameters_.ly);
    }


    S_qw.resize(3*Parameters_.ns);
    for(int i=0;i<S_qw.size();i++){
        S_qw[i].clear();
    }




    for(int curLine = 0; getline(Sqw_in, line_temp); curLine++) {

        if( curLine>0 && ( !line_temp.empty() ) ){

            stringstream line_temp_ss(line_temp, stringstream::in);

            line_temp_ss >> nx >> ny >> k_ind >> omega_val >> omega_ind;

            for(int type=0;type<3;type++){
                line_temp_ss >> real_temp >> imag_temp;
                S_qw[k_ind + (type*Parameters_.ns)].push_back( complex<double> (real_temp, imag_temp) );
            }


            w_max = omega_val;
            n_wpoints = omega_ind + 1;

        }

    }

    dw = w_max/(1.0*(n_wpoints - 1));

    cout<<"dw = "<<dw<<endl;
    cout<<"w_max = "<<w_max<<endl;
    cout<<"n_wpoints = "<<n_wpoints<<endl;




    for(int nx_=0;nx_<Parameters_.lx;nx_++){
        for(int ny_=0;ny_<Parameters_.ly;ny_++){
            k_ind = Coordinates_.Nc(nx_,ny_);

            Sq_dyn[nx_][ny_]=zero_complex;
            for(int wi=0;wi<S_qw[0].size();wi++){
                for(int type=0;type<3;type++){
                    Sq_dyn[nx_][ny_] += (S_qw[k_ind + (type*Parameters_.ns)][wi]*dw);
                }

            }
        }

    }


    //Calculate Sq_static from A(q)'s ;;HERE
    Aq.resize(3*Parameters_.ns);
    delAq.resize(3*Parameters_.ns);
    Aq_avg.resize(3*Parameters_.ns);
    Aq2_avg.resize(3*Parameters_.ns);

    for(int i=0;i<3*Parameters_.ns;i++){
        Aq_avg[i]=complex<double>(0,0);
        Aq2_avg[i]=complex<double>(0,0);
    }


    int nx_temp, ny_temp, k_index_temp, k_index, wi_temp;
    double double_temp, temp1, temp2;


    for(int ms=0;ms<No_Of_Inputs;ms++){

        ifstream file_Aq_in(Aq_inputs[ms].c_str());
        getline(file_Aq_in, line_temp);

        //file_out_Aq<<"#nx   ny   k_ind     Aq.real()    Aq.imag()"<<endl;

        for(int nx=0;nx<Parameters_.lx;nx++){
            for(int ny=0;ny<Parameters_.ly;ny++){
                k_index=Coordinates_.Nc(nx,ny);

                getline(file_Aq_in, line_temp);
                stringstream line_temp_ss(line_temp, stringstream::in);
                line_temp_ss>>nx_temp>>ny_temp>>k_index_temp;
                assert((nx_temp==nx) && (ny_temp==ny) && (k_index==k_index_temp));
                // file_out_Aq<<nx<<"   "<<ny<<"   "<<k_index<<"   ";


                for(int type=0;type<3;type++){
                    //file_out_Aq<<Aq[k_index + (type*Parameters_.ns)].real()<<"   "<<Aq[k_index + (type*Parameters_.ns)].imag()<<"    ";

                    line_temp_ss>>temp1>>temp2;
                    //Aq[k_index + (type*Parameters_.ns)]=complex<double> (temp1, temp2);
                    Aq2_avg[k_index + (type*Parameters_.ns)] +=complex<double> (temp1, temp2)*conj(complex<double> (temp1, temp2));
                    Aq_avg[k_index + (type*Parameters_.ns)] +=complex<double> (temp1, temp2);

                }
                //file_out_Aq<<endl;

            }
            //file_out_Aq<<endl;
            getline(file_Aq_in, line_temp);
        }
    }


    for(int i=0;i<3*Parameters_.ns;i++){
        Aq_avg[i]=Aq_avg[i]*(1.0/(No_Of_Inputs*1.0));
        Aq2_avg[i]=Aq2_avg[i]*(1.0/(No_Of_Inputs*1.0));
    }


    for(int nx_=0;nx_<Parameters_.lx;nx_++){
        for(int ny_=0;ny_<Parameters_.ly;ny_++){
            k_index=Coordinates_.Nc(nx_,ny_);

            Sq_static[nx_][ny_]=zero_complex;
            for(int type=0;type<3;type++){

                Sq_static[nx_][ny_] += Aq2_avg[k_index + (type*Parameters_.ns)]
                        - (Aq_avg[k_index + (type*Parameters_.ns)]*conj(Aq_avg[k_index + (type*Parameters_.ns)]));

            }
        }
    }



    ofstream DynamicalSq_out(Sq_dynamical_out.c_str());
    ofstream StaticSq_out(Sq_static_out.c_str());
    DynamicalSq_out<<"#nx    ny    k_ind   S(q).real  Sq.imag"<<endl;
    StaticSq_out<<"#nx    ny    k_ind   S(q).real  Sq.imag"<<endl;

    for(int nx_=0;nx_<Parameters_.lx;nx_++){
        for(int ny_=0;ny_<Parameters_.ly;ny_++){
            k_index=Coordinates_.Nc(nx_,ny_);

            DynamicalSq_out << nx_<<"     "<<ny_<<"     "<<k_index<<"     "<<Sq_dyn[nx_][ny_].real()<<"    "<<Sq_dyn[nx_][ny_].imag()<<endl;
            StaticSq_out << nx_<<"     "<<ny_<<"     "<<k_index<<"     "<<Sq_static[nx_][ny_].real()<<"    "<<Sq_static[nx_][ny_].imag()<<endl;

        }
        DynamicalSq_out<<endl;
        StaticSq_out<<endl;
    }



}























#endif


