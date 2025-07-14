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
#include "../../Matrix.h"
#include "ParametersEngine_GenericHamilSpinsO3.h"
#include "Spin_dynamics_VNE_GenericHamilSpinsO3.h"
#include "../../tensor_type.h"
#include "../../FftComplex.h"
#define PI acos(-1.0)

using namespace std;

#ifndef ST_Fourier_GenericHamilSpinsO3_
#define ST_Fourier_GenericHamilSpinsO3_

class ST_Fourier_GenericHamilSpinsO3{

public:
    ST_Fourier_GenericHamilSpinsO3(Parameters_GenericHamilSpinsO3& Parameters__, SC_SW_ENGINE_VNE_GenericHamilSpinsO3& SC_SW_ENGINE_VNE_GenericHamilSpinsO3__)
        : Parameters_(Parameters__), SC_SW_ENGINE_VNE_GenericHamilSpinsO3_(SC_SW_ENGINE_VNE_GenericHamilSpinsO3__)
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


    int lx, ly;
    int SitesOffset;

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

    Mat_3_Complex_doub S_tr_2point;
    Mat_2_Complex_doub S_tr_1point;
    Mat_2_doub S_tr;
    Mat_2_doub C_tr_; //Space-time dispaced correlation, avg over Ensemble
    Mat_2_doub C_Quantum_tr_;
    Mat_2_doub C_Classical_tr_;
    Mat_3_doub C_tr; //Space-time dispaced correlation, avg over Ensemble
    Mat_3_doub C_Quantum_tr;
    Mat_3_doub C_Classical_tr;

    Mat_3_Complex_doub F_rw_2point;
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

    Parameters_GenericHamilSpinsO3& Parameters_;
    SC_SW_ENGINE_VNE_GenericHamilSpinsO3& SC_SW_ENGINE_VNE_GenericHamilSpinsO3_;



    void Initialize_engine(string inputfile_);
    void Read_parameters();

    void Calculate_Fw_and_Aq(string fileout, string fileout_Aq);

    double matchstring(string file, string match);
    string matchstring2(string file, string match);

};





double ST_Fourier_GenericHamilSpinsO3::matchstring(string file, string match)
{
    string test;
    string line;
    ifstream readFile(file);
    double amount;
    bool pass = false;
    while (std::getline(readFile, line))
    {
        std::istringstream iss(line);
        if (std::getline(iss, test, '=') && pass == false)
        {
            // ---------------------------------
            if (iss >> amount && test == match)
            {
                // cout << amount << endl;
                pass = true;
            }
            else
            {
                pass = false;
            }
            // ---------------------------------
            if (pass)
                break;
        }
    }
    if (pass == false)
    {
        string errorout = match;
        errorout += "= argument is missing in the input file!";
        throw std::invalid_argument(errorout);
    }
    cout << match << " = " << amount << endl;
    return amount;
}

string ST_Fourier_GenericHamilSpinsO3::matchstring2(string file, string match)
{

    string line;
    ifstream readFile(file);
    string amount;
    int offset;

    if (readFile.is_open())
    {
        while (!readFile.eof())
        {
            getline(readFile, line);

            if ((offset = line.find(match, 0)) != string::npos)
            {
                amount = line.substr(offset + match.length() + 1);
            }
        }
        readFile.close();
    }
    else
    {
        cout << "Unable to open input file while in the Parameters class." << endl;
    }

    cout << match << " = " << amount << endl;
    return amount;
}





void ST_Fourier_GenericHamilSpinsO3::Initialize_engine(string inputfile_){



    lx = int(matchstring(inputfile_, "lx"));
    cout << "Lx = " << lx << endl;

    ly = int(matchstring(inputfile_, "ly"));
    cout << "Ly = " << ly << endl;

    SitesOffset=int(matchstring(inputfile_,"SitesOffset"));
    cout<<"SitesOffset = "<<SitesOffset<<endl;

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


void ST_Fourier_GenericHamilSpinsO3::Calculate_Fw_and_Aq(string fileout, string fileout_Aq){

    if(!Use_FFT){
        cout<<"FFT is NOT used"<<endl;
    }
    else{
        cout<<"FFT is used"<<endl;
    }

    int GaussianCenteredAtTmaxby2=1;
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



        for(int r=0;r<3*SitesOffset;r++){
            line_temp_ss>>double_temp;
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
    for(int pos_ix=0;pos_ix<lx;pos_ix++){

        Mat_1_Complex_doub Vec_1;
        Vec_1.resize(time_steps);

        for(int pos_iy=0;pos_iy<ly;pos_iy++){

            for(int type=0;type<3;type++){

                pos_i =pos_ix + lx*pos_iy;
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
                        F_rw[index][wi] = Vec_1[wi].real();
                        //F_rw[index][wi] = Vec_1[wi];

                    }
                }

            }

        }

        vector< complex<double> >().swap( Vec_1 );

        cout<<"FT for S[t][rx="   << pos_ix <<  ", ry] is completed for all ry, and all 6 types of S"<<endl;

    }



    double kx, ky;
    int k_index;
    complex<double> temp;

    for(int type=0;type<3;type++){
        for(int nx=0;nx<lx;nx++){
            for(int ny=0;ny<ly;ny++){

                kx = (2*nx*PI)/(1.0*lx);
                ky = (2*ny*PI)/(1.0*ly);
                k_index = nx + lx*ny + (type*Parameters_.ns);

                temp = complex<double>(0,0);
                for(int x_i=0;x_i<lx;x_i++){
                    for(int y_i=0;y_i<ly;y_i++){

                        pos_i = x_i + lx*y_i;
                        index  = pos_i + (type*Parameters_.ns);
                        temp += S_tr[(GCATm2*time_steps)/2][index]*exp(iota*(1.0*( (x_i)*kx +  (y_i)*ky ) ) );
                        //temp += F_rw[index][10].real()*exp(iota*(1.0*( (x_i)*kx +  (y_i)*ky ) ) );

                    }
                }
                Aq[k_index] = temp;

            }
        }
    }




    for(int type=0;type<3;type++){
        for(int nx=0;nx<lx;nx++){
            cout<<"FT r --> k="<<nx<<", ny"<<endl;
            for(int ny=0;ny<ly;ny++){

                kx = (2*nx*PI)/(1.0*lx);
                ky = (2*ny*PI)/(1.0*ly);
                k_index = nx + lx*ny + (type*Parameters_.ns);

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(pos_i, index, temp)
#endif

                for(int wi=0;wi<n_wpoints;wi++){

                    temp = complex<double> (0,0);
                    for(int x_i=0;x_i<lx;x_i++){
                        for(int y_i=0;y_i<ly;y_i++){

                            pos_i = x_i + lx*y_i;
                            index  = pos_i + (type*Parameters_.ns);

                            temp += F_rw[index][wi]*exp(iota*(-1.0*( (x_i)*(kx) +  (y_i)*ky ) ) );

                        }
                    }
                    //file_out2<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<temp.real()<<"   "<<temp.imag()<<"    "<<temp2.real()<<"   "<<temp2.imag()<<"    "<<temp3.real()<<"   "<<temp3.imag()<<endl;
                    //file_out_full<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<wi<<"   "<<temp.real()<<"   "<<temp.imag()<<"    "<<"none"<<"   "<<"none"<<"    "<<"none"<<"   "<<"none"<<endl;
                    F_qw[k_index][wi]=temp*(1.0/(lx*ly));
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
    for(int nx=0;nx<lx;nx++){
        for(int ny=0;ny<ly;ny++){
            k_index=nx + lx*ny;

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

    for(int nx=0;nx<lx;nx++){
        for(int ny=0;ny<ly;ny++){
            k_index=nx + lx*ny;

            file_out_Aq<<nx<<"   "<<ny<<"   "<<k_index<<"   ";

            for(int type=0;type<3;type++){
                file_out_Aq<<Aq[k_index + (type*Parameters_.ns)].real()<<"   "<<Aq[k_index + (type*Parameters_.ns)].imag()<<"    ";

            }
            file_out_Aq<<endl;

        }
        //file_out_Aq<<endl;
    }


    cout<<"Calculate_Fw_and_Aq completed"<<endl;




    if(false) { //For triangular lattices (example Kagome)

        assert(lx==ly);
        //Create Path Gamma--> M---->K--->Gamma
        int n1, n2;
        Mat_1_intpair k_path;
        pair_int temp_pair;


        // ---k_path---------

        //--------\Gamma to M----------------
        n1=0;
        n2=0;
        while (n1<=int(lx/2))
        {
            temp_pair.first = n1;
            temp_pair.second = n2;
            k_path.push_back(temp_pair);
            n2++;
            n1++;
        }
        //----------------------------------

        //--------\M to K----------------
        n1=int(lx/2)+1;
        n2=int(ly/2)-1;
        while (n1<=int((2*lx)/3))
        {
            temp_pair.first = n1;
            temp_pair.second = n2;
            k_path.push_back(temp_pair);
            n2--;
            n1++;
        }
        //----------------------------------

        //--------K to \Gamma----------------
        n1=int((2*lx)/3)-2;
        n2=int((ly)/3)-1;
        while (n1>=0)
        {
            temp_pair.first = n1;
            temp_pair.second = n2;
            k_path.push_back(temp_pair);
            n2--;
            n1--;
            n1--;
        }
        //----------------------------------


	//--------\Gamma to K'----------
	n1=1;
	n2=2;
    while(n1<=int((lx)/3)){
	temp_pair.first = n1;
        temp_pair.second = n2;
        k_path.push_back(temp_pair);
	n1=n1+1;
        n2=n2+2;
	}


        temp_pair.first = 0;
        temp_pair.second = 0;
        k_path.push_back(temp_pair);
        temp_pair.first = 0;
        temp_pair.second = 0;
        k_path.push_back(temp_pair);

        //----------------------------------
        cout<<"PRINTING PATH"<<endl;
        for (int k_point = 0; k_point < k_path.size(); k_point++)
        {
            cout<<k_path[k_point].first<< "   "<<k_path[k_point].second<<endl;
        }


        //----k_path done-------


        string fileout_TL = "TL_" + fileout;
        ofstream file_out_TL(fileout_TL.c_str());

        file_out_TL<<"#w_details:  "<<"   "<<w_max<<"  "<<w_min<<"  "<<dw<<"  "<<n_wpoints<<endl;
        file_out_TL<<"#nx   ny   k_ind   wi*dw   wi   F_qw[k_ind][wi].real()   F_qw[k_ind][wi].imag()"<<endl;
        for (int k_point = 0; k_point < k_path.size(); k_point++)
        {
                int nx = k_path[k_point].first;
                int ny = k_path[k_point].second;
                k_index=nx + lx*ny;

                for(int wi=0;wi<n_wpoints;wi++){

                    //file_out2<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<temp.real()<<"   "<<temp.imag()<<"    "<<temp2.real()<<"   "<<temp2.imag()<<"    "<<temp3.real()<<"   "<<temp3.imag()<<endl;
                    file_out_TL<<nx<<"   "<<ny<<"   "<<k_point<<"   "<<wi*dw<<"   "<<wi<<"   ";

                    for(int type=0;type<3;type++){

                        file_out_TL<<F_qw[k_index + (type*Parameters_.ns)][wi].real()<<"   "<<F_qw[k_index + (type*Parameters_.ns)][wi].imag()<<"   ";
                    }

                    file_out_TL<<endl;

                }

                file_out_TL<<endl;


        }



    }



    if(true) { //For square lattices (example Kagome)

        assert(lx==ly);
        //Create Path Gamma--> M---->X--->G-->Y-->M-->G
        int n1, n2;
        Mat_1_intpair k_path;
        pair_int temp_pair;


        // ---k_path---------

        //--------[\Gamma to M]----------------
        n1=0;
        n2=0;
        while (n1<=int(lx/2))
        {
            temp_pair.first = n1;
            temp_pair.second = n2;
            k_path.push_back(temp_pair);
            n2++;
            n1++;
        }
        //----------------------------------

        //--------(M to X]----------------
        n1=int(lx/2);
        n2=int(ly/2)-1;
        while (n2>=0)
        {
            temp_pair.first = n1;
            temp_pair.second = n2;
            k_path.push_back(temp_pair);
            n2--;
        }
        //----------------------------------

        //--------(X to \Gamma]----------------
        n1=(lx/2)-1;
        n2=0;
        while (n1>=0)
        {
            temp_pair.first = n1;
            temp_pair.second = n2;
            k_path.push_back(temp_pair);
            n1--;
        }
        //----------------------------------

        //--------(G to Y]----------------
        n1=0;
        n2=1;
        while (n2<=ly/2)
        {
            temp_pair.first = n1;
            temp_pair.second = n2;
            k_path.push_back(temp_pair);
            n2++;
        }
        //----------------------------------


        //--------(Y to M]----------------
        n1=1;
        n2=ly/2;
        while (n1<=lx/2)
        {
            temp_pair.first = n1;
            temp_pair.second = n2;
            k_path.push_back(temp_pair);
            n1++;
        }
        //----------------------------------

        //--------(M to G]----------------
        n1=lx/2-1;
        n2=ly/2-1;
        while (n1>=0)
        {
            temp_pair.first = n1;
            temp_pair.second = n2;
            k_path.push_back(temp_pair);
            n1--;
            n2--;
        }
        //----------------------------------

        temp_pair.first = 0;
        temp_pair.second = 0;
        k_path.push_back(temp_pair);
        temp_pair.first = 0;
        temp_pair.second = 0;
        k_path.push_back(temp_pair);

        //----------------------------------
        cout<<"PRINTING PATH"<<endl;
        for (int k_point = 0; k_point < k_path.size(); k_point++)
        {
            cout<<k_path[k_point].first<< "   "<<k_path[k_point].second<<endl;
        }


        //----k_path done-------


        string fileout_TL = "SQL_" + fileout;
        ofstream file_out_TL(fileout_TL.c_str());

        file_out_TL<<"#w_details:  "<<"   "<<w_max<<"  "<<w_min<<"  "<<dw<<"  "<<n_wpoints<<endl;
        file_out_TL<<"#nx   ny   k_ind   wi*dw   wi   F_qw[k_ind][wi].real()   F_qw[k_ind][wi].imag()"<<endl;
        for (int k_point = 0; k_point < k_path.size(); k_point++)
        {
            int nx = k_path[k_point].first;
            int ny = k_path[k_point].second;
            k_index=nx + lx*ny;

            for(int wi=0;wi<n_wpoints;wi++){

                //file_out2<<nx<<"   "<<ny<<"   "<<k_ind<<"   "<<wi*dw<<"   "<<temp.real()<<"   "<<temp.imag()<<"    "<<temp2.real()<<"   "<<temp2.imag()<<"    "<<temp3.real()<<"   "<<temp3.imag()<<endl;
                file_out_TL<<nx<<"   "<<ny<<"   "<<k_point<<"   "<<wi*dw<<"   "<<wi<<"   ";

                for(int type=0;type<3;type++){

                    file_out_TL<<F_qw[k_index + (type*Parameters_.ns)][wi].real()<<"   "<<F_qw[k_index + (type*Parameters_.ns)][wi].imag()<<"   ";
                }

                file_out_TL<<endl;

            }

            file_out_TL<<endl;


        }



    }






}

void ST_Fourier_GenericHamilSpinsO3::Read_parameters(){


    w_min=SC_SW_ENGINE_VNE_GenericHamilSpinsO3_.w_min;
    w_max=SC_SW_ENGINE_VNE_GenericHamilSpinsO3_.w_max;

    dw=SC_SW_ENGINE_VNE_GenericHamilSpinsO3_.dw;
    w_conv=SC_SW_ENGINE_VNE_GenericHamilSpinsO3_.w_conv;
    Use_FFT=SC_SW_ENGINE_VNE_GenericHamilSpinsO3_.Use_FFT;

    no_of_processors=SC_SW_ENGINE_VNE_GenericHamilSpinsO3_.no_of_processors;


}





















#endif


