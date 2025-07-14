#include <stdlib.h>
#include <bits/stdc++.h>
#include "Spin_dynamics_VNE_GenericHamilSpinsO3.h"
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;





void SC_SW_ENGINE_VNE_GenericHamilSpinsO3::Read_parameters(string filename){

    //Here
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
    string Conf_Initialize_ = "conf_initialize = ";
    string AnsatzType_= "AnsatzType = ";
    string conf_seed_, Conf_Seed_ = "conf_seed = ";


    string restart_, Restart_ = "Restart = ";
    string restart_time_, Restart_Time_ = "Restart_time = ";
    string restart_classical_file_, Restart_Classical_File_ = "Restart_Classical_angles = ";


    string predictor_corrector_ ,Predictor_Corrector_ = "Predictor_Corrector = ";
    string use_fft_, Use_FFT_ = "Use_FFT = ";

    string save_final_time_, Save_Final_Time_ = "Save_Final_Time_Results = ";
    string Save_Classical_File_ = "Save_Classical_angles = ";


    string insitu_spacetimefourier_, Insitu_SpaceTimeFourier_ = "Insitu_SpaceTimeFourier = ";



    int offset;
    string line;
    ifstream inputfile(filename.c_str());


    if(inputfile.is_open())
    {
        while(!inputfile.eof())
        {
            getline(inputfile,line);


            if ((offset = line.find(Insitu_SpaceTimeFourier_, 0)) != string::npos) {
                insitu_spacetimefourier_ = line.substr (offset + Insitu_SpaceTimeFourier_.length());		}

            if ((offset = line.find(Restart_Time_, 0)) != string::npos) {
                restart_time_ = line.substr (offset + Restart_Time_.length());		}


            if ((offset = line.find(Restart_Classical_File_, 0)) != string::npos) {
                restart_classical_file_ = line.substr (offset + Restart_Classical_File_.length());		}

            if ((offset = line.find(Save_Final_Time_, 0)) != string::npos) {
                save_final_time_= line.substr (offset + Save_Final_Time_.length());		}


            if ((offset = line.find(Predictor_Corrector_, 0)) != string::npos) {
                predictor_corrector_= line.substr (offset + Predictor_Corrector_.length());		}


            if ((offset = line.find(Use_FFT_, 0)) != string::npos) {
                use_fft_= line.substr (offset + Use_FFT_.length());		}


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

            if ((offset = line.find(Conf_Initialize_, 0)) != string::npos) {
                conf_initialize = line.substr (offset + Conf_Initialize_.length());		}

            if ((offset = line.find(AnsatzType_, 0)) != string::npos) {
                AnsatzType = line.substr (offset + AnsatzType_.length());		}

            if ((offset = line.find(Spins_R_T_Out_, 0)) != string::npos) {
                spins_r_t_out = line.substr (offset + Spins_R_T_Out_.length());		}

            if ((offset = line.find(SKW_Out_, 0)) != string::npos) {
                Skw_out = line.substr (offset + SKW_Out_.length());		}

            if ((offset = line.find(SKW_Out_Full_, 0)) != string::npos) {
                Skw_out_full = line.substr (offset + SKW_Out_Full_.length());		}

            if ((offset = line.find(No_Of_Processors_, 0)) != string::npos) {
                no_of_processors_ = line.substr (offset + No_Of_Processors_.length());		}

            if ((offset = line.find(Conf_Seed_, 0)) != string::npos) {
                conf_seed_ = line.substr (offset + Conf_Seed_.length());		}

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
    conf_seed = atoi(conf_seed_.c_str());


    ns =  Parameters_.ns;
    n_Spins_=1;


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




    if(RESTART){
        cout<<"This is a Restart from time = "<<Restart_Time<<endl;
        cout<<"Files "<< Restart_DM_file<<", "<<Restart_Classical_file<<" are used"<<endl;
    }
    {
        cout <<"This starts from equibrium configuration 0 time using file "<<conf_input<<endl;
    }

    cout<<"PARAMETERS::::::::"<<endl;

    cout<<"No. of sites = "<<ns<<endl;
    cout<<"w_min = "<<w_min<<endl;
    cout<<"w_max = "<<w_max<<endl;
    cout<<"dw = "<<dw<<endl;
    cout<<"w_convolution = "<<w_conv<<endl;
    cout<<"time_max = "<<time_max<<endl;
    cout<<"dt = "<<dt_<<endl;
    cout<<"Runge_Kutta_order = "<<Runge_Kutta_order<<endl;

    if(SAVE){
        cout<<"Final Results for time = "<<time_max<<endl;
        cout<<"will be written in "<< Save_DM_file<<", "<<Save_Classical_file<<endl;
    }



}




void SC_SW_ENGINE_VNE_GenericHamilSpinsO3::Create_Connections_wrt_site(){

    LocalConnectionSize.resize(ns);
    LocalConnectionSites.resize(ns);
    LocalConnectionOprs.resize(ns);
    LocalConnectionValue.resize(ns);


    string temp_opr_str;
    Mat_1_string oprs_list;
    Mat_1_int oprs_site;

    for(int FileNo=0;FileNo<Parameters_.ConnectionFiles.size();FileNo++){
        for(int connection_no=0;connection_no<Parameters_.Connections[FileNo].size();connection_no++){
            // cout<<"here 2"<<endl;
            // connection_stream.str("");

            stringstream connection_stream;
            connection_stream<<Parameters_.Connections[FileNo][connection_no];

            int n_oprs;
            oprs_list.clear();
            oprs_site.clear();


            double connection_val;
            connection_stream>>n_oprs;
            oprs_list.resize(n_oprs);
            oprs_site.resize(n_oprs);

            for(int opr_no=(n_oprs-1);opr_no>=0;opr_no--){
                connection_stream>>temp_opr_str;
                oprs_list[opr_no]=temp_opr_str;
            }
            for(int opr_no=(n_oprs-1);opr_no>=0;opr_no--){
                connection_stream>>oprs_site[opr_no];
            }
            connection_stream>>connection_val;


            for(int site=0;site<ns;site++){

                bool site_present=false;
                for(int site_ind=0;site_ind<oprs_site.size();site_ind++){
                    if(site==oprs_site[site_ind]){
                        site_present=true;
                    }
                }

                if(site_present){
                    LocalConnectionSize[site].push_back(oprs_site.size());
                    LocalConnectionValue[site].push_back(connection_val);
                    LocalConnectionSites[site].push_back(oprs_site);
                    LocalConnectionOprs[site].push_back(oprs_list);
                }

            }



        }}



}



void SC_SW_ENGINE_VNE_GenericHamilSpinsO3::Initialize_engine(){


    Create_Connections_wrt_site();

    TIME_STEP_GLOBAL=0;

    n_wpoints=(int) ((w_max-w_min)/dw + 0.5);
    if(!RESTART){
        Restart_Time=0.0;
    }


    time_steps=(int) ( (time_max - Restart_Time)/(fabs(dt_)) + 0.5 );


    //assert(false);

    Theta.resize(ns);
    Phi.resize(ns);
    Moment_Size.resize(ns);

    Theta_time.resize(ns);
    Phi_time.resize(ns);


    S_kw.resize(Parameters_.ns);
    for(int i=0;i<Parameters_.ns;i++){
        S_kw[i].resize(n_wpoints);
    }



    Aux_S_x_eq.resize(Parameters_.ns);
    Aux_S_y_eq.resize(Parameters_.ns);
    Aux_S_z_eq.resize(Parameters_.ns);
    Aux_S_x.resize(Parameters_.ns);
    Aux_S_y.resize(Parameters_.ns);
    Aux_S_z.resize(Parameters_.ns);




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

double SC_SW_ENGINE_VNE_GenericHamilSpinsO3::random1()
{

    return dis1_(Generator1_);
}



void SC_SW_ENGINE_VNE_GenericHamilSpinsO3::Set_Initial_configuration_using_Ansatz(){

    //AFM on square lattice
    if(AnsatzType=="S1_AFM"){
    for(int sitex=0;sitex<Parameters_.lx;sitex++){
        for(int sitey=0;sitey<Parameters_.ly;sitey++){
            int site = sitex + sitey*Parameters_.lx;
            Theta[site]=PI*0.5*(1.0 +  pow(-1.0,sitex+sitey)) + 0.01*PI*random1();
            Phi[site]=0.01*random1()*PI;
            Moment_Size[site]=1.0;
        }}
    }

    if(AnsatzType=="S1_AFM_Tau1by2_AFO"){
        for(int sitex=0;sitex<Parameters_.lx;sitex++){
            for(int sitey=0;sitey<Parameters_.ly;sitey++){
                int site = sitex + sitey*Parameters_.lx;
                Theta[site]=PI*0.5*(1.0 +  pow(-1.0,sitex+sitey)) + 0.05*PI*random1();
                Phi[site]=0.05*random1()*PI;
                Moment_Size[site]=1.0;
            }}
        for(int sitex=0;sitex<Parameters_.lx;sitex++){
            for(int sitey=0;sitey<Parameters_.ly;sitey++){
                int site = sitex + sitey*Parameters_.lx;
                Theta[site+Parameters_.lx*Parameters_.ly]=PI*0.5*(1.0 +  pow(-1.0,sitex+sitey)) + 0.05*PI*random1();
                Phi[site+Parameters_.lx*Parameters_.ly]=0.05*random1()*PI;
                Moment_Size[site+Parameters_.lx*Parameters_.ly]=0.5;
            }}


    }


    if(AnsatzType=="S1_AFM_Tau1by2_OD"){
        for(int sitex=0;sitex<Parameters_.lx;sitex++){
            for(int sitey=0;sitey<Parameters_.ly;sitey++){
                int site = sitex + sitey*Parameters_.lx;
                Theta[site]=PI*0.5*(1.0 +  pow(-1.0,sitex+sitey)) + 0.05*PI*random1();
                Phi[site]=0.05*random1()*PI;
                Moment_Size[site]=1.0;
            }}
        for(int sitex=0;sitex<Parameters_.lx;sitex++){
            for(int sitey=0;sitey<Parameters_.ly;sitey++){
                int site = sitex + sitey*Parameters_.lx;
                // Theta[site+Parameters_.lx*Parameters_.ly]=PI*0.5*(1.0 +  pow(-1.0,sitex+sitey)) + 1.0*PI*random1();
                // Phi[site+Parameters_.lx*Parameters_.ly]=1.0*random1()*PI;
                Theta[site+Parameters_.lx*Parameters_.ly]=1.0*PI*random1();
                Phi[site+Parameters_.lx*Parameters_.ly]=2.0*random1()*PI;
                Moment_Size[site+Parameters_.lx*Parameters_.ly]=0.5;
            }}


    }

    cout<<"t=0 is Ansatz + small_noXrandom configuration"<<endl;

}


void SC_SW_ENGINE_VNE_GenericHamilSpinsO3::Set_Initial_configuration(){

    //AFM on square lattice
    for(int sitex=0;sitex<Parameters_.lx;sitex++){
        for(int sitey=0;sitey<Parameters_.ly;sitey++){
            int site = sitex + sitey*Parameters_.lx;
        Theta[site]=PI*0.5*(1.0 +  pow(-1.0,sitex+sitey)) + 0.01*PI*random1();
        Phi[site]=0.01*random1()*PI;
        Moment_Size[site]=1.0;
        }}

    cout<<"t=0 is random configuration"<<endl;

}

void SC_SW_ENGINE_VNE_GenericHamilSpinsO3::Read_equilibrium_configuration(){


    ifstream file_in;
    file_in.open(conf_input.c_str());

    double temp_theta,temp_phi, temp_moment_size;
    int temp_site;
    string temp_string;
    getline(file_in, temp_string);


    for(int site=0;site<ns;site++){
        file_in>>temp_site;
        file_in>>temp_theta;
        file_in>>temp_phi;
        file_in>>temp_moment_size;
        Theta[temp_site]=temp_theta;
        Phi[temp_site]=temp_phi;
        Moment_Size[temp_site]=temp_moment_size;
        cout<<temp_site<<"  "<<temp_theta<<"   "<<temp_phi<<"   "<<temp_moment_size<<endl;
    }


    cout<<"Inital classical spin configuration is read from given input file "<<conf_input<<endl;


}


void SC_SW_ENGINE_VNE_GenericHamilSpinsO3::Map_Variables_to_Y(Mat_1_doub & AuxSx, Mat_1_doub & AuxSy, Mat_1_doub & AuxSz,  Mat_1_Complex_doub & Y_)
{


    int Y_size;
    Y_size = 3*n_Spins_*Parameters_.ns;


    Y_.resize(Y_size);
    //Convention Aux_Sx--Aux_Sy---Aux_Sz--

    int index=0;
    for(int i=0;i<Parameters_.ns;i++){
        Y_[index]=complex<double>(AuxSx[i],0.0);
        //cout<<Y_[index]<<endl;
        index++;
    }
    for(int i=0;i<Parameters_.ns;i++){
        Y_[index]=complex<double>(AuxSy[i],0.0);
        //cout<<Y_[index]<<endl;
        index++;
    }
    for(int i=0;i<Parameters_.ns;i++){
        Y_[index]=complex<double>(AuxSz[i],0.0);
        //cout<<Y_[index]<<endl;
        index++;
    }

    assert(index==Y_size);

}


void SC_SW_ENGINE_VNE_GenericHamilSpinsO3::Map_Y_to_Variables(Mat_1_Complex_doub & Y_, Mat_1_doub & AuxSx,
                                                              Mat_1_doub & AuxSy, Mat_1_doub & AuxSz)
{

    int Y_size;
    Y_size = 3*n_Spins_*Parameters_.ns;

    AuxSx.resize(Parameters_.ns);
    AuxSy.resize(Parameters_.ns);
    AuxSz.resize(Parameters_.ns);

    //Convention Aux_Sx--Aux_Sy---Aux_Sz--


    int index=0;
    for(int i=0;i<Parameters_.ns;i++){
        AuxSx[i]=Y_[index].real();
        index++;
    }
    for(int i=0;i<Parameters_.ns;i++){
        for(int spin_no=0;spin_no<n_Spins_;spin_no++){
            AuxSy[i]=Y_[index].real();
            index++;
        }
    }
    for(int i=0;i<Parameters_.ns;i++){
        AuxSz[i]=Y_[index].real();
        index++;
    }

    assert(index==Y_size);

}




void SC_SW_ENGINE_VNE_GenericHamilSpinsO3::IndexMapping_bw_Y_and_Variables()
{

    int Y_size;
    // = 3*n_Spins_*Parameters_.ns + 4*Parameters_.n_orbs*Parameters_.n_orbs*Parameters_.ns*Parameters_.ns;

    Y_size = 3*n_Spins_*Parameters_.ns;



    index_to_AuxSx.resize(Y_size);
    index_to_AuxSy.resize(Y_size);
    index_to_AuxSz.resize(Y_size);

    AuxSx_to_index.resize(Parameters_.ns);
    AuxSy_to_index.resize(Parameters_.ns);
    AuxSz_to_index.resize(Parameters_.ns);


    //Convention Aux_Sx--Aux_Sy---Aux_Sz-- Red_Den_mat[][][][]

    int index=0;
    for(int i=0;i<Parameters_.ns;i++){
        AuxSx_to_index[i] = index;
        index_to_AuxSx[index]=i;
        index++;
    }
    for(int i=0;i<Parameters_.ns;i++){
        AuxSy_to_index[i] = index;
        index_to_AuxSy[index]=i;
        index++;
    }
    for(int i=0;i<Parameters_.ns;i++){
        AuxSz_to_index[i] = index;
        index_to_AuxSz[index]=i;
        index++;
    }

    assert(index==Y_size);

}




complex<double> SC_SW_ENGINE_VNE_GenericHamilSpinsO3::GetLocalOprExp(string opr_str, int opr_site, Mat_1_Complex_doub & Y_){

    complex<double> value;

        if(opr_str == "Sz"){
            value = Y_[AuxSz_to_index[opr_site]].real();
        }
        else if(opr_str == "Sx"){
            value = Y_[AuxSx_to_index[opr_site]].real();
        }
        else if(opr_str == "Sy"){
            value = Y_[AuxSy_to_index[opr_site]].real();
        }
        else{
            cout<<"used opr = \""<<opr_str<<"\""<<endl;
            cout<<"Only Sz, Sx, Sy oprs allowed right now"<<endl;
            assert(false);
        }

    return value;

}

double SC_SW_ENGINE_VNE_GenericHamilSpinsO3::Calculate_TotalE(Mat_1_Complex_doub & Y_){

    // string outfiletempstr="Ebond.txt";
    // ofstream outfiletemp(outfiletempstr.c_str());

    complex<double> temp_TotalE_=0.0;
    double TotalE_=0.0;

    for(int FileNo=0;FileNo<Parameters_.ConnectionFiles.size();FileNo++){

        Mat_1_Complex_doub E_array;
        int Np_=1;

#ifdef _OPENMP
        Np_= omp_get_max_threads();
#endif

        E_array.resize(Np_);

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
        for(int connection_no=0;connection_no<Parameters_.Connections[FileNo].size();connection_no++){

            int thread_id=0;
#ifdef _OPENMP
            thread_id = omp_get_thread_num();
#endif


            stringstream connection_stream;
            connection_stream<<Parameters_.Connections[FileNo][connection_no];

            int n_oprs;
            string temp_opr_str;
            Mat_1_string oprs_list;
            Mat_1_int oprs_site;
            oprs_list.clear();
            oprs_site.clear();

            double connection_val;
            connection_stream>>n_oprs;
            oprs_list.resize(n_oprs);
            oprs_site.resize(n_oprs);

            for(int opr_no=(n_oprs-1);opr_no>=0;opr_no--){
                connection_stream>>temp_opr_str;
                oprs_list[opr_no]=temp_opr_str;
            }
            for(int opr_no=(n_oprs-1);opr_no>=0;opr_no--){
                connection_stream>>oprs_site[opr_no];
            }
            connection_stream>>connection_val;


            complex<double> E_conn=1.0;
            for(int opr_no=0;opr_no<oprs_list.size();opr_no++){
                int opr_site = oprs_site[opr_no];
                string opr_str = oprs_list[opr_no];
                E_conn = E_conn*GetLocalOprExp(opr_str, opr_site, Y_);
                //outfiletemp<<"("<<opr_str<<","<<opr_site<<")="<<GetLocalOprExp(opr_str, opr_site, Y_)<<"  ";
            }

           // outfiletemp<<E_conn.real()<<"  "<<connection_val<<endl;

            E_conn = E_conn*connection_val;

            E_array[thread_id] += E_conn;
        }

        for(int th_id=0;th_id<Np_;th_id++){
            temp_TotalE_ +=E_array[th_id];
        }

    }



    if(temp_TotalE_.imag()>0.000001){
        cout<<"temp_TotalE_.real() : "<<temp_TotalE_.real()<<endl;
        cout<<"temp_TotalE_.imag() : "<<temp_TotalE_.imag()<<endl;
        assert(temp_TotalE_.imag()<0.000001);
    }

    TotalE_ =  temp_TotalE_.real();

    return TotalE_;

}




void SC_SW_ENGINE_VNE_GenericHamilSpinsO3::Start_Engine(){

#ifdef _OPENMP
    omp_set_num_threads(no_of_processors);
    int N_p = omp_get_max_threads();
    cout<<"Max threads which can be used parallely = "<<N_p<<endl;
    cout<<"No. of threads used parallely = "<<no_of_processors<<endl;
#endif



    string Energy_out = "ClassicalEnergy_vs_t.txt";
    ofstream Energy_file_out(Energy_out.c_str());
    Energy_file_out<<"#time   ClassicalEnergy"<<endl;


    ofstream file_out;
    string list_outfiles;
    list_outfiles = "Spin_"+spins_r_t_out;
    file_out.open(list_outfiles);
    file_out<<"# Sz------, Sx-----,Sy-----"<<endl;
    file_out<<scientific<<setprecision(15);


    for(int ts=0;ts<=time_steps;ts++){

        TIME_STEP_GLOBAL=ts;

        if(SAVE && (ts==time_steps)){
            Write_final_time_result();
        }


        //create Hamiltonian for time_step=ts
        if(ts==0){

            if(!RESTART){
                cout<<"Starting from Equilibrium condition"<<endl;
            }
            else{
                cout<<"Restart is done"<<endl;
                //Read_Restart_Data();
            }


            for(int pos=0;pos<Parameters_.ns;pos++){
                Aux_S_x[pos] = Moment_Size[pos]*sin(Theta[pos])*cos(Phi[pos]);
                Aux_S_y[pos] = Moment_Size[pos]*sin(Theta[pos])*sin(Phi[pos]);
                Aux_S_z[pos] = Moment_Size[pos]*cos(Theta[pos]);
            }

            Map_Variables_to_Y(Aux_S_x, Aux_S_y, Aux_S_z, YVec0);
        }



        Parameters_.PrintingNoOfTimeSlices=time_steps-1;
        int ts_gap = (time_steps-1)/Parameters_.PrintingNoOfTimeSlices;
        if((ts==0) || (ts==time_steps) ||
            ((ts%ts_gap)==0)
            ){
            for(int spin_no=0;spin_no<n_Spins_;spin_no++){
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
            }
        }


        Energy_file_out<<scientific<<setprecision(15)<<ts*dt_<<"   "<<Calculate_TotalE(YVec0)<<endl;
        Evolve_classical_spins_Runge_Kutta(0);
        YVec0=YVec1;
        //Update_Hamiltonian_Classical_dof();
    }



}



void SC_SW_ENGINE_VNE_GenericHamilSpinsO3::Evolve_classical_spins_Runge_Kutta(int ts){


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
    else if(Runge_Kutta_order==8){
        RungeKuttaEight(YVec0, YVec1);
    }


}



void SC_SW_ENGINE_VNE_GenericHamilSpinsO3::Derivative(Mat_1_Complex_doub & Y_, Mat_1_Complex_doub & dYbydt){



    dYbydt.resize(Y_.size());
    for(int i=0;i<Y_.size();i++){
        dYbydt[i]=zero_complex;
    }


#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
    for(int pos=0;pos<Parameters_.ns;pos++){

        for(int term_no=0;term_no<LocalConnectionSize[pos].size();term_no++){

            complex<double> coeff_temp=LocalConnectionValue[pos][term_no];
            complex<double> coeff_x, coeff_y, coeff_z;
            int pos_times_z=0;
            int pos_times_x=0;
            int pos_times_y=0;

            for(int opr_no=0;opr_no<LocalConnectionSize[pos][term_no];opr_no++){
                if(LocalConnectionSites[pos][term_no][opr_no]==pos){
                    if(LocalConnectionOprs[pos][term_no][opr_no]=="Sx"){
                    pos_times_x +=1;
                    }
                    if(LocalConnectionOprs[pos][term_no][opr_no]=="Sy"){
                        pos_times_y +=1;
                    }
                    if(LocalConnectionOprs[pos][term_no][opr_no]=="Sz"){
                        pos_times_z +=1;
                    }
                }
                else{
                    complex<double> opr_value = GetLocalOprExp(LocalConnectionOprs[pos][term_no][opr_no], LocalConnectionSites[pos][term_no][opr_no], Y_);
                    coeff_temp = coeff_temp*opr_value;
                }
            }

            double Sx_pos, Sy_pos, Sz_pos;
            Sx_pos = (GetLocalOprExp("Sx", pos, Y_)).real();
            Sy_pos = (GetLocalOprExp("Sy", pos, Y_)).real();
            Sz_pos = (GetLocalOprExp("Sz", pos, Y_)).real();

            if(pos_times_x>0){
                coeff_x =  coeff_temp*(1.0*pos_times_x)*pow(Sx_pos,pos_times_x-1);
                coeff_x = coeff_x*pow(Sy_pos,pos_times_y)*pow(Sz_pos,pos_times_z);

                //zxy - yxz
                dYbydt[AuxSz_to_index[pos]] += 1.0*coeff_x*Sy_pos;
                dYbydt[AuxSy_to_index[pos]] -= 1.0*coeff_x*Sz_pos;
            }
            if(pos_times_y>0){
                coeff_y =  coeff_temp*(1.0*pos_times_y)*pow(Sy_pos,pos_times_y-1);
                coeff_y = coeff_y*pow(Sx_pos,pos_times_x)*pow(Sz_pos,pos_times_z);

                //xyz - zyx
                dYbydt[AuxSx_to_index[pos]] += 1.0*coeff_y*Sz_pos;
                dYbydt[AuxSz_to_index[pos]] -= 1.0*coeff_y*Sx_pos;
            }
            if(pos_times_z>0){
                coeff_z =  coeff_temp*(1.0*pos_times_z)*pow(Sz_pos,pos_times_z-1);
                coeff_z = coeff_z*pow(Sx_pos,pos_times_x)*pow(Sy_pos,pos_times_y);

                //yzx - xzy
                dYbydt[AuxSy_to_index[pos]] += 1.0*coeff_z*Sx_pos;
                dYbydt[AuxSx_to_index[pos]] -= 1.0*coeff_z*Sy_pos;
            }


        }

    }



}



void SC_SW_ENGINE_VNE_GenericHamilSpinsO3::RungeKuttaOne(Mat_1_Complex_doub & Yn, Mat_1_Complex_doub & Ynp1){


    Ynp1.resize(Yn.size());
    Mat_1_Complex_doub K1;

    //step_no==0
    Derivative(Yn,K1);

    for(int i=0;i<Yn.size();i++){
        Ynp1[i] = Yn[i] + (dt_)*(K1[i]);

    }


}

void SC_SW_ENGINE_VNE_GenericHamilSpinsO3::RungeKuttaFour(Mat_1_Complex_doub & Yn, Mat_1_Complex_doub & Ynp1){



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


void SC_SW_ENGINE_VNE_GenericHamilSpinsO3::RungeKuttaSix(Mat_1_Complex_doub & Yn, Mat_1_Complex_doub & Ynp1){



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




void SC_SW_ENGINE_VNE_GenericHamilSpinsO3::RungeKuttaEight(Mat_1_Complex_doub & Yn, Mat_1_Complex_doub & Ynp1){


    //Shanks Eighth order Runge Kutta

    Mat_1_Complex_doub Y_temp;
    Y_temp.resize(Yn.size());

    Ynp1.resize(Yn.size());
    Mat_1_Complex_doub K0, K1, K2, K3, K4, K5, K6, K7, K8, K9;


    //step_no==0
    Derivative(Yn,K0);

    //PARALELLIZE THIS
    for(int i=0;i<Yn.size();i++){
        Y_temp[i] = Yn[i]  + ((4.0/27.0)*K0[i]*dt_);
    }

    //step_no==1
    Derivative(Y_temp,K1);
    for(int i=0;i<Yn.size();i++){
        Y_temp[i] = Yn[i]  + (K0[i]*dt_*(1.0/18.0)) + (K1[i]*dt_*(3.0/18.0));
    }

    //step_no==2
    Derivative(Y_temp,K2);
    for(int i=0;i<Yn.size();i++){
        Y_temp[i] = Yn[i]  + (K0[i]*dt_*(1.0/12.0)) + (K1[i]*dt_*(0.0)) + (K2[i]*dt_*(3.0/12.0));
    }

    //step_no==3
    Derivative(Y_temp,K3);
    for(int i=0;i<Yn.size();i++){
        Y_temp[i] = Yn[i]  + (K0[i]*dt_*(1.0/8.0)) + (K1[i]*dt_*(0.0)) + (K2[i]*dt_*(0.0)) + (K3[i]*dt_*(3.0/8.0)) ;
    }


    //step_no==4
    Derivative(Y_temp,K4);
    for(int i=0;i<Yn.size();i++){
        Y_temp[i] = Yn[i]  + (K0[i]*dt_*(13.0/54.0)) + (K1[i]*dt_*(0.0)) + (K2[i]*dt_*(-27.0/54.0)) + (K3[i]*dt_*(42.0/54.0)) + (K4[i]*dt_*(8.0/54.0)) ;
    }

    //step_no==5
    Derivative(Y_temp,K5);
    for(int i=0;i<Yn.size();i++){
        Y_temp[i] = Yn[i]  + (K0[i]*dt_*(389.0/4320.0)) + (K1[i]*dt_*(0.0)) + (K2[i]*dt_*(-54.0/4320.0)) + (K3[i]*dt_*(966.0/4320.0)) + (K4[i]*dt_*(-824.0/4320.0)) + (K5[i]*dt_*(243.0/4320.0)) ;
    }


    //step_no==6
    Derivative(Y_temp,K6);
    for(int i=0;i<Yn.size();i++){
        Y_temp[i] = Yn[i]  + (K0[i]*dt_*(-231.0/20.0)) + (K1[i]*dt_*(0.0)) + (K2[i]*dt_*(81.0/20.0)) + (K3[i]*dt_*(-1164.0/20.0)) + (K4[i]*dt_*(656.0/20.0)) + (K5[i]*dt_*(-122.0/20.0)) + (K6[i]*dt_*(800.0/20.0)) ;
    }


    //step_no==7
    Derivative(Y_temp,K7);
    for(int i=0;i<Yn.size();i++){
        Y_temp[i] = Yn[i]  + (K0[i]*dt_*(-127.0/288.0)) + (K1[i]*dt_*(0.0)) + (K2[i]*dt_*(18.0/288.0)) + (K3[i]*dt_*(-678.0/288.0)) + (K4[i]*dt_*(456.0/288.0)) + (K5[i]*dt_*(-9.0/288.0)) + (K6[i]*dt_*(576.0/288.0)) + (K7[i]*dt_*(4.0/288.0) ) ;
    }


    //step_no==8
    Derivative(Y_temp,K8);
    for(int i=0;i<Yn.size();i++){
        Y_temp[i] = Yn[i]  + (K0[i]*dt_*(1481.0/820.0)) + (K1[i]*dt_*(0.0)) + (K2[i]*dt_*(-81.0/820.0)) + (K3[i]*dt_*(7104.0/820.0)) + (K4[i]*dt_*(-3376.0/820.0)) + (K5[i]*dt_*(72.0/820.0)) + (K6[i]*dt_*(-5040.0/820.0)) + (K7[i]*dt_*(-60.0/820.0)) + (K8[i]*dt_*(720.0/820.0) ) ;
    }

    //step_no==9
    Derivative(Y_temp,K9);


    for(int i=0;i<Yn.size();i++){
        Ynp1[i] = Yn[i] + (dt_*(1.0/840))*( (41.0)*K0[i] +  (0.0)*K1[i] + (0.0)*K2[i] + (27.0)*K3[i] + (272.0)*K4[i] + (27.0)*K5[i] + (216.0)*K6[i] + (0.0)*K7[i] + (216.0)*K8[i] + (41.0)*K9[i]);

    }


}


void SC_SW_ENGINE_VNE_GenericHamilSpinsO3::Write_final_time_result(){


    ofstream file_out_Classical(Save_Classical_file.c_str());

    for(int p=0;p<ns;p++){
        file_out_Classical<<p<<"   "<<Theta[p]<<"    "<<Phi[p]<<endl;
    }


}



void SC_SW_ENGINE_VNE_GenericHamilSpinsO3::Read_Restart_Data(){

    complex<double> one(1,0);
    complex<double> iota(0,1);


    ifstream file_in_Classical(Restart_Classical_file.c_str());

    int temp_p;
    double temp_theta, temp_phi;
    for(int p=0;p<Parameters_.ns;p++){
                file_in_Classical>>temp_p>>temp_theta>>temp_phi;
                Theta[temp_p]=temp_theta;
                Phi[temp_p]=temp_phi;
    }


}









