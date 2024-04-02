#include <stdlib.h>
#include <bits/stdc++.h>
#include "Spin_dynamics_Rotor.h"
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;


void SC_SW_ENGINE_Rotor::Create_Scheduler(){


Time_.resize(time_steps+1);
Gamma_.resize(time_steps+1);
Js_.resize(time_steps+1);

double time_normalized;
double gamma_temp, js_temp;
for(int time_ind=0;time_ind<=time_steps;time_ind++){

time_normalized=(time_ind*dt_)/(time_max);

//assert(time_normalized>=0 && time_normalized<=1.0);

for(int time_ind2=0;time_ind2<Time_bare.size()-1;time_ind2++){
    if(time_normalized>=Time_bare[time_ind2] && time_normalized<=Time_bare[time_ind2+1] ){
    gamma_temp = ((Time_bare[time_ind2+1]-time_normalized)*Gamma_bare[time_ind2]
                 -(Time_bare[time_ind2]-time_normalized)*Gamma_bare[time_ind2+1])*
                 (1.0/(Time_bare[time_ind2+1]-Time_bare[time_ind2]));

    js_temp = ((Time_bare[time_ind2+1]-time_normalized)*Js_bare[time_ind2]
                 -(Time_bare[time_ind2]-time_normalized)*Js_bare[time_ind2+1])*
                 (1.0/(Time_bare[time_ind2+1]-Time_bare[time_ind2]));

    break;
    }
}

Time_[time_ind]=time_ind*dt_;
Gamma_[time_ind]=gamma_temp;
Js_[time_ind]=js_temp;
}


string created_schd_str = "Created_Scheduler.txt";
ofstream created_schd_stream(created_schd_str.c_str());

created_schd_stream<<"# time   Gamma    Js"<<endl;


int ts_gap = (time_steps-1)/Parameters_.PrintingNoOfTimeSlices;

for(int time_ind=0;time_ind<Time_.size();time_ind++){
    if((time_ind==0) || (time_ind==time_steps) ||
        ((time_ind%ts_gap)==0)
            ){
created_schd_stream<<Time_[time_ind]<<"  "<<Gamma_[time_ind]<<"  "<<Js_[time_ind]<<endl;
    }
}

}

void SC_SW_ENGINE_Rotor::Initialize_engine(){


    TIME_STEP_GLOBAL=0;
    ns = Parameters_.ns;

    time_steps=(int) ( (time_max)/(fabs(dt_)) + 0.5 );

    if(Use_Scheduler){
    Create_Scheduler();
    }


    Theta.resize(ns);
    Momentum.resize(ns);


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


   Generator_.seed(Parameters_.RandomDisorderSeed);

   GeneratorNoise_.seed(Parameters_.RandomNoiseSeed);

   StdDev=sqrt(2.0*Parameters_.DampingConst*(1.0/Parameters_.beta));
   normal_distribution<double> GaussianDistribution_temp(0.0,1.0);
   GaussianDistribution=GaussianDistribution_temp;


   EPSILON=0.000001;

   double temp_val;
   temp_val=  GaussianDistribution(GeneratorNoise_);
   cout<<"temp_val : "<<temp_val<<endl;

}

double SC_SW_ENGINE_Rotor::random1()
{

    return dis1_(Generator_);
}

void SC_SW_ENGINE_Rotor::Initial_condition(){


    //Sx=cos(theta), Sz=sin(theta)
    for(int i=0;i<ns;i++){
     Momentum[i]=0.0;
     Theta[i]=0.0; //in +x direction
    }


}


void SC_SW_ENGINE_Rotor::Start_Engine(){
#ifdef _OPENMP
    omp_set_num_threads(no_of_processors);
    int N_p = omp_get_max_threads();
    cout<<"Max threads which can be used parallely = "<<N_p<<endl;
    cout<<"No. of threads used parallely = "<<no_of_processors<<endl;
#endif


    string Energy_out = "ClassicalEnergy_vs_t.txt";
    ofstream Energy_file_out(Energy_out.c_str());
    Energy_file_out<<"#time   ClassicalEnergy"<<endl;
    Energy_file_out<<scientific<<setprecision(10);


    ofstream file_out;
    string list_outfiles;
    list_outfiles = spins_r_t_out;
    file_out.open(list_outfiles);
    file_out<<"# theta------, p (momentum)----"<<endl;
    file_out<<scientific<<setprecision(15);


    for(int ts=0;ts<=time_steps;ts++){

        TIME_STEP_GLOBAL=ts;

        //create Hamiltonian for time_step=ts
        if(ts==0){
            Initial_condition();
            Map_Variables_to_Y(Theta, Momentum, YVec0);
        }

        int ts_gap = (time_steps-1)/Parameters_.PrintingNoOfTimeSlices;
        if((ts==0) || (ts==time_steps) ||
            ((ts%ts_gap)==0)
                ){
        file_out<<(ts*dt_)<<"  ";
        for(int pos=0;pos<Parameters_.ns;pos++){
                file_out<< YVec0[theta_to_index[pos]]<<"  ";
            }
        for(int pos=0;pos<Parameters_.ns;pos++){
                file_out<< YVec0[momentum_to_index[pos]]<<"  ";
        }
        file_out<<endl;
        }


        Energy_file_out<<ts*dt_<<"   "<<Get_Classical_Energy(YVec0)<<endl;
        Evolve_classical_spins_Runge_Kutta(0);
        YVec0=YVec1;
    }


    Rotor_KinkDen_type1=Get_Kink_Density_1d();
    Rotor_KinkDen_type2= Get_Kink_Density_1d_type2();
    cout<<"kink density : "<<Rotor_KinkDen_type1<<endl;
    cout<<"kink density type 2: "<<Rotor_KinkDen_type2<<endl;


}





double SC_SW_ENGINE_Rotor::Get_Kink_Density_1d_type2(){
double kink_den=0.0;
int pos_neigh;

double Sz_pos, Sz_pos_neigh;
for(int pos=0;pos<Parameters_.ns;pos++){

pos_neigh=(pos+1+ns)%(Parameters_.ns);

Sz_pos = sin(YVec0[theta_to_index[pos]]);
Sz_pos_neigh = sin(YVec0[theta_to_index[pos_neigh]]);

kink_den += (0.5/Parameters_.ns)*(1.0 - (Sz_pos*Sz_pos_neigh) );

}

return kink_den;
}

double SC_SW_ENGINE_Rotor::Get_Kink_Density_1d(){

double kink_den=0.0;
int pos_neigh;

double Sz_pos, Sz_pos_neigh;
double sign_temp;
for(int pos=0;pos<ns;pos++){

pos_neigh=(pos+1+ns)%(Parameters_.ns);

Sz_pos = sin(YVec0[theta_to_index[pos]]);
Sz_pos_neigh = sin(YVec0[theta_to_index[pos_neigh]]);

sign_temp = (sign(Sz_pos)*sign(Sz_pos_neigh));
//sign_temp = -1;
kink_den += (0.5/(1.0*Parameters_.ns))*(1.0 - sign_temp);

}

return kink_den;
}

double SC_SW_ENGINE_Rotor::sign(double x){
    double sign_temp;
    if(x<0){
        sign_temp=-1.0;
    }
    else{
        sign_temp=1.0;
    }
    return sign_temp;
}

void SC_SW_ENGINE_Rotor::Evolve_classical_spins_Runge_Kutta(int ts){


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




void SC_SW_ENGINE_Rotor::Map_Variables_to_Y(Mat_1_doub & Theta_, Mat_1_doub & Momentum_,  Mat_1_doub & Y_)
{

    int Y_size;
     Y_size = 2*ns;

    Y_.resize(Y_size);
    //Convention Aux_Sx--Aux_Sy---Aux_Sz-- Red_Den_mat[][][][]

    int index=0;
    for(int i=0;i<ns;i++){
            Y_[index]=Theta_[i];
            index++;
    }
    for(int i=0;i<ns;i++){
            Y_[index]=Momentum_[i];
            index++;
    }

    assert(index==Y_size);

}



void SC_SW_ENGINE_Rotor::Map_Y_to_Variables(Mat_1_doub & Y_, Mat_1_doub & Theta_, Mat_1_doub & Momentum_)
{

    int Y_size;
    Y_size = 2*ns;

    Theta_.resize(ns);
    Momentum_.resize(ns);

    //Convention Aux_Sx--Aux_Sy---Aux_Sz-- Red_Den_mat[][][][]


    int index=0;
    for(int i=0;i<ns;i++){
            Theta_[i]=Y_[index];
            index++;
    }
    for(int i=0;i<ns;i++){
            Momentum_[i]=Y_[index];
            index++;
    }

    assert(index==Y_size);

}


void SC_SW_ENGINE_Rotor::IndexMapping_bw_Y_and_Variables()
{

    int Y_size;

    Y_size = 2*ns;

    index_to_theta.resize(Y_size);
    index_to_momentum.resize(Y_size);

    theta_to_index.resize(ns);
    momentum_to_index.resize(ns);


    //Convention Aux_Sx--Aux_Sy---Aux_Sz-- Red_Den_mat[][][][]

    int index=0;
    for(int i=0;i<Parameters_.ns;i++){
            theta_to_index[i] = index;
            index_to_theta[index] = i;
            index++;
    }
    for(int i=0;i<Parameters_.ns;i++){
            momentum_to_index[i] = index;
            index_to_momentum[index] = i;
            index++;
    }

    assert(index==Y_size);


}




double SC_SW_ENGINE_Rotor::Get_Classical_Energy(Mat_1_doub & Y_){

    double Class_Energy;
    Class_Energy=0.0;

    for(int pos=0;pos<Parameters_.ns;pos++){
            //Class_Energy += -1.0*Parameters_.hz_mag*sin(Y_[theta_to_index[pos]]);
            Class_Energy += -1.0*Gamma_[TIME_STEP_GLOBAL]*Parameters_.hx_mag*cos(Y_[theta_to_index[pos]]);


    for(int pos_neigh=0;pos_neigh<ns;pos_neigh++){
        if(abs(Parameters_.Jzz_longrange(pos,pos_neigh))>EPSILON){
        Class_Energy += 0.5*Js_[TIME_STEP_GLOBAL]*Parameters_.Jzz_longrange(pos,pos_neigh)*
                        sin(Y_[theta_to_index[pos]])*sin(Y_[theta_to_index[pos_neigh]]);
        }
    }

    }

    return Class_Energy;

}




void SC_SW_ENGINE_Rotor::Derivative(Mat_1_doub & Y_, Mat_1_doub & dYbydt){


    assert(Use_Scheduler);

    dYbydt.resize(Y_.size());
    for(int i=0;i<Y_.size();i++){
        dYbydt[i]=0.0;
    }

 //HERE

//#ifdef _OPENMP
//#pragma omp parallel for default(shared) private(sx, sy, sz, pos_neigh, Spin_no)
//#endif
    for(int pos=0;pos<Parameters_.ns;pos++){

    dYbydt[theta_to_index[pos]] += 1.0*Y_[momentum_to_index[pos]]/(Rotor_mass);

    //Noise
    dYbydt[momentum_to_index[pos]] += GaussianDistribution(GeneratorNoise_)*(StdDev/sqrt(dt_)) +
                                      -1.0*((Parameters_.DampingConst/Rotor_mass)*Y_[momentum_to_index[pos]])
                                      -1.0*(Gamma_[TIME_STEP_GLOBAL]*Parameters_.hx_mag*sin(Y_[theta_to_index[pos]]));


    for(int pos_neigh=0;pos_neigh<ns;pos_neigh++){
        if(abs(Parameters_.Jzz_longrange(pos,pos_neigh))>EPSILON){
    dYbydt[momentum_to_index[pos]] += -1.0*Js_[TIME_STEP_GLOBAL]*Parameters_.Jzz_longrange(pos,pos_neigh)*
                            sin(Y_[theta_to_index[pos_neigh]])*cos(Y_[theta_to_index[pos]]);
        }
    }


}

}



void SC_SW_ENGINE_Rotor::RungeKuttaOne(Mat_1_doub & Yn, Mat_1_doub & Ynp1){


    Ynp1.resize(Yn.size());
    Mat_1_doub K1;

    //step_no==0
    Derivative(Yn,K1);

    for(int i=0;i<Yn.size();i++){
        Ynp1[i] = Yn[i] + (dt_)*(K1[i]);

    }


}
void SC_SW_ENGINE_Rotor::RungeKuttaFour(Mat_1_doub & Yn, Mat_1_doub & Ynp1){



    Mat_1_doub Y_temp;
    Y_temp.resize(Yn.size());

    Ynp1.resize(Yn.size());
    Mat_1_doub K1, K2, K3, K4;


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


void SC_SW_ENGINE_Rotor::RungeKuttaSix(Mat_1_doub & Yn, Mat_1_doub & Ynp1){



    //Sixth Order Runge-Kutta Method is used
    //From "https://www.ams.org/journals/mcom/1968-22-102/S0025-5718-68-99876-1/S0025-5718-68-99876-1.pdf"
    //or google "An Explicit Sixth-Order Runge-Kutta Formula By H. A. Luther"



    double nu, surd;

    nu=0.9;
    surd=-1.0*sqrt(21.0);

    Mat_1_doub Y_temp;
    Y_temp.resize(Yn.size());

    Ynp1.resize(Yn.size());
    Mat_1_doub K1, K2, K3, K4, K5, K6, K7;


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




void SC_SW_ENGINE_Rotor::RungeKuttaEight(Mat_1_doub & Yn, Mat_1_doub & Ynp1){


    //Shanks Eighth order Runge Kutta

    Mat_1_doub Y_temp;
    Y_temp.resize(Yn.size());

    Ynp1.resize(Yn.size());
    Mat_1_doub K0, K1, K2, K3, K4, K5, K6, K7, K8, K9;


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



void SC_SW_ENGINE_Rotor::Read_parameters(string filename){


    string rotor_mass_, Rotor_Mass_ = "Rotor_mass = ";

    string time_max_, Time_Max_ = "time_max = ";
    string dt__, dT_ = "dt = ";
    string runge_kutta_order_, Runge_Kutta_Order_ = "Runge_Kutta_order = ";

    string Spins_R_T_Out_ = "spins_r_t_out = ";

    string no_of_processors_, No_Of_Processors_ = "no_of_threads = ";

    string Scheduler_File_ = "Scheduler_File = ";

    string use_scheduler_, Use_Scheduler_ = "Use_Scheduler = ";

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

            if ((offset = line.find(Use_Scheduler_, 0)) != string::npos) {
                use_scheduler_ = line.substr (offset + Use_Scheduler_.length());		}


            if ((offset = line.find(Time_Max_, 0)) != string::npos) {
                time_max_ = line.substr (offset + Time_Max_.length());		}


            if ((offset = line.find(Rotor_Mass_, 0)) != string::npos) {
                rotor_mass_ = line.substr (offset + Rotor_Mass_.length());		}

            if ((offset = line.find(dT_, 0)) != string::npos) {
                dt__= line.substr (offset + dT_.length());		}

            if ((offset = line.find(Runge_Kutta_Order_, 0)) != string::npos) {
                runge_kutta_order_ = line.substr (offset + Runge_Kutta_Order_.length());		}


             if ((offset = line.find(Scheduler_File_, 0)) != string::npos) {
                 Scheduler_File = line.substr (offset + Scheduler_File_.length());		}


            if ((offset = line.find(Spins_R_T_Out_, 0)) != string::npos) {
                spins_r_t_out = line.substr (offset + Spins_R_T_Out_.length());		}


            if ((offset = line.find(No_Of_Processors_, 0)) != string::npos) {
                no_of_processors_ = line.substr (offset + No_Of_Processors_.length());		}

        }
        inputfile.close();
    }
    else
    {cout<<"Unable to open input file while in the model class."<<endl;}




    Rotor_mass=atof(rotor_mass_.c_str());
    time_max=atof(time_max_.c_str());
    dt_ =atof(dt__.c_str());
    Runge_Kutta_order=atoi(runge_kutta_order_.c_str());

    no_of_processors=atoi(no_of_processors_.c_str());


    if(use_scheduler_ == "true"){
        Use_Scheduler=true;
        cout<<"Scheduler is used i.e. time dependent hamiltonian"<<endl;
    }
    else{
        Use_Scheduler=false;
    }


    cout<<"PARAMETERS::::::::"<<endl;

    cout<<"time_max = "<<time_max<<endl;
    cout<<"dt = "<<dt_<<endl;
    cout<<"Runge_Kutta_order = "<<Runge_Kutta_order<<endl;
    cout<<"Temperature = "<<Parameters_.temp<<endl;
    cout<<"TotalSites = "<<Parameters_.ns<<endl;


    string line2;
    double temp_t, temp_h, temp_J;
    if(Use_Scheduler){
    ifstream scheduler_stream(Scheduler_File.c_str());
    while(getline(scheduler_stream,line2)){
    stringstream line_ss(line2);
    line_ss>>temp_t>>temp_h>>temp_J;
    Time_bare.push_back(temp_t);
    Gamma_bare.push_back(temp_h);
    Js_bare.push_back(temp_J);
    }
    }




}







