#include <iostream>  //for cin and cout
#include <math.h>  // for pow
#include <stdlib.h>  //for div(q,n).rem(quot),abs(int n)
#include <time.h>
#include <fstream>
#include <sstream>
#include "Matrix.h"

#include "src/Model_1orbHubbard_old/Basis_1orb_SF.h"
#include "src/Model_1orbHubbard_old/Model_1orb_SF.h"
#include "src/Model_1orbHubbard_old/Spin_dynamics_VNE_1orbHubbard_engine.h"
#include "src/Model_1orbHubbard_old/Space_Time_FourierTransform.h"
#include "src/Model_3orbPnictides/Coordinates.h"
#include "src/Model_3orbPnictides/Hamiltonian.h"

#include "src/Model_3orbPnictides/MFParams.h"
#include "src/Model_3orbPnictides/Observables.h"
#include "src/Model_3orbPnictides/ParametersEngine.h"
#include "src/Model_3orbPnictides/Space_Time_FourierTransform_3orb.h"
#include "src/Model_3orbPnictides/Spin_dynamics_VNE_3orbPnictides_engine.h"

#include "src/MCMF_1orbHubbard/Coordinates_MCMF.h"
#include "src/MCMF_1orbHubbard/Hamiltonian_MCMF.h"
#include "src/MCMF_1orbHubbard/MFParams_MCMF.h"
#include "src/MCMF_1orbHubbard/Observables_MCMF.h"
#include "src/MCMF_1orbHubbard/ParametersEngine_MCMF.h"
#include "src/MCMF_1orbHubbard/Space_Time_FourierTransform_1orb_MCMF.h"
#include "src/MCMF_1orbHubbard/Spin_dynamics_VNE_1orb_engine_MCMF.h"

#include "src/MultiOrb_SpinFermion/Coordinates_MultiOrbSF.h"
#include "src/MultiOrb_SpinFermion/Hamiltonian_MultiOrbSF.h"
#include "src/MultiOrb_SpinFermion/MFParams_MultiOrbSF.h"
#include "src/MultiOrb_SpinFermion/Observables_MultiOrbSF.h"
#include "src/MultiOrb_SpinFermion/ParametersEngine_MultiOrbSF.h"
#include "src/MultiOrb_SpinFermion/Space_Time_FourierTransform_MultiOrbSF.h"
#include "src/MultiOrb_SpinFermion/Spin_dynamics_VNE_MultiOrbSF.h"


#include "src/RotorDynamics/ParametersEngine_Rotor.h"
#include "src/RotorDynamics/Spin_dynamics_Rotor.h"


int main(int argc, char** argv){


    string ex_string_original =argv[0];
    string model_=argv[1];

    string input;
    input=argv[2];

    if( (model_!="1orbHubard")
            &&
            (model_!="3orbPnictides")
            &&
            (model_!="1orb_MCMF")
             &&
             (model_!="MultiOrbSF")
            &&
            (model_!="Rotor")
            ){
        cout<<"You are using "<<argv[1]<<" , this model is not present"<<endl;
        cout<<"Please something from following :"<<endl;
        cout<<"1orbHubard, 3orbPnictides"<<endl;
        assert(model_=="1orbHubard" || model_=="3orbPnictides" || model_=="1orb_MCMF" || model_=="MultiOrbSF");
    }
    cout<<ex_string_original<<endl;

    string ex_string;

    //ex_string.substr(ex_string_original.length()-5);

    ex_string=ex_string_original.substr (ex_string_original.length() - 5);
    cout<<ex_string<<endl;

    if(ex_string == "amics"){


        bool S_kw_using_Von_Nuemann=true;


        /*
        if(S_kw_using_ED){

            BASIS_1_orb_SF Basis;
            MODEL_1_orb_SF Model;
            SC_SW_ENGINE_ED Skw_Engine;


            Skw_Engine.Read_parameters(input, Model, Basis);


            Basis.Initialize_Basis();
            Skw_Engine.Initialize_engine(Basis, Model);

            Skw_Engine.Read_equilibrium_configuration(Basis, Model);

            Skw_Engine.Start_Engine(Basis, Model);

        }
        */

        if(S_kw_using_Von_Nuemann){

            if(model_=="1orbHubardXX"){
                BASIS_1_orb_SF Basis;
                MODEL_1_orb_SF Model;
                SC_SW_ENGINE_VNE_1orbHubbard Skw_Engine;

                Skw_Engine.Read_parameters(input, Model, Basis);


                Basis.Initialize_Basis();
                Skw_Engine.Initialize_engine(Model, Basis);

                if(!Skw_Engine.RESTART){
                    Skw_Engine.Read_equilibrium_configuration(Model, Basis);
                }

                Skw_Engine.Start_Engine(Model, Basis);
            }


            if(model_=="3orbPnictides"){

                Parameters Parameters_;
                Parameters_.Initialize(input);

                Coordinates Coordinates_(Parameters_.lx, Parameters_.ly);

                mt19937_64 Generator_(Parameters_.RandomSeed);
                MFParams MFParams_(Parameters_,Coordinates_,Generator_);

                Hamiltonian Hamiltonian_(Parameters_,Coordinates_,MFParams_);
                Observables Observables_(Parameters_,Coordinates_,MFParams_,Hamiltonian_);


                SC_SW_ENGINE_VNE_3orbPnictides Skw_Engine_(Parameters_,Coordinates_,MFParams_,Hamiltonian_,Observables_);
                Skw_Engine_.Read_parameters(input);
                Skw_Engine_.Initialize_engine();

                if(!Skw_Engine_.RESTART){
                    Skw_Engine_.Read_equilibrium_configuration();
                }

                Skw_Engine_.Start_Engine();

            }

            if(model_=="1orb_MCMF"){

                Parameters_MCMF Parameters_;
                Parameters_.Initialize(input);

                Coordinates_MCMF Coordinates_(Parameters_.lx, Parameters_.ly);

                mt19937_64 Generator_(Parameters_.RandomSeed);
                MFParams_MCMF MFParams_(Parameters_,Coordinates_,Generator_);

                Hamiltonian_MCMF Hamiltonian_(Parameters_,Coordinates_,MFParams_);
                Observables_MCMF Observables_(Parameters_,Coordinates_,MFParams_,Hamiltonian_);


                SC_SW_ENGINE_VNE_1orb_MCMF Skw_Engine_(Parameters_,Coordinates_,MFParams_,Hamiltonian_,Observables_);
                Skw_Engine_.Read_parameters(input);
                Skw_Engine_.Initialize_engine();
                Skw_Engine_.IndexMapping_bw_Y_and_Variables();

                if(!Skw_Engine_.RESTART){
                    Skw_Engine_.Read_equilibrium_configuration();
                }

                Skw_Engine_.Start_Engine();

            }

            if(model_=="MultiOrbSF"){

                Parameters_MultiOrbSF Parameters_;
                Parameters_.Initialize(input);


                double Avg_KinkDen_type1, Avg_KinkDen_type2;
                Avg_KinkDen_type1=0;Avg_KinkDen_type2=0;
                for(int noiseseed_no=0;noiseseed_no<Parameters_.RandomNoiseSeed_array.size();noiseseed_no++){
                 Parameters_.RandomNoiseSeed = Parameters_.RandomNoiseSeed_array[noiseseed_no];

                 cout<<"================================================="<<endl;
                 cout<<"FOR NOISE SEED = "<<Parameters_.RandomNoiseSeed<<endl;
                 cout<<"================================================="<<endl;



                Coordinates_MultiOrbSF Coordinates_(Parameters_.lx, Parameters_.ly, Parameters_.n_orbs);
                Coordinates_MultiOrbSF Coordinatestemp_(Parameters_.lx, Parameters_.ly, Parameters_.n_orbs);


                mt19937_64 Generator_(Parameters_.RandomSeed);
                mt19937_64 Generator2_(Parameters_.RandomDisorderSeed); //for random disorder

                MFParams_MultiOrbSF MFParams_(Parameters_,Coordinates_,Generator_, Generator2_);

                Hamiltonian_MultiOrbSF Hamiltonian_(Parameters_,Coordinates_, Coordinatestemp_, MFParams_);
                Observables_MultiOrbSF Observables_(Parameters_,Coordinates_,MFParams_,Hamiltonian_);


                SC_SW_ENGINE_VNE_MultiOrbSF Skw_Engine_(Parameters_,Coordinates_,MFParams_,Hamiltonian_,Observables_);
                Skw_Engine_.Read_parameters(input);
                Skw_Engine_.Initialize_engine();
                Skw_Engine_.IndexMapping_bw_Y_and_Variables();

                if(!Skw_Engine_.RESTART){
                    if(Skw_Engine_.conf_initialize=="Read"){
                    Skw_Engine_.Read_equilibrium_configuration();
                    }
                    else if(Skw_Engine_.conf_initialize=="Random"){
                    Skw_Engine_.Set_Initial_configuration();
                    }
                    else{
                     cout<< "ERROR : conf_initialize should be Read, Random"<<endl;
                     assert(false);
                    }
                }

                Skw_Engine_.Start_Engine();

            Avg_KinkDen_type1 += (1.0/(Parameters_.RandomNoiseSeed_array.size()))*Skw_Engine_.Rotor_KinkDen_type1;
            Avg_KinkDen_type2 += (1.0/(Parameters_.RandomNoiseSeed_array.size()))*Skw_Engine_.Rotor_KinkDen_type2;
            }//Random Noise Seeds

            cout<<"Avg_KinkDen_type1 : "<<Avg_KinkDen_type1<<endl;
            cout<<"Avg_KinkDen_type2 : "<<Avg_KinkDen_type2<<endl;
            }


            if(model_=="Rotor"){


                bool Print_Cij_for_every_noise=true;
                bool Print_Cij_Noise_Avgd=true;

                Parameters_Rotor Parameters_;
                Parameters_.Initialize(input);

                double Avg_KinkDen_type1, Avg_KinkDen_type2;
                Avg_KinkDen_type1=0;Avg_KinkDen_type2=0;


                Matrix<double> C_ij_type1_Noise_Avgd, C_ij_type2_Noise_Avgd;
                C_ij_type1_Noise_Avgd.resize(Parameters_.ns,Parameters_.ns);
                C_ij_type2_Noise_Avgd.resize(Parameters_.ns,Parameters_.ns);

                double q2_type1, q2_type2;
                double q4_type1, q4_type2;
                int ns = Parameters_.ns;
                for(int noiseseed_no=0;noiseseed_no<Parameters_.RandomNoiseSeed_array.size();noiseseed_no++){
                 Parameters_.RandomNoiseSeed = Parameters_.RandomNoiseSeed_array[noiseseed_no];

                 cout<<"================================================="<<endl;
                 cout<<"FOR NOISE SEED = "<<Parameters_.RandomNoiseSeed<<endl;
                 cout<<"================================================="<<endl;


                SC_SW_ENGINE_Rotor Skw_Engine_(Parameters_);

                Skw_Engine_.Read_parameters(input);
                Skw_Engine_.Initialize_engine();
                Skw_Engine_.IndexMapping_bw_Y_and_Variables();
                Skw_Engine_.Start_Engine();

            Avg_KinkDen_type1 += (1.0/(Parameters_.RandomNoiseSeed_array.size()))*Skw_Engine_.Rotor_KinkDen_type1;
            Avg_KinkDen_type2 += (1.0/(Parameters_.RandomNoiseSeed_array.size()))*Skw_Engine_.Rotor_KinkDen_type2;

            for(int pos_i=0;pos_i<Parameters_.ns;pos_i++){
                for(int pos_j=0;pos_j<Parameters_.ns;pos_j++){
               C_ij_type1_Noise_Avgd(pos_i, pos_j) += (1.0/(Parameters_.RandomNoiseSeed_array.size()))*
                                                        Skw_Engine_.C_ij_type1[pos_i][pos_j];
               C_ij_type2_Noise_Avgd(pos_i, pos_j) += (1.0/(Parameters_.RandomNoiseSeed_array.size()))*
                                                        Skw_Engine_.C_ij_type2[pos_i][pos_j];
                }
            }


            if(Print_Cij_for_every_noise){
            string Cij_1_file_String= "Cij_type_1_NoiseSeed_"+to_string(Parameters_.RandomNoiseSeed)+".txt";
            ofstream Cij_1_file(Cij_1_file_String.c_str());

            string Cij_2_file_String= "Cij_type_2_NoiseSeed_"+to_string(Parameters_.RandomNoiseSeed)+".txt";
            ofstream Cij_2_file(Cij_2_file_String.c_str());

            for(int pos_i=0;pos_i<Parameters_.ns;pos_i++){
                for(int pos_j=0;pos_j<Parameters_.ns;pos_j++){
              Cij_1_file<<pos_i<<"  "<<pos_j<<"  "<<  Skw_Engine_.C_ij_type1[pos_i][pos_j]<<endl;
              Cij_2_file<<pos_i<<"  "<<pos_j<<"  "<<  Skw_Engine_.C_ij_type2[pos_i][pos_j]<<endl;
              }
            }
            }



                }//Random Noise Seeds

                if(Print_Cij_Noise_Avgd){
                string Cij_1_avg_file_String= "Cij_type_1_NoiseAvgd.txt";
                ofstream Cij_1_avg_file(Cij_1_avg_file_String.c_str());

                string Cij_2_avg_file_String= "Cij_type_2_NoiseAvgd.txt";
                ofstream Cij_2_avg_file(Cij_2_avg_file_String.c_str());

                for(int pos_i=0;pos_i<Parameters_.ns;pos_i++){
                    for(int pos_j=0;pos_j<Parameters_.ns;pos_j++){
                  Cij_1_avg_file<<pos_i<<"  "<<pos_j<<"  "<<  C_ij_type1_Noise_Avgd(pos_i, pos_j)<<endl;
                  Cij_2_avg_file<<pos_i<<"  "<<pos_j<<"  "<<  C_ij_type2_Noise_Avgd(pos_i, pos_j)<<endl;
                  }
                }
                }

            //q2 calculation using noise Avgd Cij
            q2_type1=0.0;
            q2_type2=0.0;
            for(int pos_i=0;pos_i<Parameters_.ns;pos_i++){
                for(int pos_j=pos_i+1;pos_j<Parameters_.ns;pos_j++){
                q2_type1 += (2.0/(ns*(ns-1)))*C_ij_type1_Noise_Avgd(pos_i, pos_j)*
                            C_ij_type1_Noise_Avgd(pos_i, pos_j);
                q2_type2 += (2.0/(ns*(ns-1)))*C_ij_type2_Noise_Avgd(pos_i, pos_j)*
                            C_ij_type2_Noise_Avgd(pos_i, pos_j);
            }}

            //cout<<"Avg_KinkDen_type1 : "<<Avg_KinkDen_type1<<endl;
            //cout<<"Avg_KinkDen_type2 : "<<Avg_KinkDen_type2<<endl;

            cout<<"q2_type1 : "<<q2_type1<<endl;
            cout<<"q2_type2 : "<<q2_type2<<endl;
            cout<<"q4_type1 : "<<q2_type1*q2_type1<<endl;
            cout<<"q4_type2 : "<<q2_type2*q2_type2<<endl;
            }


        }

    }


    if(ex_string == "I_Skw"){

        cout<<"S(k,w) for Tight-binding model[Non-interacting fermions] is calculated"<<endl;


        if(model_=="1orbHubard"){
            string NISkw_filename_full = "Skw_NI_full.txt";
            string NISkw_filename_specific_kpath = "Skw_NI_specific_kpath.txt";
            Get_NI_Skw(NISkw_filename_full , NISkw_filename_specific_kpath, 2, 1.0);
        }

        if(model_=="3orbPnictides"){


            string NISkw_filename_full = "Skw_NI_full.txt";
            Parameters Parameters_;
            Parameters_.Initialize(input);

            Coordinates Coordinates_(Parameters_.lx, Parameters_.ly);

            mt19937_64 Generator_(Parameters_.RandomSeed);
            MFParams MFParams_(Parameters_,Coordinates_,Generator_);

            Hamiltonian Hamiltonian_(Parameters_,Coordinates_,MFParams_);
            Observables Observables_(Parameters_,Coordinates_,MFParams_,Hamiltonian_);

            Hamiltonian_.InteractionsCreate();
            Hamiltonian_.Check_Hermiticity();
            //char flag='V';
            Hamiltonian_.Diagonalize('V');
            double mu_ = Hamiltonian_.chemicalpotential(0.5,Parameters_.Fill);
            cout<<"mu = "<<mu_<<endl;
            cout<<Parameters_.Fill<<endl;
            //Observables_.Calculate_Nw();
            //Observables_.Calculate_Akw();
            Observables_.Calculate_Skw(mu_);

        }



    }


    if(ex_string == "urier"){

        if(model_=="1orbHubard"){

            cout<<"Space time Fourier Transform is performed on time-displaced correlation function, Ensemble Averaging is done."<<endl;
            cout<<"---------------[ <m_{i,t}m_{j,t}> - <m_{i,t}><m_{j,t}> ] is used, '< >' means Ensemble Average.--------------"<<endl;



            BASIS_1_orb_SF Basis;

            ST_Fourier SpaceTime_Fourier;

            string No_of_inputs= argv[3];
            stringstream No_of_inputs_ss(No_of_inputs, stringstream::in);
            No_of_inputs_ss>>SpaceTime_Fourier.No_Of_Inputs;

            SpaceTime_Fourier.conf_inputs.resize(SpaceTime_Fourier.No_Of_Inputs);

            for(int i=0;i<SpaceTime_Fourier.No_Of_Inputs;i++){
                SpaceTime_Fourier.conf_inputs[i]=argv[i+4];
            }


            SpaceTime_Fourier.Read_parameters(input, Basis);
            Basis.Initialize_Basis();


            SpaceTime_Fourier.Initialize_engine(Basis);
            SpaceTime_Fourier.Perform_Smarter_Averaging_on_one_point(Basis);


            for(int i=0;i<SpaceTime_Fourier.No_Of_Inputs;i++){

                string str_int;

                stringstream ss;
                ss << i;
                ss >> str_int;

                string Skw_conf = "Skw_conf_" + str_int + ".txt";
                cout<<Skw_conf<<endl;
                SpaceTime_Fourier.Calculate_Skw_from_Srt_file(Basis, SpaceTime_Fourier.conf_inputs[i],Skw_conf);

            }
            SpaceTime_Fourier.Calculate_Skw_from_Srt_file(Basis, "Temp_file_Srt.txt", "Skw_OnAveragedConfs.txt");


        }

        if(model_=="3orbPnictides"){

            Parameters Parameters_;
            Parameters_.Initialize(input);

            Coordinates Coordinates_(Parameters_.lx, Parameters_.ly);

            mt19937_64 Generator_(Parameters_.RandomSeed);
            MFParams MFParams_(Parameters_,Coordinates_,Generator_);

            Hamiltonian Hamiltonian_(Parameters_,Coordinates_,MFParams_);
            Observables Observables_(Parameters_,Coordinates_,MFParams_,Hamiltonian_);


            SC_SW_ENGINE_VNE_3orbPnictides Skw_Engine_(Parameters_,Coordinates_,MFParams_,Hamiltonian_,Observables_);
            Skw_Engine_.Read_parameters(input);
            Skw_Engine_.Initialize_engine();



            ST_Fourier_3orb SpaceTime_Fourier(Parameters_,Coordinates_,MFParams_,Hamiltonian_,Observables_, Skw_Engine_);

            string No_of_inputs= argv[3];
            stringstream No_of_inputs_ss(No_of_inputs, stringstream::in);
            No_of_inputs_ss>>SpaceTime_Fourier.No_Of_Inputs;

            SpaceTime_Fourier.conf_inputs.resize(SpaceTime_Fourier.No_Of_Inputs);

            for(int i=0;i<SpaceTime_Fourier.No_Of_Inputs;i++){
                SpaceTime_Fourier.conf_inputs[i]=argv[i+4];
            }


            SpaceTime_Fourier.Read_parameters();
            SpaceTime_Fourier.Initialize_engine();
            SpaceTime_Fourier.Perform_Smarter_Averaging_on_one_point();

            ostringstream ostr_w_conv;
            ostr_w_conv << SpaceTime_Fourier.w_conv;
            string string_w_conv = ostr_w_conv.str();


            string SpaceTimeDisplaced_Crt_file = "SpaceTimeDisplaced_Crt_w_conv" + string_w_conv + ".txt";
            SpaceTime_Fourier.Calculate_SpaceTimeDisplacedCorrelations(SpaceTimeDisplaced_Crt_file);

            string Skw_using_Crt_file = "Skw_using_Crt_w_conv" + string_w_conv + ".txt";
            SpaceTime_Fourier.Calculate_Skw_from_Crt(Skw_using_Crt_file);



            /*
            for(int i=0;i<SpaceTime_Fourier.No_Of_Inputs;i++){

                string str_int;

                stringstream ss;
                ss << i;
                ss >> str_int;

                string Skw_conf = "Skw_conf_" + str_int + ".txt";
                cout<<Skw_conf<<endl;
                SpaceTime_Fourier.Calculate_Skw_from_Srt_file(SpaceTime_Fourier.conf_inputs[i],Skw_conf);

            }


            SpaceTime_Fourier.Calculate_Skw_from_Srt_file( "Average_Srt.txt", "Skw_OnAveragedConfs.txt");
            */


        }


        if(model_=="1orb_MCMF"){

            Parameters_MCMF Parameters_;
            Parameters_.Initialize(input);

            Coordinates_MCMF Coordinates_(Parameters_.lx, Parameters_.ly);

            mt19937_64 Generator_(Parameters_.RandomSeed);
            MFParams_MCMF MFParams_(Parameters_,Coordinates_,Generator_);

            Hamiltonian_MCMF Hamiltonian_(Parameters_,Coordinates_,MFParams_);
            Observables_MCMF Observables_(Parameters_,Coordinates_,MFParams_,Hamiltonian_);


            SC_SW_ENGINE_VNE_1orb_MCMF Skw_Engine_(Parameters_,Coordinates_,MFParams_,Hamiltonian_,Observables_);
            Skw_Engine_.Read_parameters(input);
            Skw_Engine_.Initialize_engine();



            ST_Fourier_1orb_MCMF SpaceTime_Fourier(Parameters_,Coordinates_,MFParams_,Hamiltonian_,Observables_, Skw_Engine_);

            string No_of_inputs= argv[3];
            stringstream No_of_inputs_ss(No_of_inputs, stringstream::in);
            No_of_inputs_ss>>SpaceTime_Fourier.No_Of_Inputs;

            SpaceTime_Fourier.conf_inputs.resize(SpaceTime_Fourier.No_Of_Inputs);

            for(int i=0;i<SpaceTime_Fourier.No_Of_Inputs;i++){
                SpaceTime_Fourier.conf_inputs[i]=argv[i+4];
            }


#ifdef _OPENMP
            cout<<"Parallel threads are used"<<endl;
#endif
#ifndef _OPENMP
            cout<<"single thread is used"<<endl;
#endif

            SpaceTime_Fourier.Read_parameters();
            SpaceTime_Fourier.Initialize_engine();
            SpaceTime_Fourier.Space_Fourier_using_single_S=false;


            SpaceTime_Fourier.Perform_Smarter_Averaging_on_one_point();

            ostringstream ostr_w_conv;
            ostr_w_conv << SpaceTime_Fourier.w_conv;
            string string_w_conv = ostr_w_conv.str();


            string SpaceTimeDisplaced_Crt_file = "SpaceTimeDisplaced_Crt_w_conv" + string_w_conv + ".txt";
            //SpaceTime_Fourier.Calculate_SpaceTimeDisplacedCorrelations(SpaceTimeDisplaced_Crt_file);
            SpaceTime_Fourier.Calculate_SpaceTimeDisplacedCorrelations_Smarter(SpaceTimeDisplaced_Crt_file);



            string Skw_using_Crt_file = "Skw_using_Crt_w_conv" + string_w_conv + ".txt";
            //SpaceTime_Fourier.Calculate_Skw_from_Crt(Skw_using_Crt_file);
            SpaceTime_Fourier.Calculate_Skw_from_Crt_(Skw_using_Crt_file);


        }



        if(model_=="MultiOrbSF"){

            Parameters_MultiOrbSF Parameters_;
            Parameters_.Initialize(input);

            Coordinates_MultiOrbSF Coordinates_(Parameters_.lx, Parameters_.ly, Parameters_.n_orbs);
            Coordinates_MultiOrbSF Coordinatestemp_(Parameters_.lx, Parameters_.ly, Parameters_.n_orbs);


            mt19937_64 Generator_(Parameters_.RandomSeed);
            mt19937_64 Generator2_(Parameters_.RandomDisorderSeed);

            MFParams_MultiOrbSF MFParams_(Parameters_,Coordinates_,Generator_, Generator2_);

            Hamiltonian_MultiOrbSF Hamiltonian_(Parameters_,Coordinates_, Coordinatestemp_, MFParams_);
            Observables_MultiOrbSF Observables_(Parameters_,Coordinates_,MFParams_,Hamiltonian_);


            SC_SW_ENGINE_VNE_MultiOrbSF Skw_Engine_(Parameters_,Coordinates_,MFParams_,Hamiltonian_,Observables_);
            Skw_Engine_.Read_parameters(input);
            Skw_Engine_.Initialize_engine();



            ST_Fourier_MultiOrbSF SpaceTime_Fourier(Parameters_,Coordinates_,MFParams_,Hamiltonian_,Observables_, Skw_Engine_);

            string No_of_inputs= argv[3];
            stringstream No_of_inputs_ss(No_of_inputs, stringstream::in);
            No_of_inputs_ss>>SpaceTime_Fourier.No_Of_Inputs;

            SpaceTime_Fourier.conf_inputs.resize(SpaceTime_Fourier.No_Of_Inputs);

            for(int i=0;i<SpaceTime_Fourier.No_Of_Inputs;i++){
                SpaceTime_Fourier.conf_inputs[i]=argv[i+4];
            }


#ifdef _OPENMP
            cout<<"Parallel threads are used"<<endl;
#endif
#ifndef _OPENMP
            cout<<"single thread is used"<<endl;
#endif

            SpaceTime_Fourier.Read_parameters();
            SpaceTime_Fourier.Initialize_engine();
            SpaceTime_Fourier.Space_Fourier_using_single_S=false;


            SpaceTime_Fourier.Perform_Smarter_Averaging_on_one_point();

            ostringstream ostr_w_conv;
            ostr_w_conv << SpaceTime_Fourier.w_conv;
            string string_w_conv = ostr_w_conv.str();


            string SpaceTimeDisplaced_Crt_file = "SpaceTimeDisplaced_Crt_w_conv" + string_w_conv + ".txt";
            //SpaceTime_Fourier.Calculate_SpaceTimeDisplacedCorrelations(SpaceTimeDisplaced_Crt_file);
            SpaceTime_Fourier.Calculate_SpaceTimeDisplacedCorrelations_Smarter(SpaceTimeDisplaced_Crt_file);



            string Skw_using_Crt_file = "Skw_using_Crt_w_conv" + string_w_conv + ".txt";
            //SpaceTime_Fourier.Calculate_Skw_from_Crt(Skw_using_Crt_file);
            SpaceTime_Fourier.Calculate_Skw_from_Crt_(Skw_using_Crt_file);


        }

    }





    if(ex_string == "e_Fqw"){

        if(model_=="1orb_MCMF"){

            Parameters_MCMF Parameters_;
            Parameters_.Initialize(input);

            Coordinates_MCMF Coordinates_(Parameters_.lx, Parameters_.ly);

            mt19937_64 Generator_(Parameters_.RandomSeed);
            MFParams_MCMF MFParams_(Parameters_,Coordinates_,Generator_);

            Hamiltonian_MCMF Hamiltonian_(Parameters_,Coordinates_,MFParams_);
            Observables_MCMF Observables_(Parameters_,Coordinates_,MFParams_,Hamiltonian_);


            SC_SW_ENGINE_VNE_1orb_MCMF Skw_Engine_(Parameters_,Coordinates_,MFParams_,Hamiltonian_,Observables_);
            Skw_Engine_.Read_parameters(input);
            Skw_Engine_.Initialize_engine();



            ST_Fourier_1orb_MCMF SpaceTime_Fourier(Parameters_,Coordinates_,MFParams_,Hamiltonian_,Observables_, Skw_Engine_);

            string No_of_inputs= argv[3];
            string MS_tag= argv[4];

            stringstream No_of_inputs_ss(No_of_inputs, stringstream::in);
            No_of_inputs_ss>>SpaceTime_Fourier.No_Of_Inputs;


            if(SpaceTime_Fourier.No_Of_Inputs!=1){
                cout << "Only 1 Microstate allowed for this run"<<endl;
            }
            assert(SpaceTime_Fourier.No_Of_Inputs==1);

            SpaceTime_Fourier.conf_inputs.resize(SpaceTime_Fourier.No_Of_Inputs);

            for(int i=0;i<SpaceTime_Fourier.No_Of_Inputs;i++){
                SpaceTime_Fourier.conf_inputs[i]=argv[i+5];
            }


#ifdef _OPENMP
            cout<<"Parallel threads are used"<<endl;
#endif
#ifndef _OPENMP
            cout<<"single thread is used"<<endl;
#endif

            SpaceTime_Fourier.Read_parameters();
            SpaceTime_Fourier.time_max=Skw_Engine_.time_max;
            SpaceTime_Fourier.dt_ = Skw_Engine_.dt_;
            SpaceTime_Fourier.Initialize_engine();

            //SpaceTime_Fourier.Perform_Smarter_Averaging_on_one_point();


            ostringstream ostr_w_conv;
            ostr_w_conv << SpaceTime_Fourier.w_conv;
            string string_w_conv = ostr_w_conv.str();

            string Fw_file = "Fw_w_conv" + string_w_conv + "_state_" + MS_tag +".txt";
            string Aq_file = "Aq_state_" + MS_tag +".txt" ;
            SpaceTime_Fourier.Calculate_Fw_and_Aq(Fw_file, Aq_file);


            cout<<"JOB DONE"<<endl;

        }

        if(model_=="MultiOrbSF"){


            Parameters_MultiOrbSF Parameters_;
            Parameters_.Initialize(input);

            Coordinates_MultiOrbSF Coordinates_(Parameters_.lx, Parameters_.ly, Parameters_.n_orbs);
            Coordinates_MultiOrbSF Coordinatestemp_(Parameters_.lx, Parameters_.ly, Parameters_.n_orbs);



            mt19937_64 Generator_(Parameters_.RandomSeed);
            mt19937_64 Generator2_(Parameters_.RandomDisorderSeed);

            MFParams_MultiOrbSF MFParams_(Parameters_,Coordinates_,Generator_, Generator2_);

cout<<"here 0"<<endl;
            Hamiltonian_MultiOrbSF Hamiltonian_(Parameters_,Coordinates_, Coordinatestemp_, MFParams_);
cout<<"here 2"<<endl;
            Observables_MultiOrbSF Observables_(Parameters_,Coordinates_,MFParams_,Hamiltonian_);



            SC_SW_ENGINE_VNE_MultiOrbSF Skw_Engine_(Parameters_,Coordinates_,MFParams_,Hamiltonian_,Observables_);
            Skw_Engine_.Read_parameters(input);
            Skw_Engine_.Initialize_engine();



            ST_Fourier_MultiOrbSF SpaceTime_Fourier(Parameters_,Coordinates_,MFParams_,Hamiltonian_,Observables_, Skw_Engine_);

            string No_of_inputs= argv[3];
            string MS_tag= argv[4];

            stringstream No_of_inputs_ss(No_of_inputs, stringstream::in);
            No_of_inputs_ss>>SpaceTime_Fourier.No_Of_Inputs;


            if(SpaceTime_Fourier.No_Of_Inputs!=1){
                cout << "Only 1 Microstate allowed for this run"<<endl;
            }
            assert(SpaceTime_Fourier.No_Of_Inputs==1);

            SpaceTime_Fourier.conf_inputs.resize(SpaceTime_Fourier.No_Of_Inputs);

            for(int i=0;i<SpaceTime_Fourier.No_Of_Inputs;i++){
                SpaceTime_Fourier.conf_inputs[i]=argv[i+5];
            }


#ifdef _OPENMP
            cout<<"Parallel threads are used"<<endl;
#endif
#ifndef _OPENMP
            cout<<"single thread is used"<<endl;
#endif

            SpaceTime_Fourier.Read_parameters();
            SpaceTime_Fourier.time_max=Skw_Engine_.time_max;
            SpaceTime_Fourier.dt_ = Skw_Engine_.dt_;
            SpaceTime_Fourier.Initialize_engine();

            //SpaceTime_Fourier.Perform_Smarter_Averaging_on_one_point();


            ostringstream ostr_w_conv;
            ostr_w_conv << SpaceTime_Fourier.w_conv;
            string string_w_conv = ostr_w_conv.str();

            string Fw_file = "Fw_w_conv" + string_w_conv + "_state_" + MS_tag +".txt";
            string Aq_file = "Aq_state_" + MS_tag +".txt" ;
            SpaceTime_Fourier.Calculate_Fw_and_Aq(Fw_file, Aq_file);


            cout<<"JOB DONE"<<endl;

        }


    }





    if(ex_string == "2pFqw"){


        if(model_=="MultiOrbSF"){


            Parameters_MultiOrbSF Parameters_;
            Parameters_.Initialize(input);

            Coordinates_MultiOrbSF Coordinates_(Parameters_.lx, Parameters_.ly, Parameters_.n_orbs);
            Coordinates_MultiOrbSF Coordinatestemp_(Parameters_.lx, Parameters_.ly, Parameters_.n_orbs);



            mt19937_64 Generator_(Parameters_.RandomSeed);
            mt19937_64 Generator2_(Parameters_.RandomDisorderSeed);

            MFParams_MultiOrbSF MFParams_(Parameters_,Coordinates_,Generator_, Generator2_);

cout<<"here 0"<<endl;
            Hamiltonian_MultiOrbSF Hamiltonian_(Parameters_,Coordinates_, Coordinatestemp_, MFParams_);
cout<<"here 2"<<endl;
            Observables_MultiOrbSF Observables_(Parameters_,Coordinates_,MFParams_,Hamiltonian_);



            SC_SW_ENGINE_VNE_MultiOrbSF Skw_Engine_(Parameters_,Coordinates_,MFParams_,Hamiltonian_,Observables_);
            Skw_Engine_.Read_parameters(input);
            Skw_Engine_.Initialize_engine();



            ST_Fourier_MultiOrbSF SpaceTime_Fourier(Parameters_,Coordinates_,MFParams_,Hamiltonian_,Observables_, Skw_Engine_);

            string No_of_inputs= argv[3];
            string MS_tag= argv[4];

            stringstream No_of_inputs_ss(No_of_inputs, stringstream::in);
            No_of_inputs_ss>>SpaceTime_Fourier.No_Of_Inputs;


            if(SpaceTime_Fourier.No_Of_Inputs!=1){
                cout << "Only 1 Microstate allowed for this run"<<endl;
            }
            assert(SpaceTime_Fourier.No_Of_Inputs==1);

            SpaceTime_Fourier.conf_inputs.resize(SpaceTime_Fourier.No_Of_Inputs);

            for(int i=0;i<SpaceTime_Fourier.No_Of_Inputs;i++){
                SpaceTime_Fourier.conf_inputs[i]=argv[i+5];
            }


#ifdef _OPENMP
            cout<<"Parallel threads are used"<<endl;
#endif
#ifndef _OPENMP
            cout<<"single thread is used"<<endl;
#endif

            SpaceTime_Fourier.Read_parameters();
            SpaceTime_Fourier.time_max=Skw_Engine_.time_max;
            SpaceTime_Fourier.dt_ = Skw_Engine_.dt_;
            SpaceTime_Fourier.Initialize_engine();

            //SpaceTime_Fourier.Perform_Smarter_Averaging_on_one_point();


            ostringstream ostr_w_conv;
            ostr_w_conv << SpaceTime_Fourier.w_conv;
            string string_w_conv = ostr_w_conv.str();

            string Fw_file = "2point_Fw_w_conv" + string_w_conv + "_state_" + MS_tag +".txt";
            SpaceTime_Fourier.Calculate_2point_Fw(Fw_file);

            cout<<"JOB DONE"<<endl;

        }


    }







    if(ex_string == "e_Sqw"){

        if(model_=="1orb_MCMF"){


            Parameters_MCMF Parameters_;
            Parameters_.Initialize(input);

            Coordinates_MCMF Coordinates_(Parameters_.lx, Parameters_.ly);

            mt19937_64 Generator_(Parameters_.RandomSeed);
            MFParams_MCMF MFParams_(Parameters_,Coordinates_,Generator_);

            Hamiltonian_MCMF Hamiltonian_(Parameters_,Coordinates_,MFParams_);
            Observables_MCMF Observables_(Parameters_,Coordinates_,MFParams_,Hamiltonian_);


            SC_SW_ENGINE_VNE_1orb_MCMF Skw_Engine_(Parameters_,Coordinates_,MFParams_,Hamiltonian_,Observables_);
            Skw_Engine_.Read_parameters(input);
            Skw_Engine_.Initialize_engine();



            ST_Fourier_1orb_MCMF SpaceTime_Fourier(Parameters_,Coordinates_,MFParams_,Hamiltonian_,Observables_, Skw_Engine_);
            SpaceTime_Fourier.Space_Fourier_using_single_S=false;

            string No_of_inputs= argv[3];

            stringstream No_of_inputs_ss(No_of_inputs, stringstream::in);
            No_of_inputs_ss>>SpaceTime_Fourier.No_Of_Inputs;


            if( (SpaceTime_Fourier.No_Of_Inputs==1) && (!SpaceTime_Fourier.Space_Fourier_using_single_S) ){
                cout << "Use SpaceTime_Fourier.Space_Fourier_using_single_S=true"<<endl;
                assert(SpaceTime_Fourier.Space_Fourier_using_single_S);
            }


            SpaceTime_Fourier.Aq_inputs.resize(SpaceTime_Fourier.No_Of_Inputs);
            SpaceTime_Fourier.F_wq_inputs.resize(SpaceTime_Fourier.No_Of_Inputs);

            for(int i=0;i<SpaceTime_Fourier.No_Of_Inputs;i++){
                SpaceTime_Fourier.Aq_inputs[i]=argv[i+4];
            }
            for(int i=0;i<SpaceTime_Fourier.No_Of_Inputs;i++){
                SpaceTime_Fourier.F_wq_inputs[i]=argv[i+4+SpaceTime_Fourier.No_Of_Inputs];
            }





#ifdef _OPENMP
            cout<<"Parallel threads are used"<<endl;
#endif
#ifndef _OPENMP
            cout<<"single thread is used"<<endl;
#endif

            SpaceTime_Fourier.Read_parameters();
            //  SpaceTime_Fourier.Initialize_engine();

            string Sqw_file = argv[4+(2*SpaceTime_Fourier.No_Of_Inputs)] ;
            string Dqw_file = argv[5+(2*SpaceTime_Fourier.No_Of_Inputs)] ;
            //SpaceTime_Fourier.Calculate_Sqw_using_Fwq_with_negativeandpositive_Time(Sqw_file, Dqw_file);
            SpaceTime_Fourier.Calculate_Sqw_using_Fwq(Sqw_file, Dqw_file);


            //SpaceTime_Fourier.Calculate_Sqw_using_Fwq(Sqw_file, Dqw_file);
            //SpaceTime_Fourier.Calculate_Sqw_using_Aq_Fwq(Sqw_file, Dqw_file);



        }

        if(model_=="MultiOrbSF"){


            Parameters_MultiOrbSF Parameters_;
            Parameters_.Initialize(input);

            Coordinates_MultiOrbSF Coordinates_(Parameters_.lx, Parameters_.ly, Parameters_.n_orbs);
            Coordinates_MultiOrbSF Coordinatestemp_(Parameters_.lx, Parameters_.ly, Parameters_.n_orbs);

            mt19937_64 Generator_(Parameters_.RandomSeed);
            mt19937_64 Generator2_(Parameters_.RandomDisorderSeed);

            MFParams_MultiOrbSF MFParams_(Parameters_,Coordinates_,Generator_, Generator2_ );

            Hamiltonian_MultiOrbSF Hamiltonian_(Parameters_,Coordinates_, Coordinatestemp_, MFParams_);
            Observables_MultiOrbSF Observables_(Parameters_,Coordinates_,MFParams_,Hamiltonian_);


            SC_SW_ENGINE_VNE_MultiOrbSF Skw_Engine_(Parameters_,Coordinates_,MFParams_,Hamiltonian_,Observables_);
            Skw_Engine_.Read_parameters(input);
            Skw_Engine_.Initialize_engine();


            ST_Fourier_MultiOrbSF SpaceTime_Fourier(Parameters_,Coordinates_,MFParams_,Hamiltonian_,Observables_, Skw_Engine_);
            SpaceTime_Fourier.Space_Fourier_using_single_S=false;

            string No_of_inputs= argv[3];

            stringstream No_of_inputs_ss(No_of_inputs, stringstream::in);
            No_of_inputs_ss>>SpaceTime_Fourier.No_Of_Inputs;


            if( (SpaceTime_Fourier.No_Of_Inputs==1) && (!SpaceTime_Fourier.Space_Fourier_using_single_S) ){
                cout << "Use SpaceTime_Fourier.Space_Fourier_using_single_S=true"<<endl;
                assert(SpaceTime_Fourier.Space_Fourier_using_single_S);
            }


            SpaceTime_Fourier.Aq_inputs.resize(SpaceTime_Fourier.No_Of_Inputs);
            SpaceTime_Fourier.F_wq_inputs.resize(SpaceTime_Fourier.No_Of_Inputs);

            for(int i=0;i<SpaceTime_Fourier.No_Of_Inputs;i++){
                SpaceTime_Fourier.Aq_inputs[i]=argv[i+4];
            }
            for(int i=0;i<SpaceTime_Fourier.No_Of_Inputs;i++){
                SpaceTime_Fourier.F_wq_inputs[i]=argv[i+4+SpaceTime_Fourier.No_Of_Inputs];
            }





#ifdef _OPENMP
            cout<<"Parallel threads are used"<<endl;
#endif
#ifndef _OPENMP
            cout<<"single thread is used"<<endl;
#endif

            SpaceTime_Fourier.Read_parameters();
            //  SpaceTime_Fourier.Initialize_engine();

            string Sqw_file = argv[4+(2*SpaceTime_Fourier.No_Of_Inputs)] ;
            string Dqw_file = argv[5+(2*SpaceTime_Fourier.No_Of_Inputs)] ;
            //SpaceTime_Fourier.Calculate_Sqw_using_Fwq_with_negativeandpositive_Time(Sqw_file, Dqw_file);
            SpaceTime_Fourier.Calculate_Sqw_using_Fwq(Sqw_file, Dqw_file);
            //SpaceTime_Fourier.Calculate_Sqw_using_Aq_Fwq(Sqw_file, Dqw_file);



        }
    }

    if(ex_string == "_rule"){

        if(model_=="1orb_MCMF"){


            Parameters_MCMF Parameters_;
            Parameters_.Initialize(input);

            Coordinates_MCMF Coordinates_(Parameters_.lx, Parameters_.ly);

            mt19937_64 Generator_(Parameters_.RandomSeed);
            MFParams_MCMF MFParams_(Parameters_,Coordinates_,Generator_);

            Hamiltonian_MCMF Hamiltonian_(Parameters_,Coordinates_,MFParams_);
            Observables_MCMF Observables_(Parameters_,Coordinates_,MFParams_,Hamiltonian_);


            SC_SW_ENGINE_VNE_1orb_MCMF Skw_Engine_(Parameters_,Coordinates_,MFParams_,Hamiltonian_,Observables_);
            Skw_Engine_.Read_parameters(input);
            Skw_Engine_.Initialize_engine();



            ST_Fourier_1orb_MCMF SpaceTime_Fourier(Parameters_,Coordinates_,MFParams_,Hamiltonian_,Observables_, Skw_Engine_);
            SpaceTime_Fourier.Space_Fourier_using_single_S=false;

            string No_of_inputs= argv[3];

            stringstream No_of_inputs_ss(No_of_inputs, stringstream::in);
            No_of_inputs_ss>>SpaceTime_Fourier.No_Of_Inputs;


            if( (SpaceTime_Fourier.No_Of_Inputs==1) && (!SpaceTime_Fourier.Space_Fourier_using_single_S) ){
                cout << "Use SpaceTime_Fourier.Space_Fourier_using_single_S=true"<<endl;
                assert(SpaceTime_Fourier.Space_Fourier_using_single_S);
            }


            SpaceTime_Fourier.Aq_inputs.resize(SpaceTime_Fourier.No_Of_Inputs);

            for(int i=0;i<SpaceTime_Fourier.No_Of_Inputs;i++){
                SpaceTime_Fourier.Aq_inputs[i]=argv[i+4];
            }


#ifdef _OPENMP
            cout<<"Parallel threads are used"<<endl;
#endif
#ifndef _OPENMP
            cout<<"single thread is used"<<endl;
#endif

            SpaceTime_Fourier.Read_parameters();
            //  SpaceTime_Fourier.Initialize_engine();

            //string Sqw_file_in, string Sq_static_out, string Sq_dynamical_out
            string Sqw_file_in = argv[4+(SpaceTime_Fourier.No_Of_Inputs)] ;
            string Sq_static_file = argv[5+(SpaceTime_Fourier.No_Of_Inputs)] ;
            string Sq_dynamical_file = argv[6+(SpaceTime_Fourier.No_Of_Inputs)] ;
            SpaceTime_Fourier.Sum_Rule_For_Way1(Sqw_file_in, Sq_static_file, Sq_dynamical_file);

        }

        if(model_=="MultiOrbSF"){


            Parameters_MultiOrbSF Parameters_;
            Parameters_.Initialize(input);

            Coordinates_MultiOrbSF Coordinates_(Parameters_.lx, Parameters_.ly, Parameters_.n_orbs);
            Coordinates_MultiOrbSF Coordinatestemp_(Parameters_.lx, Parameters_.ly, Parameters_.n_orbs);

            mt19937_64 Generator_(Parameters_.RandomSeed);
            mt19937_64 Generator2_(Parameters_.RandomDisorderSeed);
            MFParams_MultiOrbSF MFParams_(Parameters_,Coordinates_,Generator_, Generator2_);

            Hamiltonian_MultiOrbSF Hamiltonian_(Parameters_,Coordinates_, Coordinatestemp_,MFParams_);
            Observables_MultiOrbSF Observables_(Parameters_,Coordinates_,MFParams_,Hamiltonian_);


            SC_SW_ENGINE_VNE_MultiOrbSF Skw_Engine_(Parameters_,Coordinates_,MFParams_,Hamiltonian_,Observables_);
            Skw_Engine_.Read_parameters(input);
            Skw_Engine_.Initialize_engine();



            ST_Fourier_MultiOrbSF SpaceTime_Fourier(Parameters_,Coordinates_,MFParams_,Hamiltonian_,Observables_, Skw_Engine_);
            SpaceTime_Fourier.Space_Fourier_using_single_S=false;

            string No_of_inputs= argv[3];

            stringstream No_of_inputs_ss(No_of_inputs, stringstream::in);
            No_of_inputs_ss>>SpaceTime_Fourier.No_Of_Inputs;


            if( (SpaceTime_Fourier.No_Of_Inputs==1) && (!SpaceTime_Fourier.Space_Fourier_using_single_S) ){
                cout << "Use SpaceTime_Fourier.Space_Fourier_using_single_S=true"<<endl;
                assert(SpaceTime_Fourier.Space_Fourier_using_single_S);
            }


            SpaceTime_Fourier.Aq_inputs.resize(SpaceTime_Fourier.No_Of_Inputs);

            for(int i=0;i<SpaceTime_Fourier.No_Of_Inputs;i++){
                SpaceTime_Fourier.Aq_inputs[i]=argv[i+4];
            }


#ifdef _OPENMP
            cout<<"Parallel threads are used"<<endl;
#endif
#ifndef _OPENMP
            cout<<"single thread is used"<<endl;
#endif

            SpaceTime_Fourier.Read_parameters();
            //  SpaceTime_Fourier.Initialize_engine();

            //string Sqw_file_in, string Sq_static_out, string Sq_dynamical_out
            string Sqw_file_in = argv[4+(SpaceTime_Fourier.No_Of_Inputs)] ;
            string Sq_static_file = argv[5+(SpaceTime_Fourier.No_Of_Inputs)] ;
            string Sq_dynamical_file = argv[6+(SpaceTime_Fourier.No_Of_Inputs)] ;
            SpaceTime_Fourier.Sum_Rule_For_Way1(Sqw_file_in, Sq_static_file, Sq_dynamical_file);

        }

    }


    if(ex_string == "ion_w"){
        if(model_=="1orb_MCMF"){


            Parameters_MCMF Parameters_;
            Parameters_.Initialize(input);

            Coordinates_MCMF Coordinates_(Parameters_.lx, Parameters_.ly);

            mt19937_64 Generator_(Parameters_.RandomSeed);
            MFParams_MCMF MFParams_(Parameters_,Coordinates_,Generator_);

            Hamiltonian_MCMF Hamiltonian_(Parameters_,Coordinates_,MFParams_);
            Observables_MCMF Observables_(Parameters_,Coordinates_,MFParams_,Hamiltonian_);


            SC_SW_ENGINE_VNE_1orb_MCMF Skw_Engine_(Parameters_,Coordinates_,MFParams_,Hamiltonian_,Observables_);
            Skw_Engine_.Read_parameters(input);
            Skw_Engine_.Initialize_engine();



            ST_Fourier_1orb_MCMF SpaceTime_Fourier(Parameters_,Coordinates_,MFParams_,Hamiltonian_,Observables_, Skw_Engine_);
            SpaceTime_Fourier.Space_Fourier_using_single_S=false;



            SpaceTime_Fourier.Read_parameters();
            //  SpaceTime_Fourier.Initialize_engine();

            string Sqw_file_in = argv[3];
            string Sqw_file_out = argv[4];


            SpaceTime_Fourier.Convolute_the_spectrum(Sqw_file_in, Sqw_file_out);




        }

        if(model_=="MultiOrbSF"){


            Parameters_MultiOrbSF Parameters_;
            Parameters_.Initialize(input);

            Coordinates_MultiOrbSF Coordinates_(Parameters_.lx, Parameters_.ly, Parameters_.n_orbs);
            Coordinates_MultiOrbSF Coordinatestemp_(Parameters_.lx, Parameters_.ly, Parameters_.n_orbs);


            mt19937_64 Generator_(Parameters_.RandomSeed);
            mt19937_64 Generator2_(Parameters_.RandomDisorderSeed);
            MFParams_MultiOrbSF MFParams_(Parameters_,Coordinates_,Generator_, Generator2_);

            Hamiltonian_MultiOrbSF Hamiltonian_(Parameters_,Coordinates_,Coordinatestemp_, MFParams_);
            Observables_MultiOrbSF Observables_(Parameters_,Coordinates_,MFParams_,Hamiltonian_);


            SC_SW_ENGINE_VNE_MultiOrbSF Skw_Engine_(Parameters_,Coordinates_,MFParams_,Hamiltonian_,Observables_);
            Skw_Engine_.Read_parameters(input);
            Skw_Engine_.Initialize_engine();



            ST_Fourier_MultiOrbSF SpaceTime_Fourier(Parameters_,Coordinates_,MFParams_,Hamiltonian_,Observables_, Skw_Engine_);
            SpaceTime_Fourier.Space_Fourier_using_single_S=false;



            SpaceTime_Fourier.Read_parameters();
            //  SpaceTime_Fourier.Initialize_engine();

            string Sqw_file_in = argv[3];
            string Sqw_file_out = argv[4];


            SpaceTime_Fourier.Convolute_the_spectrum(Sqw_file_in, Sqw_file_out);


        }
    }





    return 0;
}
