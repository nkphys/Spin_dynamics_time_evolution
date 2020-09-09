#include <iostream>  //for cin and cout
#include <math.h>  // for pow
#include <stdlib.h>  //for div(q,n).rem(quot),abs(int n)
#include <time.h>
#include <fstream>
#include <sstream>
//#include "Spin_dynamics_ED_engine.h"
#include "src/Model_1orbHubbard_old/Basis_1orb_SF.h"
#include "src/Model_1orbHubbard_old/Model_1orb_SF.h"
#include "src/Model_1orbHubbard_old/Spin_dynamics_VNE_1orbHubbard_engine.h"
#include "src/Model_1orbHubbard_old/Space_Time_FourierTransform.h"
#include "src/Model_3orbPnictides/Coordinates.h"
#include "src/Model_3orbPnictides/Hamiltonian.h"
#include "Matrix.h"
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
            ){
        cout<<"You are using "<<argv[1]<<" , this model is not present"<<endl;
        cout<<"Please something from following :"<<endl;
        cout<<"1orbHubard, 3orbPnictides"<<endl;
        assert(model_=="1orbHubard" || model_=="3orbPnictides" || model_=="1orb_MCMF");
    }
    cout<<ex_string_original<<endl;

    string ex_string;

    //ex_string.substr(ex_string_original.length()-5);

    ex_string=ex_string_original.substr (ex_string_original.length() - 5);
    cout<<ex_string<<endl;

    if(ex_string == "amics"){


        bool S_kw_using_ED= false;
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

                if(!Skw_Engine_.RESTART){
                    Skw_Engine_.Read_equilibrium_configuration();
                }

                Skw_Engine_.Start_Engine();

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
            SpaceTime_Fourier.Initialize_engine();

            //SpaceTime_Fourier.Perform_Smarter_Averaging_on_one_point();


            ostringstream ostr_w_conv;
            ostr_w_conv << SpaceTime_Fourier.w_conv;
            string string_w_conv = ostr_w_conv.str();

            string Fw_file = "Fw_w_conv" + string_w_conv + "_state_" + MS_tag +".txt";
            string Aq_file = "Aq_state_" + MS_tag +".txt" ;
            SpaceTime_Fourier.Calculate_Fw_and_Aq(Fw_file, Aq_file);



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
            SpaceTime_Fourier.Calculate_Sqw_using_Aq_Fwq(Sqw_file, Dqw_file);




        }

    }






    return 0;
}
