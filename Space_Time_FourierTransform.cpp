#include "Space_Time_FourierTransform.h"
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;



void ST_Fourier::Initialize_engine(BASIS_1_orb_SF & Basis){

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



    /*

    SS_rt.resize(Basis.Max_pos +1);
    sqsq_rt.resize(Basis.Max_pos +1);
    TT_rt.resize(Basis.Max_pos +1);
    for(int pos_i=0;pos_i<Basis.Max_pos +1;pos_i++){
        SS_rt[pos_i].resize(Basis.Max_pos +1);
        sqsq_rt[pos_i].resize(Basis.Max_pos +1);
        TT_rt[pos_i].resize(Basis.Max_pos +1);
        for(int pos_j=0;pos_j<Basis.Max_pos +1;pos_j++){
            SS_rt[pos_i][pos_j].resize(time_steps);
            sqsq_rt[pos_i][pos_j].resize(time_steps);
            TT_rt[pos_i][pos_j].resize(time_steps);
            for(int ti=0;ti<time_steps;ti++){
                SS_rt[pos_i][pos_j][ti]=0;
                sqsq_rt[pos_i][pos_j][ti]=0;
                TT_rt[pos_i][pos_j][ti]=0;
            }
        }
    }




    S_rt.resize(Basis.Max_pos +1);
    sq_rt.resize(Basis.Max_pos +1);
    T_rt.resize(Basis.Max_pos +1);
    for(int pos_i=0;pos_i<Basis.Max_pos +1;pos_i++){
            S_rt[pos_i].resize(time_steps);
            sq_rt[pos_i].resize(time_steps);
            T_rt[pos_i].resize(time_steps);
            for(int ti=0;ti<time_steps;ti++){
                S_rt[pos_i][ti]=0;
                sq_rt[pos_i][ti]=0;
                T_rt[pos_i][ti]=0;
            }

    }
*/




    /*
    Model.Theta.resize(Runge_Kutta_order+1);
    Model.Phi.resize(Runge_Kutta_order+1);


    Model.Theta_eq.resize(Basis.Max_pos + 1);
    Model.Phi_eq.resize(Basis.Max_pos + 1);

    for(int i=0;i<Model.Theta.size();i++){

        Model.Theta[i].resize(Basis.Max_pos + 1);
        Model.Phi[i].resize(Basis.Max_pos + 1);
    }

    S_kw.resize(Basis.Max_pos + 1 );
    for(int i=0;i<=Basis.Max_pos;i++){
        S_kw[i].resize(n_wpoints);

    }



    Model.quant_s_x.resize(1);
    Model.quant_s_y.resize(1);
    Model.quant_s_z.resize(1);

    Model.Red_Den_mat.resize(Basis.Max_pos + 1);
    Model.quant_s_x_eq.resize(Basis.Max_pos + 1);
    Model.quant_s_y_eq.resize(Basis.Max_pos + 1);
    Model.quant_s_z_eq.resize(Basis.Max_pos + 1);
    Model.quant_s_x[0].resize(Basis.Max_pos + 1);
    Model.quant_s_y[0].resize(Basis.Max_pos + 1);
    Model.quant_s_z[0].resize(Basis.Max_pos + 1);
    for(int i=0;i<=Basis.Max_pos;i++){
        Model.Red_Den_mat[i].resize(2);
        for(int si=0;si<2;si++){
            Model.Red_Den_mat[i][si].resize(Basis.Max_pos + 1);
            for(int j=0;j<=Basis.Max_pos;j++){
                Model.Red_Den_mat[i][si][j].resize(2);

            }

        }

    }





    Model.center = int((Basis.Max_pos + 1)*0.5 + 0.5) - 1;



    Model.Jval_array.resize(Basis.Max_pos + 1);
    Model.Bval_array.resize(Basis.Max_pos + 1);


    if(Model.Single_classical_S_in_B==true){
        cout<<"only one classical spin is coupled with quantum bath"<<endl;
        for(int pos=0;pos<Basis.Max_pos + 1;pos++){

            if(pos== Model.center){
                Model.Jval_array[pos]=Model.Jval;
                Model.Bval_array[pos]=Model.B_mag;
            }
            else{
                Model.Jval_array[pos]=0;
                Model.Bval_array[pos]=0;
            }
        }
        cout<<"Magnetic Field at site = "<<Model.center<<endl;
    }
    else{
        for(int pos=0;pos<Basis.Max_pos + 1;pos++){
            Model.Jval_array[pos]=Model.Jval;
            Model.Bval_array[pos]=0.0;
        }

    }


*/



}
void ST_Fourier::Perform_Averaging_on_one_point(BASIS_1_orb_SF & Basis){

    string Temp_file_Srt="Temp_file_Srt.txt";
    ofstream Temp_file_Srt_out(Temp_file_Srt.c_str());


    Temp_file_Srt_out<<"#time   Sz   Sx  Sy  sz  sx  sy"<<endl;

    string line_temp;
    double row_time, double_temp;
    int row_ts;

    Mat_1_doub Sr;
    Sr.resize(6*Basis.Lx*Basis.Ly);


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

            for(int r=0;r<6*Basis.Lx*Basis.Ly;r++){
                line_temp_ss>>double_temp;

                Sr[r] += double_temp;
            }

            Specific_conf_in.clear();
            Specific_conf_in.seekg(0, ios::beg);


        }

        for(int r=0;r<6*Basis.Lx*Basis.Ly;r++){
            Temp_file_Srt_out<<Sr[r]/(No_Of_Inputs)<<"  ";
        }
        Temp_file_Srt_out<<endl;

    }


}


void ST_Fourier::Perform_Smarter_Averaging_on_one_point(BASIS_1_orb_SF & Basis){

    string Temp_file_Srt="Temp_file_Srt.txt";
    ofstream Temp_file_Srt_out(Temp_file_Srt.c_str());


    Temp_file_Srt_out<<"#time   Sz   Sx  Sy  sz  sx  sy"<<endl;

    string line_temp;
    double row_time, double_temp;
    int row_ts;


    S_tr.resize(time_steps+1);

    for(int ts=0;ts<=time_steps;ts++){
        S_tr[ts].resize(6*Basis.Lx*Basis.Ly);
        for(int r=0;r<6*Basis.Lx*Basis.Ly;r++){
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

            for(int r=0;r<6*Basis.Lx*Basis.Ly;r++){
                line_temp_ss>>double_temp;

                S_tr[ts][r]+= double_temp;
            }


        }

    }



    for(int ts=0;ts<=time_steps;ts++){
        Temp_file_Srt_out<<ts*dt_<<"  ";
        for(int r=0;r<6*Basis.Lx*Basis.Ly;r++){
            Temp_file_Srt_out<<S_tr[ts][r]/(No_Of_Inputs)<<"  ";
        }
        Temp_file_Srt_out<<endl;
    }


}

void ST_Fourier::Calculate_Skw_from_Srt_file(BASIS_1_orb_SF & Basis, string filename, string fileout){


    S_rw.resize(Basis.Max_pos +1);
    s_quantum_rw.resize(Basis.Max_pos +1);
    T_rw.resize(Basis.Max_pos +1);

    for(int pos_i=0;pos_i<Basis.Max_pos +1;pos_i++){
        S_rw[pos_i].resize(Basis.Max_pos +1);
        s_quantum_rw[pos_i].resize(Basis.Max_pos +1);
        T_rw[pos_i].resize(Basis.Max_pos +1);

        for(int pos_j=0;pos_j<Basis.Max_pos +1;pos_j++){
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

    Sz_eq.resize(Basis.Max_pos+1);Sx_eq.resize(Basis.Max_pos+1);Sy_eq.resize(Basis.Max_pos+1);
    sz_eq.resize(Basis.Max_pos+1);sx_eq.resize(Basis.Max_pos+1);sy_eq.resize(Basis.Max_pos+1);

    Sz_t.resize(Basis.Max_pos+1);Sx_t.resize(Basis.Max_pos+1);Sy_t.resize(Basis.Max_pos+1);
    sz_t.resize(Basis.Max_pos+1);sx_t.resize(Basis.Max_pos+1);sy_t.resize(Basis.Max_pos+1);


    string line_temp;
    double temp_waste;
    ifstream file_Srt_in;
    file_Srt_in.open(filename.c_str());

    getline(file_Srt_in, line_temp);//First line is a waste
    getline(file_Srt_in, line_temp);
    stringstream line_temp_ss(line_temp, stringstream::in);

    line_temp_ss>>temp_waste;

    for(int pos_i=0;pos_i<Basis.Max_pos +1;pos_i++){
        line_temp_ss>>Sz_eq[pos_i];
    }
    for(int pos_i=0;pos_i<Basis.Max_pos +1;pos_i++){
        line_temp_ss>>Sx_eq[pos_i];
    }
    for(int pos_i=0;pos_i<Basis.Max_pos +1;pos_i++){
        line_temp_ss>>Sy_eq[pos_i];
    }
    for(int pos_i=0;pos_i<Basis.Max_pos +1;pos_i++){
        line_temp_ss>>sz_eq[pos_i];
    }
    for(int pos_i=0;pos_i<Basis.Max_pos +1;pos_i++){
        line_temp_ss>>sx_eq[pos_i];
    }
    for(int pos_i=0;pos_i<Basis.Max_pos +1;pos_i++){
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
        for(int pos_i=0;pos_i<Basis.Max_pos +1;pos_i++){
            line_temp2_ss>>Sz_t[pos_i];
        }
        for(int pos_i=0;pos_i<Basis.Max_pos +1;pos_i++){
            line_temp2_ss>>Sx_t[pos_i];
        }
        for(int pos_i=0;pos_i<Basis.Max_pos +1;pos_i++){
            line_temp2_ss>>Sy_t[pos_i];
        }
        for(int pos_i=0;pos_i<Basis.Max_pos +1;pos_i++){
            line_temp2_ss>>sz_t[pos_i];
        }
        for(int pos_i=0;pos_i<Basis.Max_pos +1;pos_i++){
            line_temp2_ss>>sx_t[pos_i];
        }
        for(int pos_i=0;pos_i<Basis.Max_pos +1;pos_i++){
            line_temp2_ss>>sy_t[pos_i];
        }






#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
        for(int wi=0;wi<n_wpoints;wi++){

            for(int pos_i=0;pos_i<Basis.Max_pos +1;pos_i++){
                for(int pos_j=0;pos_j<Basis.Max_pos +1;pos_j++){



                    S_rw[pos_i][pos_j][wi] += cos((wi * dw) * (ts* dt_))*exp(-0.5*(ts*dt_*w_conv*ts*dt_*w_conv))*dt_*(
                                ( (Sz_eq[pos_i]*Sz_t[pos_j]  )  )
                                +
                                ( (Sx_eq[pos_i]*Sx_t[pos_j]  )  )
                                +
                                ( (Sy_eq[pos_i]*Sy_t[pos_j]  )  )
                                );


                    s_quantum_rw[pos_i][pos_j][wi] += cos((wi * dw) * (ts* dt_))*exp(-0.5*(ts*dt_*w_conv*ts*dt_*w_conv))*dt_*(
                                ( (sz_eq[pos_i]*( sz_t[pos_j] - 1.0*sz_eq[pos_j])  )  )
                                +
                                ( (sx_eq[pos_i]*( sx_t[pos_j] - 1.0*sx_eq[pos_j])  )  )
                                +
                                ( (sy_eq[pos_i]*( sy_t[pos_j] - 1.0*sy_eq[pos_j])  )  )
                                );




                    /*    T_rw[pos_i][pos_j][wi] +=exp(iota * (-wi * dw) * (ts* dt_))*exp(-0.5*(ts*dt_*w_conv*ts*dt_*w_conv))*dt_* (


                            );*/






                }
            }
        }

    }







    //Now "rw - space" to "kw - space"

    S_kw.resize(Basis.Lx);
    s_quantum_kw.resize(Basis.Lx);
    T_kw.resize(Basis.Lx);


    for(int x_i=0;x_i<Basis.Lx;x_i++){

        S_kw[x_i].resize(Basis.Ly);
        s_quantum_kw[x_i].resize(Basis.Ly);
        T_kw[x_i].resize(Basis.Ly);

        for(int y_j=0;y_j<Basis.Ly;y_j++){

            S_kw[x_i][y_j].resize(n_wpoints);
            s_quantum_kw[x_i][y_j].resize(n_wpoints);
            T_kw[x_i][y_j].resize(n_wpoints);

        }
    }





    double kx, ky;
    int pos_i, pos_j;
    for(int nx=0;nx<Basis.Lx;nx++){
        for(int ny=0;ny<Basis.Ly;ny++){

            kx = (2*nx*PI)/(1.0*Basis.Lx);
            ky = (2*ny*PI)/(1.0*Basis.Ly);


#ifdef _OPENMP
#pragma omp parallel for default(shared) private(pos_i,pos_j)
#endif

            for(int wi=0;wi<n_wpoints;wi++){
                double temp3=0;
                double temp2=0;
                double temp=0;
                for(int x_i=0;x_i<Basis.Lx;x_i++){
                    for(int y_i=0;y_i<Basis.Ly;y_i++){

                        pos_i = y_i*Basis.Lx + x_i;

                        for(int x_j=0;x_j<Basis.Lx;x_j++){
                            for(int y_j=0;y_j<Basis.Ly;y_j++){

                                pos_j = y_j*Basis.Lx + x_j;

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

    int k_ind=0;
    for(int nx=0;nx<Basis.Lx;nx++){
        for(int ny=0;ny<Basis.Ly;ny++){

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

void ST_Fourier::Read_parameters(string filename, BASIS_1_orb_SF & Basis){


    string Geometry_ = "Geometry = ";
    string boundary_conditions_, Boundary_Conditions_ = "Boundary conditions = ";

    string lx_, Lx_ = "Lx = ";
    string ly_, Ly_ = "Ly = ";
    //string no_of_inputs_, No_Of_Inputs_ = "Number_Of_Configurations = "




    string w_min_, W_Min_ = "w_min = ";
    string w_max_, W_Max_ = "w_max = ";
    string dw_, dW_ = "dw = ";
    string w_conv_, W_conv_ = "w_convolution = ";

    //string Conf_Input_ = "conf_input = ";

    string SKW_Out_ = "Skw_out = ";
    string SKW_Out_Full_ = "Skw_out_Full = ";

    string no_of_processors_, No_Of_Processors_ = "no_of_threads = ";


    int offset;
    string line;
    ifstream inputfile(filename.c_str());


    if(inputfile.is_open())
    {
        while(!inputfile.eof())
        {
            getline(inputfile,line);


            if ((offset = line.find(W_conv_, 0)) != string::npos) {
                w_conv_ = line.substr (offset + W_conv_.length());		}


            if ((offset = line.find(Geometry_, 0)) != string::npos) {
                Basis.Geometry = line.substr (offset + Geometry_.length());		}

            if ((offset = line.find(Boundary_Conditions_, 0)) != string::npos) {
                boundary_conditions_ = line.substr (offset + Boundary_Conditions_.length());		}

            if ((offset = line.find(Lx_, 0)) != string::npos) {
                lx_ = line.substr (offset + Lx_.length());		}

            if ((offset = line.find(Ly_, 0)) != string::npos) {
                ly_ = line.substr (offset + Ly_.length());		}

            //if ((offset = line.find(No_Of_Inputs_, 0)) != string::npos) {
            //  no_of_inputs_ = line.substr (offset + No_Of_Inputs_.length());		}

            if ((offset = line.find(W_Min_, 0)) != string::npos) {
                w_min_ = line.substr (offset + W_Min_.length());		}

            if ((offset = line.find(W_Max_, 0)) != string::npos) {
                w_max_ = line.substr (offset + W_Max_.length());		}

            if ((offset = line.find(dW_, 0)) != string::npos) {
                dw_= line.substr (offset + dW_.length());		}

            //  if ((offset = line.find(Conf_Input_, 0)) != string::npos) {
            //    conf_input = line.substr (offset + Conf_Input_.length());		}

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
    {cout<<"Unable to open input file while in the Model class."<<endl;}


    w_min=atof(w_min_.c_str());
    w_max=atof(w_max_.c_str());

    dw=atof(dw_.c_str());
    w_conv=atof(w_conv_.c_str());

    if (boundary_conditions_ == "PBC"){
        Basis.PBC =true;
    }
    else{
        Basis.PBC = false;
    }

    no_of_processors=atoi(no_of_processors_.c_str());


    Basis.Lx=atoi(lx_.c_str());
    Basis.Ly=atoi(ly_.c_str());

    //No_Of_Inputs=atoi(no_of_inputs_.c_str());



    cout<<"PARAMETERS::::::::"<<endl;
    cout<<"Boundary condition = "<<boundary_conditions_<<endl;
    cout<<"Lx = "<<Basis.Lx<<endl;
    cout<<"Ly = "<<Basis.Ly<<endl;
    cout<<"w_min = "<<w_min<<endl;
    cout<<"w_max = "<<w_max<<endl;
    cout<<"dw = "<<dw<<endl;
    cout<<"w_convolution = "<<w_conv<<endl;
    // cout<<"No of Configurations = "<< No_Of_Inputs<<endl;



}



































