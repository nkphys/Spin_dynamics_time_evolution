#include "Model_1orb_SF.h"
#include <stdlib.h>

using namespace std;
#define PI acos(-1.0)



void MODEL_1_orb_SF::Create_Hamil(BASIS_1_orb_SF &Basis, int time_step_){


    Hamil.clear();

    Jval_array_in_H=Jval_array;


   /*
    for(int i=0;i<Jval_array_in_H.size();i++){
        Jval_array_in_H[i]=0.05;
    }*/

    complex<double> one(1.0, 0.0);
    complex<double> iota(0.0, 1.0);

    Hamil.resize(2*(Basis.Max_pos + 1));
    for(int i=0;i<Hamil.size();i++){
        Hamil[i].resize(Hamil.size());
    }

    //tight_binding
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
    for(int pos=0;pos<=Basis.Max_pos;pos++){

        for(int neigh_no=0;neigh_no<Basis.Neighbour[pos].size();neigh_no++){

            //hopping for spin_up
            Hamil[pos][Basis.Neighbour[pos][neigh_no]].real( (-1.0)*t_hop );
            Hamil[pos][Basis.Neighbour[pos][neigh_no]].imag(0);



            //hopping for spin_dn
            Hamil[pos + Basis.Max_pos + 1][Basis.Neighbour[pos][neigh_no] +  Basis.Max_pos + 1].real( (-1.0)*t_hop );
            Hamil[pos + Basis.Max_pos + 1][Basis.Neighbour[pos][neigh_no] + Basis.Max_pos + 1].imag(0);



        }

    }


    //sz(i)Sz(i)
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
    for(int pos=0;pos<=Basis.Max_pos;pos++){


        Hamil[pos][pos] = Jval_array_in_H[pos]*one*S_mag*cos(Theta[time_step_][pos]);

        Hamil[pos + Basis.Max_pos + 1][pos + Basis.Max_pos + 1] = Jval_array_in_H[pos]*one*(-1.0)*S_mag*cos(Theta[time_step_][pos]);

    }



    //s_{plus}(i)S_{minus}(i) + s_{minus}(i)S_{plus}(i)
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
    for(int pos=0;pos<=Basis.Max_pos;pos++){


        Hamil[pos][pos + Basis.Max_pos + 1] = sin(Theta[time_step_][pos])*exp(iota * (-1.0) * Phi[time_step_][pos])*S_mag*Jval_array_in_H[pos];
        Hamil[pos + Basis.Max_pos + 1][pos] = sin(Theta[time_step_][pos])*exp(iota * (1.0) * Phi[time_step_][pos])*S_mag*Jval_array_in_H[pos];


    }


    //Svec[i].Svec[j]
    /*
    complex<double> value(0,0);
    int pos_neigh;
    for(int pos=0;pos<=Basis.Max_pos;pos++){
        for(int neigh=0;neigh<Basis.Neighbour[pos].size();neigh++){
            pos_neigh = Basis.Neighbour[pos][neigh];

            value +=   0.5*Jval_p*one*S_mag*S_mag* ( cos(Theta[time_step_][pos])*cos(Theta[time_step_][pos_neigh]) +
                                                     (sin(Theta[time_step_][pos])*sin(Theta[time_step_][pos_neigh])*
                                                      cos( (Theta[time_step_][pos] - Theta[time_step_][pos_neigh]))
                                                      )
                                                     );

        }

    }

    for(int pos=0;pos<=Basis.Max_pos;pos++){
        Hamil[pos][pos] += value;
        Hamil[pos + Basis.Max_pos + 1][pos + Basis.Max_pos + 1] +=value;
    }
*/




}

void MODEL_1_orb_SF::Create_Hamil(BASIS_1_orb_SF &Basis, int time_step_, string str){



    assert (str=="biased_quantum_spins");

    double Jval_fake=2.0;

    Hamil.clear();


    complex<double> one(1.0, 0.0);
    complex<double> iota(0.0, 1.0);

    Hamil.resize(2*(Basis.Max_pos + 1));
    for(int i=0;i<Hamil.size();i++){
        Hamil[i].resize(Hamil.size());
    }

    //tight_binding
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
    for(int pos=0;pos<=Basis.Max_pos;pos++){

        for(int neigh_no=0;neigh_no<Basis.Neighbour[pos].size();neigh_no++){

            //hopping for spin_up
            Hamil[pos][Basis.Neighbour[pos][neigh_no]].real( (-1.0)*t_hop );
            Hamil[pos][Basis.Neighbour[pos][neigh_no]].imag(0);



            //hopping for spin_dn
            Hamil[pos + Basis.Max_pos + 1][Basis.Neighbour[pos][neigh_no] +  Basis.Max_pos + 1].real( (-1.0)*t_hop );
            Hamil[pos + Basis.Max_pos + 1][Basis.Neighbour[pos][neigh_no] + Basis.Max_pos + 1].imag(0);



        }

    }


    //sz(i)Sz(i)
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
    for(int pos=0;pos<=Basis.Max_pos;pos++){

        Hamil[pos][pos] = Jval_fake*one*S_mag*cos(Theta[time_step_][pos]);

        Hamil[pos + Basis.Max_pos + 1][pos + Basis.Max_pos + 1] = Jval_fake*one*(-1.0)*S_mag*cos(Theta[time_step_][pos]);

    }



    //s_{plus}(i)S_{minus}(i) + s_{minus}(i)S_{plus}(i)
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
    for(int pos=0;pos<=Basis.Max_pos;pos++){

        Hamil[pos][pos + Basis.Max_pos + 1] = sin(Theta[time_step_][pos])*exp(iota * (-1.0) * Phi[time_step_][pos])*S_mag*Jval;
        Hamil[pos + Basis.Max_pos + 1][pos] = sin(Theta[time_step_][pos])*exp(iota * (1.0) * Phi[time_step_][pos])*S_mag*Jval;

    }


    //Svec[i].Svec[j]
    /*
    complex<double> value(0,0);
    int pos_neigh;
    for(int pos=0;pos<=Basis.Max_pos;pos++){
        for(int neigh=0;neigh<Basis.Neighbour[pos].size();neigh++){
            pos_neigh = Basis.Neighbour[pos][neigh];

            value +=   0.5*Jval_p*one*S_mag*S_mag* ( cos(Theta[time_step_][pos])*cos(Theta[time_step_][pos_neigh]) +
                                                     (sin(Theta[time_step_][pos])*sin(Theta[time_step_][pos_neigh])*
                                                      cos( (Theta[time_step_][pos] - Theta[time_step_][pos_neigh]))
                                                      )
                                                     );

        }

    }

    for(int pos=0;pos<=Basis.Max_pos;pos++){
        Hamil[pos][pos] += value;
        Hamil[pos + Basis.Max_pos + 1][pos + Basis.Max_pos + 1] +=value;
    }
*/




}


void MODEL_1_orb_SF::Diagonalize_Hamil(){

    Diagonalize(Hamil, Evals, Eigvecs);


}

void MODEL_1_orb_SF::Get_quantum_Spins(int ts){

    quant_s_x.resize(ts +1);
    quant_s_y.resize(ts +1);
    quant_s_z.resize(ts +1);

    int size = (int) ( (Hamil.size()*0.5) + 0.5 );

    quant_s_x[ts].resize(size);
    quant_s_y[ts].resize(size);
    quant_s_z[ts].resize(size);


    complex<double> zero(0.0,0.0);
    complex<double> iota(0.0,1.0);
    complex<double> splus, sminus, sz;
    complex<double> temp;

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(splus,sminus,sz,temp)
#endif
    for(int pos=0;pos<size;pos++){

        splus=zero;
        sminus=zero;
        sz=zero;
        for(int n=0;n<Evals.size();n++){


            splus += conj(Eigvecs[n][pos])*Eigvecs[n][pos+size]*fermi(Evals,n,Mu_);
            sz +=  ((conj(Eigvecs[n][pos])*Eigvecs[n][pos]) - (conj(Eigvecs[n][pos+size])*Eigvecs[n][pos+size]))*fermi(Evals,n,Mu_);
        }

        sminus = conj(splus);

        temp=(splus + sminus);
        quant_s_x[ts][pos] = real(temp);

        temp =iota*(-1.0)*(splus - sminus);
        quant_s_y[ts][pos] = real(temp);

        quant_s_z[ts][pos] = real(sz);

    }


}


void MODEL_1_orb_SF::Get_red_den_mat(int ts){

    //Red_Den_mat


    int size = (int) ( (Hamil.size()*0.5) + 0.5 ); //total no. of sites




    complex<double> zero(0.0,0.0);
    complex<double> iota(0.0,1.0);
    complex<double> temp;

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
    for(int pos_i=0;pos_i<size;pos_i++){
        for(int pos_j=0;pos_j<size;pos_j++){
            for(int s_i=0;s_i<2;s_i++){
                for(int s_j=0;s_j<2;s_j++){

                    Red_Den_mat[pos_i][s_i][pos_j][s_j]=zero;
                    for(int n=0;n<Evals.size();n++){
                        Red_Den_mat[pos_i][s_i][pos_j][s_j]+=
                                conj(Eigvecs[n][pos_j + s_j*(size)])*Eigvecs[n][pos_i + s_i*(size)]*fermi(Evals,n,Mu_);
                    }
                }
            }
        }

    }
}



double MODEL_1_orb_SF::Calculate_mu(Mat_1_doub &Evals_Temp){

    double mu_old,mu_new;

    double N_diff;
    double N_eps;
    N_eps=10e-5;


    double N_total_temp;
    double d_mu_by_dn;
    d_mu_by_dn=0.05*(t_hop*4.0)/(1.0*Evals_Temp.size());


    mu_new = (Evals_Temp[N_total - 1] +  Evals_Temp[N_total])*0.5; //start with mu for Temp=0
    N_diff=100.0;

    while(abs(N_diff)>N_eps)
    {
        mu_old=mu_new;
        N_total_temp=0.0;
        for(int i=0;i<Evals_Temp.size();i++){
            N_total_temp += fermi(Evals_Temp, i, mu_old);
        }

        N_diff = (1.0*N_total) - N_total_temp;


        mu_new = mu_old + d_mu_by_dn*(N_diff);

    }

    return mu_old;
}




double MODEL_1_orb_SF::fermi(Mat_1_doub &Evals_Temp,int n, double mu){


    double val;
    val = 1.0/( exp(Beta*(Evals_Temp[n] - mu))  + 1.0);

    return val;

}


void MODEL_1_orb_SF::Create_Hamil(BASIS_1_orb_SF &Basis, Mat_1_doub &Phi_vec,
                                  Mat_1_doub &delta_phi_vec, Mat_1_doub &Theta_vec,
                                  Mat_1_doub &delta_theta_vec, Mat_2_Complex_doub &Hamil_Temp, double factor){



    //    cout<<"Model created for "<<t_hop<<" "<<Jval<<"  "<<Jval_p << "  "<<S_mag<<endl;

    Hamil_Temp.clear();


    complex<double> one(1.0, 0.0);
    complex<double> iota(0.0, 1.0);

    Hamil_Temp.resize(2*(Basis.Max_pos + 1));
    for(int i=0;i<Hamil.size();i++){
        Hamil_Temp[i].resize(Hamil_Temp.size());
    }

    //tight_binding
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
    for(int pos=0;pos<=Basis.Max_pos;pos++){

        for(int neigh_no=0;neigh_no<Basis.Neighbour[pos].size();neigh_no++){

            //hopping for spin_up
            Hamil_Temp[pos][Basis.Neighbour[pos][neigh_no]].real( (-1.0)*t_hop );
            Hamil_Temp[pos][Basis.Neighbour[pos][neigh_no]].imag(0);



            //hopping for spin_dn
            Hamil_Temp[pos + Basis.Max_pos + 1][Basis.Neighbour[pos][neigh_no] +  Basis.Max_pos + 1].real( (-1.0)*t_hop );
            Hamil_Temp[pos + Basis.Max_pos + 1][Basis.Neighbour[pos][neigh_no] + Basis.Max_pos + 1].imag(0);



        }

    }


    //sz(i)Sz(i)
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
    for(int pos=0;pos<=Basis.Max_pos;pos++){

        Hamil_Temp[pos][pos] = Jval*one*S_mag*cos(Theta_vec[pos]+(factor*delta_theta_vec[pos]));

        Hamil_Temp[pos + Basis.Max_pos + 1][pos + Basis.Max_pos + 1] = Jval*one*(-1.0)*S_mag*cos(Theta_vec[pos]+(factor*delta_theta_vec[pos]));

    }



    //s_{plus}(i)S_{minus}(i) + s_{minus}(i)S_{plus}(i)
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
    for(int pos=0;pos<=Basis.Max_pos;pos++){

        Hamil_Temp[pos][pos + Basis.Max_pos + 1] =
                sin(Theta_vec[pos]+(factor*delta_theta_vec[pos]))*exp(iota * (-1.0) * (Phi_vec[pos]+(factor*delta_phi_vec[pos])))*S_mag*Jval;
        Hamil_Temp[pos + Basis.Max_pos + 1][pos] =
                sin(Theta_vec[pos]+(factor*delta_theta_vec[pos]))*exp(iota * (1.0) * (Phi_vec[pos]+(factor*delta_phi_vec[pos])))*S_mag*Jval;

    }




}


void MODEL_1_orb_SF::Diagonalize_Hamil(Mat_2_Complex_doub &Hamil_Temp, Mat_1_doub &Evals_Temp, Mat_2_Complex_doub &Eigvecs_Temp){

    Diagonalize(Hamil_Temp, Evals_Temp, Eigvecs_Temp);


}

void MODEL_1_orb_SF::Get_quantum_Spins(Mat_1_doub &quant_s_y_new,Mat_1_doub &quant_s_x_new,Mat_1_doub &quant_s_z_new,
                                       Mat_2_Complex_doub &Hamil_Temp, Mat_2_Complex_doub &Eigvecs_Temp,Mat_1_doub &Evals_Temp, double mu){



    int size = (int) ( (Hamil_Temp.size()*0.5) + 0.5 );

    quant_s_y_new.resize(size);
    quant_s_x_new.resize(size);
    quant_s_z_new.resize(size);


    complex<double> zero(0.0,0.0);
    complex<double> iota(0.0,1.0);
    complex<double> splus, sminus, sz;
    complex<double> temp;

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(splus,sminus,sz,temp)
#endif
    for(int pos=0;pos<size;pos++){

        splus=zero;
        sminus=zero;
        sz=zero;
        for(int n=0;n<Evals_Temp.size();n++){


            splus += conj(Eigvecs_Temp[n][pos])*Eigvecs_Temp[n][pos+size]*fermi(Evals_Temp, n,mu);
            sz +=  ((conj(Eigvecs_Temp[n][pos])*Eigvecs_Temp[n][pos]) - (conj(Eigvecs_Temp[n][pos+size])*Eigvecs_Temp[n][pos+size]))*fermi(Evals_Temp,n,mu);
        }


        sminus = conj(splus);

        temp=(splus + sminus);
        quant_s_x_new[pos] = real(temp);

        temp =iota*(-1.0)*(splus - sminus);
        quant_s_y_new[pos] = real(temp);

        quant_s_z_new[pos] = real(sz);

    }


}
