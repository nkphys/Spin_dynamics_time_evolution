#include "functions.h"

extern "C" void   zheev_(char *,char *,int *,std::complex<double> *, int *, double *,
                         std::complex<double> *,int *, double *, int *);

extern "C" void   dsyev_(char *,char *,int *, double *, int *, double *,
                         double *,int *, int *);


void Print_vector_in_file(Mat_1_doub vec, string filename){

    ofstream outfile(filename.c_str());

    for(int j=0;j<vec.size();j++){
        outfile<<vec[j]<<endl;
    }

}



bool present_before(Mat_1_int nup_2, Mat_1_int ndn_2, Mat_2_int nup_2_group, Mat_2_int ndn_2_group, int &pos){


    bool check;
    check=false;


    if(nup_2_group.size()==0){
        pos=0;
        check=false;
    }

    else{
        bool c_up, c_dn, c_both;
        c_up=false;c_dn=false;

        for(int i=0;i<nup_2_group.size();i++){

            if(nup_2_group[i]==nup_2){
                c_up=true;
            }
            else{
                c_up=false;
            }
            if(ndn_2_group[i]==ndn_2){
                c_dn=true;
            }
            else{
                c_dn=false;
            }


            c_both= (c_up && c_dn) ;
            if(c_both==true){
                check=true;
                pos=i;
            }



        }
    }
    return check;


}



string NumberToString ( int Number )
{
    ostringstream ss;
    ss << Number;
    return ss.str();
}

static bool sort_using_greater_than(double u, double v)
{
    return u > v;
}

bool comp_greater(double i, double j){
    return i > j;
}


bool comp_greater_pair_double_int(pair_double_int i, pair_double_int j){
    return i.first > j.first;
}


int Find_int_in_intarray(int num, Mat_1_int &array){

    int pos;
    bool not_found=true;
    int ind=0;
    while(not_found){
        assert(ind<array.size());
        if(num==array[ind]){
            pos=ind;
            not_found = false;
        }
        ind++;
    }

    return pos;
}

void Print_Matrix_COO(Matrix_COO &A){

    Mat_2_doub B;
    B.resize(A.nrows);
    for(int i=0;i<B.size();i++){
        B[i].resize(A.nrows);
        for(int j=0;j<A.nrows;j++){
            B[i][j]=0.0;
        }
    }

    for(int i=0;i<A.value.size();i++){
        B[A.rows[i]][A.columns[i]]=A.value[i];
    }

    cout<<"--------------------PRINTING THE MATRIX:--------------------"<<endl;
    for(int i=0;i<B.size();i++){
        for(int j=0;j<B.size();j++){

            cout<<B[i][j]<<"  ";
        }
        cout<<endl;
    }
    cout<<"-------------------------------------------------------------"<<endl;


}


void Diagonalize(Matrix_COO &X, double & EG, Mat_1_doub & vecG){

    int LDA=X.nrows, info;

    /* Local arrays */
    double* eval = new double[X.nrows];

    double* mat = new double[X.nrows*X.nrows];


    for(int i=0;i<X.nrows;i++){
        for(int j=i;j<X.nrows;j++){
            mat[i*(X.nrows)+j] = 0;
        }
    }

    for(int i=0;i<X.value.size();i++){
        int r=X.rows[i];
        int c=X.columns[i];
        mat[r*(X.nrows)+c] = X.value[i];
        mat[c*(X.nrows)+r] = X.value[i];

    }

    char jobz='V';
    char uplo='U';
    int n=LDA;
    vector<double> work(3);
    int lwork= -1;

    //info=LAPACKE_dsyev(LAPACK_ROW_MAJOR,  'V', 'U', X.nrows, mat , LDA, eval );
    dsyev_(&jobz,&uplo,&n,&mat[0],&LDA, &eval[0],&(work[0]),&lwork,&info);
    lwork = int((work[0]));
    work.resize(lwork);
     dsyev_(&jobz,&uplo,&n,&mat[0],&LDA, &eval[0],&(work[0]),&lwork,&info);
    /* Check for convergence */
    if( info > 0 ) {
        cout<< "The LAPACKE_dsyev failed to diagonalize."<<endl;

    }



    EG=eval[0];
    free(eval);
    vecG.clear();
    vecG.resize(X.nrows);
    for(int i=0;i<X.nrows;i++){
        vecG[i] =mat[i*X.nrows];
    }
    free(mat);


}


double dot_product(Mat_1_doub &vec1, Mat_1_doub &vec2){
    //This dot_product is parallelized always, create another one with _PARALLELIZE_AT_MATRICES_LEVEL
    double temp;
    double temp1=0;
    assert(vec1.size()==vec2.size());
    // bool Parallelize_dot_product;
    // Parallelize_dot_product=true;

    // if(!Parallelize_dot_product){goto skiploop_143;}
    //#pragma omp parallel for default(shared) reduction(+:temp1)
    //skiploop_143:
    for(int i=0;i<vec1.size();i++){
        temp1 = temp1 + (vec1[i])*(vec2[i]);
    }

    temp = temp1;

    return temp;

}

void Matrix_COO_vector_multiplication(string COO_type, Matrix_COO & A,Mat_1_doub &u,Mat_1_doub &v){

    v.clear();
    v.resize(u.size());


    for (int i=0;i<v.size();i++){
        v[i]=0;}

    if(COO_type =="U"){
        for (int n=0;n<A.value.size();n++){

            if(A.rows[n]==A.columns[n]){
                v[A.rows[n]] = v[A.rows[n]] + u[A.columns[n]]*A.value[n];
            }
            else{
                v[A.rows[n]] = v[A.rows[n]] + u[A.columns[n]]*A.value[n];
                v[A.columns[n]] = v[A.columns[n]] + u[A.rows[n]]*A.value[n];
            }

        }
    }
    else{
        for (int n=0;n<A.value.size();n++){

            v[A.rows[n]] = v[A.rows[n]] + u[A.columns[n]]*A.value[n];



        }
    }


}



void Subtract( Mat_1_doub &temp1, double x, Mat_1_doub &temp2){

    assert(temp1.size()==temp2.size());
    for(int k=0;k<temp1.size();k++){

        temp1[k]=temp1[k] - x*temp2[k];

    }

}


void Diagonalize(Mat_2_Complex_doub &Matrix, Mat_1_doub & Evals, Mat_2_Complex_doub &Eigvecs){



    int size = Matrix.size();
    int LDA=size, info;
    /* Local arrays */
    double* evals = new double[size];

    //complex<double>* mat = new complex<double>[Matrix.size()*Matrix.size()];

    //MKL_Complex16 mat[size*size]; //HERE
    //complex<double>* mat = new complex<double>[size*size];

    complex<double>* mat = new complex<double>[size*size];


    //MKL_Complex16* mat[size*size]
    for(int i=0;i<size;i++){
        for(int j=0;j<=i;j++){

           // mat[i*(size)+j].real=real(Matrix[i][j]);
           // mat[i*(size)+j].imag=imag(Matrix[i][j]);
            mat[i*(size)+j] = Matrix[i][j];


        }
    }

   // info=LAPACKE_zheev(LAPACK_ROW_MAJOR,  'V', 'L', size, mat , LDA, evals);
//-------------
    char jobz='V';
    char uplo='L'; //WHY ONLY 'L' WORKS?
    int n=Matrix.size();
    int lda=Matrix.size();
    vector<complex<double>> work(3);
    vector<double> rwork(3*n -2);
    int lwork= -1;

    // query:
    zheev_(&jobz,&uplo,&n,&mat[0],&lda,&(evals[0]),&(work[0]),&lwork,&(rwork[0]),&info);
    //lwork = int(real(work[0]))+1;
    lwork = int((work[0].real()));
    work.resize(lwork);
    // real work:
    zheev_(&jobz,&uplo,&n,&mat[0],&lda,&(evals[0]),&(work[0]),&lwork,&(rwork[0]),&info);
    if (info!=0) {
        std::cerr<<"info="<<info<<"\n";
        perror("diag: zheev: failed with info!=0.\n");
    }
    //---------------



    /* Check for convergence */
    if( info > 0 ) {

        cout<< "The LAPACKE_dsyev failed to diagonalize."<<endl;

    }




    Evals.resize(size);
    for(int i=0;i<size;i++){
        Evals[i]=evals[i];
    }

    free(evals);

    Eigvecs.resize(size);
    for(int i=0;i<size;i++){
        Eigvecs[i].resize(size);
    }

    for(int i=0;i<size;i++){
        for(int j=0;j<size;j++){
            //Eigvecs[i][j].real(mat[j*size+i].real);
            //Eigvecs[i][j].imag(mat[j*size+i].imag);
            Eigvecs[i][j] = mat[j*size+i];
        }
    }
    //mkl_free(mat);


}




void Sum(Matrix_COO A, Matrix_COO B, Matrix_COO & C, double value1, double value2){

    int a_i=0;
    int b_j=0;
    int row_a_i,col_a_i, row_b_j,col_b_j;
    Matrix_COO temp;

    if((A.nrows == B.nrows) && (A.ncols == B.ncols)){

        temp.nrows = A.nrows;
        temp.ncols = A.ncols;

        temp.value.clear();
        temp.rows.clear();
        temp.columns.clear();

        while(a_i<A.value.size() || b_j<B.value.size()){

            if(a_i<A.value.size()){
                row_a_i=A.rows[a_i];
                col_a_i=A.columns[a_i];
            }
            else{
                row_a_i=A.nrows+1;
                col_a_i=A.ncols+1;
            }

            if(b_j<B.value.size()){
                row_b_j=B.rows[b_j];
                col_b_j=B.columns[b_j];
            }
            else{
                row_b_j=B.nrows+1;
                col_b_j=B.ncols+1;
            }

            //if element of B comes before element of A
            if( (row_b_j < row_a_i)  ||
                    (row_b_j == row_a_i && col_b_j < col_a_i)
                    )
            {

                temp.value.push_back(value2*B.value[b_j]);

                temp.rows.push_back(row_b_j);
                temp.columns.push_back(col_b_j);
                b_j=b_j+1;

            }
            //if element of A comes before element of B
            if( (row_b_j > row_a_i) ||
                    (row_b_j == row_a_i && col_b_j > col_a_i)
                    )

            {
                temp.value.push_back(value1*A.value[a_i]);
                temp.rows.push_back(row_a_i);
                temp.columns.push_back(col_a_i);
                a_i=a_i+1;

            }
            //if elements of B, A comes at same place
            if(row_b_j==row_a_i && col_b_j==col_a_i){
                temp.value.push_back(value1*A.value[a_i] + value2*B.value[b_j]);
                temp.rows.push_back(row_a_i);
                temp.columns.push_back(col_a_i);
                a_i=a_i+1;
                b_j=b_j+1;
            }
        }

    }
    else{cout<<"Error in doing Sum"<<endl;}

    C=temp;
    temp.value.clear();
    temp.rows.clear();
    temp.columns.clear();
}

void Calculate_recursive_GF(Mat_1_doub A, Mat_1_doub B2, complex<double> &Recursive_GF, double omega,
                            double eta, double GS_energy){

    complex<double> temp_n;
    temp_n.imag(0);temp_n.real(0);

    double val1,val2;
    omega = omega + GS_energy;



    for (int i=A.size()-1;i>=0;i--){
        if(i==(A.size()-1)){
            val1 = omega - A[i];
            val2 = eta;

            temp_n.real(val1/( (val1*val1) + (val2*val2) ) );
            temp_n.imag(-val2/( (val1*val1) + (val2*val2) ) );

        }
        else{
            val1= omega -A[i] - (B2[i+1]*temp_n.real());
            val2= eta - (B2[i+1]*temp_n.imag());

            temp_n.real(val1/( (val1*val1) + (val2*val2) ) );
            temp_n.imag(-val2/( (val1*val1) + (val2*val2) ) );

        }

    }


    Recursive_GF = temp_n;


}


void Get_NI_Skw(string filename_full , string filename_specific_k , int dim, double thop){


    ofstream fileout_full(filename_full.c_str());

    complex<double> zero(0,0);
    complex<double> one(1,0);
    complex<double> iota(0,1);
    double PI_ = acos(-1.0);


    double mu=0;//for half-filling


    double dw=0.01;
    double w_min=-0.02;
    double w_max=4.1;
    int n_wpoints;
    n_wpoints=(int) ((w_max-w_min)/dw + 0.5);
    double BETA=500;
    int Lx=1;
    int Ly=32;
    Mat_3_Complex_doub Srw;
    Srw.resize(n_wpoints);

    for(int nw =0;nw<n_wpoints;nw++){
        Srw[nw].resize(Lx*Ly);
        for(int i =0;i<Lx*Ly;i++){
            Srw[nw][i].resize(Lx*Ly);
        }
    }

    //pos_i = y_i*Basis.Lx + x_i;


    Mat_2_Complex_doub Pauli_x,Pauli_y,Pauli_z;
    Pauli_x.resize(2);Pauli_y.resize(2);Pauli_z.resize(2);

    for(int i=0;i<2;i++){
        Pauli_x[i].resize(2); Pauli_y[i].resize(2);Pauli_z[i].resize(2);
        for(int j=0;j<2;j++){
            Pauli_x[i][j]=0;Pauli_y[i][j]=0;Pauli_z[i][j]=0;
        }
    }

    Pauli_x[0][1]=1.0;Pauli_x[1][0]=1.0;
    Pauli_y[0][1]=-1.0*iota;Pauli_y[1][0]=1.0*iota;
    Pauli_z[0][0]=1.0;Pauli_z[1][1]=-1.0;



    if(dim==2){

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
        for(int nw =0;nw<n_wpoints;nw++){
            double w_val=nw*dw;

            //pos_j
            for(int xj=0;xj<Lx;xj++){
                for(int yj=0;yj<Ly;yj++){
                    int pos_j = yj*Lx + xj;

                    //pos_l
                    for(int xl=0;xl<Lx;xl++){
                        for(int yl=0;yl<Ly;yl++){
                            int pos_l = yl*Lx + xl;

                            Srw[nw][pos_j][pos_l]=zero;


                            //qvec
                            for(int nx=0;nx<Lx;nx++){
                                double qx= (2.0*PI_*nx)/(1.0*Lx);
                                for(int ny=0;ny<Ly;ny++){
                                double qy= (2.0*PI_*ny)/(1.0*Ly);

                                    //qvec_p
                                    for(int nx_p=0;nx_p<Lx;nx_p++){
                                        double qx_p= (2.0*PI_*nx_p)/(1.0*Lx);
                                        for(int ny_p=0;ny_p<Ly;ny_p++){
                                            double qy_p= (2.0*PI_*ny_p)/(1.0*Ly);

                                            if( (nx !=nx_p) || (ny !=ny_p)  ){

                                                for(int alpha=0;alpha<2;alpha++){
                                                    for(int beta=0;beta<2;beta++){

                                                        for(int alpha_p=0;alpha_p<2;alpha_p++){
                                                            for(int beta_p=0;beta_p<2;beta_p++){

                                                                Srw[nw][pos_j][pos_l] += one*( ( Pauli_z[alpha][beta]*Pauli_z[alpha_p][beta_p]  +  Pauli_x[alpha][beta]*Pauli_x[alpha_p][beta_p] + Pauli_y[alpha][beta]*Pauli_y[alpha_p][beta_p])*
                                                                                               exp(iota*(  (xj-xl)*(qx-qx_p) +  (yj-yl)*(qy-qy_p) ))*Lorentzian(w_val + E_NI(qx,qy,thop)  - E_NI(qx_p,qy_p,thop) )*
                                                                                               (1.0/(1.0 + exp(BETA*(E_NI(qx,qy,thop) - mu)) ))*
                                                                                               (1.0/(1.0 + exp(-BETA*(E_NI(qx_p,qy_p,thop) - mu)) ))
                                                                                               );
                                                                //cout<<Pauli_z[alpha][beta]*Pauli_z[alpha_p][beta_p]  +  Pauli_x[alpha][beta]*Pauli_x[alpha_p][beta_p] + Pauli_y[alpha][beta]*Pauli_y[alpha_p][beta_p]<<endl;



                                                            }
                                                        }


                                                    }

                                                }
                                            }

                                        }
                                    }


                                }

                            }



                        }
                    }
                }
            }
        }

    }


    //kvec
    for(int nx=0;nx<Lx;nx++){
        double kx= (2.0*PI_*nx)/(1.0*Lx);
        for(int ny=0;ny<Ly;ny++){
        double ky= (2.0*PI_*ny)/(1.0*Ly);

        for(int nw=0;nw<n_wpoints;nw++){
            double w_val=nw*dw;


            complex<double> val = zero;
        //pos_j
        for(int xj=0;xj<Lx;xj++){
            for(int yj=0;yj<Ly;yj++){
                int pos_j = yj*Lx + xj;

                //pos_l
                for(int xl=0;xl<Lx;xl++){
                    for(int yl=0;yl<Ly;yl++){
                        int pos_l = yl*Lx + xl;


                    val += one*( exp( iota*(-kx*(xj-xl) - ky*(yj-yl)) )*
                                 Srw[nw][pos_j][pos_l]
                                );




                    }}
            }}


        fileout_full<< nx << "    "<<ny <<"    "<< nw <<"    "<<w_val << "   "<<val.real()<<"   "<<val.imag()<<endl;
        }
        fileout_full<<endl;

        }
    }




}


double E_NI(double kx, double ky, double thop){

    double temp;
    temp = -2*thop*(cos(ky));
    return temp;
}


double Lorentzian(double x){

    double PI_ = acos(-1.0);
    double eta=0.2;
    double temp;
    temp = (1.0/PI_)*(0.5*eta)/(  ((x)*(x))  + ((0.5*eta)*(0.5*eta))  );

    return temp;
}



