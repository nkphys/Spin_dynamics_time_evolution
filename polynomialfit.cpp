/******************************************************
*************Chi-square fitting**************
Polynomial Fitting
******************************************************/
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <stdlib.h>
#include "tensor_type.h"
using namespace std;


/*******
 Function that performs Gauss-Elimination and returns the Upper triangular matrix and solution of equations:
There are two options to do this in C.
1. Pass the augmented matrix (a) as the parameter, and calculate and store the upperTriangular(Gauss-Eliminated Matrix) in it.
2. Use malloc and make the function of pointer type and return the pointer.
This program uses the first option.
********/
void gaussEliminationLS(int m, int n, Mat_2_doub & a, Mat_1_doub & x){
    int i,j,k;
    for(i=0;i<m-1;i++){
        //Partial Pivoting
        for(k=i+1;k<m;k++){
            //If diagonal element(absolute vallue) is smaller than any of the terms below it
            if(fabs(a[i][i])<fabs(a[k][i])){
                //Swap the rows
                for(j=0;j<n;j++){                
                    double temp;
                    temp=a[i][j];
                    a[i][j]=a[k][j];
                    a[k][j]=temp;
                }
            }
        }
        //Begin Gauss Elimination
        for(k=i+1;k<m;k++){
            double  term=a[k][i]/ a[i][i];
            for(j=0;j<n;j++){
                a[k][j]=a[k][j]-term*a[i][j];
            }
        }
         
    }
    //Begin Back-substitution
    for(i=m-1;i>=0;i--){
        x[i]=a[i][n-1];
        for(j=i+1;j<n-1;j++){
            x[i]=x[i]-a[i][j]*x[j];
        }
        x[i]=x[i]/a[i][i];
    }
             
}
/*******
Function that prints the elements of a matrix row-wise
Parameters: rows(m),columns(n),matrix[m][n] 
*******/
void printMatrix(int m, int n, Mat_2_doub matrix){
    int i,j;
    for(i=0;i<m;i++){
        for(j=0;j<n;j++){
            printf("%lf\t",matrix[i][j]);
        }
        printf("\n");
    } 
}
int main(){

    //no. of data-points
    int N;  
    //degree of polynomial
    int n;  
    //printf("Enter the no. of data-points:\n");
    //cin>>N;
    //arrays to store the c and y-axis data-points
    vector<double> x, y;
    
    /*
    printf("Enter the x-axis values:\n");
    int i,j;
    for(i=0;i<N;i++){
    cin>>x[i];
    }
    printf("Enter the y-axis values:\n");
    for(i=0;i<N;i++){
    cin>>y[i];
    }


    printf("Enter the degree of polynomial to be used:\n");
    cin>>n;
   */
	
    n=50;
    string line;
    string Scheduler_File="dmrg_sched_nofilter.txt";
    double temp_t, temp_h, temp_J;
    Mat_1_doub Time_bare, Gamma_bare, Js_bare; 
    ifstream scheduler_stream(Scheduler_File.c_str());
    while(getline(scheduler_stream,line)){
    stringstream line_ss(line);
    line_ss>>temp_t>>temp_h>>temp_J;
    Time_bare.push_back(temp_t);
    Gamma_bare.push_back(temp_h);
    Js_bare.push_back(temp_J);
    }
    x=Time_bare;
    y=Gamma_bare;

    N=x.size();



    // an array of size 2*n+1 for storing N, Sig xi, Sig xi^2, ...., etc. which are the independent components of the normal matrix
    vector<double> X;
    X.resize(2*n+1);  
    for(int i=0;i<=2*n;i++){
        X[i]=0;
        for(int j=0;j<N;j++){
            X[i]=X[i]+pow(x[j],i);
        }
    }
    //the normal augmented matrix
    Mat_2_doub B;
    B.resize(n+1);
    for(int i=0;i<B.size();i++){
    B[i].resize(n+2);
    }

    // rhs
    Mat_1_doub Y;
    Y.resize(n+1);
    for(int i=0;i<=n;i++){
        Y[i]=0;
        for(int j=0;j<N;j++){
            Y[i]=Y[i]+pow(x[j],i)*y[j];
        }
    }
    for(int i=0;i<=n;i++){
        for(int j=0;j<=n;j++){
            B[i][j]=X[i+j];
        }
    }
    for(int i=0;i<=n;i++){
        B[i][n+1]=Y[i];
    }
    Mat_1_doub A;
    A.resize(n+1);
    printf("The polynomial fit is given by the equation:\n");
    printMatrix(n+1,n+2,B);
    gaussEliminationLS(n+1,n+2,B,A);
    for(int i=0;i<=n;i++){
if(i==0){
        printf("%lf*x**%d",A[i],i);
}
else{
if(A[i]>0){
	printf("+%lf*x**%d",A[i],i);
}
else{
	printf("%lf*x**%d",A[i],i);
}
}
    }
    cout<<endl;
return 0;
}
