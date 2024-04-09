#include <iostream>
#include <math.h>  //fabs(double x) =|x|
#include <algorithm>
#include <stdlib.h>  //for div(q,n).rem(quot),rand
#include <time.h>
#include <fstream>
#include <sstream>
#include <limits>
#include <iomanip>
#include <stdio.h>
#include <random>
#define PI cos(-1.0)

using namespace std;

int main(){

int Lx=8;
int Ly=8;

int L=Lx*Ly;
double J_val=-1.0;

string file_out_str = "Jzz_2d.txt";
ofstream file_out(file_out_str.c_str());

vector< vector<double>> Jzz_conn;
Jzz_conn.resize(L);

for (int i=0;i<L;i++){
Jzz_conn[i].resize(L);
for (int j=0;j<L;j++){
Jzz_conn[i][j]=0.0;
}
}

int i_neigh, j_neigh, neigh;
int site;
for (int ix=0;ix<Lx;ix++){
for (int iy=0;iy<Ly;iy++){
site=ix + iy*Lx;

//+x
i_neigh = ( ix + 1 + Lx )%Lx;
j_neigh = iy;
neigh = i_neigh + j_neigh*Lx;
Jzz_conn[site][neigh]=J_val;
Jzz_conn[neigh][site]=J_val;


//+y
j_neigh = ( iy + 1 + Ly )%Ly;
i_neigh = ix;
neigh = i_neigh + j_neigh*Lx;
Jzz_conn[site][neigh]=J_val;
Jzz_conn[neigh][site]=J_val;


}
}

for (int i=0;i<L;i++){
for (int j=0;j<L;j++){
file_out<<Jzz_conn[i][j]<<" ";
}
file_out<<endl;
}

return 0;
}
