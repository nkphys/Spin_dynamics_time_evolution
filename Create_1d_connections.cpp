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

int L=120;
double J_val=-1.0;
string file_out_str = "Jzz_1d.txt";
ofstream file_out(file_out_str.c_str());

vector< vector<double>> Jzz_conn;
Jzz_conn.resize(L);

for (int i=0;i<L;i++){
Jzz_conn[i].resize(L);
for (int j=0;j<L;j++){
Jzz_conn[i][j]=0.0;
}
}

int neigh;
for (int i=0;i<L;i++){
neigh = ( i + 1 + L )%L;
Jzz_conn[i][neigh]=J_val;
Jzz_conn[neigh][i]=J_val;
}

for (int i=0;i<L;i++){
for (int j=0;j<L;j++){
file_out<<Jzz_conn[i][j]<<" ";
}
file_out<<endl;
}

return 0;
}
