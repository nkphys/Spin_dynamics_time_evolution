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

double Lorentzian(double x, double brd)
{
    double temp;

    temp = (1.0 / PI) * ((brd / 2.0) / ((x * x) + ((brd * brd) / 4.0)));

    return temp;
}


int main(){


mt19937_64 Generator_;
Generator_.seed(11);
normal_distribution<double> GaussianDistribution2(2.0,1.0);
normal_distribution<double> GaussianDistribution;

GaussianDistribution=GaussianDistribution2;

int Length=1000000;
vector<double> gaussarray;
gaussarray.resize(Length);
for(int i=0;i<Length;i++){
gaussarray[i]=GaussianDistribution(Generator_);
}

double x=-10.0;
double px=0.0;
while(x<10){
px=0.0;
for(int i=0;i<Length;i++){
px += Lorentzian((x-gaussarray[i]), 0.1);
}
cout<<x<<"  "<<px/Length<<endl;
x=0.1+x;
}



return 0;
}
