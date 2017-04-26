#include <iostream>
#include <cmath>
#include "zeta0.h"

using std::pow;
using std::sqrt;
clock_t begin = clock();
double zeta0(int n, int s, double *time){
	double output = 0.0;
	for (int i=1; i<=n;i++){
	output+= 1.0/pow(i,s);
	}
	double pi;
	pi=sqrt(output*6);
	clock_t end = clock();
	*time=double(end - begin) / CLOCKS_PER_SEC;
	return pi;

}

