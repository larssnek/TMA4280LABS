#include <iostream>
#include <cmath>
#include "zeta0.h"
#include <stdlib.h>
#include <stdio.h>
using namespace std;

int main(int argc, char** argv) {
	
	if(argc<2){
	cout << "error, # of input arguments should be one" << endl;
	}
	else{
		int n=atof(argv[1]);
		int s=2;
		double pi;
		double time;
		pi=zeta0(n,s, &time);
		cout <<scientific<<"For n=" << n <<" we get pi="<< pi << endl;
		double error=abs(M_PI-pi);
		printf("Error for n=%d is %.16e\n",n,error );
		printf("Run time for n=%d is %.3e\n",n,time );
	}
}
