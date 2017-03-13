#include <iostream>
#include "mach0.h"
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
using namespace std;

int main(int argc, char** argv) {
	
	if(argc<2){
	cout << "error, # of input arguments should be one" << endl;
	}
	else{
		int n=atof(argv[1]);
		double pi;
		pi=mach0(n);
		cout <<scientific<<"For n=" << n <<" we get pi="<< pi << endl;
		double error=abs(M_PI-pi);
		printf("Error for n=%d is %.16e\n\n",n,error );
	}
	return 0;
}
