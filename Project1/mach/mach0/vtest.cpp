#include <iostream>
#include <cmath>
#include <fstream>
#include "vtest.h"
#include "mach0.h"
using namespace std;

int main()
{
	ofstream myfile;
	myfile.open("vtest_mach.txt");
	myfile << "Verification test: \n";
	double error, pi_approx;
	int n;
	int k=24;
	for(int i = 1; i <=k; i++){
		n=pow(2,i);
		pi_approx=mach0(n);
		error=abs(M_PI-pi_approx);
		myfile  <<"n=" << n << ", |pi-pi_approx|="<< error << endl;
	}
	myfile.close();
	return 0;
}
