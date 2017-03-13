#include <iostream>
#include <cmath>
#include "mach0.h"
#include <stdio.h>
#include <stdlib.h>

using namespace std;
int main()
	{
	int ret=0;
	double expected_value=3.141621029;
	int n=3;
	double pi_approx=mach0(n);
	//cout<<scientific << "expected value of pi with Machin formula" << endl << "using calculator: " << expected_value << endl << "vs computed value: " << pi_approx << endl;
	double rel_err=1-expected_value/pi_approx;
	 printf("Relative error for n=3 rel_err = %.16e\n\n", rel_err);
  	return 0;
}
