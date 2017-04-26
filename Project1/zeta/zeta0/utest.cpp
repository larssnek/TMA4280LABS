#include <iostream>
#include <cmath>
#include "zeta0.h"
#include <stdio.h>
#include <stdlib.h>
int main(){
	double corr_value =2.857738033;
	int n = 3, s = 2;

	double Sz=zeta0(n,s);

	double rel_err=1-corr_value/Sz;
	printf("Relative error for n=3 rel_err = %.16e\n\n", rel_err);

	return 0;
}



