#include "mach2.h"
#include "mpi.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <omp.h>
#include <sstream>
using namespace std;

int main(int argc,char **argv)
{
	int rank, size;
	MPI_Init(&argc, &argv);	
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	// Print format will be data << "n: " << "error: " << "time: " << endl;
	int n;
	int const max_n=16;
	double time;
	double error;
	ofstream data;
        std::stringstream ss;
        ss << "Data_" << size << ".txt";
	if(rank==0){
		data.open(ss.str() , ios_base::app);
	}
	int start=3;
	for(int i =start; i<=max_n; i++){
		n=pow(2,i);
		mach2(n, &error, &time);
		if(rank==0){
			data << n << " " << error << " " << time << " " << endl;
		}
	}

	if(rank==0){	
		data.close();
	}
		
	MPI_Finalize();
	
	return 0;
}
