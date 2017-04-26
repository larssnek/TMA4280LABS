#include "mach1.h"
#include "mpi.h"
#include <cmath>
#include <iostream>
#include <fstream>
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
	double time;
	double error;
	///////////////////////////////////////////////////////////////////////////////////
	/* Used for printing errors and time for different n */	
	
	
	int const max_n=16;
	ofstream data;
	std::stringstream ss;
        ss << "Data_" << size << ".txt";
	if(rank==0){
		//ofstream data;
		data.open(ss.str() , ios_base::app);
	}
	 //I get a problem when number of processes > 2*n, so this must be accounted for
	int start=3;
	for(int i =start; i<=max_n; i++){
		n=pow(2,i);
		mach1(n, &error, &time);
		if(rank==0){
			data << n << " " << error << " " << time << " " << endl;
		}
	}


	if(rank==0){	
		data.close();
	}
	

///////////////////////////////////////////////////////////////////////////
	/* Used when not needing printing */
	/*
	n=pow(2,10);	
	mach1(n, &error, &time);
	*/	
	MPI_Finalize();
	
	return 0;
}
