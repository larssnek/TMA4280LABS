#include <iostream>
#include <fstream>
#include <cmath>
#include "mpi.h"
#include <ctime>
#include <stdlib.h>
#include <ctime>
using namespace std;

double starttime=MPI_Wtime();

int main(int argc, char** argv)
{
	ofstream data;
	MPI_Init(&argc,&argv);
	if(argc<3){
		cout << "Error, program needs two arguments" <<endl;
	}
	else{
		int s=2;
		int n=atof(argv[1]);
		int boolic=atof(argv[2]);	// 1 for MPI_Allreduce and 0 for recursive doubling sum
		double const corr_value = (M_PI*M_PI)/6.;
		MPI_Status status;
		int rank, size, tag;
		MPI_Comm_size(MPI_COMM_WORLD,&size);
		MPI_Comm_rank(MPI_COMM_WORLD,&rank);
		tag=100;
		double *vec=new double[n];
		int r=(n)%size; //rest to balance partitioning of array between processes
		int m=(n-r)/size; // least number of elements to each process
		int p = m + (rank < r ? 1 : 0);
		double part_sum;
		int ind=p; //indicates how many elements already sent to processes	

		// Make vector elements v_i
		if(rank==0){
			double pi=0., a=1./5., b=1./239.;
			for (int i=1; i<=n; ++i) {
				vec[i-1]=1.0/std::pow(i,s);
			}
			// set the size for root
			int const p0 = p;
			// update offset for rank 1
			for(int i=1;i< size;i++){
				if(r>0){
					r-=1;
					p=m+1;
				}
				else{
					p=m;
				}
				MPI_Send(vec+ind,p,MPI_DOUBLE,i,tag,MPI_COMM_WORLD);
				ind+=p;
			}
			// sum for root
			p=p0;
		
		}
		else{
			MPI_Recv(vec, p, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
		}
		for(int i=0; i<p;i++){
			part_sum+=vec[i];
		}
		double Sn;
		if(boolic){
		
			MPI_Allreduce(&part_sum, &Sn,1,MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD);
		}
		else{
			 Sn=part_sum;
			double sigma_q;
			int p(rank), P(size);
			int q;
			for (int d = 0; d <std::log2(P); d++){
				q=(p^(1<<d));
				MPI_Send(&Sn,1,MPI_DOUBLE,q ,tag,MPI_COMM_WORLD);
				MPI_Recv(&sigma_q,1, MPI_DOUBLE,q ,tag,MPI_COMM_WORLD, &status);
				Sn +=sigma_q;
			}
		}
		if(rank==0){
			double pi=std::sqrt(6.0*Sn);
			printf("The global sum is pi=%.16e\n",pi );
			double time=MPI_Wtime()-starttime;
			double error=abs(corr_value -Sn);
			printf("Error for n=%d and p=%d is %.16e\n",n,size,error );			
			std::cout << "Run time: " << time << std::endl;

		}
		delete [] vec;
		MPI_Finalize();
	}

}
