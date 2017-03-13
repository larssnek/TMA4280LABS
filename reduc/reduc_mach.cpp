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
	
	MPI_Init(&argc,&argv);
	if(argc<3){
		cout << "Error, program needs two arguments" <<endl;
	}
	else{
		int n=atof(argv[1]);
		int boolic=atof(argv[2]);	// 1 for MPI_Allreduce and 0 for recursive doubling sum
		double *vec=new double[2*n];
		double const corr_value =M_PI/4.;
		double *vec1=vec;
		MPI_Status status;
		int rank, size, tag;
		MPI_Comm_size(MPI_COMM_WORLD,&size);
		MPI_Comm_rank(MPI_COMM_WORLD,&rank);
		tag=100;
		int r=(2*n)%size; //rest to balance partitioning of array between processes
		int m=(2*n-r)/size; // least number of elements to each process
		double* recvec=new double[m+1]; //Array used to receive
		int p = m + (rank < r ? 1 : 0);
		double part_sum;
		int ind=p; //indicates how many elements already sent to processes	

		// Make vector elements v_i
		if(rank==0){
			double pi=0., a=1./5., b=1./239.;
			for (int i=1; i<=n; ++i) {
				vec[2*(i-1)]=4.0*pow(-1.0,i-1)/(2.0*i-1)*pow(a,2*i-1);
				vec[2*(i-1)+1]=pow(-1.0,i)/(2.0*i-1)*pow(b,2*i-1); //multiply by -1 so that we can just sum up for each process
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

				vec1=vec+ind;
				MPI_Send(vec1,p,MPI_DOUBLE,i,tag,MPI_COMM_WORLD);
				ind+=p;
			}
			// sum for root
			p=p0;
			recvec=vec;
		}
		else{
			MPI_Recv(recvec, p, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
		}
		for(int i=0; i<p;i++){
			part_sum+=recvec[i];
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
			
			cout<< scientific << "The global sum: " <<4.0*Sn << endl;
			cout <<scientific  << "Error for n=" << n<<" and the number of processes=" <<size << ": " << abs(corr_value-Sn) << endl;
			
			double duration = MPI_Wtime()-starttime;	
			cout << "Run time: " << duration << endl;
		}
		MPI_Finalize();
	}
}
