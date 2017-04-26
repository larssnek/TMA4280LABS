#include <iostream>
#include <cmath>
#include "mpi.h"
#include "omp.h"
#include <ctime>
#include "zeta2.h"
#include <fstream>

double starttime=MPI_Wtime();

using namespace std;

void zeta2(int s, int n, double *error, double *time)
{
	double *vec=new double[n];
	clock_t start;
	double const corr_value = (M_PI*M_PI)/6.;
	double *vec1=vec;
	MPI_Status status;
	int rank, size, tag;
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	tag=100;

	int r=(n)%size; //rest to balance partitioning of array between processes
	int m=(n-r)/size; // least number of elements to each process
	double* recvec=new double[m+1]; //Array used to receive
	int p = m + (rank < r ? 1 : 0);
	double summer;
	int ind=p; //indicates how many elements already sent to processes	
	// to test wich threads
	int nthreads;
	int t_id;
	
	
	if(rank==0){
		// Make elements v_i
		double pi=0., a=1./5., b=1./239.;
		#pragma omp parallel for schedule(static)
		for (int i=1; i<=n; ++i) {
			 vec[i-1]=1.0/std::pow(i,s);			
		}
                // set the size for root
                int const p0 = p;
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
/////////////////////////////////////////////////////////////////////////////////
			t_id = omp_get_thread_num();
			//cout << "loop 1 " << "and thread #" << t_id << endl;
			if(t_id==0){
			nthreads=omp_get_num_threads();
			//cout << "Number of threads: " << nthreads << endl;
			}
		}
                // sum for root
		p=p0;
		recvec=vec;
	}
	else{
		MPI_Recv(recvec, p, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
	}
	summer=0.0;
	#pragma omp parallel for schedule(static) reduction(+:summer) private(t_id, nthreads)
	for(int i=0; i<p;i++){
		summer+=recvec[i];
	/////////////////////////////////////////////////////////////////////////////////		
		t_id = omp_get_thread_num();
		//cout << "loop 2 " << "and thread #" << t_id << endl;
		if(t_id==0){
			nthreads=omp_get_num_threads();
		//	cout << "Number of threads: " << nthreads << endl;
		}
	}
	double Sn;
	MPI_Reduce(&summer, &Sn,1,MPI_DOUBLE, MPI_SUM,0,MPI_COMM_WORLD);
	
	if(rank==0){
		// Prints global sum, error, number of processes and running time
		*error=abs(corr_value-Sn);
		cout<< scientific << "The global sum: " <<sqrt(6.0*Sn) << endl;
		cout <<scientific  << "Error for n=" << n<<" and the number of processes=" <<size << ": " << *error << endl;
		*time=MPI_Wtime()-starttime; 	
		cout << "Run time: " << *time << endl;
	}



}
