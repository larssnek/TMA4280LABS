
/**
 * C program to solve the two-dimensional Poisson equation on
 * a unit square using one-dimensional eigenvalue decompositions
 * and fast sine transforms.
 *
 * Einar M. Rønquist
 * NTNU, October 2000
 * Revised, October 2001
 * Revised by Eivind Fonn, February 2015
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "omp.h"
#include <string.h>

#define PI 3.14159265358979323846
#define true 1
#define false 0

typedef double real;
typedef int bool;

// Function prototypes
real *mk_1D_array(size_t n, bool zero);
real **mk_2D_array(size_t n1, size_t n2, bool zero);
void transpose(real **bt, real **b, real *bvalues, real *btvalues, int *bsize, int *displ,int *nrows, int size, int rank, int m);
real rhs(real x, real y);
void fst_(real *v, int *n, real *w, int *nn);
void fstinv_(real *v, int *n, real *w, int *nn);

int main(int argc, char **argv)
{
	if (argc < 2) {
		printf("Usage:\n");
		printf("  poisson n\n\n");
		printf("Arguments:\n");
		printf("  n: the problem size (must be a power of 2)\n");
	}
	
	// The number of grid points in each direction is n+1
	// The number of degrees of freedom in each direction is n-1
	int n = atoi(argv[1]);

	int m = n - 1;
	int nn = 4 * n;
	real h = 1.0 / n;
	//Used for MPI partitioning
	int rank, size;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Status status;

	//start time
	double start;
	if (rank == 0){
	        start = omp_get_wtime();
	}
	int index=0;
	int tag=100;

	// Distribute rows evenly between the MPI-processes
	int rows=m/size;
	int *nrows = calloc(size, sizeof(int));
	
	
	for(size_t i=0; i<size; i++){
		nrows[i]=rows;
	}
	//Distribute restrows
	int rest=m%size;
	
	for(size_t i =1; i<=rest; i++){
		nrows[size-i]++;
	}

	int *bsize = calloc(size, sizeof(int));
	#pragma omp parallel for schedule(static)	
	for(size_t i=0; i<size; i++){
		bsize[i]=nrows[rank]*nrows[i];
	}
	
	//Make a displacement vector to keep track for each rank
	int *displacement=calloc(size+1, sizeof(int));
	int *local_displacement=calloc(size+1,sizeof(int));
	displacement[0]=0;	
	local_displacement[0]=0;
	for(size_t i = 1; i<=size; i++){
		displacement[i]=displacement[i-1]+bsize[i-1];
		local_displacement[i]=local_displacement[i-1]+nrows[i-1];
	}
	//int threads=0;	
	// Grid points
	real *grid = mk_1D_array(n+1, false);
	#pragma omp parallel for schedule(static)
	for (size_t i = 0; i < n+1; i++) {
		grid[i] = i * h;
	//	threads=omp_get_num_threads();
	}
	
	

	 // The diagonal of the eigenvalue matrix of T
	real *diag = mk_1D_array(m, false);
	#pragma omp parallel for schedule(static)
	for (size_t i = 0; i < m; i++) {
		diag[i] = 2.0 * (1.0 - cos((i+1) * PI / n)); 
	}


	// Initialize the right hand side data
	real **b = mk_2D_array(nrows[rank], m, false);
	real **bt = mk_2D_array(nrows[rank], m, false);
	real **z = mk_2D_array(nrows[rank],nn, true);             
	real * bvalues = mk_1D_array(nrows[rank]*m,false);
	real * btvalues = mk_1D_array(nrows[rank]*m,false);

	#pragma omp parallel for schedule(static)
	for (size_t i = 0; i < nrows[rank]; i++) {
		for (size_t j = 0; j < m; j++) {
			b[i][j]=h*h*rhs(grid[local_displacement[rank]+i+1], grid[j+1]); 
			
		}
	}
	

	// Calculate Btilde^T = S^-1 * (S * B)^T      
	#pragma omp parallel for schedule(static)
	for (size_t i = 0; i < nrows[rank]; i++) {
		fst_(b[i], &n, z[i], &nn);
	}
	
	
	transpose(bt, b,bvalues,btvalues, bsize,displacement,nrows,size,rank,m);
	
	
	#pragma omp parallel for schedule(static)	
	for (size_t i = 0; i < nrows[rank]; i++) {
		fstinv_(bt[i], &n, z[i], &nn);
	}

	// Solve Lambda * Xtilde = Btilde
	#pragma omp parallel for schedule(static)
	for (size_t i = 0; i < nrows[rank]; i++) {
		for (size_t j = 0; j < m; j++) {
			bt[i][j] = bt[i][j] / (diag[j] + diag[local_displacement[rank]+i]); 
		}
	}

	// Calculate X = S^-1 * (S * Xtilde^T)
	#pragma omp parallel for schedule(static)
	for (size_t i = 0; i < nrows[rank]; i++) {
		fst_(bt[i], &n, z[i], &nn);
	}
	
	transpose(b, bt,bvalues,btvalues, bsize,displacement,nrows,size,rank,m);
	
	#pragma omp parallel for schedule(static)
	for (size_t i = 0; i < nrows[rank]; i++) {
		fstinv_(b[i], &n, z[i], &nn);
	}

	// Calculate maximal value of solution
	double u_max = 0.0;
	#pragma omp parallel for schedule(static) reduction(max:u_max)
	for (size_t i = 0; i < nrows[rank]; i++) {
		for (size_t j = 0; j < m; j++) {
			u_max = u_max > b[i][j] ? u_max : b[i][j];
		}
	}

	// All to find reduce to find maximum global sum
	double global_u_max=0.0;
	MPI_Reduce(&u_max, &global_u_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	
	if(rank==0){
		double duration = (omp_get_wtime() - start);
            	printf("run time: %e\n", duration);	
		printf("global_u_max = %e\n", global_u_max);
		FILE *f = fopen("umax.txt", "a");
		if (f == NULL)
		{
			printf("Error opening file!\n");
    			exit(1);
		}

		/* print umax to file */
		fprintf(f, "%d %f %f \n",n, global_u_max, duration);
		fclose(f);

	}
	
/*	
	// Free memory
	free(btvalues);
	free(bvalues);
	for(int i =0; i<nrows[rank]; i++){
		free(z[i]);
		free(bt[i]);
		free(b[i]);
	}
	free(z);
	free(bt);
	free(b);
	free(diag);
	free(grid);
	free(local_displacement);
	free(displacement);
	free(bsize);
	free(nrows);
*/
	MPI_Finalize();	
	return 0;
}

real rhs(real x, real y) {
    //return 2 * (y - y*y + x - x*x);
    return 5*PI*PI*sin(PI*x)*sin(2*PI*y);
}




void transpose(real **bt, real **b, real *bvalues, real *btvalues, int *bsize, int *displ,int *nrows, int size, int rank,int m){

	int np = nrows[rank];

  // Copy sub-blocks with packing 
	real *Bpck = btvalues;
	for (int p = 0, off_rp = 0; p < size; off_rp+=nrows[p], ++p){
		for (int i = 0; i < np; ++i, Bpck+=nrows[p]){
			memcpy(Bpck, b[i] + off_rp, nrows[p]*sizeof(real)); 
		}
	}


	Bpck = btvalues;

  // Exchange blocks 
	MPI_Alltoallv(btvalues, bsize, displ, MPI_DOUBLE, bvalues, bsize, displ, MPI_DOUBLE, MPI_COMM_WORLD);
  
	
  // Transpose blocks
	int count=0; 
	Bpck = bvalues;
	for (int p = 0, off_rp = 0; p < size; off_rp+=nrows[p], ++p){
		for (int j = 0; j < nrows[p]; ++j){
			for (int i = 0; i < np; ++i){
				bt[i][off_rp + j] = Bpck[count];
				count++;
			}
		}
	}

}





real *mk_1D_array(size_t n, bool zero)
{
    if (zero) {
        return (real *)calloc(n, sizeof(real));
    }
    return (real *)malloc(n * sizeof(real));
}

real **mk_2D_array(size_t n1, size_t n2, bool zero)
{
    real **ret = (real **)malloc(n1 * sizeof(real *));

    if (zero) {
        ret[0] = (real *)calloc(n1 * n2, sizeof(real));
    }
    else {
        ret[0] = (real *)malloc(n1 * n2 * sizeof(real));
    }

    for (size_t i = 1; i < n1; i++) {
        ret[i] = ret[i-1] + n2;
    }
    return ret;
}
