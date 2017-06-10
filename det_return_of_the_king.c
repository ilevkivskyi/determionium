#include <stdio.h>
#include "mpi.h"
#include <cblas.h>
#include <lapacke.h>
#include <math.h>

#define REQUEST 5
#define DATA 8
#define RESULT 40

#define N 100000
#define SIZE (t+(int)floor(DEC*t_RC)+1)
#define DEC 8

#define PI 3.141592653589793

typedef struct{ double re; double im; } dcmplx;
typedef struct{ double t; double t_RC; double la; double mu; double T; } task_t;

double decay(int tp, int t, double t_RC){ 
	if (tp<t) return 1-exp(-tp/t_RC);
	return exp(-(tp-t)/t_RC) - exp(-tp/t_RC);
}

double re_F(int tp, int tpp, int t, double mu, double T){

	int d;

	d = tp-tpp;
	if (tp == tpp) return 0.5+mu*T/(2*PI);
	return ((1+T*(cos(mu*d) - 1)-cos(PI*d))*(cos(2*PI*d/N)-1) + T*sin(mu*d)*sin(2*PI*d/N))
	        /(2*N-2*N*cos(2*PI*d/N));
}

double im_F(int tp, int tpp, int t, double mu, double T){

	int d;

	d = tp-tpp;
	if (tp == tpp) return 0;
	return (T*sin(mu*d)*(cos(2*PI*d/N)-1) - sin(2*PI*d/N)*(1+T*(cos(mu*d) - 1)-cos(PI*d)))
	        /(2*N-2*N*cos(2*PI*d/N));
}

dcmplx det(int t, double t_RC, double la, double mu, double T){

	int i, j, info=0, pivot[SIZE];
	dcmplx *mat, elem, result;
	
	for(i=0; i<SIZE; i++){
        pivot[i]=0;
    }
	mat = (dcmplx *)malloc(sizeof(dcmplx)*SIZE*SIZE);
    for (i=0; i<SIZE; i++){
		for (j = 0; j<SIZE; j++){
			mat[i+SIZE*j].re = -re_F(i,j,t,mu,T)+re_F(i,j,t,mu,T)*cos(la*decay(i,t,t_RC))
			                   -im_F(i,j,t,mu,T)*sin(la*decay(i,t,t_RC));
			mat[i+SIZE*j].im = -im_F(i,j,t,mu,T)+re_F(i,j,t,mu,T)*sin(la*decay(i,t,t_RC))
			                   +im_F(i,j,t,mu,T)*cos(la*decay(i,t,t_RC));
			if (i==j) mat[i+SIZE*j].re += 1;
		}
	}
	info = LAPACKE_zgetrf(LAPACK_COL_MAJOR, SIZE, SIZE, mat, SIZE, pivot);
	result.re = 0;
	result.im = 0;
	for (i=0;i<SIZE;i++){
		elem = mat[i+SIZE*i];
		result.re += log(elem.im*elem.im + elem.re*elem.re)/2; 
        result.im += atan2(elem.im, elem.re);
        if(pivot[i] != i+1) result.im += PI;
	} 
	free(mat); 
	result.im -= la*t/2;		
	return result;
}

int main(int argc, char *argv[]){
    int numtasks, rank, tag, dest, src, rc, tasks_left=1, working;
    dcmplx res;
    task_t task;
    double val;
    FILE *in, *out;
    MPI_Status status;
     
    // Initialize and test MPI
    
    rc = MPI_Init(&argc, &argv);
    if (rc != MPI_SUCCESS){
        printf ("Error starting MPI program. Terminating.\n");
        MPI_Abort(MPI_COMM_WORLD, rc);
    }
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    printf ("Process number %d says hi!\n", rank);
    
    if(rank == 0){ // Server node
		in = fopen(argv[1],"r");
		out = fopen(argv[2], "w");
		working = numtasks-1;
		while(working){
			MPI_Recv(&task, 5, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			tag = status.MPI_TAG;
			src = status.MPI_SOURCE;
			if(tag==REQUEST){
				dest = src;
				rc = fscanf(in, "%le %le %le %le %le", &(task.t), &(task.t_RC), &(task.la), &(task.mu), &(task.T));
				if(rc==EOF){
					task.t = -1;
					working -= 1;
				}
				MPI_Send(&task, 5, MPI_DOUBLE, dest, DATA, MPI_COMM_WORLD);
				printf("Processing time = %d\n", (int)(task.t));
			} else {
				MPI_Recv(&val, 1, MPI_DOUBLE, src, RESULT, MPI_COMM_WORLD, &status);
				res.re = val;
				MPI_Recv(&val, 1, MPI_DOUBLE, src, RESULT, MPI_COMM_WORLD, &status);
				res.im = val;
				fprintf(out, "%.10e %.10e %.10e %.10e %.10e %.10e %.10e\n", 
				               task.t, task.t_RC, task.la, task.mu, task.T, res.re, res.im);
			}
		}
		printf("Process number 0 says bye!\n");
	}
	else {   // working nodes
		while(tasks_left){
			val=0;
			MPI_Send(&val, 1, MPI_DOUBLE, 0, REQUEST, MPI_COMM_WORLD);   // send request for new task
			MPI_Recv(&task, 5, MPI_DOUBLE, 0, DATA, MPI_COMM_WORLD, &status); 
			if(task.t >= 0){    // if received a task: process it
				res = det((int)(task.t), task.t_RC, task.la, task.mu, task.T);
				MPI_Send(&task, 5, MPI_DOUBLE, 0, RESULT, MPI_COMM_WORLD);
				val = res.re;
				MPI_Send(&val, 1, MPI_DOUBLE, 0, RESULT, MPI_COMM_WORLD);
				val = res.im;
				MPI_Send(&val, 1, MPI_DOUBLE, 0, RESULT, MPI_COMM_WORLD);
			} else tasks_left = 0;
		}
		printf("Process number %d says bye!\n", rank); 
	}
    
    // Finalize MPI
    rc = MPI_Finalize();
	return rc;
}
