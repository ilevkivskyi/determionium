#include <stdio.h>
#include "mpi.h"
#include "/bgsys/include/essl.h"
#include <math.h>
#define cset(x,a,b) (RE(x)=a, IM(x)=b)
#define MU_N ((int)(MU*N/(2*PI)))
#define mtime rank

// Complexity eqs. to approx. (T_NUM+Q+S)N^3

#define N 4000  // Number of electrons
#define T_NUM 14 // Number of transparencies

// Precission

#define Q 9 // Number of terms in approximation of exponent by Taylor series
#define S 5 // Number of squarings to make the norm < approx. 0.2

// Physical parameters of the problem

#define MU (0.015*PI) // Voltage bias
#define LA PI // Coupling constant

// Constants

#define PI 3.141592653589793

// Auxiliary functions

double im_mat(int t, int n){
    if(n == 0){
        return LA*t*1.0/N;
    }
    return LA*(((cos(2*PI*n*t/N) - 1) * (cos(2*PI*n/N) - 1) + sin(2*PI*n*t/N) * sin(2*PI*n/N)) 
    / (2 * N - 2 * N * cos(2*PI*n/N)));
}

double re_mat(int t, int n){
    if(n == 0){
        return 0.0;
    }
    return -LA*(((cos(2*PI*n/N) - 1) * sin(2*PI*n*t/N) - (cos(2*PI*n*t/N) - 1) * sin(2*PI*n/N)) 
    / (2 * N - 2 * N * cos(2*PI*n/N)));
}

int main(int argc, char *argv[]){
    
    // Initialize variables
    
    int numtasks, rank, rc, i, j, k, index, info;
    int pivot[N+MU_N];
    for(i=0; i<N+MU_N; i++){
        pivot[i]=0;
    }
    dcmplx *a, *b, *c, *temp, alpha, beta, results[T_NUM+1];
    double T[T_NUM] = {0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.475, 0.5, 0.7, 1.0};
    a = (dcmplx *)malloc(sizeof(dcmplx)*N*N);
    b = (dcmplx *)malloc(sizeof(dcmplx)*N*N);
    c = (dcmplx *)malloc(sizeof(dcmplx)*N*N);
    FILE *out;
    
    // Initialize and test MPI
    
    rc = MPI_Init(&argc, &argv);
    if (rc != MPI_SUCCESS){
        // printf ("Error starting MPI program. Terminating.\n");
        MPI_Abort(MPI_COMM_WORLD, rc);
    }
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // printf ("Number of tasks= %d My rank= %d, \n", numtasks, rank);

    // Initially, fill each matrix with same elemets divided by 2^S
    
    for (i=0; i<N; i++){
       for (j = 0; j<N; j++){
           if (abs(i - j) < N/2){
               cset (a[N*j+i], re_mat(mtime, i-j)/pow(2, S), im_mat(mtime, i-j)/pow(2, S));
           }
           else {
               cset (a[N*j+i], 0.0, 0.0);
           }
       }
    }
    zcopy(N*N, a, 1, b, 1);
    zcopy(N*N, a, 1, c, 1);
    
    // Main loop to calculate the exponent
    
    cset(beta, 1.0, 0.0);
    for(i=Q; i>1; i--){
        cset(alpha, 1.0/i, 0.0);
        zgemm("n", "n", N, N, N, alpha, a, N, b, N, beta, c, N);
        temp = c;
        c = a;
        a = temp;
        zcopy(N*N, b, 1, c, 1);
    }
   
    // Add unit matrix and perform squarings
    
    for(i=0; i<N; i++){
        RE(a[N*i+i]) += 1.0;
    }
    
    cset(alpha, 1.0, 0.0);
    cset(beta, 0.0, 0.0);
    for(i=0; i<S; i++){
        zcopy(N*N, a, 1, b, 1);
        zcopy(N*N, a, 1, c, 1);
        zgemm("n", "n", N, N, N, alpha, b, N, c, N, beta, a, N);
    }
    
    for(i=0; i<N; i++){
        RE(a[N*i+i]) += -1.0;
    }
    
    // Now 'a' contains the exponent of the initial matrix minux one; perform some checks

    //if (rank==0 || rank==5 || rank==27 || rank==55 || rank==100 || rank==127){
        //printf("My time= %d : ", mtime);
        //printf("My intermediate result = %e + i* %e\n", RE(a[10*N+7]), IM(a[10*N+7]));
    //}

    // Calculating the determinant for each transparency

    free(b);
    free(c);
    c = (dcmplx *)malloc(sizeof(dcmplx)*(N+MU_N)*(N+MU_N));
    
        
    for(k=0; k<T_NUM; k++){
        
        // Create the full matrix of four blocks
        
        for(i=0; i<N/2; i++){
            for(j=0; j<N/2; j++){
                index = N*j+i;
                cset(c[(N+MU_N)*j+i], RE(a[index]) * (1.0-T[k]), IM(a[index]) * (1.0-T[k]));
            }
        }
        for(i=N/2; i<N+MU_N; i++){
            for(j=N/2; j<N+MU_N; j++){
                index = N*(j-N/2)+i-N/2;
                cset(c[(N+MU_N)*j+i], RE(a[index]) * T[k], IM(a[index]) * T[k]);
            }
        }
        for(i=0; i<N/2; i++){
            for(j=N/2; j<N+MU_N; j++){
                index = N*(j-N/2)+i;
                cset(c[(N+MU_N)*j+i], RE(a[index])*sqrt(T[k]*(1.0-T[k])), IM(a[index])*sqrt(T[k]*(1.0-T[k])));
            }
        }
        for(i=N/2; i<N+MU_N; i++){
            for(j=0; j<N/2; j++){
                index = N*j+i-N/2;
                cset(c[(N+MU_N)*j+i], RE(a[index])*sqrt(T[k]*(1.0-T[k])), IM(a[index])*sqrt(T[k]*(1.0-T[k])));
            }
        }
        for(i=0; i<N+MU_N; i++){
        RE(c[(N+MU_N)*i+i]) += 1.0;
        }
        // Perform the LU decomposition and sum logs of diagonal elements to get the log(det)
        
        info = 0;
        zgetrf(N+MU_N, N+MU_N, c, N+MU_N, pivot, &info);
        cset(results[k], 0.0, 0.0);
        for(i=0; i<N+MU_N; i++){
            RE(results[k]) += log(IM(c[i*(N+MU_N)+i])*IM(c[i*(N+MU_N)+i]) + RE(c[i*(N+MU_N)+i])*RE(c[i*(N+MU_N)+i]))/2; 
            IM(results[k]) += atan2(IM(c[i*(N+MU_N)+i]), RE(c[i*(N+MU_N)+i]));
            if(pivot[i] != i){
                IM(results[k]) -= PI;
            }
        }
        IM(results[k]) += -LA*mtime/2;
    }
    cset(results[T_NUM], mtime, 0.0);
    free(c);
    free(a);
    
    // Here we wait and gather results for writing to disk
    
    //rc = MPI_Barrier(MPI_COMM_WORLD);
    //if (rank == 0){
        //c = (cmplx *)malloc(sizeof(cmplx)*numtasks*(T_NUM+1));
    //}
    //else{
        //c = (cmplx *)malloc(sizeof(cmplx)*1);
    //}
    //rc = MPI_Gather(results, T_NUM+1, MPI_COMPLEX,
                       //c, T_NUM+1, MPI_COMPLEX,
                       //0, MPI_COMM_WORLD);

    rc = MPI_Barrier(MPI_COMM_WORLD);
    for(i=0; i<numtasks; i++){
        if(rank == i){
            out = fopen("/bgscratch/levkivsk/out_4000_pi_0015_d_1.txt", "a");
            for(k=0; k<T_NUM+1; k++){
                fprintf(out, "%.10e+%.10ej, ", RE(results[k]), IM(results[k]));
            }
            fprintf(out, "\n");
            rc = fclose(out);
        }
        
        rc = MPI_Barrier(MPI_COMM_WORLD);
    }
    
    // Free the rest and finalize
    
    MPI_Finalize();
    return 0;
}
