/*
C library for equilibrium chemistry calculations: linalg module

@author: Nick Gibbons
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "linalg.h"

int solve_matrix(double* A, double* B, double *X, int N){
    /*
    Bare bones Gaussian matrix solver for A*X = B 
       A : pointer to NxN matrix of doubles
       B : pointer to N vector of doubles
       X : pointer to N solution vector of doubles
       N : problem size (int)
    */
    int i,j,k,p;
    double pval,a;
    k=0; // working diagonal for row reduction

    while (k<N) {

    // Find the partial-pivot index
    pval = -1.0;
    p = k;
    for (i=k; i<N; i++) {
        a = fabs(A[i*N+k]);
        if (a>pval){
           pval=a;
           p=i;
           }
    }

    // Row swap
    for (j=0; j<N; j++) {
        a = A[k*N+j];
        A[k*N+j] = A[p*N+j];
        A[p*N+j] = a;
    }
    a = B[k];
    B[k] = B[p];
    B[p] = a;

    if (fabs(A[k*N+k])<1e-16) return 1; // Check for singular matrix error
    // Multiply Pivot Row
    a = 1.0/A[k*N+k];
    for (j=0; j<N; j++) {
        A[k*N+j]*=a;
    }
    B[k]*=a;

    // Add multiples of pivot row to descending rows
    for (i=k+1; i<N; i++) { // should skip if k==N
        a = A[i*N+k];
        for (j=0; j<N; j++){
            A[i*N+j] -= a*A[k*N+j];
        }
        B[i] -= a*B[k];
    }
    k++;
    //for (i=0; i<N; i++) {
    //    printf("        [");
    //    for (j=0; j<N; j++){
    //        printf("%f ", A[i*N+j]);
    //    }
    //    printf("]  [ %f ]\n", B[i]);
    //}
    //printf("\n");
    } // end while loop for row reduction


    // Now compute X using back substitution
    for (k=N-1; k>=0; k--){
        a = 0.0;
        for (j=k+1; j<N; j++) a+= A[k*N+j]*X[j]; // Should skip when k=N
        X[k] = B[k] - a;
    }
    return 0;
}

int iterate_solve_matrix(double* A, double* B, double *X, double* X2, int N){
    /*
    Bare bones Gauss Seidel matrix solver for A*X = B 
       A : pointer to NxN matrix of doubles
       B : pointer to N vector of doubles
       X : pointer to N solution vector of doubles
       X2: pointer to N sized workspace for the algorithm
       N : problem size (int)
    */
    const int ITERATION_LIMIT=20;
    const double tol=1e-8;

    int n,i,j,timeout;
    double s1, s2,errorL2;

    // Reset x to prevent NaNs or memory garbage from polluting the iterations
    for (i=0; i<N; i++) X[i] = 0.0; 

    n=0;
    timeout=0;
    while (n<ITERATION_LIMIT){
        for (i=0; i<N; i++) X2[i] = 0.0; 

        for (i=0; i<N; i++){
            s1 = 0.0;
            s2 = 0.0;

            for (j=0; j<i; j++) s1 += A[i*N+j]*X2[j];
            for (j=i+1; j<N; j++) s2 += A[i*N+j]*X[j];

            X2[i] = (B[i] - s1 - s2)/A[i*N+i];
        }

        errorL2 = 0.0;
        for (i=0; i<N; i++) errorL2 += (X2[i]-X[i])*(X2[i]-X[i]);
        errorL2 = sqrt(errorL2); 
        for (i=0; i<N; i++) X[i] = X2[i];

        printf("Iter: %d [", n);
        for (i=0; i<N; i++) printf(" %f ",X[i]);
        printf("] error: %e \n", errorL2);

        if (errorL2<tol) {
            timeout=1;
            break;
        }
        n++;
    }
    return timeout;
}

void test_solve_matrix(){
    double XX[4];
    int N=4;
    int i;

    double A[16] = {0, 3, 5, 4, 6, 1, 7, 5, 1, 3, 2, 6, 6, 0, 2, 0};
    double X[4] = {6.0, 1.0, 4.0, 1.0};
    double B[4] = {27.0, 70.0, 23.0, 44.0};

    printf("Calling solve_matrix\n");
    solve_matrix(A, B, XX, N);

    for (i=0; i<N; i++){
        printf("%d X: %f XX: %f\n",i , X[i], XX[i]);
    }
    return;
}

void test_iterate_solve_matrix(){
    double XX[4];
    int N=4;
    int i;

    // TODO: This method does not converge if there is a zero on the centerline.
    // Can this ever happen during an equilibrium problem?

    double A[16] = {10., -1., 2., 0., -1., 11., -1., 3., 2., -1., 10., -1., 0., 3., -1., 8.};
    double X[4] = {1.0, 2.0,-1.0, 1.0};
    double X2[4]= {0.0, 0.0, 0.0, 0.0};
    double B[4] = { 6.0, 25.0,-11.0, 15.0};

    printf("Calling iterate_solve_matrix\n");
    iterate_solve_matrix(A, B, XX, X2, N);

    for (i=0; i<N; i++){
        printf("%d X: %f XX: %f\n",i , X[i], XX[i]);
    }
    return;
}

#ifdef TEST
int main(){
    test_solve_matrix();
    test_iterate_solve_matrix();
    return 0;
}
#endif
