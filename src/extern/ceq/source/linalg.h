#ifndef linalg_h
#define linalg_h

int solve_matrix(double* A, double* B, double *X, int N);
int iterate_solve_matrix(double* A, double* B, double *X, double *X2, int N);

#endif
