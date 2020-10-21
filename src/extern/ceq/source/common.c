/*
C library for equilibrium chemistry calculations
 - This file contains useful utility functions shared between pt.c and rhou.c

@author: Nick Gibbons
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

const double TRACELIMIT=1e-6;   // Trace species limiter (for ns/n)
const double tol=1e-9;
const int attempts=50;

double update_limit_factor(double x, double dx, double fac){
    /*
    Compute lambda for variable updates, checking for dx==0.0 hardware exception
    Inputs:
        x      : reference quantity
        dx     : proposed raw change in quantity
        fac    : limit parameter  
    Outputs:
        lambda : update limit factor, generally x' = x + lambda*dx
    */
    double lambda;
    if (dx==0.0) return 1.0;
    lambda = fmin(1.0, fac*fabs(x)/fabs(dx));
    return lambda;
}



// FIXME: Should bi0 be increased if the species get unlocked? Maybe you need to keep a bi00
void handle_trace_species_locking(double* a, double n, int nsp, int nel, double* ns, double* bi0, double* dlnns, int verbose){
    /*
    Check for small species and lock appropriately
    Inputs:
        a      : elemental composition array [nel,nsp]
        n      : total moles/mixture kg  
        nsp    : total number of species
        nel    : total  number of elements 
    Outputs:
        ns     : species moles/mixture kg [nsp]
        bi0    : Initial Nuclear moles/mixture [nel]
        dlnns : vector of species mole/mixture corrections [nsp]
    */ 
    int s,i;
    double bi;

    for (s=0; s<nsp; s++){
        if (ns[s]/n<TRACELIMIT){
            if (verbose>1) printf("    Locking species: %d (%f)\n", s, dlnns[s]);
            ns[s] = 0.0;
            dlnns[s] = 0.0; // This species is considered converged now
        }
    }

    // Sometimes the initial bi0 arrays can contain trace amounts of species in them,
    // if the solver has gone and eliminated all of these species, we need to set bi0
    // to zero to make sure that the convergence check works.
    // This routine does that, albeit a bit crudely. Watch out for possible problems
    // if some of those elements get unlocked after this goes off!
    for (i=0; i<nel; i++){
        bi = 0.0;
        for (s=0; s<nsp; s++){
            bi += a[i*nsp + s]*ns[s];
        }
        if (bi<1e-16) {
            if (verbose>1) printf("        element bi[%d]: %f, zeroing b0\n",i,bi);
            bi0[i] = 0.0;
        }
    }
    return;
}

void composition_guess(double* a,double* M,double* X0,int nsp,int nel,double* ns,double* np,double* bi0){
    /*
    Unified setup of initial composition from mole fractions
    Inputs:
        a      : elemental composition array [nel,nsp]
        X0    : Intial Mole fractions [nsp]
        M      : species molecular masses [nsp]
        nsp    : total number of species
        nel    : total  number of elements 
    Outputs:
        ns     : species moles/mixture kg [nsp]
        n      : total moles/mixture kg  [1]
        bi0    : Initial Nuclear moles/mixture [nel]
    */ 
    int i,s;
    double M0,n;

    M0 = 0.0;
    for (s=0; s<nsp; s++) M0 += M[s]*X0[s];
    for (s=0; s<nsp; s++) ns[s] = X0[s]/M0;

    for (i=0; i<nel; i++){
        bi0[i] = 0.0;
        for (s=0; s<nsp; s++){
            bi0[i] += a[i*nsp + s]*X0[s]/M0;
            }
    }


    n = 0.0;
    for (s=0; s<nsp; s++) n += ns[s];
    for (s=0; s<nsp; s++) ns[s] = fmax(ns[s], n*TRACELIMIT*100.0);
    *np = n;

    // Auto lock species with missing elements
    // Code removed (16/10/2020)
    return;
}

int all_but_one_species_are_trace(int nsp, double* ns){
    /*
    If only one species is left in the calculation assume we have found the answer
    Inputs:
        nsp : number of species 
        ns  : species moles/mixture kg [nsp]
    Returns:
        -1 if false, equal to the index of the nontrace species otherwise
    */
    int i,s,ntrace;

    ntrace=0;
    for (s=0; s<nsp; s++) if (ns[s]==0.0) ntrace++;

    i=-1;
    if (ntrace==nsp-1) { // Pseudo convergence criteria, all the species but one are trace
        for (s=0; s<nsp; s++) if (ns[s]!=0.0) i=s;
        return i;
    }
    else {
        return -1;
    }
}

double constraint_errors(double* S,double* a,double* bi0,double* ns,int nsp,int nel,int neq,double* dlnns){
    /*
    Unified computation of current error to determine whether to break loop
    Inputs:
        S      : Corrections  [neq]
        a      : elemental composition array [nel,nsp]
        bi0    : Initial Nuclear moles/mixture [nel]
        ns     : species moles/mixture kg [nsp]
        nsp    : total number of species
        nel    : total  number of elements 
        neq    : total number of equations 
        dllns  : raw change in log(ns) [nsp]

    Outputs:
        unified error number (method as determined by this function)
    */ 
    int s,i,n,nS;
    double bi,errorrms,error;

    // Compute the change in current variables (note this is the unlimited dlnns! not the real change)
    errorrms = 0.0;
    n = 0;
    nS = neq-nel;

    for (i=0; i<nS; i++){
        errorrms += S[i]*S[i];
        n += 1;
    }

    for (s=0; s<nsp; s++){
        errorrms += dlnns[s]*dlnns[s];
        n += 1;
    }

    for (i=0; i<nel; i++) {
        bi = 0.0;
        for (s=0; s<nsp; s++) bi += a[i*nsp + s]*ns[s];
        error = bi - bi0[i];
        errorrms += error*error;
        n += 1;
    }
    errorrms /= n; // Implicit typecast, is this a good idea?
    errorrms = sqrt(errorrms);
    return errorrms;
}

void check_ill_posed_matrix_row(double* A, double* B, int neq, int k, int offset){
    /*
    Check for singular rows within the solve matrix. These may arise if all of the species containing
    one element have been zeroed out of the calculation, for example.
    Inputs:
        A      : Linear Solve Matrix [neq,neq]
        B      : Linear Solve RHS [neq]
        neq    : number of equations in the matrix
        k      : the row in question
        offset : how many entries before the pi_i's.

        Notes:
         - For example, in the pt problem the unknowns are [dln n, pi_1, pi_2, ...], so offset=1
    */
    double row_is_bad;
    int i;

    row_is_bad = 0.0;
    for (i=0; i<neq; i++) row_is_bad += fabs(A[k*neq+i]);

    // If a singular row entry is detected, set that row to a trivial equation.
    // This equation is designed to force pi_k to be zero.
    if (row_is_bad<1e-16) {
        for (i=0; i<neq; i++) A[k*neq+i] = 0.0;
        A[k*neq + k+offset] = 1.0;
        B[k] = 0.0;
    }
}
