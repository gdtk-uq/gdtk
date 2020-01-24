/*
C library for equilibrium chemistry calculations
    References:
        "Computer Program for Calculation of Complex Equilibrium Compositions and Applications"
        Nasa Reference Publication 1311, October 1995
        Sanford Gordon and Bonnie J. McBride

        "NASA Glenn Coefficients for Calculating Thermodynamic Properties of Individual Species"
        NASA/TP - 2002-211556, September 2002
        Bonnie J. McBride, Michael J. Zehe, and Sanford Gordon

@author: Nick Gibbons
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "thermo.h"
#include "linalg.h"
#include "common.h"
#include "pt.h"

static void Assemble_Matrices(double* a,double* bi0,double* G0_RTs,double p,double* ns,double n,
                          int nsp,int nel,double* A, double* B){
    /*
    Construct Iteration Matrix for reduced Newton Rhapson step, (eqn 2.24 and 2.26 from cea_I)
    */
    double lnns, lnn, lnp, akjaijnj, akjnjmuj, mus_RTj, bk;
    double nss, nsmus, coeffsum;
    int k,neq,i,j,s;
    neq = nel+1;
    lnn = log(n);
    lnp = log(p/1e5); // Standard pressure for the tables is one BAR

    // Equation 2.24: k-> equation index, i-> variable index
    for (k=0; k<nel; k++){

        if (bi0[k]<1e-16) { // Check for missing missing element equations and Lock
            for (i=0; i<neq; i++) A[k*neq+i] = 0.0;
            A[k*neq + k+1] = 1.0;   
            B[k] = 0.0;
            continue;
        }

        bk = 0.0; for (s=0; s<nsp; s++) bk += a[k*nsp + s]*ns[s];
        A[k*neq + 0] = bk;

        for (i=0; i<nel; i++){
            akjaijnj = 0.0;
            for (j=0; j<nsp; j++){
                akjaijnj += a[k*nsp+j]*a[i*nsp+j]*ns[j];
            }
            A[k*neq + i+1] = akjaijnj;
        }
        akjnjmuj = 0.0;
        for (j=0; j<nsp; j++){
            if (ns[j]==0.0) continue;
            mus_RTj = G0_RTs[j] + log(ns[j]) - lnn + lnp;
            akjnjmuj += a[k*nsp+j]*ns[j]*mus_RTj;

        }
        B[k] = bi0[k] - bk + akjnjmuj;
    }

    // Equation 2.26 - > (only the pii entries, we're highjacking k here to go across the last row)
    for (k=0; k<nel; k++){
        bk = 0.0; for (s=0; s<nsp; s++) bk += a[k*nsp + s]*ns[s];
        A[nel*neq + k+1] = bk;
    }

    // Equation 2.26 - > (now the rest)
    nss = 0.0;
    nsmus = 0.0;
    for (j=0; j<nsp; j++){
        if (ns[j]==0.0) continue;
        mus_RTj = G0_RTs[j] + log(ns[j]) - lnn + lnp; // I guess its okay to compute this again
        nss += ns[j];
        nsmus += ns[j]*mus_RTj;
    }
    A[nel*neq + 0]  = nss - n;
    B[nel] = n - nss + nsmus;
    
    //for (i=0; i<neq; i++){
    //    printf("    [");
    //    for (j=0; j<neq; j++){
    //        printf("%f, ", A[i*neq+j]);
    //    }
    //    printf("  %f]", B[i]);
    //    printf("\n");
    //}
    return;
}

static void species_corrections(double* S,double* a,double* G0_RTs,double p,double n,double* ns,
                        int nsp, int nel, double* dlnns, int verbose){
    /*
    Compute delta_log(ns) from the reduced iteration equations from 
    equation 2.18m using the other deltas in S
    Inputs:
        S      : Corrections (pi1, pi2, pi3 ... dlog(n) [nel+1]
        a      : elemental composition array [nel,nsp]
        G0_RTs : Gibbs free energy of each species, divided by RT [nsp]
        p      : pressure 
        n      : total moles/mixture kg 
        ns     : species moles/mixture kg [nsp]
        nsp    : total number of species
        nel    : total  number of elements 

    Outputs:
        dllns : change in log(ns) [nsp]
    */
    double dlnn,aispii,mu_RTs,lnn,lnp;
    int s,i;
    dlnn = S[0];
    lnn = log(n);
    lnp = log(p/1e5);

    for (s=0; s<nsp; s++) {
        if (ns[s]==0.0) { dlnns[s] = 0.0; continue;}
        mu_RTs = G0_RTs[s] + log(ns[s]) - lnn + lnp;

        aispii = 0.0;
        for (i=0; i<nel; i++){
            aispii += a[i*nsp+s]*S[i+1]; // S[i+1] = pi_i, the lagrange multiplier
        }
        dlnns[s] = -mu_RTs + dlnn + aispii;
        //printf("    dlnns[%d] = %f (%f %f %f)\n", s, dlnns[s], -mu_RTs, dlnn, aispii);
    }
    return; 
}

static void handle_singularity(double* S,double* a,double* G0_RTs,double p,double n,double* ns,
                        int nsp, int nel, double* dlnns, int verbose){
    /*
    Handle crash by possibly resetting species compositions to fix
    Inputs:
        S      : Corrections array (pi1, pi2, pi3 ... dlog(n) [nel+1]
        a      : elemental composition array [nel,nsp]
        G0_RTs : Gibbs free energy of each species, divided by RT [nsp]
        p      : pressure 
        n      : total moles/mixture kg 
        ns     : species moles/mixture kg [nsp]
        nsp    : total number of species
        nel    : total  number of elements 
        verbose: flag to print debugging  information
    */
    const double RESET=4.0;
    int s;

    if (verbose>0) printf("    ### Singular Matrix!: Unlocking to %f\n",log(RESET*n*TRACELIMIT));

    for (s=0; s<nsp; s++){
        if (ns[s]!=0.0) continue;  // Ignore non trace species
        ns[s] = fmax(RESET*n*TRACELIMIT, ns[s]); // Reset trace trace species to a small but finite number
        species_corrections(S, a, G0_RTs, p, n, ns, nsp, nel, dlnns, 0); // approximately predict dlnns
        if (dlnns[s]<0.0) ns[s] = 0.0; // Re-zero any species with negative predicted dlnns
        if (verbose>1) printf("   faux dlnns: %f changed to: %e \n", dlnns[s], ns[s]);
    }
    return; 
}

static void update_unknowns(double* S,double* dlnns,int nsp,int nel,double* ns,double* n,int verbose){
    /*
    Add corrections to unknown values (ignoring lagrange multipliers)
    Inputs:
        S : vector of corrections from matrix step [nel+1]
        dlnns : vector of species mole/mixture corrections [nsp]
        nsp : number of species
    Outputs:
        ns : vector of species mole/mixtures [nsp]
        n  : pointer to total moles/mixture (passed by reference!) [1]
    */
    int s,i;
    double lnns,lnn,n_copy,lambda,newns,rdlnns,bi;
    const char pstring[] = "  s: %d lnns: % f rdlnns: % f dlnns: %f TR: % e lambda: % f\n"; 

    lnn = log(*n); // compute the log of the thing n is pointing to
    lambda = fmin(1.0, 0.5*fabs(lnn)/fabs(S[0]));
    n_copy = exp(lnn + lambda*S[0]); 
    *n = n_copy;   // thing pointed to by n set to exp(lnn + S[0]);

    for (s=0; s<nsp; s++){
        if (ns[s]==0.0) {
            if (verbose>1) printf(pstring, s, -1.0/0.0, 0.0, dlnns[s], 0.0, 0.0);
            dlnns[s] = 0.0;
            continue;
        }
        lnns = log(ns[s]);
        lambda = fmin(1.0, 0.5*fabs(lnn)/fabs(dlnns[s]));
        newns = exp(lnns + lambda*dlnns[s]);
        rdlnns = log(newns) - lnns;
        ns[s] = newns;
        lnns = log(newns);

        if (verbose>1) printf(pstring, s, lnns, rdlnns, dlnns[s], ns[s]/n_copy/TRACELIMIT, lambda);
    }

    return;
}

int solve_pt(double p,double T,double* X0,int nsp,int nel,double* lewis,double* M,double* a,double* X1,int verbose){
    /*
    Compute the equilibrium composition X1 at a fixed temperature and pressure
    Inputs:
        p     : Pressure (Pa)
        T     : Temperature (K)
        X0    : Intial Mole fractions [nsp]
        nsp   : number of species 
        nel   : number of elements 
        lewis : Nasa Lewis Thermodynamic Database Data [nsp*3*9]
        M     : Molar Mass of each species (kg/mol) [nsp]
        a     : elemental composition array [nel,nsp]
        verbose: print debugging information

    Output:
        X1 : Equilibrium Mole Fraction [nsp]  
    */
    double *A, *B, *S, *G0_RTs, *ns, *bi0, *dlnns; // Dynamic arrays
    double *lp;
    int neq,s,i,k,errorcode;
    double n,M1,errorrms;

    errorcode=0;
    neq= nel+1;
    A     = (double*) malloc(sizeof(double)*neq*neq); // Iteration Jacobian
    B     = (double*) malloc(sizeof(double)*neq);     // Iteration RHS
    S     = (double*) malloc(sizeof(double)*neq);     // Iteration unknown vector
    G0_RTs= (double*) malloc(sizeof(double)*nsp);     // Species Gibbs Free Energy
    ns    = (double*) malloc(sizeof(double)*nsp);     // Species moles/mixture mass
    bi0   = (double*) malloc(sizeof(double)*nel);     // starting composition coefficients
    dlnns = (double*) malloc(sizeof(double)*nsp);     // starting composition coefficients

    composition_guess(a, M, X0, nsp, nel, ns, &n, bi0);
    n*=1.1;

    for (s=0; s<nsp; s++){
        lp = lewis + 9*3*s;
        G0_RTs[s] = compute_G0_RT(T, lp);
    }

    // Main Iteration Loop: 
    for (k=0; k<=attempts; k++){
        // 1: Perform an update of the equations
        Assemble_Matrices(a, bi0, G0_RTs, p, ns, n, nsp, nel, A, B);
        errorcode = solve_matrix(A, B, S, neq);
        if (errorcode!=0) {
            handle_singularity(S, a, G0_RTs, p, n, ns, nsp, nel, dlnns, verbose);
            continue;
        }
        species_corrections(S, a, G0_RTs, p, n, ns, nsp, nel, dlnns, verbose);
        update_unknowns(S, dlnns, nsp, nel, ns, &n, verbose);
        handle_trace_species_locking(a, n, nsp, nel, ns, bi0, dlnns, verbose);
        errorrms = constraint_errors(S, a, bi0, ns, nsp, nel, neq, dlnns);

        if (verbose>0){
            printf("iter %2d: [%f]",k,n);
            for (s=0; s<nsp; s++) printf(" %f",ns[s]);
            printf("  (%e)\n", errorrms);
        }

        // Exit loop if all but one species are trace 
        i = all_but_one_species_are_trace(nsp, ns);
        if (i>=0){
            if (verbose>0) printf("Pseudo convergence! Remaining species: (%d)\n", i);
            n = 1.0/M[i];
            ns[i] = n;
            errorrms = 0.0;
            break;
        }

        // Exit loop if convergence is achieved
        if (errorrms<tol) break;

        // Exit loop if nans appearing in dlnns
        if (isnan(errorrms)) {
            printf("Solver nan'd, exiting!\n");
            errorcode=1;
            break;
        }

        // Exit loop if too many attempts are undertaken
        if (k==attempts) {
            printf("Solver not converged, exiting!\n");
            errorcode=1;
            break;
        }
    }
    
    if ((verbose>0)&&(errorcode==0)) printf("Converged in %d iter, error: %e\n", k, errorrms);
    if ((verbose>0)&&(errorcode>0)) printf("Convergence failure in %d iter, error: %e\n", k, errorrms);
    // Compute output composition
    M1 = 1.0/n;
    for (s=0; s<nsp; s++) X1[s] = M1*ns[s];

    free(A);
    free(B);
    free(S);
    free(G0_RTs);
    free(ns);
    free(bi0);
    free(dlnns);
    return errorcode;
}

#ifdef TEST
int main(){
    printf("Called main in pt.c!\n");
    return 0;
}
#endif
