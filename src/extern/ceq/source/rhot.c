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
#include "rhou.h"

static void Assemble_Matrices(double* a,double* bi0, double rho, double T, double* ns, int nsp,
                              int nel, double* A, double* B, double* G0_RTs){
    /*
    Construct Iteration Matrix for reduced Newton Rhapson rhoT step, (eqn 2.45 from cea_I)
    */
    double mus_RTj, bk, akjaijnj, akjnjmuj;
    int k,neq,i,j,s;
    neq = nel;

    // Equation 2.45: k-> equation index, i-> variable index
    for (k=0; k<nel; k++){
        bk = 0.0; for (s=0; s<nsp; s++) bk += a[k*nsp + s]*ns[s];

        for (i=0; i<nel; i++){
            akjaijnj = 0.0;
            for (j=0; j<nsp; j++){
                akjaijnj += a[k*nsp+j]*a[i*nsp+j]*ns[j];
            }
            A[k*neq + i] = akjaijnj; // Lagrange multiplier matrix entry
        }
        akjnjmuj = 0.0;
        for (j=0; j<nsp; j++){
            if (ns[j]==0.0) continue;
            mus_RTj = G0_RTs[j] + log(rho*ns[j]*Ru*T/1e5);
            akjnjmuj += a[k*nsp+j]*ns[j]*mus_RTj;
        }
        B[k] = bi0[k] - bk + akjnjmuj; // RHS of kth equation 2.45
        check_ill_posed_matrix_row(A, B, neq, k, 0);
    }

    //for (i=0; i<neq; i++){
    //    printf("   |");
    //    for (j=0; j<neq; j++){
    //        printf("%+08.6f ", A[i*neq+j]);
    //    }
    //    printf("| %+08.6f\n", B[i]);
    //}
    //printf("______________________________________________________\n");
    return;
}

static void species_corrections(double* S, double* a, double* G0_RTs, double rho, double T,
                         double* ns, int nsp, int nel, double* dlnns, int verbose){
    /*
    Compute delta_log(ns) from the reduced iteration equations from 
    equation 2.18m using the other deltas in S
    Inputs:
        S      : Corrections (pi1, pi2, pi3 ...) [nel]
        a      : elemental composition array [nel,nsp]
        G0_RTs : Gibbs free energy of each species, divided by RT [nsp]
        rho    : goal density (kg/m3)
        T      : goal temperature (K)
        ns     : species moles/mixture kg [nsp]
        nsp    : total number of species
        nel    : total  number of elements 

    Outputs:
        dllns : change in log(ns) [nsp]
    */
    double mus_RTs,aispii;
    int s,i;

    for (s=0; s<nsp; s++) {
        if (ns[s]==0.0) { dlnns[s] = 0.0; continue;}
        mus_RTs = G0_RTs[s] + log(rho*ns[s]*Ru*T/1e5);

        aispii = 0.0;
        for (i=0; i<nel; i++){
            aispii += a[i*nsp+s]*S[i]; // S[i] = pi_i, the lagrange multiplier
        }
        dlnns[s] = -mus_RTs + aispii;
    }
    return; 
}

static void handle_singularity(double* S, double* a, double* G0_RTs, double rho, double T,
                         double* ns, double n, int nsp, int nel, double* dlnns, int verbose){
    /*
    Compute delta_log(ns) from the reduced iteration equations from 
    equation 2.18m using the other deltas in S
    Inputs:
        S      : Corrections (pi1, pi2, pi3 ... ) [nel]
        a      : elemental composition array [nel,nsp]
        G0_RTs : Gibbs free energy of each species, divided by RT [nsp]
        rho    : goal density (kg/m3)
        T      : goal temperature (K)
        ns     : species moles/mixture kg [nsp]
        n      : mixture moles/mixture kg
        nsp    : total number of species
        nel    : total  number of elements 
        dllns : change in log(ns) [nsp]
    */
    const double RESET=4.0;
    int s;

    if (verbose>1) printf("    ### Singular Matrix!: Unlocking to %f\n",log(RESET*n*TRACELIMIT));

    for (s=0; s<nsp; s++){
        if (ns[s]!=0.0) continue;  // Ignore non trace species

        ns[s] = fmax(RESET*n*TRACELIMIT, ns[s]); // Reset trace trace species to a small but finite number
        species_corrections(S,a,G0_RTs,rho,T,ns,nsp,nel,dlnns,0);
        if (dlnns[s]<0.0) ns[s] = 0.0; // Re-zero any species with negative predicted dlnns
        if (verbose>1) printf("   faux dlnns: %f changed to: %e \n", dlnns[s], ns[s]);
    }
    return; 
}

static void update_unknowns(double* S, double* dlnns, int nsp, double* ns, double* np, int verbose){
    /*
    Add corrections to unknown values (ignoring lagrange multipliers)
    Inputs:
        S : vector of corrections from matrix step [nel]
        dlnns : vector of species mole/mixture corrections [nsp]
        nsp : number of species
    Outputs:
        ns : vector of species mole/mixtures [nsp]
        np : pointer to n, total moles/mixture (passed by reference!) [1]
    */
    int s;
    double lnns,n,lnn,lambda;
    const char pstring[] = "  s: %d lnns: % f rdlnns: % f dlnns: %f TR: % e lambda: % f\n"; 
    n = *np;
    lnn=log(n);

    for (s=0; s<nsp; s++){
        if (ns[s]==0.0) {
            if (verbose>1) printf(pstring, s, 0.0, 0.0, dlnns[s], 0.0, 0.0);
            dlnns[s] = 0.0;
            continue;
        }
        lnns = log(ns[s]);
        lambda = update_limit_factor(lnn, dlnns[s], 1.0);
        ns[s] = exp(lnns + lambda*dlnns[s]);
        if (verbose>1) printf(pstring, s, lnns, lambda*dlnns[s], dlnns[s], 0.0, lambda);
    }
    n = 0.0; for (s=0; s<nsp; s++) n+=ns[s];
    *np = n;
    return;
}

int solve_rhot(double rho, double T, double* X0, int nsp, int nel, double* lewis, double* M, double* a,
               double* X1, int verbose){
    /*
    Compute the equilibrium composition X1 at a fixed volume (density) and temperature
    Inputs:
        rho   : target Density (kg/m3)
        T     : target Temperature (K)
        X0    : Intial Mole fractions [nsp]
        nsp   : number of species 
        nel   : number of elements 
        lewis : Nasa Lewis Thermodynamic Database Data [nsp*3*9]
        M     : Molar Mass of each species (kg/mol) [nsp]
        a     : elemental composition array [nel,nsp]
        verbose: print debugging information

    Output:
        X1  : Equilibrium Mole Fraction [nsp]  
    */
    double *A, *B, *S, *G0_RTs, *ns, *bi0, *dlnns; // Dynamic arrays
    double *lp;
    int neq,s,i,k,errorcode;
    double n,M1,errorrms;

    errorcode=0;
    neq= nel;
    errorrms=1e99;
    A     = (double*) malloc(sizeof(double)*neq*neq); // Iteration Jacobian
    B     = (double*) malloc(sizeof(double)*neq);     // Iteration RHS
    S     = (double*) malloc(sizeof(double)*neq);     // Iteration unknown vector
    G0_RTs= (double*) malloc(sizeof(double)*nsp);     // Species Gibbs Free Energy
    ns    = (double*) malloc(sizeof(double)*nsp);     // Species moles/mixture mass
    bi0   = (double*) malloc(sizeof(double)*nel);     // starting composition coefficients
    dlnns = (double*) malloc(sizeof(double)*nsp);     // raw change in log(ns)

    composition_guess(a, M, X0, nsp, nel, ns, &n, bi0);

    // Gibbs free energy is known to begin with because we have T
    for (s=0; s<nsp; s++){
        lp = lewis + 9*3*s;
        G0_RTs[s] = compute_G0_RT(T, lp);
    }

    // Begin Iterations
    for (k=0; k<attempts; k++){
        Assemble_Matrices(a,bi0,rho,T,ns,nsp,nel,A,B,G0_RTs);
        errorcode = solve_matrix(A, B, S, neq);
        if (errorcode!=0) {
            handle_singularity(S,a,G0_RTs,rho,T,ns,n,nsp,nel,dlnns,verbose);
            continue;
        }
        species_corrections(S,a,G0_RTs,rho,T,ns,nsp,nel,dlnns,verbose);
        update_unknowns(S, dlnns, nsp, ns, &n, verbose);
        handle_trace_species_locking(a, n, nsp, nel, ns, bi0, dlnns, verbose);
        errorrms = constraint_errors(S, a, bi0, ns, nsp, nel, neq, dlnns);

        if (verbose>0){
            printf("iter %2d: [%f]",k,T);
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

        if (k>=attempts) {
            printf("Solver not converged, exiting!\n");
            errorcode=1;
            break;
        }
    }
    
    if ((verbose>0)&&(errorcode==0)) printf("Converged in %d iter, error: %e\n", k, errorrms);
    if ((verbose>0)&&(errorcode>0)) printf("Convergence failure in %d iter, error: %e\n", k, errorrms);
    // Compute output composition
    n = 0.0;
    for (s=0; s<nsp; s++) n += ns[s];
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
    printf("Called main in rhou.c!\n");
    return 0;
}
#endif
