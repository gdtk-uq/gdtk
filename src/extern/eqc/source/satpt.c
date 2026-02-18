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
#include "satpt.h"

static void Assemble_Matrices(double* a, double* bi0, double* G0_RTs, double p, double* ns, double* lnns,
                              double n, int nsp, int nel, int ic, double Gc_RT, double* A, double* B){
    /*
    Construct Iteration Matrix for reduced Newton Rhapson step, (eqn 2.24 and 2.26 from cea_I)
    Order of rows is (35) from eqc first, then 36, then the condensed species eqn.
    Order of columns is dln(n), then the piis, then the dbc0 correction.
    */
    double lnn, lnp, akjaijnj, akjnjmuj, mus_RTj, bk;
    double nss, nsmus;
    int k,neq,i,j,s;
    neq = nel+2;
    lnn = log(n);
    lnp = log(p/1e5); // Standard pressure for the tables is one BAR

    // Equation 2.24: k-> equation index, i-> variable index
    for (k=0; k<nel; k++){
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
            mus_RTj = G0_RTs[j] + lnns[j] - lnn + lnp;
            akjnjmuj += a[k*nsp+j]*ns[j]*mus_RTj;

        }
        // Add unknown bc0 term if this k is the unknown elemental constituent
        if (k==ic) {
            A[k*neq + nel+1] = -bi0[ic];
        } else {
            A[k*neq + nel+1] = 0.0; // Need to zero the extra matrix column
        }

        B[k] = bi0[k] - bk + akjnjmuj;// CEA equations
        check_ill_posed_matrix_row(A, B, neq, k, 1);
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
        mus_RTj = G0_RTs[j] + lnns[j] - lnn + lnp; // I guess its okay to compute this again
        nss += ns[j];
        nsmus += ns[j]*mus_RTj;
    }
    A[nel*neq + 0]  = nss - n;
    B[nel] = n - nss + nsmus;
    A[nel*neq + nel+1] = 0.0; // Need to zero the extra matrix column

    // New equation for the condensed carbon Gibbs energy
    // See derivation from 04/02/26
    int iic = nel+1; // this eqn row in the matrix
    for (i=0; i<neq; i++) {
        A[iic*neq + i] = 0.0;
    }
    A[iic*neq + ic+1] = 1.0;
    B[iic] = Gc_RT;
    
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

static double compute_residual(double* pij, double* a, double* G0_RTs, double p, double T, double n,
                               double* ns, double* lnns, double* bi0, int nsp, int nel, int verbose){
    /*
    Compute the L2 of the PT Lagrangian derivatives. Note as per the paper, we actually compute
    ns[s] times equation (44). This ensures that the equations are nonsingular for ns[s] -> 0.0
    but still have the same root as the original equations.

    Inputs:
        pij    : Corrections (dlog(n), pi1, pi2, pi3 ...)  [nel+1]
        a      : elemental composition array [nel,nsp]
        G0_RTs : Gibbs free energy of each species, divided by RT [nsp]
        p      : pressure  (Pa)
        T      : Temperature (K)
        n      : total moles/mixture kg
        ns     : species moles/mixture kg [nsp]
        lnns   : natural log of the species moles/mixture kg [nsp]
        nsp    : total number of species
        nel    : total number of elements

    Outputs:
        double : L2 norm of the residual of the nonlinear equations being solved.
    */

    double residual = 0.0;

    for (int s=0; s<nsp; s++){
        // ns*log(ns) is badly behaved at ns==0.0, but should be fine at any other nonnegative number
        double nss = ns[s];
        double nslogns = nss*lnns[s];

        // Equation ?? from the eqc paper
        double Rs = nss*G0_RTs[s] + nslogns - nss*log(n) + nss*log(p/1e5);         // CEA equations
        double pijj = 0.0;
        for (int j=0; j<nel; j++){
            pijj += pij[1+j]*a[j*nsp + s];
        }
        Rs -= nss*pijj;
        if (verbose>1) printf("    Fs[%d]=%e\n", s, Rs);
        residual += Rs*Rs;
    }

    // Equation ?? from the eqc paper
    for (int j=0; j<nel; j++){
        double Rs = 0.0;
        for (int s=0; s<nsp; s++){
            Rs += a[j*nsp + s]*ns[s];
            if (verbose>1) printf("        a[%d,%d]=%e n[%d]=%e\n", j, s, a[j*nsp+s], s, ns[s]);
        }
        Rs -= bi0[j];
        if (verbose>1) printf("    Fj[%d]=%e bi0[%d]=%e\n", j, Rs, j, bi0[j]);
        residual += Rs*Rs;
    }

    // Equation ?? from the eqc paper
    double Rs = n;
    for (int s=0; s<nsp; s++){
        Rs -= ns[s];
    }
    if (verbose>1) printf("    Fn=%e\n", Rs);
    residual += Rs*Rs;

    return sqrt(residual);
}

static void species_corrections(double* S,double* a,double* G0_RTs,double p,double n,double* ns,double* lnns,
                        int nsp, int nel, double* dlnns, int verbose){
    /*
    Compute delta_log(ns) from the reduced iteration equations from 
    equation 2.18m using the other deltas in S
    Inputs:
        S      : Corrections (dlog(n), pi1, pi2, pi3 ...)  [nel+1]
        a      : elemental composition array [nel,nsp]
        G0_RTs : Gibbs free energy of each species, divided by RT [nsp]
        p      : pressure 
        n      : total moles/mixture kg 
        ns     : species moles/mixture kg [nsp]
        lnns   : natural log of the species moles/mixture kg [nsp]
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
        mu_RTs = G0_RTs[s] + lnns[s] - lnn + lnp;

        aispii = 0.0;
        for (i=0; i<nel; i++){
            aispii += a[i*nsp+s]*S[i+1]; // S[i+1] = pi_i, the lagrange multiplier
        }
        dlnns[s] = -mu_RTs + dlnn + aispii;
    }
    return; 
}

static void handle_singularity(double* S,double* a,double* G0_RTs,double p,double n,double* ns,double* lnns,
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
        lnns   : natural log of the species moles/mixture kg [nsp]
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
        species_corrections(S, a, G0_RTs, p, n, ns, lnns, nsp, nel, dlnns, 0); // approximately predict dlnns
        if (dlnns[s]<0.0) ns[s] = 0.0; // Re-zero any species with negative predicted dlnns
        if (verbose>1) printf("   faux dlnns: %f changed to: %e \n", dlnns[s], ns[s]);
    }
    return; 
}

static void update_unknowns(double* S,double* dlnns,int nsp,int nel,int ic,double* ns,double* lnns,double* n,double* b0,double* lnb,int verbose){
    /*
    Add corrections to unknown values (ignoring lagrange multipliers)
    Inputs:
        S     : vector of corrections from matrix step [nel+2]
        dlnns : vector of species mole/mixture corrections [nsp]
        nsp   : number of species
        nsp   : number of elements
        ic    : index of unknown saturating element
    Outputs:
        ns    : vector of species mole/mixtures [nsp]
        lnns  : natural log of the species moles/mixture kg [nsp]
        n     : pointer to total moles/mixture (passed by reference!) [1]
        b0    : pointer to bi0 constraint vector [nel]
        lnb   : log(bi0[ic]), the log of the unknown saturating quantity
    */
    int s;
    double newlnns,lnn,n_copy,lambda,newns,rdlnns,lnbcopy,newlnb;
    const char pstring[] = "  s: %d lnns: % f rdlnns: % f dlnns: %f TR: % e lambda: % f\n";

    lnn = log(*n); // compute the log of the thing n is pointing to
    lambda = update_limit_factor(lnn, S[0], relaxation_limit);
    n_copy = exp(lnn + lambda*S[0]); 
    *n = n_copy;   // thing pointed to by n set to exp(lnn + S[0]);

    lnbcopy = *lnb;
    lambda = update_limit_factor(lnbcopy, S[nel+1], relaxation_limit);
    newlnb = lnbcopy + lambda*S[nel+1];
    *lnb = newlnb;
    b0[ic] = exp(newlnb);

    for (s=0; s<nsp; s++){
        if (ns[s]==0.0) {
            if (verbose>1) printf(pstring, s, 0.0, 0.0, dlnns[s], 0.0, 0.0);
            dlnns[s] = 0.0;
            continue;
        }
        lambda = update_limit_factor(lnn, dlnns[s], relaxation_limit);
        newlnns = lnns[s] + lambda*dlnns[s];
        newns = exp(newlnns);
        rdlnns = newlnns - lnns[s];
        ns[s] = newns;
        lnns[s] = newlnns;

        if (verbose>1) printf(pstring, s, lnns[s], rdlnns, dlnns[s], ns[s]/n_copy/TRACELIMIT, lambda);
    }

    return;
}

int solve_satpt(double p,double T,double Gc,int ic,double* X0,int nsp,int nel,double* lewis,double* M,double* a,double* X1, int verbose){
    /*
    Compute the equilibrium composition X1 at a fixed temperature and pressure
    Inputs:
        p     : Pressure (Pa)
        T     : Temperature (K)
        Gc    : Gibbs energy of saturating condensed species
        ic    : index of saturating element
        X0    : Intial Mole fractions [nsp]
        nsp   : number of species 
        nel   : number of elements 
        lewis : Nasa Lewis Thermodynamic Database Data [nsp*3*9]
        M     : Molar Mass of each species (kg/mol) [nsp]
        a     : elemental composition array [nel,nsp]
        verbose: print debugging information

    Output:
        X1    : Equilibrium Mole Fraction [nsp]
    */
    double *A, *B, *S, *G0_RTs, *ns, *lnns, *bi0, *dlnns; // Dynamic arrays
    double *lp;
    int neq,s,i,k,errorcode,solvecode;
    double n,M1,errorrms;
    double Gc_RT, lnbc0;

    errorcode=0;
    neq= nel+2;
    errorrms=1e99;
    A     = (double*) malloc(sizeof(double)*neq*neq); // Iteration Jacobian
    B     = (double*) malloc(sizeof(double)*neq);     // Iteration RHS
    S     = (double*) malloc(sizeof(double)*neq);     // Iteration unknown vector
    G0_RTs= (double*) malloc(sizeof(double)*nsp);     // Species Gibbs Free Energy
    ns    = (double*) malloc(sizeof(double)*nsp);     // Species moles/mixture mass
    lnns  = (double*) malloc(sizeof(double)*nsp);     // Log of the species specific molarities
    bi0   = (double*) malloc(sizeof(double)*nel);     // starting composition coefficients
    dlnns = (double*) malloc(sizeof(double)*nsp);     // starting composition coefficients

    Gc_RT = Gc/Ru/T;

    composition_guess(a, M, X0, nsp, nel, ns, &n, bi0);
    for (s=0; s<nsp; s++) lnns[s] = log(ns[s]);
    n*=1.1;
    lnbc0 = log(bi0[ic]);

    for (s=0; s<nsp; s++){
        lp = lewis + 9*3*s;
        G0_RTs[s] = compute_G0_RT(T, lp);
    }

    // Main Iteration Loop: 
    for (k=0; k<=attempts; k++){
        // 1: Perform an update of the equations
        Assemble_Matrices(a, bi0, G0_RTs, p, ns, lnns, n, nsp, nel, ic, Gc_RT, A, B);
        solvecode = solve_matrix(A, B, S, neq);
        if (solvecode!=0) {
            handle_singularity(S, a, G0_RTs, p, n, ns, lnns, nsp, nel, dlnns, verbose);
            continue;
        }
        species_corrections(S, a, G0_RTs, p, n, ns, lnns, nsp, nel, dlnns, verbose);
        update_unknowns(S, dlnns, nsp, nel, ic, ns, lnns, &n, bi0, &lnbc0, verbose);
        //handle_trace_species_locking(a, n, nsp, nel, ns, bi0, dlnns, verbose);
        errorrms = compute_residual(S, a, G0_RTs, p, T, n, ns, lnns, bi0, nsp, nel, verbose);


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
    free(lnns);
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
