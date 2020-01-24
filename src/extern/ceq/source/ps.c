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
#include "ps.h"

static void Assemble_Matrices(double* a,double* bi0, double pt,double st,double T,double n,double* ns,
                              int nsp, int nel,double* A, double* B, double* mu_RTs, double* H_RTs,
                              double* Cp_Rs, double* S_Rs,  double* lewis){
    /*
    Construct Iteration Matrix for reduced Newton Rhapson step, (eqn 2.24 2.28 2.26 from cea_I)
    Inputs:
        a      : elemental composition array [nel,nsp]
        bi0    : Initial nuclear composition vector [nel]
        pt     : target pressure (Pa)
        st     : target specific entropy (J/kg)
        T      : current guess of temperature (K)
        n      : current guess of mixture weight (mol/kg)
        ns     : current guess of species mixture weights (mol/kg) [nsp]
        nsp    : total number of species
        nel    : total  number of elements 
        mu_RTs : Chemical potential of each species, divided by RT [nsp]
        H_RTs  : Molar enthalpy of each species, divided by RT [nsp]
        Cp_Rs  : Molar specific heat @ constant pressure of each species, divided by Ru [nsp]
        S_Rs   : Molar entropy at of each species, divided by Ru [nsp]
        lewis  : Lewis thermodynamic data of each species
    Outputs:
        A      : Jacobian Inversion Matrix [neq,neq]
        B      : RHS Jacobian Inversion Vector [neq]

    Notes:
        Fast index varies over variables: dlnn, dlnT, pi1, pi2, ..., pi_nel
        Slow index varies over equations 2.24_1, 2.24_2, ..., 2.24_nel, 2.25, 2.28
    */
    double lnns, lnn, lnp, akjaijnj, akjnjmuj, mus_RTj, bk;
    double akjnjHj, nsHs, aijnjSj, njSj, njCpj, njHjSj, njSjmuj;
    double nss, nsmus, coeffsum, sk, G0_RTs, S0_Rs;
    int k,neq,i,j,s,nep;
    double *lp;
    neq = nel+2;
    lnn = log(n);
    lnp = log(pt/1e5); // Standard pressure for the tables is one BAR

    sk = 0.0;
    for (s=0; s<nsp; s++){
        lp = lewis + 9*3*s;
        H_RTs[s] = compute_H0_RT(T, lp);
        Cp_Rs[s] = compute_Cp0_R(T, lp);
        G0_RTs = compute_G0_RT(T, lp);
        S0_Rs  = compute_S0_R(T, lp);                    // entropy at one BAR

        if (ns[s]!=0.0){
            mu_RTs[s] = G0_RTs + log(ns[s]) - lnn + lnp;
            S_Rs[s]   = S0_Rs  - log(ns[s]) + lnn - lnp;     // entropy at current pressure
        }
        else {
            // Triggered if ns[s]==0. Note that in this situation mu_RTs and S_Rs are actually
            // equal to log(0.0) which is negative infinity. However, mu_RTs and S_Rs are always
            // multiplied by ns before being used, and x*log(x) -> 0 as x->0, so we can set them
            // to whatever to keep the floating point math happy
            mu_RTs[s] = 0.0; 
            S_Rs[s]   = 0.0;
        }

        sk += ns[s]*S_Rs[s]*Ru; 
    }

    // FIXME: Consider having one function for each of these rather than a single huge one
    // Equation 2.24: k-> equation index, i-> variable index
    for (k=0; k<nel; k++){

        if (bi0[k]<1e-16) { // Check for missing missing element equations and Lock
            for (i=0; i<neq; i++) A[k*neq+i] = 0.0;
            A[k*neq + k+2] = 1.0;  
            B[k] = 0.0;
            continue;
        }

        bk = 0.0; for (s=0; s<nsp; s++) bk += a[k*nsp + s]*ns[s];
        A[k*neq + 0] = bk;       // dlnn entry

        akjnjHj = 0.0;
        for (j=0; j<nsp; j++) akjnjHj += a[k*nsp+j]*ns[j]*H_RTs[j];
        A[k*neq + 1] = akjnjHj;  // dlnT entry

        for (i=0; i<nel; i++){
            akjaijnj = 0.0;
            for (j=0; j<nsp; j++){
                akjaijnj += a[k*nsp+j]*a[i*nsp+j]*ns[j];
            }
            A[k*neq + i+2] = akjaijnj; // pii entry, i+2 because dlnn, dlnT come first
        }

        akjnjmuj = 0.0;
        for (j=0; j<nsp; j++){
            //if (ns[s]==0.0) continue;
            //mus_RTj = G0_RTs[j] + log(ns[j]) - lnn + lnp;
            akjnjmuj += a[k*nsp+j]*ns[j]*mu_RTs[j];
        }
        B[k] = bi0[k] - bk + akjnjmuj; // rhs
    }

    // Equation 2.26 - > (only the pii entries)
    for (k=0; k<nel; k++){
        bk = 0.0; for (s=0; s<nsp; s++) bk += a[k*nsp + s]*ns[s];
        A[nel*neq + k+2] = bk;
    }

    // Equation 2.26 - > (now the rest)
    nss = 0.0;
    nsmus = 0.0;
    nsHs = 0.0;
    for (s=0; s<nsp; s++){
        //if (ns[s]==0.0) continue;
        //mus_RTj = G0_RTs[s] + log(ns[s]) - lnn + lnp;
        nss += ns[s];
        nsmus += ns[s]*mu_RTs[s];
        nsHs += ns[s]*H_RTs[s];
    }
    A[nel*neq + 0] = nss - n; // dlnn entry
    A[nel*neq + 1] = nsHs;     // dlnT entry
    B[nel] = n - nss + nsmus;

    // Equation 2.28 (pii entries)
    nep = nel+1;
    for (i=0; i<nel; i++){
        aijnjSj = 0.0;
        for (j=0; j<nsp; j++){
            aijnjSj += a[i*nsp+j]*ns[j]*S_Rs[j];
        }
        A[nep*neq + i+2] = aijnjSj; // pii entry, i+2 because dlnn, dlnT come first
    }

    // Now the rest of 2.28
    njSj = 0.0;
    njCpj = 0.0;
    njHjSj = 0.0;
    njSjmuj = 0.0;
    for (j=0; j<nsp; j++){
        njSj    += ns[j]*S_Rs[j];
        njCpj   += ns[j]*Cp_Rs[j];
        njHjSj  += ns[j]*H_RTs[j]*S_Rs[j];
        njSjmuj += ns[j]*S_Rs[j]*mu_RTs[j];
    }
    A[nep*neq + 0] = njSj; // dlnn entry
    A[nep*neq + 1] = njCpj + njHjSj; // dnlT entry
    B[nep] = (st - sk)/Ru + n - nss + njSjmuj; // rhs
    
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

static void species_corrections(double* S,double* a,double* mu_RTs,double* H_RTs, double p,double n,
                                double* ns, int nsp, int nel, double* dlnns, int verbose){
    /*
    Compute delta_log(ns) from the reduced iteration equations from 
    equation 2.18 using the other deltas in S
    Inputs:
        S      : Corrections (dlog(n), pi1, pi2, pi3, ... , dlog(T)) [nel+1]
        a      : elemental composition array [nel,nsp]
        mu_RTs : Chemical potential of each species, divided by RT [nsp]
        H_RTs  : Molar enthalpy of each species, divided by RT [nsp]
        p      : pressure 
        n      : total moles/mixture kg 
        ns     : species moles/mixture kg [nsp]
        nsp    : total number of species
        nel    : total  number of elements 

    Outputs:
        dllns : change in log(ns) [nsp]
    */
    double dlnn,aispii,lnn,lnp,dlnT;
    int s,i;
    dlnn = S[0];
    dlnT = S[1];
    lnn = log(n);
    lnp = log(p/1e5);

    for (s=0; s<nsp; s++) {
        if (ns[s]==0.0) { dlnns[s] = 0.0; continue;}
        //mu_RTs = G0_RTs[s] + log(ns[s]) - lnn + lnp;

        aispii = 0.0;
        for (i=0; i<nel; i++){
            aispii += a[i*nsp+s]*S[i+2]; // S[i+2] = pi_i, the lagrange multiplier
        }
        dlnns[s] = -mu_RTs[s] + dlnn + aispii + H_RTs[s]*dlnT;
        //printf("    dlnns[%d] = %f (%f %f %f)\n", s, dlnns[s], -mu_RTs, dlnn, aispii);
    }
    return; 
}

static void handle_singularity(double* S,double* a,double* mu_RTs,double* H_RTs,double p,double n,
                               double* ns, int nsp, int nel, double* dlnns, int verbose){
    /*
    Handle crash by possibly resetting species compositions to fix
    Inputs:
        S      : Corrections array (pi1, pi2, pi3 ... dlog(n) [nel+1]
        a      : elemental composition array [nel,nsp]
        mu_RTs : Chemical potential of each species, divided by RT [nsp]
        H_RTs  : Molar enthalpy of each species, divided by RT [nsp]
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

        ns[s] = fmax(RESET*n*TRACELIMIT, ns[s]); // Reset trace trace species to small but finite number
        species_corrections(S,a,mu_RTs,H_RTs,p,n,ns,nsp,nel,dlnns,0); // approximately predict dlnns

        if (dlnns[s]<0.0) ns[s] = 0.0; // Re-zero any species with negative predicted dlnns
        if (verbose>1) printf("   faux dlnns: %f changed to: %e \n", dlnns[s], ns[s]);
    }
    return; 
}

static void update_unknowns(double* S,double* dlnns,int nsp,int nel,double* ns,double* n,double* T,
                            int verbose){
    /*
    Add corrections to unknown values (ignoring lagrange multipliers)
    Inputs:
        S : vector of corrections from matrix step [nel+1]
        dlnns : vector of species mole/mixture corrections [nsp]
        nsp : number of species
    Outputs:
        ns : vector of species mole/mixtures [nsp]
        n  : pointer to total moles/mixture (passed by reference!) [1]
        T  : pointer to current Temperature (passed by reference!) [1]
    */
    int s,i;
    double lnns,lnn,n_copy,lambda,newns,rdlnns,bi,lnT,T_copy;
    const char pstring[] = "  s: %d lnns: % f rdlnns: % f dlnns: %f TR: % e lambda: % f\n"; 

    lnn = log(*n); // compute the log of the thing n is pointing to
    lambda = fmin(1.0, 0.5*fabs(lnn)/fabs(S[0]));
    n_copy = exp(lnn + lambda*S[0]); 
    *n = n_copy;   // thing n in pointing to set to n_copy

    lnT = log(*T); // compute the log of the thing T is pointing to
    lambda = fmin(1.0, 0.5*fabs(lnT)/fabs(S[1]));
    T_copy = exp(lnT + lambda*S[1]); 
    *T = T_copy;   // thing T is pointing to set to T_copy;

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

static double temperature_guess(int nsp, double st, double pt, double n, double* ns, double* lewis){
    /*
    Guess a first iteration temperature assuming constant Cp from 298 K and target p and s
    Inputs:
        nsp   : Number of species
        st    : Target mixture specific entropy (J/kg)
        pt    : Target pressure (Pa)
        n     : Mixture mixture weight [nsp]
        ns    : Intiial composition mixture weights [nsp]
        lewis : Nasa Lewis Thermodynamic Database Data [nsp*3*9]

    Output:
        T : Temperature Guess (K)
    */
    int s;
    double* lp;
    double cp,s0,Cps298,S0s298,nss,T;

    cp = 0.0;
    s0 = 0.0;
    for (s=0; s<nsp; s++){
        lp = lewis + 9*3*s;
        Cps298 = compute_Cp0_R(298.15, lp)*Ru;
        S0s298 = compute_S0_R(298.15, lp)*Ru;

        nss = ns[s];
        cp += nss*Cps298;
        s0 += nss*S0s298;
    }

    // Guess with fixed cp at current composition
    T = 298.15*exp((st - s0 + n*Ru*log(pt/1e5))/cp)/5.0; // See 06/08/2019
    //printf("s0: %e\n", s0);
    T = fmin(fmax(T, 200.0),20000.0); // Limit in case of bad initial state.
    return T;
}

int solve_ps(double pt,double st,double* X0,int nsp,int nel,double* lewis,double* M,double* a,
               double* X1, double* Teq, int verbose){
    /*
    Compute the equilibrium composition X1 at a fixed pressure and specific entropy
    Inputs:
        pt    : target pressre (Pa)
        st    : target specific entropy (J/kg)
        X0    : Intial Mole fractions [nsp]
        nsp   : number of species 
        nel   : number of elements 
        lewis : Nasa Lewis Thermodynamic Database Data [nsp*3*9]
        M     : Molar Mass of each species (kg/mol) [nsp]
        a     : elemental composition array [nel,nsp]
        verbose: print debugging information

    Output:
        X1  : Equilibrium Mole Fraction [nsp]  
        Teq: Equilibrium Temperature 
    */
    double *A, *B, *S, *mu_RTs, *H_RTs, *S_Rs, *Cp_Rs, *ns, *bi0, *dlnns; // Dynamic arrays
    int neq,s,i,k,errorcode,ntrace;
    double M0,n,M1,errorL2,errorL22,thing,T,errorrms;

    errorcode=0;
    neq= nel+2;
    A     = (double*) malloc(sizeof(double)*neq*neq); // Iteration Jacobian
    B     = (double*) malloc(sizeof(double)*neq);     // Iteration RHS
    S     = (double*) malloc(sizeof(double)*neq);     // Iteration unknown vector
    mu_RTs= (double*) malloc(sizeof(double)*nsp);     // Species chemical potential
    H_RTs = (double*) malloc(sizeof(double)*nsp);     // Species Molar Enthalpy 
    S_Rs  = (double*) malloc(sizeof(double)*nsp);     // Species Molar Entropy
    Cp_Rs = (double*) malloc(sizeof(double)*nsp);     // Species Specific Heat @ Const. pressure 
    ns    = (double*) malloc(sizeof(double)*nsp);     // Species moles/mixture mass
    bi0   = (double*) malloc(sizeof(double)*nel);     // starting composition coefficients
    dlnns = (double*) malloc(sizeof(double)*nsp);     // raw change in log(ns)

    composition_guess(a, M, X0, nsp, nel, ns, &n, bi0);
    T = temperature_guess(nsp, st, pt, n, ns, lewis);
    if (verbose>0) printf("Guess T from ps: %f\n", T);

    // Begin Iterations
    for (k=0; k<attempts; k++){
        Assemble_Matrices(a,bi0,pt,st,T,n,ns,nsp,nel,A,B,mu_RTs,H_RTs,Cp_Rs,S_Rs,lewis);
        errorcode = solve_matrix(A, B, S, neq);
        if (errorcode!=0) {
            handle_singularity(S,a,mu_RTs,H_RTs,pt,n,ns,nsp,nel,dlnns,verbose);
            continue;
        }
        species_corrections(S,a,mu_RTs,H_RTs,pt,n,ns,nsp,nel,dlnns,verbose); // FIXME? Using old n,T
        update_unknowns(S,dlnns,nsp,nel,ns,&n,&T,verbose);
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

        if (k>=attempts-1) {
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
    *Teq = T;

    free(A);
    free(B);
    free(S);     
    free(mu_RTs);
    free(H_RTs);
    free(S_Rs);
    free(Cp_Rs);
    free(ns);
    free(bi0);
    free(dlnns);
    return errorcode;
}

#ifdef TEST
int main(){
    printf("Called main in ps.c! Bad!\n");
    return 0;
}
#endif
