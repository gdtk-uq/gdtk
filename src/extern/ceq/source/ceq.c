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
#include "pt.h"
#include "rhou.h"
#include "ps.h"
#include "ceq.h"


int pt(double p,double T,double* X0,int nsp,int nel,double* lewis,double* M,double* a,double* X1, int verbose){
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
    return solve_pt(p, T, X0, nsp, nel, lewis, M, a, X1, verbose);
}

int rhou(double rho,double u,double* X0,int nsp,int nel,double* lewis,double* M,double* a,
         double* X1, double* T, int verbose){
    /*
    Compute the equilibrium composition X1 at a fixed volume and internal energy 
    Inputs:
        rho   : target Density (kg/m3)
        u     : target internal energy (J/kg)
        X0    : Intial Mole fractions [nsp]
        nsp   : number of species 
        nel   : number of elements 
        lewis : Nasa Lewis Thermodynamic Database Data [nsp*3*9]
        M     : Molar Mass of each species (kg/mol) [nsp]
        a     : elemental composition array [nel,nsp]
        verbose: print debugging information

    Output:
        X1 : Equilibrium Mole Fraction [nsp]  
        T  : Equilibrium Temperature 
    */
    return solve_rhou(rho, u, X0, nsp, nel, lewis, M, a, X1, T, verbose);
}

int ps(double pt,double st,double* X0,int nsp,int nel,double* lewis,double* M,double* a,
         double* X1, double* T, int verbose){
    /*
    Compute the equilibrium composition X1 at a fixed pressure and specific entropy 
    Inputs:
        pt    : target pressure (Pa)
        st    : target specific entropy (J/kg)
        X0    : Intial Mole fractions [nsp]
        nsp   : number of species 
        nel   : number of elements 
        lewis : Nasa Lewis Thermodynamic Database Data [nsp*3*9]
        M     : Molar Mass of each species (kg/mol) [nsp]
        a     : elemental composition array [nel,nsp]
        verbose: print debugging information

    Output:
        X1 : Equilibrium Mole Fraction [nsp]  
        T  : Equilibrium Temperature 
    */
    return solve_ps(pt, st, X0, nsp, nel, lewis, M, a, X1, T, verbose);
}

double get_u(double T, double* X, int nsp, double* lewis, double* M){
    /*
    Compute thermal equilibrium u from known composition and primitives
    Inputs:
        T     : Temperature (K)
        X     : Composition [nsp]
        nsp   : number of species 
        lewis : Nasa Lewis Thermodynamic Database Data [nsp*3*9]
        M     : Molar Mass of each species (kg/mol) [nsp]
        verbose: print debugging information

    Output:
        u : internal energy per unit mass
    */
    int s;
    double Mmix, u, ns, U0_RTs;
    double* lp;

    Mmix = 0.0; for (s=0; s<nsp; s++) Mmix+=X[s]*M[s];
    
    u = 0.0;
    for (s=0; s<nsp; s++){
        ns = X[s]/Mmix;
        lp = lewis + 9*3*s;
        U0_RTs = compute_H0_RT(T, lp) - 1.0;
        u += ns*U0_RTs*Ru*T;
    }
    return u;
}

double get_h(double T, double* X, int nsp, double* lewis, double* M){
    /*
    Compute thermal equilibrium h from known composition and primitives
    Inputs:
        T     : Temperature (K)
        X     : Composition [nsp]
        nsp   : number of species 
        lewis : Nasa Lewis Thermodynamic Database Data [nsp*3*9]
        M     : Molar Mass of each species (kg/mol) [nsp]
        verbose: print debugging information

    Output:
        h : enthalpy per unit mass
    */
    int s;
    double Mmix, h, ns, H0_RTs;
    double* lp;

    Mmix = 0.0; for (s=0; s<nsp; s++) Mmix+=X[s]*M[s];
    
    h = 0.0;
    for (s=0; s<nsp; s++){
        ns = X[s]/Mmix;
        lp = lewis + 9*3*s;
        H0_RTs = compute_H0_RT(T, lp);
        h += ns*H0_RTs*Ru*T;
    }
    return h;
}

double get_cp(double T, double* X, int nsp, double* lewis, double* M){
    /*
    Compute thermal equilibrium cp from known composition and primitives
    Inputs:
        T     : Temperature (K)
        X     : Composition [nsp]
        nsp   : number of species 
        lewis : Nasa Lewis Thermodynamic Database Data [nsp*3*9]
        M     : Molar Mass of each species (kg/mol) [nsp]
        verbose: print debugging information

    Output:
        cp : specific heat at constant pressure per unit mass
    */
    int s;
    double Mmix, ns, cp,Cp0_Rs;
    double* lp;

    Mmix = 0.0; for (s=0; s<nsp; s++) Mmix+=X[s]*M[s];
    
    cp = 0.0;
    for (s=0; s<nsp; s++){
        ns = X[s]/Mmix;
        lp = lewis + 9*3*s;
        Cp0_Rs = compute_Cp0_R(T, lp);
        cp += ns*Cp0_Rs*Ru;
    }
    return cp;
}

double get_s0(double T, double* X, int nsp, double* lewis, double* M){
    /*
    Compute specific entropy @ 1 BAR from known composition and primitives
    Inputs:
        T     : Temperature (K)
        X     : Composition [nsp]
        nsp   : number of species 
        lewis : Nasa Lewis Thermodynamic Database Data [nsp*3*9]
        M     : Molar Mass of each species (kg/mol) [nsp]
        verbose: print debugging information

    Output:
        s0 : specific entropy at one Bar and temperature T (J/kg/K)
    */
    int s;
    double Mmix, ns, s0, S0_Rs;
    double* lp;

    Mmix = 0.0; for (s=0; s<nsp; s++) Mmix+=X[s]*M[s];
    
    s0 = 0.0;
    for (s=0; s<nsp; s++){
        ns = X[s]/Mmix;
        lp = lewis + 9*3*s;
        S0_Rs = compute_S0_R(T, lp);
        s0 += ns*S0_Rs*Ru;
    }
    return s0;
}

double get_s(double T, double p, double* X, int nsp, double* lewis, double* M){
    /*
    Compute specific entropy from known composition and primitives
    Inputs:
        T     : Temperature (K)
        p     : Pressure (Pa)
        X     : Composition Mole Fractions [nsp]
        nsp   : number of species 
        lewis : Nasa Lewis Thermodynamic Database Data [nsp*3*9]
        M     : Molar Mass of each species (kg/mol) [nsp]
        verbose: print debugging information

    Output:
        smix : specific entropy of mixture (J/kg/K)
    */
    int s;
    double Mmix, ns, n, S0_Rs, S_Rs, lnn, lnp, smix;
    double* lp;

    Mmix = 0.0; for (s=0; s<nsp; s++) Mmix += X[s]*M[s];
    n = 0.0;  for (s=0; s<nsp; s++) n += X[s]/Mmix;
    lnn = log(n);
    lnp = log(p/1e5);
    
    smix = 0.0;
    for (s=0; s<nsp; s++){
        ns = X[s]/Mmix;
        if (ns==0.0) continue;

        lp = lewis + 9*3*s;
        S0_Rs = compute_S0_R(T, lp);             // entropy @ one BAR divided by Ru
        S_Rs = S0_Rs  - log(ns) + lnn - lnp;     // entropy at current pressure divided by Ru
        smix += ns*S_Rs*Ru;
    }
    return smix;
}

int batch_pt(int N, double* p,double* T,double* X0,int nsp,int nel,double* lewis,double* M,double* a,
             double* X1, int verbose){
    /*
    Compute the equilibrium composition X1 at an array of fixed temperatures and pressures
    Inputs:
        N     : Number of T and p's
        p     : Pressure (Pa) [N]
        T     : Temperature (K) [N]
        X0    : Intial Mole fractions [N,nsp]
        nsp   : number of species 
        nel   : number of elements 
        lewis : Nasa Lewis Thermodynamic Database Data [nsp*3*9]
        M     : Molar Mass of each species (kg/mol) [nsp]
        a     : elemental composition array [nel,nsp]
        verbose: print debugging information

    Output:
        X1 : Equilibrium Mole Fraction [N,nsp]  
    */
    double pi, Ti, *X0i, *X1i;
    int i,result,s;

    for (i=0; i<N; i++){
        Ti = T[i];
        pi = p[i];
        X0i = X0 + i*nsp;
        X1i = X1 + i*nsp;

        result = solve_pt(pi, Ti, X0i, nsp, nel, lewis, M, a, X1i, verbose);
        if (result!=0){
            printf("Error in batch_pt @ position: %d\n", i);
            printf("pi: %f Ti: %f \n", pi, Ti);
            for (s=0; s<nsp; s++) printf("Xs[%d]=%e\n",s,X0i[s]);
            printf("Retrying with debugging on\n");
            solve_pt(pi, Ti, X0i, nsp, nel, lewis, M, a, X1i, 2);
            return result;
        }
    }
    return 0;
}

int batch_rhou(int N, double* rho,double* u,double* X0,int nsp,int nel,double* lewis,double* M,
               double* a, double* X1, double* T, int verbose){
    /*
    Compute the equilibrium compositions X1 at a number of fixed densities and internal energies 
    Inputs:
        N     : number of rho's and u's
        rho   : target Densities (kg/m3) [N]
        u     : target internal energies (J/kg) [N]
        X0    : Intial Mole fractions [N,nsp]
        nsp   : number of species 
        nel   : number of elements 
        lewis : Nasa Lewis Thermodynamic Database Data [nsp*3*9]
        M     : Molar Mass of each species (kg/mol) [nsp]
        a     : elemental composition array [nel,nsp]
        verbose: print debugging information

    Output:
        X1 : Equilibrium Mole Fractions [N,nsp]  
        T  : Equilibrium Temperatures  [N]
    */
    double rhoi, ui, *X0i, *X1i, *Ti;
    int i,result;

    for (i=0; i<N; i++){
        rhoi = rho[i];
        ui = u[i];
        X0i = X0 + i*nsp;
        X1i = X1 + i*nsp;
        Ti = T + i;

        result = solve_rhou(rhoi, ui, X0i,nsp,nel, lewis, M, a, X1i, Ti, verbose);
        if (result!=0){
            printf("Error in batch_rhou @ position: %d\n", i);
            return result;
        }
    }
    return 0;
}

int batch_u(int N, double* T, double* X, int nsp, double* lewis, double* M, double* u){
    /*
    Compute thermal equilibrium u from an array of known composition and primitives
    Inputs:
        N     : number of points in input arrays
        p     : Pressure (Pa) [N]
        T     : Temperature (K) [N]
        X     : Composition [N,nsp]
        nsp   : number of species 
        lewis : Nasa Lewis Thermodynamic Database Data [nsp*3*9]
        M     : Molar Mass of each species (kg/mol) [nsp]
        verbose: print debugging information

    Output:
        u : internal energy per unit mass [N]
    */
    int i;
    double Ti,*Xi;

    for (i=0; i<N; i++){
        Ti = T[i];
        Xi = X + i*nsp;
        u[i] = get_u(Ti, Xi, nsp, lewis, M);
    }
    return 0;
}

/*
int main(){
    printf("Called main in ceq.c!\n");
    return 0;
}
*/
