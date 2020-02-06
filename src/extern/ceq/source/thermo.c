/*
C library for equilibrium chemistry calculations: Module for thermodynamic table evaluation

References:
    "NASA Glenn Coefficients for Calculating Thermodynamic Properties of Individual Species"
    NASA/TP - 2002-211556, September 2002
    Bonnie J. McBride, Michael J. Zehe, and Sanford Gordon

Notes: This module contains routines for evaluating standard state thermodynamic variables
       using the NASA Glenn thermodynamic database (thermo.inp), previously known as the 
       NASA Lewis thermodynamic database. It assumes the 9 component curve fit coefficients
       are arranged into an array called "lewis", ordered into three segements of nine,
       one for each range 200-1000K, 1000-6000K, and 6000-20,000K.

@author: Nick Gibbons
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "thermo.h"

const double Ru=8.3144621;      // Universal Gas Constant (J/mol/K)

double compute_Cp0_R(double Tin, double* lewis){
    /*
    Molar Specific Heat at constant Pressure at 1 BAR divided by Ru for a single species (unitless)
        Tin   : Temperature (K)
        lewis : thermo data for current species [3,9] 
    */
    double* lp;
    double Cp0_R,T;
    int iT;

    T = fmax(fmin(Tin,20000.0),200.0);
    iT = T <= 1000.0 ? 0 : (T<=6000 ? 1 : 2);
    lp = lewis + iT*9;

    Cp0_R =  lp[0]/T/T
           + lp[1]/T
           + lp[2]
           + lp[3]*T
           + lp[4]*T*T
           + lp[5]*T*T*T
           + lp[6]*T*T*T*T;
      return Cp0_R;
}

double compute_H0_RT(double Tin, double* lewis){
    /*
    Molar enthalpy at 1 BAR divided by Ru*T for a single species (unitless)
        Tin   : Temperature (K)
        lewis : thermo data for current species [3,9] 
    */
    double* lp;
    double H0_RT,T;
    int iT;

    T = fmax(fmin(Tin,20000.0),200.0);
    iT = T <= 1000.0 ? 0 : (T<=6000 ? 1 : 2);
    lp = lewis + iT*9;

    H0_RT =    -1.0*lp[0]/T/T
           + log(T)*lp[1]/T
           +    1.0*lp[2]
           +    0.5*lp[3]*T
           +  1/3.0*lp[4]*T*T
           +   0.25*lp[5]*T*T*T
           +    0.2*lp[6]*T*T*T*T + lp[7]/T;
      return H0_RT;
}

double compute_S0_R(double Tin, double* lewis){
    /*
    Molar entropy at 1 BAR divided by Ru for a single species (unitless)
        Tin   : Temperature (K)
        lewis : thermo data for current species [3,9] 
    */
    double* lp;
    double S0_R,T;
    int iT;

    T = fmax(fmin(Tin,20000.0),200.0);
    iT = T <= 1000.0 ? 0 : (T<=6000 ? 1 : 2);
    lp = lewis + iT*9;

    S0_R =     -0.5*lp[0]/T/T
           +   -1.0*lp[1]/T
           + log(T)*lp[2]
           +    1.0*lp[3]*T
           +    0.5*lp[4]*T*T
           +  1/3.0*lp[5]*T*T*T
           +   0.25*lp[6]*T*T*T*T + lp[8];
    return S0_R;
}

double compute_G0_RT(double T, double* lewis){
    /*
    Compute the Molar Gibss free energy divided by Ru*T for a single  species (unitless)
        Tin   : Temperature (K)
        lewis : thermo data for current species [3,9] 
    */
    return compute_H0_RT(T, lewis) - compute_S0_R(T, lewis);
}

void test_thermo(){
   // Lewis Data for molecular hydrogen, comparing to python code lewis_thermo.py
   double l[27]={ 4.078323210e+04,-8.009186040e+02, 8.214702010e+00,-1.269714457e-02, 1.753605076e-05,
                 -1.202860270e-08, 3.368093490e-12,                  2.682484665e+03,-3.043788844e+01,
                  5.608128010e+05,-8.371504740e+02, 2.975364532e+00, 1.252249124e-03,-3.740716190e-07,
                  5.936625200e-11,-3.606994100e-15,                  5.339824410e+03,-2.202774769e+00,
                  4.966884120e+08,-3.147547149e+05, 7.984121880e+01,-8.414789210e-03, 4.753248350e-07,
                 -1.371873492e-11, 1.605461756e-16,                  2.488433516e+06,-6.695728110e+02};
   int Ntests = 7;
   double Ttests[7] = {201.0, 500.0, 1000.0, 3000.0, 6000.0, 15000.0, 19999.0};
   double H0_RTs[7] = {-1.6442077188332007, 1.4150066814925424, 2.487099719417329, 3.5572630593766785,
                        4.1746893324131973, 4.3433993080128346, 4.0915125910525205};
   double S0_Rs[7] = {14.378045733318185, 17.528404682885633, 19.991159671959334, 24.401660937194425,
                      27.699098528354956, 31.895255444424265, 32.858161650483794};
   int i;
   for (i=0; i<Ntests; i++){
       double T = Ttests[i];
       printf("%f %f %f %f %f\n",T, H0_RTs[i], compute_H0_RT(T, l), S0_Rs[i], compute_S0_R(T, l)); 
   }
    
}

#ifdef TEST
int main(){
    test_thermo();
    return 0;
}
#endif
