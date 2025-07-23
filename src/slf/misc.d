// misc.d - miscelleneous functions and support code for slf
// @author: NNG


module misc;

import std.math;

import io;

import nm.bbla;
import ntypes.complex;
import nm.number;


/* This function borrowed from zero_rk:
   Approximate inverse error function implemented from:
   "A handy approximation for the error function and its inverse"
   by Sergei Winitzki
*/
@nogc pure
double erfc_inv(double q) {
    if(q <= 0.0) return float.infinity;
    if(q >= 2.0) return -float.infinity;
    double p = 1 - q;
    double sgn = p < 0 ? -1 : 1;

    double logpt = log((1 - p)*(1 + p));
    double pt1 = 2/(PI*0.147) + 0.5*logpt;
    double pt2 = 1/(0.147) * logpt;

    return(sgn*sqrt(-pt1 + sqrt(pt1*pt1 - pt2)));
}


@nogc
void gaussian_initial_condition(ref const Parameters pm, number[] U){
    double sigma = 0.1;

    foreach(i; 0 .. pm.N){
        size_t idx = i*pm.neq;
        double factor = 0.5*tanh(2*6.0*pm.Z[i] - 6.0) + 0.5;
        foreach(isp; 0 .. pm.nsp) U[idx+isp] = (1.0 - factor)*pm.Y0[isp] + factor*pm.Y1[isp];
        U[idx+pm.nsp].re = 1500.0*exp(-(pm.Z[i]-0.5)*(pm.Z[i]-0.5)/2.0/sigma/sigma) + 300.0;
        version(complex_numbers) { U[idx+pm.nsp].im = 0.0; }
    }
}

@nogc pure
number evaluate_scalar_dissipation(number D, double Z){
    /*
        Using equation (30) from Pitsch and Peters, 1998
    */
    double arg = erfc_inv(2.0*Z);
    number chi = D*exp(-2.0*arg*arg);
    return chi;
}

@nogc pure
number evaluate_scalar_dissipation_2(number rho, double rho_ox, number D, double Z){
    /*
        Using equation (31) from Pitsch and Peters, 1998
    */
    double arg = erfc_inv(2.0*Z);
    number chi = D*exp(-2.0*arg*arg);
    number stau = sqrt(rho_ox/rho);
    chi *= 3.0*(stau+1.0)*(stau+1.0)/(2.0*stau+1.0)/4.0;
    return chi;
}
