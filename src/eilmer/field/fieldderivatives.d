/*
    fieldderivatives.d
    Machinery for approximating high order derivatives in 2D. (Strap yourself in because yeah...)

    Author: Nick Gibbons
    Version: 2022-08-01: Machine generated code from Maxima
*/

import std.math;
import std.stdio;

@nogc pure
double denominator(double dxN, double dxS, double dxE, double dxW,
                   double dyN, double dyS, double dyE, double dyW){
    /*
        All of the finite difference formulas share the same denominator, a multi-line
        abomination that is computed here.
    */
    double D = (((dxE*dxS-dxE*dxN)*dyN+(dxE*dxN-dxN*dxS)*dyE)*dyS+(dxN-dxE)*dxS*dyE*dyN)
        *dyW^^2
        +(((dxE*dxN-dxE*dxW)*dyN+(dxN*dxW-dxE*dxN)*dyE)*dyS^^2
         +((dxE*dxW-dxE*dxS)*dyN^^2+(dxN*dxS-dxN*dxW)*dyE^^2+(dxE^^2*dxN-dxE*dxN^^2)*dxW
                                  +(dxE*dxN^^2-dxE^^2*dxN)*dxS)
          *dyS+(dxE*dxS-dxS*dxW)*dyE*dyN^^2
         +((dxS*dxW-dxN*dxS)*dyE^^2+(dxE*dxS^^2-dxE^^2*dxS)*dxW-dxE*dxN*dxS^^2
                                  +dxE^^2*dxN*dxS)
          *dyN+((dxN^^2*dxS-dxN*dxS^^2)*dxW+dxE*dxN*dxS^^2-dxE*dxN^^2*dxS)*dyE)
         *dyW+(dxE-dxN)*dxW*dyE*dyN*dyS^^2
        +((dxS-dxE)*dxW*dyE*dyN^^2+((dxN-dxS)*dxW*dyE^^2+(dxE*dxN-dxE*dxS)*dxW^^2
                                                      +(dxE^^2*dxS-dxE^^2*dxN)*dxW)
                                  *dyN
                                 +((dxN*dxS-dxE*dxN)*dxW^^2+(dxE*dxN^^2-dxN^^2*dxS)*dxW)
                                  *dyE)
         *dyS+((dxE-dxN)*dxS*dxW^^2+(dxN-dxE)*dxS^^2*dxW)*dyE*dyN;
    return D;
}

// Subcomponents needed for the x derivative

const string Nx = "-((dxS-dxE)*dyE*dyS*dyW^^2+((dxE-dxW)*dyE*dyS^^2
                          +((dxW-dxS)*dyE^^2-dxE^^2*dxW+dxE^^2*dxS)*dyS
                          +(dxS^^2*dxW-dxE*dxS^^2)*dyE)
                          *dyW+(dxE-dxS)*dxW^^2*dyE*dyS)";
 
const string Sx = "((dxN-dxE)*dyE*dyN*dyW^^2+((dxE-dxW)*dyE*dyN^^2
                         +((dxW-dxN)*dyE^^2-dxE^^2*dxW+dxE^^2*dxN)*dyN
                         +(dxN^^2*dxW-dxE*dxN^^2)*dyE)
                         *dyW+(dxE-dxN)*dxW^^2*dyE*dyN)";

const string Ex = "((dxS-dxN)*dyN*dyS*dyW^^2+((dxN-dxW)*dyN*dyS^^2
                         +((dxW-dxS)*dyN^^2-dxN^^2*dxW+dxN^^2*dxS)*dyS
                         +(dxS^^2*dxW-dxN*dxS^^2)*dyN)
                         *dyW+(dxN-dxS)*dxW^^2*dyN*dyS)";

const string Wx = "-((dxN-dxE)*dyE*dyN*dyS^^2+((dxE-dxS)*dyE*dyN^^2
                          +((dxS-dxN)*dyE^^2-dxE^^2*dxS+dxE^^2*dxN)*dyN
                          +(dxN^^2*dxS-dxE*dxN^^2)*dyE)
                          *dyS+(dxE-dxN)*dxS^^2*dyE*dyN)";

const string Ix = "-((((dxS-dxN)*dyN+(dxE-dxS)*dyE)*dyS+(dxN-dxE)*dyE*dyN)*dyW^^2
 +(((dxN-dxW)*dyN+(dxW-dxE)*dyE)*dyS^^2+((dxW-dxS)*dyN^^2+(dxS-dxW)*dyE^^2
                                                       +(dxE^^2-dxN^^2)*dxW
                                                       +(dxN^^2-dxE^^2)*dxS)
                                       *dyS+(dxE-dxW)*dyE*dyN^^2
                                      +((dxW-dxN)*dyE^^2+(dxS^^2-dxE^^2)*dxW
                                                       -dxN*dxS^^2+dxE^^2*dxN)
                                       *dyN
                                      +((dxN^^2-dxS^^2)*dxW+dxE*dxS^^2-dxE*dxN^^2)
                                       *dyE)
  *dyW+(dxE-dxN)*dyE*dyN*dyS^^2
 +((dxS-dxE)*dyE*dyN^^2+((dxN-dxS)*dyE^^2+(dxN-dxS)*dxW^^2+dxE^^2*dxS-dxE^^2*dxN)
                       *dyN+((dxS-dxE)*dxW^^2-dxN^^2*dxS+dxE*dxN^^2)*dyE)
  *dyS+((dxE-dxN)*dxW^^2+(dxN-dxE)*dxS^^2)*dyE*dyN)";

@nogc pure 
double dfdx_2d(double xi, double xN, double xS, double xE, double xW,
               double yi, double yN, double yS, double yE, double yW,
               double Pi, double PN, double PS, double PE, double PW){
    /*
        Compute dfdx using a sequence of five points in space, and the taylor series
        expansion dP = dpdx*dx + dpdy*dy + dp2dx2(dx^2 - dy^2)/2 + dp2dxdy*dx*dy.
         - We can get away with this because we assume dp2dx2 == dp2dy2 and 
           dp2dxdy = dp2dydx, which seems reasonable.
        
        Computed by Maxima file fd.max
    */

    double dxN = fmax(1e-12, xN - xi); double dyN = yN - yi;
    double dxS = xS - xi;              double dyS = yS - yi;
    double dxE = xE - xi;              double dyE = fmax(1e-12, yE - yi);
    double dxW = xW - xi;              double dyW = yW - yi;

    double D = denominator(dxN, dxS, dxE, dxW, dyN, dyS, dyE, dyW);

    double _Nx = mixin(Nx);
    double _Sx = mixin(Sx);
    double _Ex = mixin(Ex);
    double _Wx = mixin(Wx);
    double _Ix = mixin(Ix);

    return (_Ix*Pi + _Nx*PN + _Sx*PS + _Ex*PE + _Wx*PW)/D;
}

const string Ny = "((dxE*dxS*dyS-dxE*dxS*dyE)*dyW^^2+((-dxE*dxW*dyS^^2)+dxS*dxW*dyE^^2
                                                  +(dxE*dxS^^2-dxE^^2*dxS)*dxW)
                                 *dyW+dxE*dxW*dyE*dyS^^2
                                +((-dxS*dxW*dyE^^2)-dxE*dxS*dxW^^2
                                                  +dxE^^2*dxS*dxW)
                                 *dyS+(dxE*dxS*dxW^^2-dxE*dxS^^2*dxW)*dyE)";

const string Sy = "-((dxE*dxN*dyN-dxE*dxN*dyE)*dyW^^2+((-dxE*dxW*dyN^^2)+dxN*dxW*dyE^^2
                                                   +(dxE*dxN^^2-dxE^^2*dxN)*dxW)
                                  *dyW+dxE*dxW*dyE*dyN^^2
                                 +((-dxN*dxW*dyE^^2)-dxE*dxN*dxW^^2
                                                   +dxE^^2*dxN*dxW)
                                  *dyN+(dxE*dxN*dxW^^2-dxE*dxN^^2*dxW)*dyE)";

const string Ey = "-((dxN*dxS*dyS-dxN*dxS*dyN)*dyW^^2+((-dxN*dxW*dyS^^2)+dxS*dxW*dyN^^2
                                                   +(dxN*dxS^^2-dxN^^2*dxS)*dxW)
                                  *dyW+dxN*dxW*dyN*dyS^^2
                                 +((-dxS*dxW*dyN^^2)-dxN*dxS*dxW^^2
                                                   +dxN^^2*dxS*dxW)
                                  *dyS+(dxN*dxS*dxW^^2-dxN*dxS^^2*dxW)*dyN)";

const string Wy = "((dxE*dxN*dyN-dxE*dxN*dyE)*dyS^^2+((-dxE*dxS*dyN^^2)+dxN*dxS*dyE^^2
                                                  +(dxE*dxN^^2-dxE^^2*dxN)*dxS)
                                 *dyS+dxE*dxS*dyE*dyN^^2
                                +((-dxN*dxS*dyE^^2)-dxE*dxN*dxS^^2
                                                  +dxE^^2*dxN*dxS)
                                 *dyN+(dxE*dxN*dxS^^2-dxE*dxN^^2*dxS)*dyE)";

const string Iy = "(((dxN-dxE)*dxS*dyS+(dxE*dxN-dxN*dxS)*dyN+(dxE*dxS-dxE*dxN)*dyE)*dyW^^2
 +((dxE-dxN)*dxW*dyS^^2+(dxS-dxE)*dxW*dyN^^2+(dxN-dxS)*dxW*dyE^^2
                      +((dxN-dxE)*dxS^^2+(dxE^^2-dxN^^2)*dxS+dxE*dxN^^2-dxE^^2*dxN)
                       *dxW)
  *dyW+((dxN*dxW-dxE*dxN)*dyN+(dxE*dxN-dxE*dxW)*dyE)*dyS^^2
 +((dxE*dxS-dxS*dxW)*dyN^^2+(dxS*dxW-dxN*dxS)*dyE^^2+(dxE-dxN)*dxS*dxW^^2
                          +(dxN^^2-dxE^^2)*dxS*dxW+(dxE^^2*dxN-dxE*dxN^^2)*dxS)
  *dyS+(dxE*dxW-dxE*dxS)*dyE*dyN^^2
 +((dxN*dxS-dxN*dxW)*dyE^^2+(dxN*dxS-dxE*dxN)*dxW^^2+(dxE^^2*dxN-dxN*dxS^^2)*dxW
                          +dxE*dxN*dxS^^2-dxE^^2*dxN*dxS)
  *dyN
 +((dxE*dxN-dxE*dxS)*dxW^^2+(dxE*dxS^^2-dxE*dxN^^2)*dxW-dxE*dxN*dxS^^2
                          +dxE*dxN^^2*dxS)
  *dyE)";

@nogc pure
double dfdy_2d(double xi, double xN, double xS, double xE, double xW,
               double yi, double yN, double yS, double yE, double yW,
               double Pi, double PN, double PS, double PE, double PW){

    double dxN = xN - xi; double dyN = yN - yi;
    double dxS = xS - xi; double dyS = yS - yi;
    double dxE = xE - xi; double dyE = yE - yi;
    double dxW = xW - xi; double dyW = yW - yi;

    double D = denominator(dxN, dxS, dxE, dxW, dyN, dyS, dyE, dyW);

    double _Ny = mixin(Ny);
    double _Sy = mixin(Sy);
    double _Ey = mixin(Ey);
    double _Wy = mixin(Wy);
    double _Iy = mixin(Iy);

    return (_Iy*Pi + _Ny*PN + _Sy*PS + _Ey*PE + _Wy*PW)/D;
}

const string Nxx = "-(((2*dxE*dxW-2*dxE*dxS)*dyS+(2*dxE*dxS-2*dxS*dxW)*dyE)*dyW
 +(2*dxS-2*dxE)*dxW*dyE*dyS)";

const string Sxx = "(((2*dxE*dxW-2*dxE*dxN)*dyN+(2*dxE*dxN-2*dxN*dxW)*dyE)*dyW
 +(2*dxN-2*dxE)*dxW*dyE*dyN)";

const string Exx = "(((2*dxN*dxW-2*dxN*dxS)*dyS+(2*dxN*dxS-2*dxS*dxW)*dyN)*dyW
 +(2*dxS-2*dxN)*dxW*dyN*dyS)";

const string Wxx = "-(((2*dxE*dxS-2*dxE*dxN)*dyN+(2*dxE*dxN-2*dxN*dxS)*dyE)*dyS
 +(2*dxN-2*dxE)*dxS*dyE*dyN)";

const string Ixx = "-((((2*dxN-2*dxE)*dxW+(2*dxE-2*dxN)*dxS)*dyS
 +((2*dxE-2*dxS)*dxW+2*dxN*dxS-2*dxE*dxN)*dyN
 +((2*dxS-2*dxN)*dxW-2*dxE*dxS+2*dxE*dxN)*dyE)
 *dyW
 +(((2*dxS-2*dxN)*dxW-2*dxE*dxS+2*dxE*dxN)*dyN
  +((2*dxE-2*dxS)*dxW+2*dxN*dxS-2*dxE*dxN)*dyE)
  *dyS+((2*dxN-2*dxE)*dxW+(2*dxE-2*dxN)*dxS)*dyE*dyN)";

double d2fdx2_2d(double xi, double xN, double xS, double xE, double xW,
                 double yi, double yN, double yS, double yE, double yW,
                 double Pi, double PN, double PS, double PE, double PW){

    double dxN = xN - xi; double dyN = yN - yi;
    double dxS = xS - xi; double dyS = yS - yi;
    double dxE = xE - xi; double dyE = yE - yi;
    double dxW = xW - xi; double dyW = yW - yi;

    double D = denominator(dxN, dxS, dxE, dxW, dyN, dyS, dyE, dyW);

    double _Nxx = mixin(Nxx);
    double _Sxx = mixin(Sxx);
    double _Exx = mixin(Exx);
    double _Wxx = mixin(Wxx);
    double _Ixx = mixin(Ixx);

    return (_Ixx*Pi + _Nxx*PN + _Sxx*PS + _Exx*PE + _Wxx*PW)/D;
}

const string Nxy = "-((dxE*dyS-dxS*dyE)*dyW^^2+((-dxE*dyS^^2)+dxS*dyE^^2+dxE*dxS^^2-dxE^^2*dxS)*dyW
                         +dxW*dyE*dyS^^2+((-dxW*dyE^^2)-dxE*dxW^^2+dxE^^2*dxW)*dyS
                         +(dxS*dxW^^2-dxS^^2*dxW)*dyE)";

const string Sxy = "((dxE*dyN-dxN*dyE)*dyW^^2+((-dxE*dyN^^2)+dxN*dyE^^2+dxE*dxN^^2-dxE^^2*dxN)*dyW
                        +dxW*dyE*dyN^^2+((-dxW*dyE^^2)-dxE*dxW^^2+dxE^^2*dxW)*dyN
                        +(dxN*dxW^^2-dxN^^2*dxW)*dyE)";

const string Exy = "((dxN*dyS-dxS*dyN)*dyW^^2+((-dxN*dyS^^2)+dxS*dyN^^2+dxN*dxS^^2-dxN^^2*dxS)*dyW
                        +dxW*dyN*dyS^^2+((-dxW*dyN^^2)-dxN*dxW^^2+dxN^^2*dxW)*dyS
                        +(dxS*dxW^^2-dxS^^2*dxW)*dyN)";

const string Wxy = "-((dxE*dyN-dxN*dyE)*dyS^^2+((-dxE*dyN^^2)+dxN*dyE^^2+dxE*dxN^^2-dxE^^2*dxN)*dyS
                         +dxS*dyE*dyN^^2+((-dxS*dyE^^2)-dxE*dxS^^2+dxE^^2*dxS)*dyN
                         +(dxN*dxS^^2-dxN^^2*dxS)*dyE)";

const string Ixy = "-(((dxN-dxE)*dyS+(dxE-dxS)*dyN+(dxS-dxN)*dyE)*dyW^^2
 +((dxE-dxN)*dyS^^2+(dxS-dxE)*dyN^^2+(dxN-dxS)*dyE^^2+(dxN-dxE)*dxS^^2
                  +(dxE^^2-dxN^^2)*dxS+dxE*dxN^^2-dxE^^2*dxN)
  *dyW+((dxW-dxE)*dyN+(dxN-dxW)*dyE)*dyS^^2
 +((dxE-dxW)*dyN^^2+(dxW-dxN)*dyE^^2+(dxE-dxN)*dxW^^2+(dxN^^2-dxE^^2)*dxW-dxE*dxN^^2
                  +dxE^^2*dxN)
  *dyS+(dxW-dxS)*dyE*dyN^^2
 +((dxS-dxW)*dyE^^2+(dxS-dxE)*dxW^^2+(dxE^^2-dxS^^2)*dxW+dxE*dxS^^2-dxE^^2*dxS)*dyN
 +((dxN-dxS)*dxW^^2+(dxS^^2-dxN^^2)*dxW-dxN*dxS^^2+dxN^^2*dxS)*dyE)";

@nogc pure
double d2fdxdy_2d(double xi, double xN, double xS, double xE, double xW,
                  double yi, double yN, double yS, double yE, double yW,
                  double Pi, double PN, double PS, double PE, double PW){

    double dxN = xN - xi; double dyN = yN - yi;
    double dxS = xS - xi; double dyS = yS - yi;
    double dxE = xE - xi; double dyE = yE - yi;
    double dxW = xW - xi; double dyW = yW - yi;

    double D = denominator(dxN, dxS, dxE, dxW, dyN, dyS, dyE, dyW);

    double _Nxy = mixin(Nxy);
    double _Sxy = mixin(Sxy);
    double _Exy = mixin(Exy);
    double _Wxy = mixin(Wxy);
    double _Ixy = mixin(Ixy);

    return (_Ixy*Pi + _Nxy*PN + _Sxy*PS + _Exy*PE + _Wxy*PW)/D;
}

/*
int main() {
    // Little bit of test code copied from Python
	double a = 1.1;
	double c = 0.8;
	double d = 1.5;
	double e = -1.4;
	double f = 3.0;
 
    double Phi(double x, double y) { return a*x^^2 - a*y^^2 + c*x*y + d*x + e*y + f; }
    double dPhidx(double x, double y) { return 2*a*x + c*y + d; }
    double dPhidy(double x, double y) { return -2*a*y + c*y + e; }
    double d2Phidx2(double x, double y) { return 2*a; }
    double d2Phidxdy(double x, double y) { return c; }

	double xi = 1.0; double yi = 1.0;
	double Pi = Phi(xi,yi);

    double xN = 1.0001; double yN = 1.01;
    double xS = 0.999;  double yS = 0.99;
	double xE = 1.01;   double yE = 1.0;
	double xW = 0.99;   double yW = 1.0;

	double xp = 1.005;  double yp = 1.005;

	double PN = Phi(xN, yN);
	double PS = Phi(xS, yS);
	double PE = Phi(xE, yE);
	double PW = Phi(xW, yW);

    double dpdx_fd = dfdx_2d(xi, xN, xS, xE, xW, yi, yN, yS, yE, yW, Pi, PN, PS, PE, PW);
    double dpdy_fd = dfdy_2d(xi, xN, xS, xE, xW, yi, yN, yS, yE, yW, Pi, PN, PS, PE, PW);
    double d2pdx2_fd = d2fdx2_2d(xi, xN, xS, xE, xW, yi, yN, yS, yE, yW, Pi, PN, PS, PE, PW);
    double d2pdxdy_fd = d2fdxdy_2d(xi, xN, xS, xE, xW, yi, yN, yS, yE, yW, Pi, PN, PS, PE, PW);

    double dpdx = dPhidx(xi, yi);
    double dpdy = dPhidy(xi, yi);
	double d2pdx2 = d2Phidx2(xi, yi);
    double d2pdxdy = d2Phidxdy(xi, yi);

    writeln("dpdx_fd: ", dpdx_fd, " should be: ", dpdx);
    writeln("dpdy_fd: ", dpdy_fd, " should be: ", dpdy);
    writeln("d2pdx2_fd: ", d2pdx2_fd, " should be: ", d2pdx2);
    writeln("d2pdxdy_fd: ", d2pdxdy_fd, " should be: ", d2pdxdy);

    return 0;
}
*/
