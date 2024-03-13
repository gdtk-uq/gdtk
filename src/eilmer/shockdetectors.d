// shockdetector.d
// Prototyping for modular shock detection
// @author: Nick Gibbons (n.gibbons@uq.edu.au)

module shockdetectors;

import std.conv;
import std.math;
import std.stdio;

import globalconfig;
import gas;
import geom;
import flowstate;
import fvinterface;
import fvcell;
import nm.number;
import nm.complex;
import std.algorithm;

// Consider passing in lower level things to make testing easier
@nogc number PJ_ShockDetector(const FVInterface iface, const number comp_tol, const number shear_tol){
        // Change in normalised velocity to indicate a shock.
        // A value of -0.05 has been found suitable to detect the levels of
        // shock compression observed in the "sod" and "cone20" test cases.
        // It may need to be tuned for other situations, especially when
        // viscous effects are important.
        number comp, shear;
        auto cL = iface.left_cell;
        auto cR = iface.right_cell;

        if (cL && cR) {
            // We have two cells interacting.
            // Compare the relative gas velocities.
            number uL = geom.dot(cL.fs.vel, iface.n);
            number uR = geom.dot(cR.fs.vel, iface.n);
            number aL = cL.fs.gas.a;
            number aR = cR.fs.gas.a;
            number a_min = (aL < aR) ? aL : aR;
            comp = ((uR - uL)/a_min);

            number vL = geom.dot(cL.fs.vel, iface.t1);
            number vR = geom.dot(cR.fs.vel, iface.t1);
            number wL = geom.dot(cL.fs.vel, iface.t2);
            number wR = geom.dot(cR.fs.vel, iface.t2);

            number sound_speed = 0.5 * (aL + aR);
            number shear_y = fabs(vL - vR) / sound_speed;
            number shear_z = fabs(wL - wR) / sound_speed;
            shear = fmax(shear_y, shear_z);

        } else if (cL) {
            // We have left-cell with a wall on the right.
            // Use the gas velocity relative to the wall.
            Vector3 vel; vel.set(cL.fs.vel);
            vel.add(iface.gvel, to!number(-1.0));
            number uL = geom.dot(vel, iface.n);
            number aL = cL.fs.gas.a;
            comp = ((-uL)/aL);


            number vL = geom.dot(cL.fs.vel, iface.t1);
            number wL = geom.dot(cL.fs.vel, iface.t2);

            number shear_y = fabs(vL) / aL;
            number shear_z = fabs(wL) / aL;
            shear = fmax(shear_y, shear_z);

        } else if (cR) {
            // We have a right-cell with a wall on the left.
            // Use the gas velocity relative to the wall.
            Vector3 vel; vel.set(cR.fs.vel);
            vel.add(iface.gvel, to!number(-1.0));
            number uR = geom.dot(vel, iface.n);
            number aR = cR.fs.gas.a;
            comp = (uR/aR);

            number vR = geom.dot(cR.fs.vel, iface.t1);
            number wR = geom.dot(cR.fs.vel, iface.t2);

            number shear_y = fabs(vR) / aR;
            number shear_z = fabs(wR) / aR;
            shear = fmax(shear_y, shear_z);

        } else {
            return to!number(0.0);
        }

        if ((shear < shear_tol) && (comp < comp_tol)) {
            return to!number(1.0);
        }

        return to!number(0.0);
        

        // // map the compression factor to the allowed limits
        // comp = (comp < comp_tol) ? comp_tol : comp;
        // comp = (comp > to!number(0.0)) ? to!number(0.0) : comp;
        // comp /= comp_tol;

        // // clip the shear factor to the allowed limits
        // shear = 1.0 - shear/shear_tol;
        // shear = (shear < to!number(0.0)) ? to!number(0.0) : shear;

        // return comp*shear;
}

@nogc number NNG_ShockDetector(GasModel gm, in FlowState L, in FlowState R, Vector3 n, const double Mx){
    /*
        Experimental smooth wshock detector by NNG. Use the HLLC wave speeds to 
        look for sharp interfaces in pressure and velocity.
        See Notes: 6th June, 2023
        @author: Nick Gibbons
    */
    number pL = L.gas.p;
    number pR = R.gas.p;
    number uL = geom.dot(L.vel, n);
    number uR = geom.dot(R.vel, n);
    number yL = gm.gamma(L.gas);
    number yR = gm.gamma(R.gas);
    number a = 0.5*(L.gas.a + R.gas.a);
    number rho = 0.5*(L.gas.rho + R.gas.rho);

    number pstar = 0.5*(pL+pR) - 0.5*(uR-uL)*rho*a;
    number ustar = 0.5*(uL+uR) - 0.5*(pR-pL)/rho/a;
    number facL = sqrt(1.0 + (yL+1.0)/(2*yL)*(pstar/pL - 1.0));
    number facR = sqrt(1.0 + (yR+1.0)/(2*yR)*(pstar/pR - 1.0));

    number qL = (pstar>pL) ? facL : to!number(1.0);
    number qR = (pstar>pR) ? facR : to!number(1.0);
    //number q = fmax(qR, qL) - 1.0 + 1e-1;
    number q = fmax(qR, qL) - 1.0;
    number S = fmin(q/Mx, 1.0);

    return S;
}

@nogc number NNG_DiscoDectector(GasModel gm, in FlowState L, in FlowState R) {
    /*
        Experimental discontinuity detector by NNG, based on formulas derived 29/01/24, see back of diary.
        @author: Nick Gibbons
    */

    number rho1 = L.gas.rho;
    number rho2 = R.gas.rho;
    number p1 = L.gas.p;
    number p2 = R.gas.p;
    number T1 = L.gas.T;
    number T2 = R.gas.T;

    number R_L = gm.gas_constant(L.gas);
    number R_R = gm.gas_constant(R.gas);
    number R_av = 0.5*(R_L+R_R); // Check in multispecies??

    number drho = fabs(rho2 - rho1);
    number dp   = fabs(p2-p1);
    number dp_rho = dp/R_av/T1;

    immutable double factor = 0.20; // TBD: This seems to be flow dependant. Maybe make it a config 
    number argument = (drho - dp_rho)/factor;
    number S = 0.5*tanh(argument - 4.0) + 0.5;

    return S;
}

@nogc number NNG_DiscoDectector1(GasModel gm, in FlowState LL, in FlowState L, in FlowState R, in FlowState RR) {
    /*
        Experimental discontinuity detector by NNG, based on formulas derived 21/02/24, see purple book.
        @author: Nick Gibbons
    */

    number rho0 =LL.gas.rho;
    number rho1 = L.gas.rho;
    number rho2 = R.gas.rho;
    number rho3 =RR.gas.rho;

    number p0 =LL.gas.p;
    number p1 = L.gas.p;
    number p2 = R.gas.p;
    number p3 =RR.gas.p;

    number y_L = gm.gamma(LL.gas);
    number y_R = gm.gamma(RR.gas);

    number pgrad = (-1.5*p0 + -0.5*p1 + 0.5*p2 + 1.5*p3)/5;
    number pavg  = (p0 + p1 + p2 + p3)/4;

    // Use the left rho and T to compute a an estimated right state
    number pright = pavg + pgrad*1.5;
    number rhoright = rho0*pow(pright/p0, 1.0/y_L);

    // Now do the same thing from the other side
    number pleft = pavg - pgrad*1.5;
    number rholeft = rho3*pow(pleft/p3, 1.0/y_R);

    number drho_p1= fabs(rholeft-rho0)/rho0;
    number drho_p2= fabs(rhoright-rho3)/rho3;
    number drho_p = drho_p1 + drho_p2;

    number width = 0.5;
    number shift = 2.0;
    number S = 0.5*tanh(3.0*(drho_p - shift)/width) + 0.5;

    return S;
}

@nogc number NNG_DiscoDectector2(GasModel gm, in FlowState LL, in FlowState L, in FlowState R, in FlowState RR) {
    /*
        Experimental discontinuity detector by NNG, based on formulas derived 21/02/24, see purple book.
        @author: Nick Gibbons
    */

    number T0 =LL.gas.T;
    number T1 = L.gas.T;
    number T2 = R.gas.T;
    number T3 =RR.gas.T;

    number p0 =LL.gas.p;
    number p1 = L.gas.p;
    number p2 = R.gas.p;
    number p3 =RR.gas.p;

    number y_L = gm.gamma(LL.gas);
    number y_R = gm.gamma(RR.gas);

    number pgrad = (-1.5*p0 + -0.5*p1 + 0.5*p2 + 1.5*p3)/5;
    number pavg  = (p0 + p1 + p2 + p3)/4;

    // Use the left rho and T to compute a an estimated right state
    number pright = pavg + pgrad*1.5;
    number Tright = T0*pow(pright/p0, (y_L-1.0)/y_L);

    // Now do the same thing from the other side
    number pleft = pavg - pgrad*1.5;
    number Tleft = T3*pow(pleft/p3, (y_R-1.0)/y_R);

    number dT_p1= fabs(Tleft-T0)/T0;
    number dT_p2= fabs(Tright-T3)/T3;
    number dT_p = dT_p1 + dT_p2;

    number width = 0.1;
    number shift = 0.5;
    number S = 0.5*tanh(3.0*(dT_p - shift)/width) + 0.5;
    S = fmax(S, 0.1);

    return S;
}
