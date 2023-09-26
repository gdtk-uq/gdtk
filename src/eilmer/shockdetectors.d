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
    Experimental continuous shock detector by NNG. Use the HLLC wave speeds from
    Toro (1994) to look for sharp interfaces in pressure and velocity.

    "Restoration of the contact surface in the HLL-Riemann solver"
    E. R. Toro and M. Spruce and W. Speares, Shock Waves (1994)

    @author: Nick Gibbons (See Notes: 6th June, 2023)
*/
    number pL = L.gas.p;
    number pR = R.gas.p;
    number uL = geom.dot(L.vel, n);
    number uR = geom.dot(R.vel, n);
    number yL = gm.gamma(L.gas);
    number yR = gm.gamma(R.gas);
    number a = 0.5*(L.gas.a + R.gas.a);
    number rho = 0.5*(L.gas.rho + R.gas.rho);

    // Toro's "Hybrid estimates", equation 15 and 16, define a factor
    // q for the wave speeds that is one for a mach wave, but greater
    // than one for a shockwave.
    number pstar = 0.5*(pL+pR) - 0.5*(uR-uL)*rho*a;
    number ustar = 0.5*(uL+uR) - 0.5*(pR-pL)/rho/a;
    number facL = sqrt(1.0 + (yL+1.0)/(2*yL)*(pstar/pL - 1.0));
    number facR = sqrt(1.0 + (yR+1.0)/(2*yR)*(pstar/pR - 1.0));

    number qL = (pstar>pL) ? facL : to!number(1.0);
    number qR = (pstar>pR) ? facR : to!number(1.0);

    // Subtracting 1 and scaling by a tuning parameter Mx, gives the shock detector
    //number q = fmax(qR, qL) - 1.0;
    number q = 0.5*(qR+qL) - 1.0;
    number S = fmin(q/Mx, 1.0);

    return S;
}
