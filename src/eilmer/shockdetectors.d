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

// Consider passing in lower level things to make testing easier
@nogc int PJ_ShockDetector(const FVInterface iface, const double tol){
        // Change in normalised velocity to indicate a shock.
        // A value of -0.05 has been found suitable to detect the levels of
        // shock compression observed in the "sod" and "cone20" test cases.
        // It may need to be tuned for other situations, especially when
        // viscous effects are important.
        int S;
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
            S = ((uR - uL)/a_min) < tol;
        } else if (cL) {
            // We have left-cell with a wall on the right.
            // Use the gas velocity relative to the wall.
            Vector3 vel; vel.set(cL.fs.vel);
            vel.add(iface.gvel, to!number(-1.0));
            number uL = geom.dot(vel, iface.n);
            number aL = cL.fs.gas.a;
            S = ((-uL)/aL) < tol;
        } else if (cR) {
            // We have a right-cell with a wall on the left.
            // Use the gas velocity relative to the wall.
            Vector3 vel; vel.set(cR.fs.vel);
            vel.add(iface.gvel, to!number(-1.0));
            number uR = geom.dot(vel, iface.n);
            number aR = cR.fs.gas.a;
            S = (uR/aR) < tol;
        } else {
            S = 0;
        }
        return S;
}

