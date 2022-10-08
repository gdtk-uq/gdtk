// face.cu
// Include file for chicken.
// PJ 2022-09-11

#ifndef FACE_INCLUDED
#define FACE_INCLUDED

#include <string>
#include <sstream>
#include <vector>
#include <array>
#include <cmath>

#include "number.cu"
#include "vector3.cu"
#include "vertex.cu"
#include "flow.cu"
#include "gas.cu"

using namespace std;

// Interpolation functions that will be used in the fluc calculator.

__host__ __device__
number van_albada_limit1(number a, number b)
// A smooth slope limiter.
{
    constexpr number eps = 1.0e-12;
    number s = (a*b + fabs(a*b) + eps)/(a*a + b*b + eps);
    return s;
}

__host__ __device__
void interp_l2r2_scalar(number qL1, number qL0, number qR0, number qR1,
                        number& qL, number& qR)
// Reconstruct values, qL,qR, at the middle interface for a stencil of 4 cell-centred values.
// Assume equal cell widths.
{
    // Set up differences and limiter values.
    number delLminus = (qL0 - qL1);
    number del = (qR0 - qL0);
    number delRplus = (qR1 - qR0);
    number sL = van_albada_limit1(delLminus, del);
    number sR = van_albada_limit1(del, delRplus);
    // The actual high-order reconstruction, possibly limited.
    qL = qL0 + sL * 0.125 * (3.0*del + delLminus);
    qR = qR0 - sR * 0.125 * (delRplus + 3.0*del);
} // end of interp_l2r2_scalar()


// Flux calculations are done in the context of a face on a cell.

struct FVFace {
    Vector3 pos; // midpoint position in space
    number area;
    Vector3 n;  // unit normal
    Vector3 t1; // unit tangent 1
    Vector3 t2; // unit tangent 2
    array<number,CQI::n> F; // flux vector for conserved quantities
    // We will keep connections to the pieces composing the face
    // as indices into global arrays.
    array<int,4> vtx{0, 0, 0, 0};
    array<int,2> left_cells{0, 0};
    array<int,2> right_cells{0, 0};

    string toString() {
        ostringstream repr;
        repr << "FVFace(pos=" << pos.toString() << ", n=" << n.toString()
             << ", t1=" << t1.toString() << ", t2=" << t2.toString() << ", area=" << area;
        repr << ", vtx=["; for (auto i : vtx) repr << i << ","; repr << "]";
        repr << ", left_cells=["; for (auto i : left_cells) repr << i << ","; repr << "]";
        repr << ", right_cells=["; for (auto i : right_cells) repr << i << ","; repr << "]";
        repr << ")";
        return repr.str();
    }

    // Specific flux calculators here...

    __host__ __device__
    void ausmdv(const FlowState& fsL, const FlowState& fsR)
    // Compute the face's flux vector from left and right flow states.
    // Wada and Liou's flux calculator, implemented from details in their AIAA paper,
    // with hints from Ian Johnston.
    // Y. Wada and M. -S. Liou (1994)
    // A flux splitting scheme with high-resolution and robustness for discontinuities.
    // AIAA-94-0083.
    {
        Vector3 velL = Vector3(fsL.vel);
        Vector3 velR = Vector3(fsR.vel);
        velL.transform_to_local_frame(n, t1, t2);
        velR.transform_to_local_frame(n, t1, t2);
        //
        number rhoL = fsL.gas.rho;
        number pL = fsL.gas.p;
        number pLrL = pL/rhoL;
        number velxL = velL.x;
        number velyL = velL.y;
        number velzL = velL.z;
        number uL = fsL.gas.e;
        number aL = fsL.gas.a;
        number keL = 0.5*(velxL*velxL + velyL*velyL + velzL*velzL);
        number HL = uL + pLrL + keL;
        //
        number rhoR = fsR.gas.rho;
        number pR = fsR.gas.p;
        number pRrR = pR/rhoR;
        number velxR = velR.x;
        number velyR = velR.y;
        number velzR = velR.z;
        number uR = fsR.gas.e;
        number aR = fsR.gas.a;
        number keR = 0.5*(velxR*velxR + velyR*velyR + velzR*velzR);
        number HR = uR + pR/rhoR + keR;
        //
        // This is the main part of the flux calculator.
        //
        // Weighting parameters (eqn 32) for velocity splitting.
        number alphaL = 2.0*pLrL/(pLrL+pRrR);
        number alphaR = 2.0*pRrR/(pLrL+pRrR);
        // Common sound speed (eqn 33) and Mach numbers.
        number am = fmax(aL, aR);
        number ML = velxL/am;
        number MR = velxR/am;
        // Left state:
        // pressure splitting (eqn 34)
        // and velocity splitting (eqn 30)
        number pLplus, velxLplus;
        number dvelxL = 0.5 * (velxL + fabs(velxL));
        if (fabs(ML) <= 1.0) {
            pLplus = pL*(ML+1.0)*(ML+1.0)*(2.0-ML)*0.25;
            velxLplus = alphaL*((velxL+am)*(velxL+am)/(4.0*am) - dvelxL) + dvelxL;
        } else {
            pLplus = pL * dvelxL / velxL;
            velxLplus = dvelxL;
        }
        // Right state:
        // pressure splitting (eqn 34)
        // and velocity splitting (eqn 31)
        number pRminus, velxRminus;
        number dvelxR = 0.5*(velxR-fabs(velxR));
        if (fabs(MR) <= 1.0) {
            pRminus = pR*(MR-1.0)*(MR-1.0)*(2.0+MR)*0.25;
            velxRminus = alphaR*(-(velxR-am)*(velxR-am)/(4.0*am) - dvelxR) + dvelxR;
        } else {
            pRminus = pR * dvelxR / velxR;
            velxRminus = dvelxR;
        }
        // The mass flux. (eqn 29)
        number massL = velxLplus*rhoL;
        number massR = velxRminus*rhoR;
        number mass_half = massL+massR;
        // Pressure flux (eqn 34)
        number p_half = pLplus + pRminus;
        // Momentum flux: normal direction
        // Compute blending parameter s (eqn 37),
        // the momentum flux for AUSMV (eqn 21) and AUSMD (eqn 21)
        // and blend (eqn 36).
        number dp = pL - pR;
        constexpr number K_SWITCH = 10.0;
        dp = K_SWITCH * fabs(dp) / fmin(pL, pR);
        number s = 0.5 * fmin(1.0, dp);
        number rvel2_AUSMV = massL*velxL + massR*velxR;
        number rvel2_AUSMD = 0.5*(mass_half*(velxL+velxR) - fabs(mass_half)*(velxR-velxL));
        number rvel2_half = (0.5+s)*rvel2_AUSMV + (0.5-s)*rvel2_AUSMD;
        // Assemble components of the flux vector (eqn 36).
        F[CQI::mass] = mass_half;
        number vely = (mass_half >= 0.0) ? velyL : velyR;
        number velz = (mass_half >= 0.0) ? velzL : velzR;
        Vector3 momentum{rvel2_half+p_half, mass_half*vely, mass_half*velz};
        momentum.transform_to_global_frame(n, t1, t2);
        F[CQI::xMom] = momentum.x;
        F[CQI::yMom] = momentum.y;
        F[CQI::zMom] = momentum.z;
        number H = (mass_half >= 0.0) ? HL : HR;
        F[CQI::totEnergy] = mass_half*H;
        // When we introduce species, get the species flux lines from the Dlang code.
        // [TODO] PJ 2022-09-11
        return;
    } // end ausmdv()

    // And one generic flux calculation function.

    __host__ __device__
    void calculate_flux(FlowState& fsL1, FlowState& fsL0, FlowState& fsR0, FlowState& fsR1, int x_order)
    // Generic fluc calculation function.
    {
        // First-order reconstruction is just a copy from the nearest cell centre.
        FlowState fsL{fsL0};
        FlowState fsR{fsR0};
        if (x_order > 1) {
            // We will interpolate only some GasState properties...
            interp_l2r2_scalar(fsL1.gas.rho, fsL0.gas.rho, fsR0.gas.rho, fsR1.gas.rho, fsL.gas.rho, fsR.gas.rho);
            interp_l2r2_scalar(fsL1.gas.e, fsL0.gas.e, fsR0.gas.e, fsR1.gas.e, fsL.gas.e, fsR.gas.e);
            // and make the rest consistent.
            fsL.gas.update_from_rhoe();
            fsR.gas.update_from_rhoe();
            // Velocity components.
            interp_l2r2_scalar(fsL1.vel.x, fsL0.vel.x, fsR0.vel.x, fsR1.vel.x, fsL.vel.x, fsR.vel.x);
            interp_l2r2_scalar(fsL1.vel.y, fsL0.vel.y, fsR0.vel.y, fsR1.vel.y, fsL.vel.y, fsR.vel.y);
            interp_l2r2_scalar(fsL1.vel.z, fsL0.vel.z, fsR0.vel.z, fsR1.vel.z, fsL.vel.z, fsR.vel.z);
        }
        // Use the reconstructed values near the face in a simple flux calculator.
        ausmdv(fsL, fsR);
        // [TODO] Allow other flux calculators.
    }

}; // end FVFace

#endif
