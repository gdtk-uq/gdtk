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
#include "gas.cu"
#include "vertex.cu"
#include "config.cu"
#include "flow.cu"

using namespace std;

namespace FluxCalc {
    // Selection of flux calculator happens at compile time,
    // setting flux_calculator below.
    array<string,2> names{"ausmdv", "sbp_asf"};
    //
    constexpr int ausmdv = 0;
    constexpr int sbp_asf = 1;
};

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


    __host__ __device__
    void sbp_asf(FlowState& fsL1, FlowState& fsL0, FlowState& fsR0, FlowState& fsR1)
    // Christine's Summation-By-Parts Alpha-Split Flux calculation function.
    {
        Vector3 velL1 = Vector3(fsL1.vel);
        Vector3 velL0 = Vector3(fsL0.vel);
        Vector3 velR0 = Vector3(fsR0.vel);
        Vector3 velR1 = Vector3(fsR1.vel);
        velL1.transform_to_local_frame(n, t1, t2);
        velL0.transform_to_local_frame(n, t1, t2);
        velR0.transform_to_local_frame(n, t1, t2);
        velR1.transform_to_local_frame(n, t1, t2);
        //
        // [TODO] some fancy calculation here.
        //
        // L1
        number rhoL1 = fsL1.gas.rho; number velxL1 = velL1.x; number velyL1 = velL1.y; number velzL1 = velL1.z; number pL1 = fsL1.gas.p; number eL1 = fsL1.gas.e;
        number rhoL1velxL1 = rhoL1*velxL1; number rhoL1velyL1 = rhoL1*velyL1; number rhoL1velzL1 = rhoL1*velzL1; number rhoL1eL1 = eL1*rhoL1;
        // L0
        number rhoL0 = fsL0.gas.rho; number velxL0 = velL0.x; number velyL0 = velL0.y; number velzL0 = velL0.z; number pL0 = fsL0.gas.p; number eL0 = fsL0.gas.e;
        number rhoL0velxL0 = rhoL0*velxL0; number rhoL0velyL0 = rhoL0*velyL0; number rhoL0velzL0 = rhoL0*velzL0; number rhoL0eL0 = eL0*rhoL0;
        // R0
        number rhoR0 = fsR0.gas.rho; number velxR0 = velR0.x; number velyR0 = velR0.y; number velzR0 = velR0.z; number pR0 = fsR0.gas.p; number eR0 = fsR0.gas.e;
        number rhoR0velxR0 = rhoR0*velxR0; number rhoR0velyR0 = rhoR0*velyR0; number rhoR0velzR0 = rhoR0*velzR0; number rhoR0eR0 = eR0*rhoR0;
        // R1
        number rhoR1 = fsR1.gas.rho; number velxR1 = velR1.x; number velyR1 = velR1.y; number velzR1 = velR1.z; number pR1 = fsR1.gas.p; number eR1 = fsR1.gas.e;
        number rhoR1velxR1 = rhoR1*velxR1; number rhoR1velyR1 = rhoR1*velyR1; number rhoR1velzR1 = rhoR1*velzR1; number rhoR1eR1 = eR1*rhoR1;
        //
        number f_c_0 = (1.0 / 12.0) * (-rhoL1*velxL1 + 7.0*rhoL0*velxL0 + 7.0*rhoR0*velxR0 - rhoR1*velxR1);
        number f_c_1 = (1.0 / 12.0) * (-rhoL1velxL1*velxL1 + 7.0*rhoL0velxL0*velxL0 + 7.0*rhoR0velxR0*velxR0 - rhoR1velxR1*velxR1);
        number f_c_2 = (1.0 / 12.0) * (-rhoL1velyL1*velxL1 + 7.0*rhoL0velyL0*velxL0 + 7.0*rhoR0velyR0*velxR0 - rhoR1velyR1*velxR1);
        number f_c_3 = (1.0 / 12.0) * (-rhoL1velzL1*velxL1 + 7.0*rhoL0velzL0*velxL0 + 7.0*rhoR0velzR0*velxR0 - rhoR1velzR1*velxR1);
        number f_c_4 = (1.0 / 12.0) * (-rhoL1eL1*velxL1 + 7.0*rhoL0eL0*velxL0 + 7.0*rhoR0eR0*velxR0 - rhoR1eR1*velxR1);
        number f_c_5 = (1.0 / 12.0) * (-rhoL1velxL1*velxL1*velxL1 + 7.0*rhoL0velxL0*velxL0*velxL0 + 7.0*rhoR0velxR0*velxR0*velxR0 - rhoR1velxR1*velxR1*velxR1);
        number f_c_6 = (1.0 / 12.0) * (-rhoL1velyL1*velyL1*velxL1 + 7.0*rhoL0velyL0*velyL0*velxL0 + 7.0*rhoR0velyR0*velyR0*velxR0 - rhoR1velyR1*velyR1*velxR1);
        number f_c_7 = (1.0 / 12.0) * (-rhoL1velzL1*velzL1*velxL1 + 7.0*rhoL0velzL0*velzL0*velxL0 + 7.0*rhoR0velzR0*velzR0*velxR0 - rhoR1velzR1*velzR1*velxR1);
        number f_c_8 = (1.0 / 12.0) * (-pL1*velxL1 + 7.0*pL0*velxL0 + 7.0*pR0*velxR0 - pR1*velxR1);
        number f_c_9 = (1.0 / 12.0) * (-pL1 + 7.0*pL0 + 7.0*pR0 - pR1);
        //
        number f_e_0 = (1.0 / 12.0) * (-rhoL1*velxR0 - rhoR0*velxL1 + 8.0*rhoL0*velxR0 + 8.0*rhoR0*velxL0 - rhoL0*velxR1 - rhoR1*velxL0);
        number f_e_1 = (1.0 / 12.0) * (-rhoL1velxL1*velxR0 - rhoR0velxR0*velxL1 + 8.0*rhoL0velxL0*velxR0 + 8.0* rhoR0velxR0*velxL0 - rhoL0velxL0*velxR1 - rhoR1velxR1*velxL0);
        number f_e_2 = (1.0 / 12.0) * (-rhoL1velyL1*velxR0 - rhoR0velyR0*velxL1 + 8.0*rhoL0velyL0*velxR0 + 8.0* rhoR0velyR0*velxL0 - rhoL0velyL0*velxR1 - rhoR1velyR1*velxL0);
        number f_e_3 = (1.0 / 12.0) * (-rhoL1velzL1*velxR0 - rhoR0velzR0*velxL1 + 8.0*rhoL0velzL0*velxR0 + 8.0* rhoR0velzR0*velxL0 - rhoL0velzL0*velxR1 - rhoR1velzR1*velxL0);
        number f_e_4 = (1.0 / 12.0) * (-rhoL1eL1*velxR0 - rhoR0eR0*velxL1 + 8.0*rhoL0eL0*velxR0 + 8.0* rhoR0eR0*velxL0 - rhoL0eL0*velxR1 - rhoR1eR1*velxL0);
        number f_e_5 = (1.0 / 12.0) * (-rhoL1velxL1*velxL1*velxR0 - rhoR0velxR0*velxR0*velxL1 + 8.0*rhoL0velxL0*velxL0*velxR0 + 8.0*rhoR0velxR0*velxR0*velxL0 - rhoL0velxL0*velxL0*velxR1 - rhoR1velxR1*velxR1*velxL0);
        number f_e_6 = (1.0 / 12.0) * (-rhoL1velyL1*velyL1*velxR0 - rhoR0velyR0*velyR0*velxL1 + 8.0*rhoL0velyL0*velyL0*velxR0 + 8.0*rhoR0velyR0*velyR0*velxL0 - rhoL0velyL0*velyL0*velxR1 - rhoR1velyR1*velyR1*velxL0);
        number f_e_7 = (1.0 / 12.0) * (-rhoL1velzL1*velzL1*velxR0 - rhoR0velzR0*velzR0*velxL1 + 8.0*rhoL0velzL0*velzL0*velxR0 + 8.0*rhoR0velzR0*velzR0*velxL0 - rhoL0velzL0*velzL0*velxR1 - rhoR1velzR1*velzR1*velxL0);
        number f_e_8 = (1.0 / 12.0) * (-pL1*velxR0 - pR0*velxL1 + 8.0*pL0*velxR0 + 8.0*pR0*velxL0 - pL0*velxR1 - pR1*velxL0);
        number f_e_9 = (1.0 / 12.0) * (-pL1 - pR0 + 8.0*pL0 + 8.0*pR0 - pL0 - pR1);
        //
        number alpha_mass = 1.0; number alpha_mom = 0.5; number alpha_ie = 0.5; number alpha_ke = 0.0; number alpha_p = 0.0;
        //
        F[CQI::mass] = alpha_mass * f_c_0 + (1.0 - alpha_mass) * f_e_0;
        number mom_x = alpha_mom * f_c_1 + (1.0 - alpha_mom) * f_e_1 + (alpha_p * f_c_9 + (1.0 - alpha_p) * f_e_9);
        number mom_y = alpha_mom * f_c_2 + (1.0 - alpha_mom) * f_e_2;
        number mom_z = alpha_mom * f_c_3 + (1.0 - alpha_mom) * f_e_3;
        //
        Vector3 momentum{mom_x, mom_y, mom_z};
        momentum.transform_to_global_frame(n, t1, t2);
        F[CQI::xMom] = momentum.x;
        F[CQI::yMom] = momentum.y;
        F[CQI::zMom] = momentum.z;
        F[CQI::totEnergy] = alpha_ie * f_c_4 + (1.0 - alpha_ie) * f_e_4 + (1.0 / 2.0) * (alpha_ke * f_c_5 + (1.0 - alpha_ke) * f_e_5 + alpha_ke * f_c_6 + (1.0 - alpha_ke) * f_e_6 + alpha_ke * f_c_7 + (1.0 - alpha_ke) * f_e_7) + alpha_p * f_c_8 + (1.0 - alpha_p) * f_e_8;
    } // end sbp_asf()


    // And one generic flux calculation function.

    __host__ __device__
    void calculate_flux(FlowState& fsL1, FlowState& fsL0, FlowState& fsR0, FlowState& fsR1, int x_order)
    // Generic flux calculation function.
    {
        // Choose the flavour of flux calculator at compile time.
        constexpr int flux_calculator = FluxCalc::sbp_asf;
        //
        if (flux_calculator ==  FluxCalc::ausmdv) {
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
        } else if (flux_calculator == FluxCalc::sbp_asf) {
            sbp_asf(fsL1, fsL0, fsR0, fsR1);
        }
    } // end calculate_flux()

}; // end FVFace

#endif
