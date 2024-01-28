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
#include "flow.cu"

using namespace std;


namespace BCCode {
    // Boundary condition codes, to decide what to do for the ghost cells.
    // Periodic boundary conditions should just work if we wrap the index in each direction.
    // There's not enough information here to have arbitrary block connections.
    constexpr int wall_with_slip = 0;
    constexpr int wall_no_slip_adiabatic = 1;
    constexpr int wall_no_slip_fixed_T = 2;
    constexpr int exchange = 3;
    constexpr int inflow = 4;
    constexpr int outflow = 5;
    constexpr int inflow_function = 6;

    array<string,7> names{"wall_with_slip", "wall_no_slip_adiabatic", "wall_no_slip_fixed_T",
            "exchange", "inflow", "outflow", "inflow_function"};
};

int BC_code_from_name(string name)
{
    if (name == "wall_with_slip") return BCCode::wall_with_slip;
    if (name == "wall_no_slip_adiabatic") return BCCode::wall_no_slip_adiabatic;
    if (name == "wall_no_slip_fixed_T") return BCCode::wall_no_slip_fixed_T;
    if (name == "wall_no_slip") return BCCode::wall_no_slip_adiabatic; // alias
    if (name == "exchange") return BCCode::exchange;
    if (name == "inflow") return BCCode::inflow;
    if (name == "outflow") return BCCode::outflow;
    if (name == "inflow_function") return BCCode::inflow_function;
    return BCCode::wall_with_slip;
}

namespace BCFunction {
    // We have a number of functions coded to use in filling the ghost cells.
    constexpr int none = 0;
    constexpr int supersonic_vortex = 1;
    constexpr int laminar_boundary_layer = 2;
    constexpr int manufactured_solution = 3;

    array<string,4> names{"none", "supersonic_vortex", "laminar_boundary_layer", "manufactured_solution"};
};

int BC_function_from_name(string name)
{
    if (name == "none") return BCFunction::none;
    if (name == "supersonic_vortex") return BCFunction::supersonic_vortex;
    if (name == "laminar_boundary_layer") return BCFunction::laminar_boundary_layer;
    if (name == "manufactured_solution") return BCFunction::manufactured_solution;
    return BCFunction::none;
}

// Interpolation functions that will be used in the fluc calculator.

__host__ __device__
number van_albada_limit1(number a, number b)
// A smooth slope limiter.
{
    constexpr number eps = (number)1.0e-12;
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
    qL = qL0 + sL * one_eighth * (three*del + delLminus);
    qR = qR0 - sR * one_eighth * (delRplus + three*del);
} // end of interp_l2r2_scalar()


// Flux calculations are done in the context of a face on a cell.

constexpr int cloud_ncmax = 2;
constexpr int cloud_nfmax = 8;
constexpr int cloud_nmax = cloud_ncmax + cloud_nfmax;

struct FVFace {
    Vector3 pos; // midpoint position in space
    number area;
    Vector3 n;  // unit normal
    Vector3 t1; // unit tangent 1
    Vector3 t2; // unit tangent 2
    array<number,CQI::n> F; // flux vector for conserved quantities
    // We will keep connections to the pieces composing the face
    // as indices into global arrays.
    array<int,4> vtx{-1,-1,-1,-1};
    array<int,2> left_cells{-1,-1};
    array<int,2> right_cells{-1,-1};
    // To apply boundary conditions at the face, we need to carry some extra information.
    // Not all faces are boundary faces, so we start with dummy values.
    int bcId{-1};
    int bcCode{-1};
    int other_blkId{-1};
    array<int,2> other_cells{-1,-1};
    int inflowId{-1};
    number TWall{(number)300.0};
    int bcFun{-1};
    // For the gradient calculations that form part of the viscous fluxes
    // we keep lists of faces and cells that form a cloud of points around
    // this face-centre.
    array<int,cloud_ncmax> cells_in_cloud{-1,-1};
    array<int,cloud_nfmax> faces_in_cloud{-1,-1,-1,-1,-1,-1,-1,-1};
    int cloud_nc = 0;
    int cloud_nf = 0;
    // Prepared least-squares solution for cloud of cell- and face-FlowStates.
    array<number,cloud_nmax> wx, wy, wz;
    // We also need the FlowState at this face-centre.  It will be set during
    // the convective-flux calculation or by the boundary-condition code for a wall.
    FlowState fs;

    string toString() const
    {
        ostringstream repr;
        repr << "FVFace(pos=" << pos.toString() << ", n=" << n.toString()
             << ", t1=" << t1.toString() << ", t2=" << t2.toString() << ", area=" << area;
        repr << ", F=["; for (auto v : F) repr << v << ","; repr << "]";
        repr << ", vtx=["; for (auto i : vtx) repr << i << ","; repr << "]";
        repr << ", left_cells=["; for (auto i : left_cells) repr << i << ","; repr << "]";
        repr << ", right_cells=["; for (auto i : right_cells) repr << i << ","; repr << "]";
        repr << ", bcId=" << bcId << ", bcCode=" << bcCode
             << ", other_blkId=" << other_blkId << ", other_cells=[" << other_cells[0]
             << ", " << other_cells[1] << "], inflowId=" << inflowId
             << ", TWall=" << TWall << ", bcFun=" << bcFun;
        repr << ", cloud_nc=" << cloud_nc << ", cloud_nf=" << cloud_nf;
        repr << ", cells_in_cloud=["; for (auto i : cells_in_cloud) repr << i <<","; repr << "]";
        repr << ", faces_in_cloud=["; for (auto i : faces_in_cloud) repr << i <<","; repr << "]";
        repr << ", fs=" << fs;
        repr << ")";
        return repr.str();
    }


    // Specific convective-flux calculators here...

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
        number keL = half*(velxL*velxL + velyL*velyL + velzL*velzL);
        number HL = uL + pLrL + keL;
        number YBL = fsL.gas.YB;
        //
        number rhoR = fsR.gas.rho;
        number pR = fsR.gas.p;
        number pRrR = pR/rhoR;
        number velxR = velR.x;
        number velyR = velR.y;
        number velzR = velR.z;
        number uR = fsR.gas.e;
        number aR = fsR.gas.a;
        number keR = half*(velxR*velxR + velyR*velyR + velzR*velzR);
        number HR = uR + pR/rhoR + keR;
        number YBR = fsR.gas.YB;
        //
        // This is the main part of the flux calculator.
        //
        // Weighting parameters (eqn 32) for velocity splitting.
        number alphaL = two*pLrL/(pLrL+pRrR);
        number alphaR = two*pRrR/(pLrL+pRrR);
        // Common sound speed (eqn 33) and Mach numbers.
        number am = fmax(aL, aR);
        number ML = velxL/am;
        number MR = velxR/am;
        // Left state:
        // pressure splitting (eqn 34)
        // and velocity splitting (eqn 30)
        number pLplus, velxLplus;
        number dvelxL = half*(velxL+fabs(velxL));
        if (fabs(ML) <= one) {
            pLplus = pL*(ML+one)*(ML+one)*(two-ML)*one_quarter;
            velxLplus = alphaL*((velxL+am)*(velxL+am)/(four*am) - dvelxL) + dvelxL;
        } else {
            pLplus = pL * dvelxL / velxL;
            velxLplus = dvelxL;
        }
        // Right state:
        // pressure splitting (eqn 34)
        // and velocity splitting (eqn 31)
        number pRminus, velxRminus;
        number dvelxR = half*(velxR-fabs(velxR));
        if (fabs(MR) <= one) {
            pRminus = pR*(MR-one)*(MR-one)*(two+MR)*one_quarter;
            velxRminus = alphaR*(-(velxR-am)*(velxR-am)/(four*am) - dvelxR) + dvelxR;
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
        constexpr number K_SWITCH = (number)10.0;
        dp = K_SWITCH * fabs(dp) / fmin(pL, pR);
        number s = half * fmin(one, dp);
        number rvel2_AUSMV = massL*velxL + massR*velxR;
        number rvel2_AUSMD = half*(mass_half*(velxL+velxR) - fabs(mass_half)*(velxR-velxL));
        number rvel2_half = (half+s)*rvel2_AUSMV + (half-s)*rvel2_AUSMD;
        // Assemble components of the flux vector (eqn 36).
        F[CQI::mass] = mass_half;
        number vely = (mass_half >= zero) ? velyL : velyR;
        number velz = (mass_half >= zero) ? velzL : velzR;
        Vector3 momentum{rvel2_half+p_half, mass_half*vely, mass_half*velz};
        momentum.transform_to_global_frame(n, t1, t2);
        F[CQI::xMom] = momentum.x;
        F[CQI::yMom] = momentum.y;
        F[CQI::zMom] = momentum.z;
        number H = (mass_half >= zero) ? HL : HR;
        F[CQI::totEnergy] = mass_half*H;
        number YB = (mass_half >= zero) ? YBL : YBR;
        F[CQI::YB] = mass_half*YB;
        return;
    } // end ausmdv()


    __host__ __device__
    void sbp_asf(const FlowState& fsL1, const FlowState& fsL0,
                 const FlowState& fsR0, const FlowState& fsR1)
    // Lachlan's and Christine's Summation-By-Parts Alpha-Split Flux calculation function.
    //
    // This flux calculator is based on the formulation in the NASA Techmical Memo
    // Travis C. Fisher, Mark H. Carpenter, Jan Nordstroem, Nail K. Yamaleev and R. Charles Swanson
    // Discretely Conservative Finite-Difference Formulations for Nonlinear Conservation Laws
    // in Split Form: Theory and Boundary Conditions
    // NASA/TM-2011-217307  November 2011
    {
        // Get local copies of the near-by flow states and transform into the frame
        // that is local to the interface, with the local x-direction aligned with
        // the face normal.
        array<GasState,4> gas{fsL1.gas, fsL0.gas, fsR0.gas, fsR1.gas};
        array<Vector3,4> vel{fsL1.vel, fsL0.vel, fsR0.vel, fsR1.vel};
        for (int i = 0; i < 4; ++i) { vel[i].transform_to_local_frame(n, t1, t2); }
        //
        // Factored terms from the conservtion equations.
        number v[10][4]; number w[10][4];
        for (int i = 0; i < 4; ++i) {
            number rho = gas[i].rho; number p = gas[i].p; number e = gas[i].e;
            number velx = vel[i].x; number vely = vel[i].y; number velz = vel[i].z;
            v[0][i] = rho;           w[0][i] = velx;  // mass
            v[1][i] = velx*rho;      w[1][i] = velx;  // x-momentum
            v[2][i] = vely*rho;      w[2][i] = velx;  // y-momentum
            v[3][i] = velz*rho;      w[3][i] = velx;  // z-momentum
            v[4][i] = e*rho;         w[4][i] = velx;  // Internal energy
            v[5][i] = velx*velx*rho; w[5][i] = velx;  // Kinetic energy, x
            v[6][i] = vely*vely*rho; w[6][i] = velx;  // Kinetic energy, y
            v[7][i] = velz*velz*rho; w[7][i] = velx;  // Kinetic energy, z
            v[8][i] = p;             w[8][i] = velx;  // Pressure work in energy equation.
            v[9][i] = p;             w[9][i] = 1;     // Pressure in momentum equation
        }
        //
        // Fluxes.
        number f_c[10]; number f_e[10];
        for (int j = 0; j < 10; ++j) {
            // Divergence-form flux (eq 3.5 in NASA/TM-2011-217307)
            f_c[j] = one_twelfth * (-v[j][0]*w[j][0] + seven*v[j][1]*w[j][1] + seven*v[j][2]*w[j][2] - v[j][3]*w[j][3]);
            // Product-rule flux (eq 3.6 in NASA/TM-2011-217307)
            f_e[j] = one_twelfth * (-v[j][0]*w[j][2] - v[j][2]*w[j][0] + eight*v[j][1]*w[j][2]
                                    + eight*v[j][2]*w[j][1] - v[j][1]*w[j][3] - v[j][3]*w[j][1]);
        }
        //
        // Define the splitting values as per White et al,
        // in the conservative skew-symmetric form of Honein and Moin.
        constexpr number alpha_mass = (number)1.0;
        constexpr number alpha_mom = (number)0.5;
        constexpr number alpha_ie = (number)0.5;
        constexpr number alpha_ke = (number)0.0;
        constexpr number alpha_p = (number)0.0;
        //
        // Assemble weighted fluxes.
        number mass_flux = alpha_mass*f_c[0] + (one-alpha_mass)*f_e[0];
        F[CQI::mass] = mass_flux;
        //
        number mom_x = alpha_mom*f_c[1] + (one-alpha_mom)*f_e[1] + (alpha_p*f_c[9] + (one-alpha_p)*f_e[9]);
        number mom_y = alpha_mom*f_c[2] + (one-alpha_mom)*f_e[2];
        number mom_z = alpha_mom*f_c[3] + (one-alpha_mom)*f_e[3];
        Vector3 momentum{mom_x, mom_y, mom_z};
        momentum.transform_to_global_frame(n, t1, t2);
        F[CQI::xMom] = momentum.x;
        F[CQI::yMom] = momentum.y;
        F[CQI::zMom] = momentum.z;
        //
        F[CQI::totEnergy] = alpha_ie*f_c[4] + (one-alpha_ie)*f_e[4] +
            half*(alpha_ke*f_c[5] + (one-alpha_ke)*f_e[5] +
                  alpha_ke*f_c[6] + (one-alpha_ke)*f_e[6] +
                  alpha_ke*f_c[7] + (one-alpha_ke)*f_e[7]) +
            alpha_p*f_c[8] + (one-alpha_p)*f_e[8];
        //
        // Choose which cell to use for the species, based on which way the wind is blowing.
        F[CQI::YB] = mass_flux * ((mass_flux >= zero) ? gas[1].YB : gas[2].YB);
    } // end sbp_asf()


    // And one generic flux calculation function.

    __host__ __device__
    void calculate_convective_flux(const FlowState& fsL1, const FlowState& fsL0,
                                   const FlowState& fsR0, const FlowState& fsR1,
                                   int flux_calc, int x_order)
    // Generic convective-flux calculation function.
    {
        if (flux_calc == 0) { // FluxCalc::ausmdv (PJ 2022-10-19 Name cannot be seen.)
            // First-order reconstruction is just a copy from the nearest cell centre.
            FlowState fsL{fsL0};
            FlowState fsR{fsR0};
            if (x_order > 1) {
                // We will interpolate only some GasState properties...
                interp_l2r2_scalar(fsL1.gas.rho, fsL0.gas.rho, fsR0.gas.rho, fsR1.gas.rho, fsL.gas.rho, fsR.gas.rho);
                interp_l2r2_scalar(fsL1.gas.e, fsL0.gas.e, fsR0.gas.e, fsR1.gas.e, fsL.gas.e, fsR.gas.e);
                interp_l2r2_scalar(fsL1.gas.YB, fsL0.gas.YB, fsR0.gas.YB, fsR1.gas.YB, fsL.gas.YB, fsR.gas.YB);
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
            // For later use in gradient calculations for viscous fluxes.
            fs.set_as_average(fsL,fsR);
        } else if (flux_calc == 1) { // FluxCalc::sbp_asf
            sbp_asf(fsL1, fsL0, fsR0, fsR1);
            fs.set_as_average(fsL0,fsR0);
        }
    } // end calculate_convective_flux()

}; // end FVFace


__host__
ostream& operator<<(ostream& os, const FVFace f)
{
    os << f.toString();
    return os;
}

#endif
