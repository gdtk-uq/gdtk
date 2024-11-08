// flux.d -- Part of the Puffin and Lorikeet flow calculators.
//
// PA Jacobs
// 2022-02-02
//
module flux;

import std.format;
import std.math;

import geom;
import gas;
import gasdyn.gasflow;
import config;
import flow;
import face;
import cell;
import interp;


@nogc
bool is_shock(Face2D f, double compression_tol=-0.01, double shear_tol=0.20)
// A compression in the normal velocity field will have
// a decrease in normal velocity in the direction of the face normal.
// If the shear velocity is too high, we will suppress the shock value.
// compression_tol: a value of -0.30 is default for Eilmer, however,
//   we expect somewhat weak shocks in our space-marching solution.
// shear_tol: a value of 0.20 is default for Eilmer.
{
    auto fsL = &(f.left_cells[0].fs);
    auto fsR = &(f.right_cells[0].fs);
    // We have two cells interacting.
    // Compare the relative gas velocities normal to the face.
    double velxL = geom.dot(fsL.vel, f.n);
    double velxR = geom.dot(fsR.vel, f.n);
    double aL = fsL.gas.a;
    double aR = fsR.gas.a;
    double a_min = (aL < aR) ? aL : aR;
    double comp = (velxR-velxL)/a_min;
    //
    double velyL = geom.dot(fsL.vel, f.t1);
    double velyR = geom.dot(fsR.vel, f.t1);
    double sound_speed = 0.5*(aL+aR);
    double shear = fabs(velyL-velyR)/sound_speed;
    //
    return ((shear < shear_tol) && (comp < compression_tol));
} // end is_shock()

@nogc
void calculate_flux(Face2D f, ref FlowState2D fsL, ref FlowState2D fsR,
                    GasModel gmodel, FluxCalcCode flux_calc, int x_order, in CQIndex cqi)
// Compute the face's flux vector from left and right flow states.
// If requested, high-order reconstruction is applied.
{
    if (x_order == 1) {
        // Low-order reconstruction is just copy the nearest cell-centre flowstate.
        fsL.copy_values_from(f.left_cells[0].fs);
        fsR.copy_values_from(f.right_cells[0].fs);
    } else {
        // High-order reconstruction from the left_cells and right_cells stencil.
        f.interp_l2r2(fsL, fsR, gmodel, false);
    }
    final switch (flux_calc) {
    case FluxCalcCode.ausmdv:
        calculate_flux_ausmdv(f, fsL, fsR, gmodel, cqi);
        break;
    case FluxCalcCode.hanel:
        calculate_flux_hanel(f, fsL, fsR, gmodel, cqi);
        break;
    case FluxCalcCode.riemann:
        calculate_flux_riemann(f, fsL, fsR, gmodel, cqi);
        break;
    case FluxCalcCode.ausmdv_plus_hanel:
        calculate_flux_ausmdv_plus_hanel(f, fsL, fsR, gmodel, cqi);
        break;
    case FluxCalcCode.riemann_plus_hanel:
        calculate_flux_riemann_plus_hanel(f, fsL, fsR, gmodel, cqi);
        break;
    }
} // end calculate_flux()

@nogc
void calculate_flux_ausmdv_plus_hanel(Face2D f, in FlowState2D fsL, in FlowState2D fsR,
                                      GasModel gmodel, in CQIndex cqi)
// Compute the face's flux vector from left and right flow states.
// We actually delegate the detailed calculation to one of the other calculators
// depending on the shock indicator.
{
    if (f.left_cells[0].shock_flag || f.right_cells[0].shock_flag) {
        calculate_flux_hanel(f, fsL, fsR, gmodel, cqi);
    } else {
        calculate_flux_ausmdv(f, fsL, fsR, gmodel, cqi);
    }
    return;
}

@nogc
void calculate_flux_riemann_plus_hanel(Face2D f, in FlowState2D fsL, in FlowState2D fsR,
                                       GasModel gmodel, in CQIndex cqi)
// Compute the face's flux vector from left and right flow states.
// We actually delegate the detailed calculation to one of the other calculators
// depending on the shock indicator.
{
    if (f.left_cells[0].shock_flag || f.right_cells[0].shock_flag) {
        calculate_flux_hanel(f, fsL, fsR, gmodel, cqi);
    } else {
        calculate_flux_riemann(f, fsL, fsR, gmodel, cqi);
    }
    return;
}

@nogc
void calculate_flux_riemann(Face2D f, in FlowState2D fsL, in FlowState2D fsR,
                            GasModel gmodel, in CQIndex cqi)
// Compute the face's flux vector from left and right flow states.
// The core of this calculation is the one-dimensional Riemann solver
// from the gasflow module.
{
    Vector3 velL = Vector3(fsL.vel);
    Vector3 velR = Vector3(fsR.vel);
    velL.transform_to_local_frame(f.n, f.t1);
    velR.transform_to_local_frame(f.n, f.t1);
    double[5] rsol = osher_riemann(fsL.gas, fsR.gas, velL.x, velR.x,
                                   f.stateLstar, f.stateRstar, f.stateX0, gmodel);
    double rho = f.stateX0.rho;
    double p = f.stateX0.p;
    double u = gmodel.internal_energy(f.stateX0);
    double velx = rsol[4];
    double vely = (velx > 0.0) ? velL.y : velR.y;
    double massFlux = rho*velx;
    Vector3 momentum = Vector3(massFlux*velx+p, massFlux*vely);
    momentum.transform_to_global_frame(f.n, f.t1);
    f.F[cqi.mass] = massFlux;
    f.F[cqi.xMom] = momentum.x;
    f.F[cqi.yMom] = momentum.y;
    f.F[cqi.totEnergy] = massFlux*(u+p/rho+0.5*(velx*velx+vely*vely));
    if (cqi.n_species > 1) {
        foreach (i; 0 .. cqi.n_species) {
            f.F[cqi.species+i] = massFlux * ((velx < 0.0) ? fsL.gas.massf[i] : fsR.gas.massf[i]);
        }
    }
    foreach (i; 0 .. cqi.n_modes) {
        f.F[cqi.modes+i] = massFlux * ((velx < 0.0) ? fsL.gas.u_modes[i] : fsR.gas.u_modes[i]);
    }
    bool allFinite = true;
    foreach (e; f.F) { if (!isFinite(e)) { allFinite = false; } }
    if (!allFinite) {
        debug { import std.stdio;  writeln("face=", f); }
        throw new Exception("At least one flux quantity is not finite.");
    }
    return;
} // end calculate_flux()

@nogc
void calculate_flux_ausmdv(Face2D f, in FlowState2D fsL, in FlowState2D fsR,
                           GasModel gmodel, in CQIndex cqi)
// Compute the face's flux vector from left and right flow states.
// Wada and Liou's flux calculator, implemented from details in their AIAA paper,
// with hints from Ian Johnston.
// Y. Wada and M. -S. Liou (1994)
// A flux splitting scheme with high-resolution and robustness for discontinuities.
// AIAA-94-0083.
{
    Vector3 velL = Vector3(fsL.vel);
    Vector3 velR = Vector3(fsR.vel);
    velL.transform_to_local_frame(f.n, f.t1);
    velR.transform_to_local_frame(f.n, f.t1);
    //
    double rhoL = fsL.gas.rho;
    double pL = fsL.gas.p;
    double pLrL = pL/rhoL;
    double velxL = velL.x;
    double velyL = velL.y;
    double uL = gmodel.internal_energy(fsL.gas);
    double aL = fsL.gas.a;
    double keL = 0.5*(velxL*velxL + velyL*velyL);
    double HL = uL + pLrL + keL;
    //
    double rhoR = fsR.gas.rho;
    double pR = fsR.gas.p;
    double pRrR = pR/rhoR;
    double velxR = velR.x;
    double velyR = velR.y;
    double uR = gmodel.internal_energy(fsR.gas);
    double aR = fsR.gas.a;
    double keR = 0.5*(velxR*velxR + velyR*velyR);
    double HR = uR + pR/rhoR + keR;
    //
    // This is the main part of the flux calculator.
    //
    // Weighting parameters (eqn 32) for velocity splitting.
    double alphaL = 2.0*pLrL/(pLrL+pRrR);
    double alphaR = 2.0*pRrR/(pLrL+pRrR);
    // Common sound speed (eqn 33) and Mach numbers.
    double am = fmax(aL, aR);
    double ML = velxL/am;
    double MR = velxR/am;
    // Left state:
    // pressure splitting (eqn 34)
    // and velocity splitting (eqn 30)
    double pLplus, velxLplus;
    double dvelxL = 0.5 * (velxL + fabs(velxL));
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
    double pRminus, velxRminus;
    double dvelxR = 0.5*(velxR-fabs(velxR));
    if (fabs(MR) <= 1.0) {
        pRminus = pR*(MR-1.0)*(MR-1.0)*(2.0+MR)*0.25;
        velxRminus = alphaR*(-(velxR-am)*(velxR-am)/(4.0*am) - dvelxR) + dvelxR;
    } else {
        pRminus = pR * dvelxR / velxR;
        velxRminus = dvelxR;
    }
    // The mass flux. (eqn 29)
    double massL = velxLplus*rhoL;
    double massR = velxRminus*rhoR;
    double mass_half = massL+massR;
    // Pressure flux (eqn 34)
    double p_half = pLplus + pRminus;
    // Momentum flux: normal direction
    // Compute blending parameter s (eqn 37),
    // the momentum flux for AUSMV (eqn 21) and AUSMD (eqn 21)
    // and blend (eqn 36).
    double dp = pL - pR;
    const double K_SWITCH = 10.0;
    dp = K_SWITCH * fabs(dp) / fmin(pL, pR);
    double s = 0.5 * fmin(1.0, dp);
    double rvel2_AUSMV = massL*velxL + massR*velxR;
    double rvel2_AUSMD = 0.5*(mass_half*(velxL+velxR) - fabs(mass_half)*(velxR-velxL));
    double rvel2_half = (0.5+s)*rvel2_AUSMV + (0.5-s)*rvel2_AUSMD;
    // Assemble components of the flux vector (eqn 36).
    f.F[cqi.mass] = mass_half;
    double vely = (mass_half >= 0.0) ? velyL : velyR;
    Vector3 momentum = Vector3(rvel2_half+p_half, mass_half*vely);
    momentum.transform_to_global_frame(f.n, f.t1);
    f.F[cqi.xMom] = momentum.x;
    f.F[cqi.yMom] = momentum.y;
    double H = (mass_half >= 0.0) ? HL : HR;
    f.F[cqi.totEnergy] = mass_half*H;
    if (cqi.n_species > 1) {
        foreach (i; 0 .. cqi.n_species) {
            double massf = (mass_half >= 0.0) ? fsL.gas.massf[i] : fsR.gas.massf[i];
            f.F[cqi.species+i] = mass_half*massf;
        }
    }
    foreach (i; 0 .. cqi.n_modes) {
        double u_mode = (mass_half >= 0.0) ? fsL.gas.u_modes[i] : fsR.gas.u_modes[i];
        f.F[cqi.modes+i] = mass_half*u_mode;
    }
    //
    bool allFinite = true;
    foreach (e; f.F) { if (!isFinite(e)) { allFinite = false; } }
    if (!allFinite) {
        debug { import std.stdio;  writeln("face=", f); }
        throw new Exception("At least one flux quantity is not finite.");
    }
    return;
} // end calculate_flux_ausmdv()

@nogc
void calculate_flux_hanel(Face2D f, in FlowState2D fsL, in FlowState2D fsR,
                          GasModel gmodel, in CQIndex cqi)
// Compute the face's flux vector from left and right flow states.
// Implemented from Y. Wada and M. S. Liou details in their AIAA paper
// Y. Wada and M. -S. Liou (1997)
// An accurate and robust flux splitting scheme for shock and contact discontinuities.
// with reference to....
// Hanel, Schwane, & Seider's 1987 paper
// On the accuracy of upwind schemes for the solution of the Navier-Stokes equations
{
    Vector3 velL = Vector3(fsL.vel);
    Vector3 velR = Vector3(fsR.vel);
    velL.transform_to_local_frame(f.n, f.t1);
    velR.transform_to_local_frame(f.n, f.t1);
    //
    double rhoL = fsL.gas.rho;
    double pL = fsL.gas.p;
    double velxL = velL.x;
    double velyL = velL.y;
    double uL = gmodel.internal_energy(fsL.gas);
    double aL = fsL.gas.a;
    double keL = 0.5*(velxL*velxL + velyL*velyL);
    double HL = uL + pL/rhoL + keL;
    //
    double rhoR = fsR.gas.rho;
    double pR = fsR.gas.p;
    double velxR = velR.x;
    double velyR = velR.y;
    double uR = gmodel.internal_energy(fsR.gas);
    double aR = fsR.gas.a;
    double keR = 0.5*(velxR*velxR + velyR*velyR);
    double HR = uR + pR/rhoR + keR;
    //
    double am = fmax(aL, aR);
    double ML = velxL/am;
    double MR = velxR/am;
    // Left state:
    // pressure splitting (eqn 7)
    // and velocity splitting (eqn 9)
    double pLplus, velxLplus;
    if (fabs(velxL) <= aL) {
        velxLplus = 1.0/(4.0*aL) * (velxL+aL)*(velxL+aL);
        pLplus = pL*velxLplus * (1.0/aL * (2.0-velxL/aL));
    } else {
        velxLplus = 0.5*(velxL+fabs(velxL));
        pLplus = pL*velxLplus * (1.0/velxL);
    }
    // Right state:
    // pressure splitting (eqn 7)
    // and velocity splitting (eqn 9)
    double pRminus, velxRminus;
    if (fabs(velxR) <= aR) {
        velxRminus = -1.0/(4.0*aR) * (velxR-aR)*(velxR-aR);
        pRminus = pR*velxRminus * (1.0/aR * (-2.0-velxR/aR));
    } else {
        velxRminus = 0.5*(velxR-fabs(velxR));
        pRminus = pR*velxRminus * (1.0/velxR);
    }
    // The mass flux.
    double massL = velxLplus * rhoL;
    double massR = velxRminus * rhoR;
    double mass_half = massL + massR;
    // Pressure flux (eqn 8)
    double p_half = pLplus + pRminus;
    // Assemble components of the flux vector (eqn 36).
    f.F[cqi.mass] = massL + massR;
    Vector3 momentum = Vector3(massL*velxL + massR*velxR + p_half, massL*velyL + massR*velyR);
    momentum.transform_to_global_frame(f.n, f.t1);
    f.F[cqi.xMom] = momentum.x;
    f.F[cqi.yMom] = momentum.y;
    f.F[cqi.totEnergy] = massL*HL + massR*HR;
    if (cqi.n_species > 1) {
        foreach (i; 0 .. cqi.n_species) {
            f.F[cqi.species+i] = massL*fsL.gas.massf[i] + massR*fsR.gas.massf[i];
        }
    }
    foreach (i; 0 .. cqi.n_modes) {
        f.F[cqi.modes+i] = massL*fsL.gas.u_modes[i] + massR*fsR.gas.u_modes[i];
    }
    //
    bool allFinite = true;
    foreach (e; f.F) { if (!isFinite(e)) { allFinite = false; } }
    if (!allFinite) {
        debug { import std.stdio;  writeln("face=", f); }
        throw new Exception("At least one flux quantity is not finite.");
    }
    return;
} // end calculate_flux_hanel()

@nogc
void simple_flux(Face2D f, in FlowState2D fs, GasModel gmodel, in CQIndex cqi)
// Computes the face's flux vector from a single flow state.
// Supersonic flow is assumed.
{
    Vector3 vel = Vector3(fs.vel);
    vel.transform_to_local_frame(f.n, f.t1);
    double rho = fs.gas.rho;
    double p = fs.gas.p;
    double u = gmodel.internal_energy(fs.gas);
    double massFlux = rho * vel.x;
    Vector3 momentum = Vector3(massFlux*vel.x+p, massFlux*vel.y);
    momentum.transform_to_global_frame(f.n, f.t1);
    f.F[cqi.mass] = massFlux;
    f.F[cqi.xMom] = momentum.x;
    f.F[cqi.yMom] = momentum.y;
    f.F[cqi.totEnergy] = massFlux * (u+p/rho+0.5*(vel.x*vel.x+vel.y*vel.y));
    if (cqi.n_species > 1) {
        foreach (i; 0 .. cqi.n_species) {
            f.F[cqi.species+i] = massFlux * fs.gas.massf[i];
        }
    }
    foreach (i; 0 .. cqi.n_modes) {
        f.F[cqi.modes+i] = massFlux * fs.gas.u_modes[i];
    }
    bool allFinite = true;
    foreach (e; f.F) { if (!isFinite(e)) { allFinite = false; } }
    if (!allFinite) {
        debug { import std.stdio;  writeln("face=", f); }
        throw new Exception("At least one flux quantity is not finite.");
    }
    return;
} // end simple_flux()
