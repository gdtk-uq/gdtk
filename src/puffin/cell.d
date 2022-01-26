// cell.d -- Part of the Puffin steady-flow calculator.
//
// PA Jacobs
// 2022-01-22
//
module cell;

import std.format;
import std.math;

import geom;
import gas;
import gasflow;
import config;
import flow;


class Cell2D {
public:
    CQIndex cqi;
    Vector3 pos;
    Vector3* p00, p10, p11, p01;
    Face2D faceN, faceE, faceS, faceW;
    double volume, xyplane_area;
    double iLen, jLen, kLen; // distances across the cell
    //
    FlowState2D fs;
    double[][3] U; // Conserved quantities at time levels.
    double[][2] dUdt; // Time-derivatives of conserved quantities.
    //
    double dt; // Local time-step size.

    this(GasModel gmodel, CQIndex cqi)
    {
        this.cqi = new CQIndex(cqi);
        pos = Vector3();
        fs = new FlowState2D(gmodel);
        foreach (i; 0 .. U.length) { U[i].length = cqi.n; }
        foreach (i; 0 .. dUdt.length) { dUdt[i].length = cqi.n; }
    }

    this(ref const(Cell2D) other)
    {
        cqi = new CQIndex(other.cqi);
        pos = Vector3(other.pos);
        fs = new FlowState2D(other.fs);
        foreach (i; 0 .. U.length) { U[i] = other.U[i].dup; }
        foreach (i; 0 .. dUdt.length) { dUdt[i] = other.dUdt[i].dup; }
        // We do not bother copying all of the other references
        // because they will be mostly useless in a newly copied cell.
        // They will need to be set up manually, later.
    }

    override
    string toString() const
    {
        string repr = "Cell2D(";
        repr ~= format("pos=%s", pos);
        repr ~= format(", volume=%g, xyplane_area=%g", volume, xyplane_area);
        repr ~= format(", iLen=%g, jLen=%g", iLen, jLen);
        repr ~= format(", fs=%s", fs);
        foreach (i; 0 .. U.length) { repr ~= format(", U[%d]=%s", i, U[i]); }
        foreach (i; 0 .. dUdt.length) { repr ~= format(", dUdt[%d]=%s", i, dUdt[i]); }
        repr ~= format(", dt=%g", dt);
        repr ~= ")";
        return repr;
    }

    @nogc
    void compute_geometry(bool axiFlag)
    // Update the geometric properties from vertex data.
    {
        xyplane_quad_cell_properties(*p00, *p10, *p11, *p01, pos, xyplane_area,
                                     iLen, jLen, kLen);
        volume = xyplane_area;
        if (axiFlag) { volume *= pos.y; }
        if (volume < 0.0 || xyplane_area < 0.0) {
            throw new Exception("volume or xyplane_area negative.");
        }
        return;
    }

    @nogc
    void encode_conserved(size_t ftl, GasModel gmodel)
    // Compute converved quantities from current flow state.
    // Input:
    //   ftl: index to the particular vector of conserved quantities to decode.
    //   gmodel: a reference to the gas model object.
    {
        auto myU = U[ftl];
        // Mass per unit volume.
        auto rho = fs.gas.rho;
        myU[cqi.mass] = rho;
        // Momentum per unit volume.
        myU[cqi.xMom] = rho*fs.vel.x;
        myU[cqi.yMom] = rho*fs.vel.y;
        // Total energy per unit volume.
        double u = gmodel.internal_energy(fs.gas);
        double ke = 0.5*(fs.vel.x*fs.vel.x + fs.vel.y*fs.vel.y);
        myU[cqi.totEnergy] = rho*(u + ke);
        if (cqi.n_species > 1) {
            foreach(i; 0 .. cqi.n_species) { myU[cqi.species+i] = rho*fs.gas.massf[i]; }
        }
        foreach(i; 0 .. cqi.n_modes) { myU[cqi.modes+i] = rho*fs.gas.u_modes[i]; }
        return;
    }

    @nogc
    void decode_conserved(size_t ftl, GasModel gmodel)
    // Compute flow state from final vector of conserved quantities.
    // Quantities are per unit volume.
    // Input:
    //   ftl: index to the particular vector of conserved quantities to decode.
    //   gmodel: a reference to the gas model object.
    // Returns nonzero if the flowstate is bad.
    {
        auto myU = U[ftl];
        bool allFinite = true;
        foreach (e; myU) { if (!isFinite(e)) { allFinite = false; } }
        if (!allFinite) {
            debug { import std.stdio;  writeln("ftl=", ftl, " cell=", this); }
            throw new Exception("At least one conserved quantity is not finite.");
        }
        // Mass per unit volume.
        auto rho = myU[cqi.mass];
        if (rho <= 0.0) {
            throw new Exception("Density zero or negative in conserved quantity vector.");
        }
        fs.gas.rho = rho;
        double dinv = 1.0 / rho;
        // Velocities from momenta.
        fs.vel.set(myU[cqi.xMom]*dinv, myU[cqi.yMom]*dinv);
        // Divide up the total energy per unit volume.
        double rE = myU[cqi.totEnergy];
        // Start with the total energy, then take out the other components.
        // Internal energy is what remains.
        double u = rE * dinv;
        double ke = 0.5*(fs.vel.x*fs.vel.x + fs.vel.y*fs.vel.y);
        u -= ke;
        // Other energies, if any.
        foreach(i; 0 .. cqi.n_modes) { fs.gas.u_modes[i] = myU[cqi.modes+i] * dinv; }
        double u_other = 0.0; foreach(ei; fs.gas.u_modes) { u_other += ei; }
        fs.gas.u = u - u_other;
        // Thermochemical species, if appropriate.
        if (cqi.n_species > 1) {
            foreach(i; 0 .. cqi.n_species) { fs.gas.massf[i] = myU[cqi.species+i] * dinv; }
            scale_mass_fractions(fs.gas.massf);
        } else {
            fs.gas.massf[0] = 1.0;
        }
        try {
            gmodel.update_thermo_from_rhou(fs.gas);
            gmodel.update_sound_speed(fs.gas);
        } catch (Exception e) {
            debug {
                import std.stdio;
                writeln("decode_conserved: cell=", this);
                writeln("rho=", rho, " dinv=", dinv, " rE=", rE, " u=", u,
                        " ke=", ke, " u_other=", u_other);
            }
            throw new Exception("decode_conserved: update thermo failed.");
        }
        return;
    } // end decode_conserved()

    @nogc
    void eval_dUdt(size_t ftl, bool axiFlag)
    // These are the spatial (RHS) terms in the semi-discrete governing equations.
    // ftl : (flow-time-level) specifies where computed derivatives are to be stored.
    //       0: predictor update.
    //       1: corrector update.
    {
        auto my_dUdt = dUdt[ftl];
        double vol_inv = 1.0 / volume;
        foreach (i; 0 .. cqi.n) {
            // Integrate the fluxes across the interfaces that bound the cell.
            double surface_integral = faceS.area*faceS.F[i] + faceW.area*faceW.F[i]
                - faceN.area*faceN.F[i] - faceE.area*faceE.F[i];
            // Then evaluate the derivatives of conserved quantity.
            // Note that conserved quantities are stored per-unit-volume.
            my_dUdt[i] = vol_inv * surface_integral;
        }
        if (axiFlag) {
            my_dUdt[cqi.yMom] += fs.gas.p * xyplane_area * vol_inv;
        }
        return;
    } // end eval_dUdt()

    @nogc
    void estimate_local_dt(double cfl)
    {
        // We assume that the cells are (roughly) aligned with the xy directions.
        dt = cfl * fmin(iLen/(fabs(fs.vel.x)+fs.gas.a), jLen/(fabs(fs.vel.y)+fs.gas.a));
        return;
    }

} // end class Cell2D


class Face2D {
public:
    CQIndex cqi;
    Vector3 pos;
    Vector3* p0, p1; // pointers to vertices at each end of face
    Vector3 n; // unit normal (to right when looking from p0 to p1)
    Vector3 t1; // unit tangent (from p0 to p1)
    double area; // per unit depth for 2D planar, per radian for axisymmetric
    //
    double[] F; // Flux vector
    FlowState2D L; // Interpolated flow state just left of face.
    FlowState2D R; // interpolated flow state just right of face.
    //
    Cell2D[2] left_cells; // References to cells on the left, starting with closest.
    Cell2D[2] right_cells; // References to cells on the right, starting with closest.
    //
    // Workspace for the Osher-type flux calculator.
    GasState stateLstar, stateRstar, stateX0;

    this(GasModel gmodel, CQIndex cqi)
    {
        this.cqi = new CQIndex(cqi);
        pos = Vector3();
        F.length = cqi.n;
        // Workspace for Osher-type flux calculator.
        stateLstar = new GasState(gmodel);
        stateRstar = new GasState(gmodel);
        stateX0 = new GasState(gmodel);
    }

    this(ref const(Face2D) other)
    {
        cqi = new CQIndex(other.cqi);
        pos = Vector3(other.pos);
        F = other.F.dup;
        stateLstar = new GasState(other.stateLstar);
        stateRstar = new GasState(other.stateRstar);
        stateX0 = new GasState(other.stateX0);
    }

    override
    string toString() const
    {
        string repr = "Face2D(";
        repr ~= format("pos=%s, n=%s, t1=%s, area=%g", pos, n, t1, area);
        repr ~= format(", F=%s", F);
        repr ~= ")";
        return repr;
    }

    @nogc
    void compute_geometry(bool axiFlag)
    // Update the geometric properties from vertex data.
    {
        t1 = *p1; t1 -= *p0; t1.normalize();
        Vector3 t2 = Vector3(0.0, 0.0, 1.0);
        cross(n, t1, t2); n.normalize();
        area = distance_between(*p1, *p0);
        pos = *p0; pos += *p1; pos *= 0.5;
        if (axiFlag) { area *= pos.y; }
        return;
    }

    @nogc
    void calculate_flux(FlowState2D fsL, FlowState2D fsR, GasModel gmodel)
    // Compute the face's flux vector from left and right flow states.
    // The core of this calculation is the one-dimensional Riemann solver
    // from the gasflow module.
    {
        Vector3 velL = Vector3(fsL.vel);
        Vector3 velR = Vector3(fsR.vel);
        velL.transform_to_local_frame(n, t1);
        velR.transform_to_local_frame(n, t1);
        double[5] rsol = osher_riemann(fsL.gas, fsR.gas, velL.x, velR.x,
                                       stateLstar, stateRstar, stateX0, gmodel);
        double rho = stateX0.rho;
        double p = stateX0.p;
        double u = gmodel.internal_energy(stateX0);
        double velx = rsol[4];
        double vely = (velx < 0.0) ? velL.y : velR.y;
        double massFlux = rho*velx;
        Vector3 momentum = Vector3(massFlux*velx+p, massFlux*vely);
        momentum.transform_to_global_frame(n, t1);
        F[cqi.mass] = massFlux;
        F[cqi.xMom] = momentum.x;
        F[cqi.yMom] = momentum.y;
        F[cqi.totEnergy] = massFlux*(u+p/rho+0.5*(velx*velx+vely*vely));
        if (cqi.n_species > 1) {
            foreach (i; 0 .. cqi.n_species) {
                F[cqi.species+i] = massFlux * ((velx < 0.0) ? fsL.gas.massf[i] : fsR.gas.massf[i]);
            }
        }
        foreach (i; 0 .. cqi.n_modes) {
            F[cqi.modes+i] = massFlux * ((velx < 0.0) ? fsL.gas.u_modes[i] : fsR.gas.u_modes[i]);
        }
        bool allFinite = true;
        foreach (e; F) { if (!isFinite(e)) { allFinite = false; } }
        if (!allFinite) {
            debug { import std.stdio;  writeln("face=", this); }
            throw new Exception("At least one flux quantity is not finite.");
        }
        return;
    } // end calculate_flux()

    @nogc
    void simple_flux(FlowState2D fs, GasModel gmodel)
    // Computes the face's flux vector from a single flow state.
    // Supersonic flow is assumed.
    {
        Vector3 vel = Vector3(fs.vel);
        vel.transform_to_local_frame(n, t1);
        double rho = fs.gas.rho;
        double p = fs.gas.p;
        double u = gmodel.internal_energy(fs.gas);
        double massFlux = rho * vel.x;
        Vector3 momentum = Vector3(massFlux*vel.x+p, massFlux*vel.y);
        momentum.transform_to_global_frame(n, t1);
        F[cqi.mass] = massFlux;
        F[cqi.xMom] = momentum.x;
        F[cqi.yMom] = momentum.y;
        F[cqi.totEnergy] = massFlux * (u+p/rho+0.5*(vel.x*vel.x+vel.y*vel.y));
        if (cqi.n_species > 1) {
            foreach (i; 0 .. cqi.n_species) {
                F[cqi.species+i] = massFlux * fs.gas.massf[i];
            }
        }
        foreach (i; 0 .. cqi.n_modes) {
            F[cqi.modes+i] = massFlux * fs.gas.u_modes[i];
        }
        bool allFinite = true;
        foreach (e; F) { if (!isFinite(e)) { allFinite = false; } }
        if (!allFinite) {
            debug { import std.stdio;  writeln("face=", this); }
            throw new Exception("At least one flux quantity is not finite.");
        }
        return;
    } // end simple_flux()

} // end class Face2D
