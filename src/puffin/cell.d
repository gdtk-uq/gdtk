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
        repr ~= ")";
        return repr;
    }

    @nogc
    void compute_geometry()
    // Update the geometric properties from vertex data.
    {
        xyplane_quad_cell_properties(*p00, *p10, *p11, *p01, pos, xyplane_area,
                                     iLen, jLen, kLen);
        volume = xyplane_area * ((Config.axisymmetric) ? pos.y : 1.0);
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
        double u_other = 0.0;
        foreach(i; 0 .. cqi.n_modes) { fs.gas.u_modes[i] = myU[cqi.modes+i] * dinv; }
        foreach(ei; fs.gas.u_modes) { u_other += ei; }
        fs.gas.u = u - u_other;
        // Thermochemical species, if appropriate.
        if (cqi.n_species > 1) {
            foreach(i; 0 .. cqi.n_species) { fs.gas.massf[i] = myU[cqi.species+i] * dinv; }
            scale_mass_fractions(fs.gas.massf);
        } else {
            fs.gas.massf[0] = 1.0;
        }
        gmodel.update_thermo_from_rhou(fs.gas);
        gmodel.update_sound_speed(fs.gas);
        return;
    } // end decode_conserved()

    @nogc
    void eval_dUdt(size_t ftl, bool axisymmetric)
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
        if (axisymmetric) {
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

    this(CQIndex cqi)
    {
        this.cqi = new CQIndex(cqi);
        pos = Vector3();
        F.length = cqi.n;
    }

    this(ref const(Face2D) other)
    {
        cqi = new CQIndex(other.cqi);
        pos = Vector3(other.pos);
        F = other.F.dup;
    }

    override
    string toString() const
    {
        string repr = "Face2D(";
        repr ~= format("pos=%s", pos);
        repr ~= ")";
        return repr;
    }

    @nogc
    void compute_geometry()
    // Update the geometric properties from vertex data.
    {
        n = *p1; n -= *p0; n.normalize();
        Vector3 t2 = Vector3(0.0, 0.0, 1.0);
        cross(t1, n, t2);
        area = distance_between(*p1, *p0);
        pos = *p0; pos += *p1; pos *= 0.5;
        if (Config.axisymmetric) { area *= pos.y; }
        return;
    }

    @nogc
    void calculate_flux(FlowState2D fsL, FlowState2D fsR, GasModel gmodel)
    // Compute the face's flux vector from left and right flow states.
    // The core of this calculation is a ond-dimensional Riemann solver.
    {
        debug { import std.stdio; writeln("[TODO] Face2D.calculate_flux"); }
        foreach (k; 0 .. cqi.n) { F[k] = 0.0; }
        return;
    }

    @nogc
    void simple_flux(FlowState2D fs, GasModel gmodel)
    // Computes the face's flux vector from a single flow state.
    // Supersonic flow is assumed.
    {
        debug { import std.stdio; writeln("[TODO] Face2D.simple_flux"); }
        foreach (k; 0 .. cqi.n) { F[k] = 0.0; }
        return;
    }

} // end class Face2D
