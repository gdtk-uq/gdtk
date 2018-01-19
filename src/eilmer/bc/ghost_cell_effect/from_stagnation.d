// from_stagnation.d

module bc.ghost_cell_effect.from_stagnation;

import std.json;
import std.string;
import std.conv;
import std.stdio;
import std.math;

import geom;
import globalconfig;
import globaldata;
import flowstate;
import fvcore;
import fvinterface;
import fvcell;
import fluidblock;
import sfluidblock;
import gas;
import bc;


class GhostCellFromStagnation : GhostCellEffect {
public:
    FlowState stagnation_condition;
    string direction_type;
    Vector3 direction_vector;
    double alpha;
    double beta;
    double mass_flux;
    double relax_factor;
private:
    FlowState inflow_condition;
    double stagnation_entropy;
    double stagnation_enthalpy;
    double p0_min;
    double p0_max;

public:    
    this(int id, int boundary, in FlowState stagnation_condition,
         string direction_type, in Vector3 vec, double alpha, double beta,
         double mass_flux, double relax_factor)
    {
        super(id, boundary, "FromStagnation");
        this.stagnation_condition = new FlowState(stagnation_condition);
        this.direction_type = direction_type;
        this.direction_vector = vec;
        this.direction_vector.normalize();
        this.alpha = alpha;
        this.beta = beta;
        this.mass_flux = mass_flux;
        this.relax_factor = relax_factor;
        auto gmodel = blk.myConfig.gmodel;
        stagnation_enthalpy = gmodel.enthalpy(stagnation_condition.gas);
        stagnation_entropy = gmodel.entropy(stagnation_condition.gas);
        inflow_condition = new FlowState(stagnation_condition);
        // The following limits are used when adjusting the stagnation pressure
        // to achieve a fixed mass flow per unit area.
        // The initially-specified stagnation condition needs to be a reasonable guess
        // because T0 will be held fixed and p0 adjusted within the following bounds.
        p0_min = 0.1 * stagnation_condition.gas.p;
        p0_max = 10.0 * stagnation_condition.gas.p;
    }

    override string toString() const 
    {
        return "FromStagnation(stagnation_condition=" ~ to!string(stagnation_condition) ~ 
            ", direction_type=" ~ direction_type ~ 
            ", direction_vector=" ~ to!string(direction_vector) ~
            ", alpha=" ~ to!string(alpha) ~ ", beta=" ~ to!string(beta) ~ 
            ", mass_flux=" ~ to!string(mass_flux) ~
            ", relax_factor=" ~ to!string(relax_factor) ~ ")";
    }

    void set_velocity_components(ref Vector3 vel, double speed, ref FVInterface face)
    {
        switch (direction_type) {
        case "uniform":
            // Given a flow direction.
            vel.set(speed*direction_vector.x, speed*direction_vector.y, speed*direction_vector.z);
            break;
        case "axial":
            // Axial-flow through a presumably circular surface.
            // [TODO] 27-Feb-2014 through 16-Oct-2015 PJ:
            // check that this fall-through is OK.
        case "radial":
            // For turbo inflow, beta sets the axial flow.
            double vz = speed * sin(beta);
            // ...and alpha sets the angle of the flow in the plane of rotation.
            double vt = speed * cos(beta) * sin(alpha);
            double vr = -speed * cos(beta) * cos(alpha);
            double x = face.pos.x;
            double y = face.pos.y;
            double rxy = sqrt(x*x + y*y);
            vel.set(vr*x/rxy - vt*y/rxy, vt*x/rxy + vr*y/rxy, vz);
            break;
        case "normal":
        default:
            // The flow direction is into the block along the local face normal.
            final switch (which_boundary) {
            case Face.north:
            case Face.east:
            case Face.top:
                // Outward-facing normal.
                vel.set(-speed*face.n.x, -speed*face.n.y, -speed*face.n.z);
                break;
            case Face.west:
            case Face.south:
            case Face.bottom:
                // Inward-facing normal.
                vel.set(speed*face.n.x, speed*face.n.y, speed*face.n.z);
            }
        }
    } // end set_velocity_components()

    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
        throw new Error("GhostCellFromStagnation.apply_unstructured_grid() not yet implemented");
    }

    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        size_t i, j, k;
        FVCell src_cell, dest_cell;
        FVInterface face;
        auto gmodel = blk.myConfig.gmodel;
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");

        double p_stag = 0.0;
        double T_stag = 0.0; // temporary

        final switch (which_boundary) {
        case Face.north:
            j = blk.jmax;
            // First, estimate the current bulk inflow condition.
            double area = 0.0;
            double rhoUA = 0.0; // current mass_flux through boundary
            double rhovxA = 0.0; // mass-weighted x-velocity
            double rhovyA = 0.0;
            double rhovzA = 0.0;
            double rhoA = 0.0;
            double pA = 0.0;
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (i = blk.imin; i <= blk.imax; ++i) {
                    auto cell = blk.get_cell(i,j,k);
                    face = cell.iface[Face.north];
                    area += face.area[0];
                    double local_rhoA = cell.fs.gas.rho * face.area[0];
                    rhoA += local_rhoA;
                    // North faces have unit normals that nominally point outward from the domain, hence '-='
                    rhoUA -= local_rhoA * dot(cell.fs.vel, face.n); // mass flux
                    rhovxA += local_rhoA * cell.fs.vel.x;
                    rhovyA += local_rhoA * cell.fs.vel.y;
                    rhovzA += local_rhoA * cell.fs.vel.z;
                    pA += cell.fs.gas.p * face.area[0];
                } // end i loop
            } // end k loop
            if (mass_flux > 0.0 && ftl == 0) {
                // Adjust the stagnation pressure to better achieve the specified mass flux.
                // Note that we only do this adjustment once, at the start of a
                // multi-level gas-dynamic update.
                double p = pA / area;
                double dp_over_p = relax_factor * 0.5 / (rhoA/area) * 
                    (mass_flux*mass_flux - rhoUA*fabs(rhoUA)/(area*area)) / p;
                double new_p0 = (1.0 + dp_over_p) * stagnation_condition.gas.p;
                new_p0 = fmin(fmax(new_p0, p0_min), p0_max);
                stagnation_condition.gas.p = new_p0;
                gmodel.update_thermo_from_pT(stagnation_condition.gas);
                stagnation_enthalpy = gmodel.enthalpy(stagnation_condition.gas);
                stagnation_entropy = gmodel.entropy(stagnation_condition.gas);
            }
            double bulk_speed = sqrt((rhovxA/rhoA)^^2 + (rhovyA/rhoA)^^2 + (rhovzA/rhoA)^^2);
            if (rhoUA < 0.0) { bulk_speed = 0.0; } // Block any outflow with stagnation condition.
            // Assume an isentropic process from a known total enthalpy.
            double enthalpy = stagnation_enthalpy - 0.5 * bulk_speed^^2;
            gmodel.update_thermo_from_hs(inflow_condition.gas, enthalpy, stagnation_entropy);
            // Now, apply the ghost-cell conditions
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (i = blk.imin; i <= blk.imax; ++i) {
                    src_cell = blk.get_cell(i,j,k);
                    face = src_cell.iface[Face.north];
                    // Velocity components may vary with position on the block face.
                    set_velocity_components(inflow_condition.vel, bulk_speed, face);
                    dest_cell = blk.get_cell(i,j+1,k);
                    dest_cell.fs.copy_values_from(inflow_condition);
                    dest_cell = blk.get_cell(i,j+2,k);
                    dest_cell.fs.copy_values_from(inflow_condition);
                } // end i loop
            } // end k loop
            break;
        case Face.east:
            i = blk.imax;
            // First, estimate the current bulk inflow condition.
            double area = 0.0;
            double rhoUA = 0.0; // current mass_flux through boundary
            double rhovxA = 0.0; // mass-weighted x-velocity
            double rhovyA = 0.0;
            double rhovzA = 0.0;
            double rhoA = 0.0;
            double pA = 0.0;
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (j = blk.jmin; j <= blk.jmax; ++j) {
                    auto cell = blk.get_cell(i,j,k);
                    face = cell.iface[Face.east];
                    area += face.area[0];
                    double local_rhoA = cell.fs.gas.rho * face.area[0];
                    rhoA += local_rhoA;
                    // East interfaces have normal that points outwards, hence '-='
                    rhoUA -= local_rhoA * dot(cell.fs.vel, face.n);
                    rhovxA += local_rhoA * cell.fs.vel.x;
                    rhovyA += local_rhoA * cell.fs.vel.y;
                    rhovzA += local_rhoA * cell.fs.vel.z;
                    pA += cell.fs.gas.p * face.area[0];
                } // end j loop
            } // end k loop
            if (mass_flux > 0.0 && ftl == 0) {
                // Adjust the stagnation pressure to better achieve the specified mass flux.
                // Note that we only do this adjustment once, at the start of a
                // multi-level gas-dynamic update.
                double p = pA / area;
                double dp_over_p = relax_factor * 0.5 / (rhoA/area) * 
                    (mass_flux*mass_flux - rhoUA*fabs(rhoUA)/(area*area)) / p;
                double new_p0 = (1.0 + dp_over_p) * stagnation_condition.gas.p;
                new_p0 = fmin(fmax(new_p0, p0_min), p0_max);
                stagnation_condition.gas.p = new_p0;
                gmodel.update_thermo_from_pT(stagnation_condition.gas);
                stagnation_enthalpy = gmodel.enthalpy(stagnation_condition.gas);
                stagnation_entropy = gmodel.entropy(stagnation_condition.gas);
            }
            double bulk_speed = sqrt((rhovxA/rhoA)^^2 + (rhovyA/rhoA)^^2 + (rhovzA/rhoA)^^2);
            if (rhoUA < 0.0) { bulk_speed = 0.0; } // Block any outflow with stagnation condition.
            // Assume an isentropic process from a known total enthalpy.
            double enthalpy = stagnation_enthalpy - 0.5 * bulk_speed^^2;
            gmodel.update_thermo_from_hs(inflow_condition.gas, enthalpy, stagnation_entropy);
            // Now, apply the ghost-cell conditions
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (j = blk.jmin; j <= blk.jmax; ++j) {
                    src_cell = blk.get_cell(i,j,k);
                    face = src_cell.iface[Face.east];
                    // Velocity components may vary with position on the block face.
                    set_velocity_components(inflow_condition.vel, bulk_speed, face);
                    dest_cell = blk.get_cell(i+1,j,k);
                    dest_cell.fs.copy_values_from(inflow_condition);
                    dest_cell = blk.get_cell(i+2,j,k);
                    dest_cell.fs.copy_values_from(inflow_condition);
                } // end j loop
            } // end k loop
            break;
        case Face.south:
            j = blk.jmin;
            // First, estimate the current bulk inflow condition.
            double area = 0.0;
            double rhoUA = 0.0; // current mass_flux through boundary
            double rhovxA = 0.0; // mass-weighted x-velocity
            double rhovyA = 0.0;
            double rhovzA = 0.0;
            double rhoA = 0.0;
            double pA = 0.0;
            for (k = blk.kmin; k <= blk.kmax; ++k) {
            for (i = blk.imin; i <= blk.imax; ++i) {
                    auto cell = blk.get_cell(i,j,k);
                    face = cell.iface[Face.south];
                    area += face.area[0];
                    double local_rhoA = cell.fs.gas.rho * face.area[0];
                    rhoA += local_rhoA;
                    rhoUA += local_rhoA * dot(cell.fs.vel, face.n);
                    rhovxA += local_rhoA * cell.fs.vel.x;
                    rhovyA += local_rhoA * cell.fs.vel.y;
                    rhovzA += local_rhoA * cell.fs.vel.z;
                    pA += cell.fs.gas.p * face.area[0];
                } // end i loop
            } // end k loop
            if ( mass_flux > 0.0 && ftl == 0 ) {
                // Adjust the stagnation pressure to better achieve the specified mass flux.
                // Note that we only do this adjustment once, at the start of a
                // multi-level gas-dynamic update.
                double p = pA / area;
                double dp_over_p = relax_factor * 0.5 / (rhoA/area) * 
                    (mass_flux*mass_flux - rhoUA*fabs(rhoUA)/(area*area)) / p;
                double new_p0 = (1.0 + dp_over_p) * stagnation_condition.gas.p;
                new_p0 = fmin(fmax(new_p0, p0_min), p0_max);
                stagnation_condition.gas.p = new_p0;
                gmodel.update_thermo_from_pT(stagnation_condition.gas);
                stagnation_enthalpy = gmodel.enthalpy(stagnation_condition.gas);
                stagnation_entropy = gmodel.entropy(stagnation_condition.gas);
            }
            double bulk_speed = sqrt((rhovxA/rhoA)^^2 + (rhovyA/rhoA)^^2 + (rhovzA/rhoA)^^2);
            if (rhoUA < 0.0) { bulk_speed = 0.0; } // Block any outflow with stagnation condition.
            // Assume an isentropic process from a known total enthalpy.
            double enthalpy = stagnation_enthalpy - 0.5 * bulk_speed^^2;
            gmodel.update_thermo_from_hs(inflow_condition.gas, enthalpy, stagnation_entropy);
            // Now, apply the ghost-cell conditions
            for (k = blk.kmin; k <= blk.kmax; ++k) {
            for (i = blk.imin; i <= blk.imax; ++i) {
                    src_cell = blk.get_cell(i,j,k);
                    face = src_cell.iface[Face.south];
                    // Velocity components may vary with position on the block face.
                    set_velocity_components(inflow_condition.vel, bulk_speed, face);
                    dest_cell = blk.get_cell(i,j-1,k);
                    dest_cell.fs.copy_values_from(inflow_condition);
                    dest_cell = blk.get_cell(i,j-2,k);
                    dest_cell.fs.copy_values_from(inflow_condition);
                } // end i loop
            } // end k loop
            break;
        case Face.west:
            i = blk.imin;
            // First, estimate the current bulk inflow condition.
            double area = 0.0;
            double rhoUA = 0.0; // current mass_flux through boundary
            double rhovxA = 0.0; // mass-weighted x-velocity
            double rhovyA = 0.0;
            double rhovzA = 0.0;
            double rhoA = 0.0;
            double pA = 0.0;
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (j = blk.jmin; j <= blk.jmax; ++j) {
                    auto cell = blk.get_cell(i,j,k);
                    face = cell.iface[Face.west];
                    area += face.area[0];
                    double local_rhoA = cell.fs.gas.rho * face.area[0];
                    rhoA += local_rhoA;
                    // West faces have unit normals that nominally point inward to the domain.
                    rhoUA += local_rhoA * dot(cell.fs.vel, face.n);
                    rhovxA += local_rhoA * cell.fs.vel.x;
                    rhovyA += local_rhoA * cell.fs.vel.y;
                    rhovzA += local_rhoA * cell.fs.vel.z;
                    pA += cell.fs.gas.p * face.area[0];
                } // end j loop
            } // end k loop
            if (mass_flux > 0.0 && ftl == 0) {
                // Adjust the stagnation pressure to better achieve the specified mass flux.
                // Note that we only do this adjustment once, at the start of a
                // multi-level gas-dynamic update.
                double p = pA / area;
                double dp_over_p = relax_factor * 0.5 / (rhoA/area) * 
                    (mass_flux*mass_flux - rhoUA*fabs(rhoUA)/(area*area)) / p;
                double new_p0 = (1.0 + dp_over_p) * stagnation_condition.gas.p;
                new_p0 = fmin(fmax(new_p0, p0_min), p0_max);
                stagnation_condition.gas.p = new_p0;
                gmodel.update_thermo_from_pT(stagnation_condition.gas);
                stagnation_enthalpy = gmodel.enthalpy(stagnation_condition.gas);
                stagnation_entropy = gmodel.entropy(stagnation_condition.gas);
            }
            double bulk_speed = sqrt((rhovxA/rhoA)^^2 + (rhovyA/rhoA)^^2 + (rhovzA/rhoA)^^2);
            if (rhoUA < 0.0) { bulk_speed = 0.0; } // Block any outflow with stagnation condition.
            // Assume an isentropic process from a known total enthalpy.
            double enthalpy = stagnation_enthalpy - 0.5 * bulk_speed^^2;
            gmodel.update_thermo_from_hs(inflow_condition.gas, enthalpy, stagnation_entropy);
            // Now, apply the ghost-cell conditions
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (j = blk.jmin; j <= blk.jmax; ++j) {
                    src_cell = blk.get_cell(i,j,k);
                    face = src_cell.iface[Face.west];
                    // Velocity components may vary with position on the block face.
                    set_velocity_components(inflow_condition.vel, bulk_speed, face);
                    dest_cell = blk.get_cell(i-1,j,k);
                    dest_cell.fs.copy_values_from(inflow_condition);
                    dest_cell = blk.get_cell(i-2,j,k);
                    dest_cell.fs.copy_values_from(inflow_condition);
                } // end j loop
            } // end k loop
            break;
        case Face.top:
            k = blk.kmax;
            // First, estimate the current bulk inflow condition.
            double area = 0.0;
            double rhoUA = 0.0; // current mass_flux through boundary
            double rhovxA = 0.0; // mass-weighted x-velocity
            double rhovyA = 0.0;
            double rhovzA = 0.0;
            double rhoA = 0.0;
            double pA = 0.0;
            for (j = blk.jmin; j <= blk.jmax; ++j) {
            for (i = blk.imin; i <= blk.imax; ++i) {
                    auto cell = blk.get_cell(i,j,k);
                    face = cell.iface[Face.top];
                    area += face.area[0];
                    double local_rhoA = cell.fs.gas.rho * face.area[0];
                    rhoA += local_rhoA;
                    // Top interfaces have normal that points outwards, hence '-='
                    rhoUA -= local_rhoA * dot(cell.fs.vel, face.n);
                    rhovxA += local_rhoA * cell.fs.vel.x;
                    rhovyA += local_rhoA * cell.fs.vel.y;
                    rhovzA += local_rhoA * cell.fs.vel.z;
                    pA += cell.fs.gas.p * face.area[0];
                } // end i loop
            } // end j loop
            if (mass_flux > 0.0 && ftl == 0) {
                // Adjust the stagnation pressure to better achieve the specified mass flux.
                // Note that we only do this adjustment once, at the start of a
                // multi-level gas-dynamic update.
                double p = pA / area;
                double dp_over_p = relax_factor * 0.5 / (rhoA/area) * 
                    (mass_flux*mass_flux - rhoUA*fabs(rhoUA)/(area*area)) / p;
                double new_p0 = (1.0 + dp_over_p) * stagnation_condition.gas.p;
                new_p0 = fmin(fmax(new_p0, p0_min), p0_max);
                stagnation_condition.gas.p = new_p0;
                gmodel.update_thermo_from_pT(stagnation_condition.gas);
                stagnation_enthalpy = gmodel.enthalpy(stagnation_condition.gas);
                stagnation_entropy = gmodel.entropy(stagnation_condition.gas);
            }
            double bulk_speed = sqrt((rhovxA/rhoA)^^2 + (rhovyA/rhoA)^^2 + (rhovzA/rhoA)^^2);
            if (rhoUA < 0.0) { bulk_speed = 0.0; } // Block any outflow with stagnation condition.
            // Assume an isentropic process from a known total enthalpy.
            double enthalpy = stagnation_enthalpy - 0.5 * bulk_speed^^2;
            gmodel.update_thermo_from_hs(inflow_condition.gas, enthalpy, stagnation_entropy);
            // Now, apply the ghost-cell conditions
            for (j = blk.jmin; j <= blk.jmax; ++j) {
            for (i = blk.imin; i <= blk.imax; ++i) {
                    src_cell = blk.get_cell(i,j,k);
                    face = src_cell.iface[Face.top];
                    // Velocity components may vary with position on the block face.
                    set_velocity_components(inflow_condition.vel, bulk_speed, face);
                    dest_cell = blk.get_cell(i,j,k+1);
                    dest_cell.fs.copy_values_from(inflow_condition);
                    dest_cell = blk.get_cell(i,j,k+2);
                    dest_cell.fs.copy_values_from(inflow_condition);
                } // end i loop
            } // end j loop
            break;
        case Face.bottom:
            k = blk.kmin;
            // First, estimate the current bulk inflow condition.
            double area = 0.0;
            double rhoUA = 0.0; // current mass_flux through boundary
            double rhovxA = 0.0; // mass-weighted x-velocity
            double rhovyA = 0.0;
            double rhovzA = 0.0;
            double rhoA = 0.0;
            double pA = 0.0;
            for (j = blk.jmin; j <= blk.jmax; ++j) {
            for (i = blk.imin; i <= blk.imax; ++i) {
                    auto cell = blk.get_cell(i,j,k);
                    face = cell.iface[Face.bottom];
                    area += face.area[0];
                    double local_rhoA = cell.fs.gas.rho * face.area[0];
                    rhoA += local_rhoA;
                    rhoUA += local_rhoA * dot(cell.fs.vel, face.n);
                    rhovxA += local_rhoA * cell.fs.vel.x;
                    rhovyA += local_rhoA * cell.fs.vel.y;
                    rhovzA += local_rhoA * cell.fs.vel.z;
                    pA += cell.fs.gas.p * face.area[0];
                } // end i loop
            } // end j loop
            if (mass_flux > 0.0 && ftl == 0) {
                // Adjust the stagnation pressure to better achieve the specified mass flux.
                // Note that we only do this adjustment once, at the start of a
                // multi-level gas-dynamic update.
                double p = pA / area;
                double dp_over_p = relax_factor * 0.5 / (rhoA/area) * 
                    (mass_flux*mass_flux - rhoUA*fabs(rhoUA)/(area*area)) / p;
                double new_p0 = (1.0 + dp_over_p) * stagnation_condition.gas.p;
                new_p0 = fmin(fmax(new_p0, p0_min), p0_max);
                stagnation_condition.gas.p = new_p0;
                gmodel.update_thermo_from_pT(stagnation_condition.gas);
                stagnation_enthalpy = gmodel.enthalpy(stagnation_condition.gas);
                stagnation_entropy = gmodel.entropy(stagnation_condition.gas);
            }
            double bulk_speed = sqrt((rhovxA/rhoA)^^2 + (rhovyA/rhoA)^^2 + (rhovzA/rhoA)^^2);
            if (rhoUA < 0.0) { bulk_speed = 0.0; } // Block any outflow with stagnation condition.
            // Assume an isentropic process from a known total enthalpy.
            double enthalpy = stagnation_enthalpy - 0.5 * bulk_speed^^2;
            gmodel.update_thermo_from_hs(inflow_condition.gas, enthalpy, stagnation_entropy);
            // Now, apply the ghost-cell conditions
            for (j = blk.jmin; j <= blk.jmax; ++j) {
            for (i = blk.imin; i <= blk.imax; ++i) {
                    src_cell = blk.get_cell(i,j,k);
                    face = src_cell.iface[Face.bottom];
                    // Velocity components may vary with position on the block face.
                    set_velocity_components(inflow_condition.vel, bulk_speed, face);
                    dest_cell = blk.get_cell(i,j,k-1);
                    dest_cell.fs.copy_values_from(inflow_condition);
                    dest_cell = blk.get_cell(i,j,k-2);
                    dest_cell.fs.copy_values_from(inflow_condition);
                } // end i loop
            } // end j loop
            break;
        } // end switch
    } // end apply()
} // end class GhostCellFixedStagnationPT
