// from_stagnation.d

module bc.ghost_cell_effect.from_stagnation;

import std.json;
import std.string;
import std.conv;
import std.stdio;
import std.math;
import ntypes.complex;
import nm.number;
import util.lua;
import util.lua_service;

import geom;
import globalconfig;
import globaldata;
import flowstate;
import fvinterface;
import lmr.fluidfvcell;
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
    string luaFileName;
    number stagnation_entropy;
    number stagnation_enthalpy;
    number p0_min;
    number p0_max;

public:
    this(int id, int boundary,
         in FlowState stagnation_condition, in string fileName,
         string direction_type, in Vector3 vec, double alpha, double beta,
         double mass_flux, double relax_factor)
    {
        super(id, boundary, "FromStagnation");
        this.stagnation_condition = FlowState(stagnation_condition);
        this.luaFileName = fileName;
        this.direction_type = direction_type;
        this.direction_vector = vec;
        this.direction_vector.normalize();
        this.alpha = alpha;
        this.beta = beta;
        this.mass_flux = mass_flux;
        this.relax_factor = relax_factor;
        auto gmodel = GlobalConfig.gmodel_master;
        stagnation_enthalpy = gmodel.enthalpy(stagnation_condition.gas);
        stagnation_entropy = gmodel.entropy(stagnation_condition.gas);
        inflow_condition = FlowState(stagnation_condition);
        // The following limits are used when adjusting the stagnation pressure
        // to achieve a fixed mass flow per unit area.
        // The initially-specified stagnation condition needs to be a reasonable guess
        // because T0 will be held fixed and p0 adjusted within the following bounds.
        p0_min = 0.1 * stagnation_condition.gas.p;
        p0_max = 10.0 * stagnation_condition.gas.p;
        version(complex_numbers) {
            this.stagnation_condition.clear_imaginary_components();
            inflow_condition.clear_imaginary_components();
        }
    }

    override void post_bc_construction()
    {
        if ((luaFileName.length > 0) && (blk.bc[which_boundary].myL == null)) {
            blk.bc[which_boundary].init_lua_State(luaFileName);
        }
    }

    override string toString() const
    {
        return "FromStagnation(stagnation_condition=" ~ to!string(stagnation_condition) ~
            ", luaFileName=\"" ~ luaFileName ~ "\"" ~
            ", direction_type=" ~ direction_type ~
            ", direction_vector=" ~ to!string(direction_vector) ~
            ", alpha=" ~ to!string(alpha) ~ ", beta=" ~ to!string(beta) ~
            ", mass_flux=" ~ to!string(mass_flux) ~
            ", relax_factor=" ~ to!string(relax_factor) ~ ")";
    }

    @nogc
    void set_velocity_components(ref Vector3 vel, number speed, ref FVInterface face, int outsign)
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
            number vz = speed * sin(beta);
            // ...and alpha sets the angle of the flow in the plane of rotation.
            number vt = speed * cos(beta) * sin(alpha);
            number vr = -speed * cos(beta) * cos(alpha);
            number x = face.pos.x;
            number y = face.pos.y;
            number rxy = sqrt(x*x + y*y);
            vel.set(vr*x/rxy - vt*y/rxy, vt*x/rxy + vr*y/rxy, vz);
            break;
        case "normal":
        default:
            // The flow direction is into the block along the local face normal.
            vel.set(-outsign*speed*face.n.x, -outsign*speed*face.n.y, -outsign*speed*face.n.z);
        }
    } // end set_velocity_components()

    // not @nogc
    void determine_stagnation_condition(double t, bool is_structured_grid)
    {
        BoundaryCondition bc = blk.bc[which_boundary];
        auto gmodel = blk.myConfig.gmodel;
        if (mass_flux > 0.0) {
            // Adjust the stagnation pressure to better achieve the specified mass flux.
            // Note that we only do this adjustment once, at the start of a
            // multi-level gas-dynamic update.
            //
            // First, estimate the current bulk inflow condition.
            number area = 0.0;
            number rhoUA = 0.0; // current mass_flux through boundary
            number rhovxA = 0.0; // mass-weighted x-velocity
            number rhovyA = 0.0;
            number rhovzA = 0.0;
            number rhoA = 0.0;
            number pA = 0.0;
            FluidFVCell cell;
            //
            foreach (i, face; bc.faces) {
                int outsign = bc.outsigns[i];
                if (is_structured_grid) {
                    cell = (outsign == 1) ? face.left_cells[0] : face.right_cells[0];
                } else {
                    cell = (outsign == 1) ? face.left_cell : face.right_cell;
                }
                area += face.area[0];
                number local_rhoA = cell.fs.gas.rho * face.area[0];
                rhoA += local_rhoA;
                rhoUA -= outsign*(local_rhoA * dot(cell.fs.vel, face.n)); // mass flux
                rhovxA += local_rhoA * cell.fs.vel.x;
                rhovyA += local_rhoA * cell.fs.vel.y;
                rhovzA += local_rhoA * cell.fs.vel.z;
                pA += cell.fs.gas.p * face.area[0];
            }
            //
            // Now, make the adjustment.
            number p = pA / area;
            number dp_over_p = relax_factor * 0.5 / (rhoA/area) *
                (mass_flux*mass_flux - rhoUA*fabs(rhoUA)/(area*area)) / p;
            number new_p0 = (1.0 + dp_over_p) * stagnation_condition.gas.p;
            new_p0 = fmin(fmax(new_p0, p0_min), p0_max);
            stagnation_condition.gas.p = new_p0;
            gmodel.update_thermo_from_pT(stagnation_condition.gas);
            stagnation_enthalpy = gmodel.enthalpy(stagnation_condition.gas);
            stagnation_entropy = gmodel.entropy(stagnation_condition.gas);
        } else if (luaFileName.length > 0) {
            // Call out to the user's function to get the current stagnation condition.
            number new_p0 = stagnation_condition.gas.p;
            number new_T0 = stagnation_condition.gas.T;
            // [FIXME] dt_global=0.0 step=0
            callUDFstagnationPT(t, 0.0, 0, 0, 0, new_p0, new_T0);
            stagnation_condition.gas.p = new_p0;
            stagnation_condition.gas.T = new_T0;
            gmodel.update_thermo_from_pT(stagnation_condition.gas);
            stagnation_enthalpy = gmodel.enthalpy(stagnation_condition.gas);
            stagnation_entropy = gmodel.entropy(stagnation_condition.gas);
        }
        version(complex_numbers) {
            stagnation_condition.clear_imaginary_components();
        }
        return;
    } // end determine_stagnation_condition()

    // not @nogc
    void callUDFstagnationPT(double t, double dt_global, int step, int gtl, int ftl,
                             ref number new_p0, ref number new_T0)
    {
        // [TODO] call stagnationPT function
        // 1. Set up for calling function
        auto L = blk.bc[which_boundary].myL;
        // 1a. Place function to call at TOS
        lua_getglobal(L, "stagnationPT");
        // 1b. Then put arguments (as single table) at TOS
        lua_newtable(L);
        lua_pushnumber(L, t); lua_setfield(L, -2, "t");
        lua_pushnumber(L, dt_global); lua_setfield(L, -2, "dt");
        lua_pushinteger(L, step); lua_setfield(L, -2, "timeStep");
        lua_pushinteger(L, gtl); lua_setfield(L, -2, "gridTimeLevel");
        lua_pushinteger(L, ftl); lua_setfield(L, -2, "flowTimeLevel");
        lua_pushinteger(L, which_boundary); lua_setfield(L, -2, "boundaryId");
        // 2. Call LuaFunction and expect two tables of ghost cell flow state
        int number_args = 1;
        int number_results = 2; // p0 and T0
        if ( lua_pcall(L, number_args, number_results, 0) != 0 ) {
            luaL_error(L, "error running user-defined b.c. ghostCell function on boundaryId %d: %s\n",
                       which_boundary, lua_tostring(L, -1));
        }
        // The returned data, p0 first then T0:
        new_p0 = to!double(luaL_checknumber(L, -2));
        new_T0 = to!double(luaL_checknumber(L, -1));
        // 4. Clear stack
        lua_settop(L, 0);
    } // end callUDFstagnationPT()

    void determine_inflow_condition_at_face(FVInterface face, FluidFVCell cell, int outsign, GasModel gmodel)
    {
        number inflow_speed = -outsign * dot(cell.fs.vel, face.n);
        // Block any outflow with stagnation condition.
        if (inflow_speed < 0.0) { inflow_speed = 0.0; }
        // Assume an isentropic process from a known total enthalpy.
        number enthalpy = stagnation_enthalpy - 0.5 * inflow_speed^^2;
        gmodel.update_thermo_from_hs(inflow_condition.gas, enthalpy, stagnation_entropy);
        // Velocity components may vary with position on the block face.
        set_velocity_components(inflow_condition.vel, inflow_speed, face, outsign);
        version(complex_numbers) {
            inflow_condition.clear_imaginary_components();
        }
        return;
    } // end determine_inflow_condition_at_face()

    override void apply_for_interface_unstructured_grid(double t, int gtl, int ftl, FVInterface face)
    {
        // TODO: debug below code. KD 13-10-2020
        // Below is the first hack at extending this boundary condition to be applied on a per interface basis.
        // Currently there is something not quite right in this implementation, it fails to work with the numerical
        // Jacobian code to form an ILU preconditioner.
        BoundaryCondition bc = blk.bc[which_boundary];
        auto gmodel = blk.myConfig.gmodel;
        // The condition is restricted to the fixed stagnation inflow condition and cannot be used with
        // the specified-mass-flux condition nor the user-defined changeable condition.
        if (mass_flux > 0.0 && ftl == 0) {
            throw new Error("GhostCellFromStagnation.apply_for_interface_unstructured_grid() with fixed mass_flux not yet implemented");
        } else if (luaFileName.length > 0) {
            throw new Error("GhostCellFromStagnation.apply_for_interface_unstructured_grid() with lua UDF not yet implemented");
        }
        int outsign = bc.outsigns[face.i_bndry];
        FluidFVCell cell = (outsign == 1) ? face.left_cell : face.right_cell;
        FluidFVCell ghost0 = (outsign == 1) ? face.right_cell : face.left_cell;
        determine_inflow_condition_at_face(face, cell, outsign, gmodel);
        ghost0.fs.copy_values_from(inflow_condition);
    }

    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
        BoundaryCondition bc = blk.bc[which_boundary];
        auto gmodel = blk.myConfig.gmodel;
        if (ftl == 0) {
            determine_stagnation_condition(t, false);
        }
        // Now, apply the ghost-cell conditions
        foreach (i, face; bc.faces) {
            int outsign = bc.outsigns[i];
            FluidFVCell cell = (outsign == 1) ? face.left_cell : face.right_cell;
            FluidFVCell ghost0 = (outsign == 1) ? face.right_cell : face.left_cell;
            determine_inflow_condition_at_face(face, cell, outsign, gmodel);
            ghost0.fs.copy_values_from(inflow_condition);
        }
    } // end apply_unstructured_grid()

    override void apply_for_interface_structured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        auto gmodel = blk.myConfig.gmodel;
        auto bc = blk.bc[which_boundary];
        int outsign = bc.outsigns[f.i_bndry];
        FluidFVCell cell = (outsign == 1) ? f.left_cells[0] : f.right_cells[0];
        determine_inflow_condition_at_face(f, cell, outsign, gmodel);
        foreach (n; 0 .. blk.n_ghost_cell_layers) {
            FluidFVCell ghost = (outsign == 1) ? f.right_cells[n] : f.left_cells[n];
            ghost.fs.copy_values_from(inflow_condition);
        }
    }

    // not @nogc
    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");
        BoundaryCondition bc = blk.bc[which_boundary];
        auto gmodel = blk.myConfig.gmodel;
        if (ftl == 0) {
            determine_stagnation_condition(t, true);
        }
        foreach (i, face; bc.faces) {
            int outsign = bc.outsigns[i];
            FluidFVCell cell = (outsign == 1) ? face.left_cells[0] : face.right_cells[0];
            determine_inflow_condition_at_face(face, cell, outsign, gmodel);
            foreach (n; 0 .. blk.n_ghost_cell_layers) {
                FluidFVCell ghost = (outsign == 1) ? face.right_cells[n] : face.left_cells[n];
                ghost.fs.copy_values_from(inflow_condition);
            }
        }
    } // end apply_structured_grid()
} // end class GhostCellFixedStagnationPT
