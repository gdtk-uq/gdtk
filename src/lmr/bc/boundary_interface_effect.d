// boundary_interface_effect.d
//
// Effects needed to compute viscous fluxes and the like.
//
// PJ and RG, 2015-04-28, initial code mainly from
//    the break-up of the Fixed_T boundary condition.
//

module bc.boundary_interface_effect;

import std.json;
import std.string;
import std.conv;
import std.stdio;
import std.math;
import ntypes.complex;
import nm.number;
import nm.brent;
import nm.bracketing;

import geom;
import json_helper;
import globalconfig;
import globaldata;
import flowstate;
import fvinterface;
import lmr.fluidfvcell;
import fluidblock;
import sfluidblock;
import gas;
import bc;
import solidfvcell;
import solidfvinterface;
import kinetics.equilibrium_update;
import mass_diffusion;

BoundaryInterfaceEffect make_BIE_from_json(JSONValue jsonData, int blk_id, int boundary)
{
    string bieType = jsonData["type"].str;
    // At the point at which we call this function, we may be inside the block-constructor.
    // Don't attempt the use the block-owned gas model.
    auto gmodel = GlobalConfig.gmodel_master;
    // If we need access to a gas model in here,
    // be sure to use GlobalConfig.gmodel_master.
    BoundaryInterfaceEffect newBIE;
    switch (bieType) {
    case "copy_cell_data":
        newBIE = new BIE_CopyCellData(blk_id, boundary);
        break;
    case "flow_state_copy_to_interface":
        auto flowstate = FlowState(jsonData["flowstate"], gmodel);
        double x0 = getJSONdouble(jsonData, "x0", 0.0);
        double y0 = getJSONdouble(jsonData, "y0", 0.0);
        double z0 = getJSONdouble(jsonData, "z0", 0.0);
        double r = getJSONdouble(jsonData, "r", 0.0);
        newBIE = new BIE_FlowStateCopy(blk_id, boundary, flowstate, x0, y0, z0, r);
        break;
    case "flow_state_copy_from_profile_to_interface":
        string fname = getJSONstring(jsonData, "filename", "");
        string match = getJSONstring(jsonData, "match", "xyz");
        newBIE = new BIE_FlowStateCopyFromProfile(blk_id, boundary, fname, match);
        break;
    case "flow_state_copy_from_history_to_interface":
        string fname = getJSONstring(jsonData, "filename", "");
        newBIE = new BIE_FlowStateCopyFromHistory(blk_id, boundary, fname);
        break;
    case "synthesise_flow_state_to_interface":
        string fname = getJSONstring(jsonData, "filename", "");
        newBIE = new BIE_SynthesiseFlowState(blk_id, boundary, fname);
        break;
    case "zero_velocity":
        newBIE = new BIE_ZeroVelocity(blk_id, boundary);
        break;
    case "zero_slip_wall_velocity":
        newBIE = new BIE_ZeroSlipWallVelocity(blk_id, boundary);
        break;
    case "translating_surface":
        Vector3 v_trans = getJSONVector3(jsonData, "v_trans", Vector3(0.0,0.0,0.0));
        newBIE = new BIE_TranslatingSurface(blk_id, boundary, v_trans);
        break;
    case "rotating_surface":
        Vector3 r_omega = getJSONVector3(jsonData, "r_omega", Vector3(0.0,0.0,0.0));
        Vector3 centre = getJSONVector3(jsonData, "centre", Vector3(0.0,0.0,0.0));
        newBIE = new BIE_RotatingSurface(blk_id, boundary, r_omega, centre);
        break;
    case "fixed_temperature":
        double Twall = getJSONdouble(jsonData, "Twall", 300.0);
        newBIE = new BIE_FixedT(blk_id, boundary, Twall);
        break;
    case "fixed_composition":
        double[] massfAtWall = getJSONdoublearray(jsonData, "wall_massf_composition", [1.0,]);
        newBIE = new BIE_FixedComposition(blk_id, boundary, massfAtWall);
        break;
    case "update_thermo_trans_coeffs":
        string thermoUpdate = getJSONstring(jsonData, "thermoUpdate", "pT");
        InterpolateOption thermoInterpOpt;
        if (thermoUpdate == "global") {
            thermoInterpOpt = GlobalConfig.thermo_interpolator;
        }
        else {
            thermoInterpOpt = thermo_interpolator_from_name(thermoUpdate);
        }
        newBIE = new BIE_UpdateThermoTransCoeffs(blk_id, boundary, thermoInterpOpt);
        break;
    case "wall_turbulent":
        newBIE = new BIE_WallTurbulent(blk_id, boundary);
        break;
    case "wall_function_interface_effect":
        newBIE = new BIE_WallFunction(blk_id, boundary);
        break;
    case "adiabatic_wall_function_interface_effect":
        newBIE = new BIE_AdiabaticWallFunction(blk_id, boundary);
        break;
    case "temperature_from_gas_solid_interface":
        int otherBlock = getJSONint(jsonData, "other_block", -1);
        string otherFaceName = getJSONstring(jsonData, "other_face", "none");
        int neighbourOrientation = getJSONint(jsonData, "neighbour_orientation", 0);
        newBIE = new BIE_TemperatureFromGasSolidInterface(blk_id, boundary,
                                                          otherBlock, face_index(otherFaceName),
                                                          neighbourOrientation);
        break;
    case "temperature_from_gas_solid_interface2":
        int otherBlock = getJSONint(jsonData, "other_block", -1);
        string otherFaceName = getJSONstring(jsonData, "other_face", "none");
        int neighbourOrientation = getJSONint(jsonData, "neighbour_orientation", 0);
        newBIE = new BIE_TemperatureFromGasSolidInterface2(blk_id, boundary,
                                                           otherBlock, face_index(otherFaceName),
                                                           neighbourOrientation);
        break;
    case "thermionic_radiative_equilibrium":
        double emissivity = getJSONdouble(jsonData, "emissivity", 0.0);
        double Ar = getJSONdouble(jsonData, "Ar", 0.0);
        double phi = getJSONdouble(jsonData, "phi", 0.0);
        int ThermionicEmissionActive = getJSONint(jsonData, "ThermionicEmissionActive", 1);
        size_t nsp = GlobalConfig.gmodel_master.n_species;
        string catalytic_type = getJSONstring(jsonData, "catalytic_type", "none");
        double[] massfAtWall = getJSONdoublearray(jsonData, "wall_massf_composition", []);
        newBIE = new BIE_ThermionicRadiativeEquilibrium(blk_id, boundary, emissivity, Ar, phi,
                             ThermionicEmissionActive, nsp, GlobalConfig.gas_model_file, catalytic_type, massfAtWall);
        break;
    case "equilibrium_composition":
        newBIE = new BIE_EquilibriumComposition(blk_id, boundary, GlobalConfig.gas_model_file);
        break;
    case "user_defined":
        string fname = getJSONstring(jsonData, "filename", "none");
        newBIE = new BIE_UserDefined(blk_id, boundary, fname);
        break;
    default:
        string errMsg = format("ERROR: The BoundaryInterfaceEffect type: '%s' is unknown.", bieType);
        throw new FlowSolverException(errMsg);
    }
    return newBIE;
} // end make_BIE_from_json()


class BoundaryInterfaceEffect {
public:
    FluidBlock blk;
    int which_boundary;
    string desc;

    this(int id, int boundary, string description)
    {
        blk = cast(FluidBlock) globalBlocks[id];
        assert(blk !is null, "Oops, this should be a FluidBlock object.");
        which_boundary = boundary;
        desc = description;
    }
    void post_bc_construction() {}
    override string toString() const
    {
        return "BoundaryInterfaceEffect()";
    }
    void apply(double t, int gtl, int ftl)
    {
        final switch (blk.grid_type) {
        case Grid_t.unstructured_grid:
            apply_unstructured_grid(t, gtl, ftl);
            break;
        case Grid_t.structured_grid:
            apply_structured_grid(t, gtl, ftl);
        }
    }
    void apply_for_interface(double t, int gtl, int ftl, FVInterface f)
    {
        final switch (blk.grid_type) {
        case Grid_t.unstructured_grid:
            apply_for_interface_unstructured_grid(t, gtl, ftl, f);
            break;
        case Grid_t.structured_grid:
            apply_for_interface_structured_grid(t, gtl, ftl, f);
            break;
        }
    }
    abstract void apply_for_interface_unstructured_grid(double t, int gtl, int ftl, FVInterface f);
    abstract void apply_unstructured_grid(double t, int gtl, int ftl);
    abstract void apply_for_interface_structured_grid(double t, int gtl, int ftl, FVInterface f);
    abstract void apply_structured_grid(double t, int gtl, int ftl);
} // end class BoundaryInterfaceEffect

class BIE_DoNothing : BoundaryInterfaceEffect {
    this(int id, int boundary) {
        super(id, boundary, "DoNothing");
    }

    override string toString() const
    {
        return "DoNothing()";
    }
    // If these methods look like they do nothing, then that's correct.
    // The intention is to do nothing.
    override void apply_for_interface_unstructured_grid(double t, int gtl, int ftl, FVInterface f) {}
    override void apply_unstructured_grid(double t, int gtl, int ftl) {}
    override void apply_for_interface_structured_grid(double t, int gtl, int ftl, FVInterface f) {}
    override void apply_structured_grid(double t, int gtl, int ftl) {}
}

class BIE_CopyCellData : BoundaryInterfaceEffect {
    this(int id, int boundary, double Twall=300.0)
    {
        super(id, boundary, "CopyCellData");
    }

    override string toString() const
    {
        return "CopyCellData()";
    }

    override void apply_for_interface_unstructured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        BoundaryCondition bc = blk.bc[which_boundary];
	if (bc.outsigns[f.i_bndry] == 1) {
	    f.fs.copy_values_from(f.left_cell.fs);
	} else {
	    f.fs.copy_values_from(f.right_cell.fs);
	}
    }

    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
        BoundaryCondition bc = blk.bc[which_boundary];
        foreach (i, f; bc.faces) {
            auto c = (bc.outsigns[i] == 1) ? f.left_cell : f.right_cell;
            f.fs.copy_values_from(c.fs);
        }
    }

    override void apply_for_interface_structured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        auto gmodel = blk.myConfig.gmodel;
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");
        BoundaryCondition bc = blk.bc[which_boundary];
        auto c = (bc.outsigns[f.i_bndry] == 1) ? f.left_cells[0] : f.right_cells[0];
        f.fs.copy_values_from(c.fs);
    }

    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        auto gmodel = blk.myConfig.gmodel;
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");
        BoundaryCondition bc = blk.bc[which_boundary];
        foreach (i, f; bc.faces) {
            auto c = (bc.outsigns[i] == 1) ? f.left_cells[0] : f.right_cells[0];
            f.fs.copy_values_from(c.fs);
        }
    } // end apply_structured_grid()
} // end class BIE_CopyCellData

class BIE_FlowStateCopy : BoundaryInterfaceEffect {
public:
    FlowState fstate;
    SourceFlow sflow;
    double x0, y0, z0, r; // conical-flow parameters

    this(int id, int boundary, in FlowState fstate, double x0, double y0, double z0, double r)
    {
        super(id, boundary, "flowStateCopy");
        // We only need to gather the freestream values once at the start of simulation.
        // Note that, at this time, the gmodel held by the block is not available.
        auto gmodel = GlobalConfig.gmodel_master;
        this.fstate = fstate.dup();
        this.x0 = x0;
        this.y0 = y0;
        this.z0 = z0;
        this.r = r;
        sflow = new SourceFlow(gmodel, fstate, r);
    }

    override string toString() const
    {
        return "flowStateCopy(fstate=" ~ to!string(fstate) ~ ")";
    }

    override void apply_for_interface_unstructured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        BoundaryCondition bc = blk.bc[which_boundary];
        f.fs.copy_values_from(fstate);
        if (r > 0.0) { compute_source_flow(f); }
    }

    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
        BoundaryCondition bc = blk.bc[which_boundary];
        foreach (i, f; bc.faces) {
            f.fs.copy_values_from(fstate);
            if (r > 0.0) { compute_source_flow(f); }
        }
    }

    override void apply_for_interface_structured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        auto gmodel = blk.myConfig.gmodel;
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");
        BoundaryCondition bc = blk.bc[which_boundary];
        f.fs.copy_values_from(fstate);
        if (r > 0.0) { compute_source_flow(f); }
    }

    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        auto gmodel = blk.myConfig.gmodel;
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");
        BoundaryCondition bc = blk.bc[which_boundary];
        foreach (i, f; bc.faces) {
            f.fs.copy_values_from(fstate);
            if (r > 0.0) { compute_source_flow(f); }
        }
    } // end apply_structured_grid()

private:
    @nogc
    void compute_source_flow(FVInterface f)
    {
        auto cqi = blk.myConfig.cqi;
        // Start by assuming uniform, parallel flow.
        number p = fstate.gas.p;
        number rho = fstate.gas.rho;
        auto gmodel = blk.myConfig.gmodel;
        number u = gmodel.internal_energy(fstate.gas);
        number velx = fstate.vel.x;
        number vely = fstate.vel.y;
        number velz = (cqi.threeD) ? fstate.vel.z : to!number(0.0);
        // (Approximate) conical inflow by adding increments.
        double dx = f.pos.x.re - x0;
        double dy = f.pos.y.re - y0;
        double dz = (cqi.threeD) ? f.pos.z.re - z0 : 0.0;
        double hypot = sqrt(dx*dx + dy*dy + dz*dz);
        double[4] deltas = sflow.get_rho_v_p_u_increments(hypot-r);
        rho += deltas[0];
        double v = fstate.vel.x.re + deltas[1];
        velx = v * dx/hypot;
        vely = v * dy/hypot;
        velz = (cqi.threeD) ? v * dz/hypot : 0.0;
        p += deltas[2];
        u += deltas[3];
        // Put into interface FlowState object.
        FlowState* fs = f.fs;
        fs.gas.p = p; // not really needed because we use rhou
        fs.gas.rho = rho;
        fs.gas.u = u;
        gmodel.update_thermo_from_rhou(fs.gas);
        fs.vel.set(velx, vely, velz);
    } // end compute_source_flow()

} // end class BIE_FlowStateCopy


class BIE_FlowStateCopyFromProfile : BoundaryInterfaceEffect {
public:
    this(int id, int boundary, string fileName, string match)
    {
        super(id, boundary, "flowStateCopyFromProfile");
        fprofile = new FlowProfile(fileName, match);
    }

    override string toString() const
    {
        return format("flowStateCopyFromProfile(filename=\"%s\", match=\"%s\")",
                      fprofile.fileName, fprofile.posMatch);
    }

    //@nogc
    override void apply_for_interface_unstructured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        BoundaryCondition bc = blk.bc[which_boundary];
        f.fs.copy_values_from(fprofile.get_flowstate(f.id, f.pos));
        fprofile.adjust_velocity(f.fs, f.pos, blk.omegaz);
    }

    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
        BoundaryCondition bc = blk.bc[which_boundary];
        foreach (i, f; bc.faces) {
            f.fs.copy_values_from(fprofile.get_flowstate(f.id, f.pos));
            fprofile.adjust_velocity(f.fs, f.pos, blk.omegaz);
        }
    }

    override void apply_for_interface_structured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        auto gmodel = blk.myConfig.gmodel;
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");
        BoundaryCondition bc = blk.bc[which_boundary];
        f.fs.copy_values_from(fprofile.get_flowstate(f.id, f.pos));
        fprofile.adjust_velocity(f.fs, f.pos, blk.omegaz);
    }

    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        auto gmodel = blk.myConfig.gmodel;
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");
        BoundaryCondition bc = blk.bc[which_boundary];
        foreach (i, f; bc.faces) {
            f.fs.copy_values_from(fprofile.get_flowstate(f.id, f.pos));
            fprofile.adjust_velocity(f.fs, f.pos, blk.omegaz);
        }
    } // end apply_structured_grid()

private:
    FlowProfile fprofile;

} // end class BIE_FlowStateCopyFromProfile


class BIE_FlowStateCopyFromHistory : BoundaryInterfaceEffect {
public:
    this(int id, int boundary, string fileName)
    {
        super(id, boundary, "flowStateCopyFromHistory");
        fhistory = new FlowHistory(fileName);
        fstate = FlowState(GlobalConfig.gmodel_master, GlobalConfig.turb_model.nturb);
    }

    override string toString() const
    {
        return format("flowStateCopyFromHistory(filename=\"%s\")", fhistory.fileName);
    }

    @nogc
    override void apply_for_interface_unstructured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        throw new Error("BIE_FlowStateCopyFromHistory.apply_for_interface_unstructured_grid() not yet implemented");
    }

    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
        BoundaryCondition bc = blk.bc[which_boundary];
        auto gmodel = blk.myConfig.gmodel;
        fhistory.set_flowstate(fstate, t, gmodel);
        foreach (i, f; bc.faces) { f.fs.copy_values_from(fstate); }
    }

    @nogc
    override void apply_for_interface_structured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        throw new Error("BIE_FlowStateCopyFromHistory.apply_for_interface_structured_grid() not yet implemented");
    }

    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");
        auto gmodel = blk.myConfig.gmodel;
        fhistory.set_flowstate(fstate, t, gmodel);
        BoundaryCondition bc = blk.bc[which_boundary];
        foreach (i, f; bc.faces) {
            f.fs.copy_values_from(fstate);
        }
    } // end apply_structured_grid()

private:
    FlowHistory fhistory;
    FlowState fstate;
} // end class BIE_FlowStateCopyFromHistory


class BIE_SynthesiseFlowState : BoundaryInterfaceEffect {
public:
    this(int id, int boundary, string fileName)
    {
        super(id, boundary, "synthesiseFlowState");
        sfs = new SyntheticFlowState(fileName);
    }

    override string toString() const
    {
        return format("synthesiseFlowState(filename=\"%s\")", sfs.fileName);
    }

    @nogc
    override void apply_for_interface_unstructured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        throw new Error("BIE_SynthesiseFlowState.apply_for_interface_unstructured_grid() not yet implemented");
    }

    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
        BoundaryCondition bc = blk.bc[which_boundary];
        auto gmodel = blk.myConfig.gmodel;
        foreach (i, f; bc.faces) {
            sfs.set_flowstate(*(f.fs), t, f.pos, gmodel);
        }
    }

    @nogc
    override void apply_for_interface_structured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        throw new Error("BIE_SynthesiseFlowState.apply_for_interface_structured_grid() not yet implemented");
    }

    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");
        auto gmodel = blk.myConfig.gmodel;
        BoundaryCondition bc = blk.bc[which_boundary];
        foreach (i, f; bc.faces) {
            sfs.set_flowstate(*(f.fs), t, f.pos, gmodel);
        }
    } // end apply_structured_grid()

private:
    SyntheticFlowState sfs;
} // end class BIE_SynthesiseFlowState


class BIE_ZeroVelocity : BoundaryInterfaceEffect {
    this(int id, int boundary)
    {
        super(id, boundary, "ZeroVelocity");
    }

    override string toString() const
    {
        return "ZeroVelocity()";
    }

    override void apply_for_interface_unstructured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        auto gmodel = blk.myConfig.gmodel;
        BoundaryCondition bc = blk.bc[which_boundary];
	f.fs.vel.set(0.0, 0.0, 0.0);
    }

    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
        auto gmodel = blk.myConfig.gmodel;
        BoundaryCondition bc = blk.bc[which_boundary];
        foreach (i, f; bc.faces) { f.fs.vel.set(0.0, 0.0, 0.0); }
    }

    override void apply_for_interface_structured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");
        BoundaryCondition bc = blk.bc[which_boundary];
        f.fs.vel.set(0.0, 0.0, 0.0);
    }

    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");
        BoundaryCondition bc = blk.bc[which_boundary];
        foreach (i, f; bc.faces) { f.fs.vel.set(0.0, 0.0, 0.0); }
    }
} // end class BIE_ZeroVelocity


class BIE_ZeroSlipWallVelocity : BoundaryInterfaceEffect {
    // This boundary interface effect should work for both moving and fixed grids
    // because the grid velocity at the face should be already set appropriately.
    this(int id, int boundary)
    {
        super(id, boundary, "ZeroSlipWallVelocity");
    }

    override string toString() const
    {
        return "ZeroSlipWallVelocity()";
    }

    override void apply_for_interface_unstructured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        auto gmodel = blk.myConfig.gmodel;
        BoundaryCondition bc = blk.bc[which_boundary];
	f.fs.vel.set(f.gvel);
    }

    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
        auto gmodel = blk.myConfig.gmodel;
        BoundaryCondition bc = blk.bc[which_boundary];
        foreach (i, f; bc.faces) {
            f.fs.vel.set(f.gvel);
        }
    }

    override void apply_for_interface_structured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");
        BoundaryCondition bc = blk.bc[which_boundary];
        f.fs.vel.set(f.gvel);
    }

    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");
        BoundaryCondition bc = blk.bc[which_boundary];
        foreach (i, f; bc.faces) {
            f.fs.vel.set(f.gvel);
        }
    } // end apply_structured_grid()
} // end class BIE_ZeroSlipWallVelocity


class BIE_TranslatingSurface : BoundaryInterfaceEffect {
    // The boundary surface is translating with fixed velocity v_trans.
    Vector3 v_trans;

    this(int id, int boundary, Vector3 v_trans)
    {
        super(id, boundary, "TranslatingSurface");
        this.v_trans = v_trans;
    }

    override string toString() const
    {
        return "TranslatingSurface(v_trans=" ~ to!string(v_trans) ~ ")";
    }

    override void apply_for_interface_unstructured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        f.fs.vel.set(v_trans);
    }

    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
        BoundaryCondition bc = blk.bc[which_boundary];
        foreach (i, f; bc.faces) {
            f.fs.vel.set(v_trans);
        }
    }

    override void apply_for_interface_structured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        f.fs.vel.set(v_trans);
    }

    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");
        BoundaryCondition bc = blk.bc[which_boundary];
        foreach (i, f; bc.faces) {
            f.fs.vel.set(v_trans);
        }
    } // end apply_structured_grid()
} // end class BIE_TranslatingSurface


class BIE_RotatingSurface : BoundaryInterfaceEffect {
    // The boundary surface is rotating with fixed angular velocity r_omega
    // about centre.
    Vector3 r_omega;
    Vector3 centre;

    this(int id, int boundary, Vector3 r_omega, Vector3 centre)
    {
        super(id, boundary, "RotatingSurface");
        this.r_omega = r_omega;
        this.centre = centre;
    }

    override string toString() const
    {
        return "RotatingSurface(r_omega=" ~ to!string(r_omega) ~
            ", centre=" ~ to!string(centre) ~ ")";
    }

    override void apply_for_interface_unstructured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        throw new Error("BIE_RotatingSurface.apply_for_interface_unstructured_grid() not implemented yet");
    }

    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
        BoundaryCondition bc = blk.bc[which_boundary];
        foreach (i, f; bc.faces) {
            f.fs.vel = cross(r_omega, f.pos-centre);
        }
    } // end apply_unstructured_grid()

    override void apply_for_interface_structured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        f.fs.vel = cross(r_omega, f.pos-centre);
    }

    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");
        BoundaryCondition bc = blk.bc[which_boundary];
        foreach (i, f; bc.faces) {
            f.fs.vel = cross(r_omega, f.pos-centre);
        }
    } // end apply_structured_grid()
} // end class BIE_RotatingSurface


class BIE_FixedT : BoundaryInterfaceEffect {
public:
    double Twall;

    this(int id, int boundary, double Twall)
    {
        super(id, boundary, "FixedT");
        this.Twall = Twall;
    }

    override string toString() const
    {
        return "FixedT(Twall=" ~ to!string(Twall) ~ ")";
    }

    override void apply_for_interface_unstructured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        auto gmodel = blk.myConfig.gmodel;
        BoundaryCondition bc = blk.bc[which_boundary];
	f.fs.gas.T = Twall;
	version(multi_T_gas) {
	    foreach(ref elem; f.fs.gas.T_modes) { elem = Twall; }
	}
    }

    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
        BoundaryCondition bc = blk.bc[which_boundary];
        foreach (i, f; bc.faces) {
            f.fs.gas.T = Twall;
            version(multi_T_gas) {
                foreach(ref elem; f.fs.gas.T_modes) { elem = Twall; }
            }
        }
    } // end apply_unstructured_grid()

    override void apply_for_interface_structured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        auto gmodel = blk.myConfig.gmodel;
        BoundaryCondition bc = blk.bc[which_boundary];
	f.fs.gas.T = Twall;
	version(multi_T_gas) {
	    foreach(ref elem; f.fs.gas.T_modes) { elem = Twall; }
	}
    }

    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        auto gmodel = blk.myConfig.gmodel;
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");
        BoundaryCondition bc = blk.bc[which_boundary];
        foreach (i, f; bc.faces) {
            f.fs.gas.T = Twall;
            version(multi_T_gas) {
                foreach(ref elem; f.fs.gas.T_modes) { elem = Twall; }
            }
        }
    } // end apply_structured_grid()
} // end class BIE_FixedT

class BIE_FixedComposition : BoundaryInterfaceEffect {
public:
    double[] massfAtWall;

    this(int id, int boundary, double[] massfAtWall)
    {
        super(id, boundary, "FixedComposition");
        this.massfAtWall = massfAtWall;
    }

    override string toString() const
    {
        return "FixedComposition(massfAtWall=" ~ to!string(massfAtWall) ~ ")";
    }

    @nogc
    override void apply_for_interface_unstructured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        uint nsp = blk.myConfig.n_species;
        BoundaryCondition bc = blk.bc[which_boundary];
        version(multi_species_gas) {
            foreach (isp; 0 .. nsp) { f.fs.gas.massf[isp] = massfAtWall[isp]; }
        }
    }

    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
        uint nsp = blk.myConfig.n_species;
        BoundaryCondition bc = blk.bc[which_boundary];
        foreach (i, f; bc.faces) {
            version(multi_species_gas) {
                foreach (isp; 0 .. nsp) { f.fs.gas.massf[isp] = massfAtWall[isp]; }
            }
        }
    } // end apply_unstructured_grid()

    @nogc
    override void apply_for_interface_structured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        uint nsp = blk.myConfig.n_species;
        BoundaryCondition bc = blk.bc[which_boundary];
        version(multi_species_gas) {
            foreach (isp; 0 .. nsp) { f.fs.gas.massf[isp] = massfAtWall[isp]; }
        }
    }

    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        uint nsp = blk.myConfig.n_species;
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");
        BoundaryCondition bc = blk.bc[which_boundary];
        foreach (i, f; bc.faces) {
            version(multi_species_gas) {
                foreach (isp; 0 .. nsp) { f.fs.gas.massf[isp] = massfAtWall[isp]; }
            }
        }
    } // end apply_structured_grid()
} // end class BIE_FixedComposition


class BIE_UpdateThermoTransCoeffs : BoundaryInterfaceEffect {
    this(int id, int boundary, InterpolateOption thermoInterpOpt)
    {
        super(id, boundary, "UpdateThermoTransCoeffs");
        mInterpOpt = thermoInterpOpt;
    }

    override string toString() const
    {
        return "UpdateThermoTransCoeffs()";
    }

    override void apply_for_interface_unstructured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        BoundaryCondition bc = blk.bc[which_boundary];
        auto gmodel = blk.myConfig.gmodel;
        final switch (mInterpOpt) {
        case InterpolateOption.pt:
            gmodel.update_thermo_from_pT(f.fs.gas);
            break;
        case InterpolateOption.rhou:
            gmodel.update_thermo_from_rhou(f.fs.gas);
            break;
        case InterpolateOption.rhop:
            gmodel.update_thermo_from_rhop(f.fs.gas);
            break;
        case InterpolateOption.rhot:
            gmodel.update_thermo_from_rhoT(f.fs.gas);
            break;
        }
        gmodel.update_trans_coeffs(f.fs.gas);
    }

    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
        BoundaryCondition bc = blk.bc[which_boundary];
        auto gmodel = blk.myConfig.gmodel;
        foreach (i, f; bc.faces) {
            final switch (mInterpOpt) {
            case InterpolateOption.pt:
                gmodel.update_thermo_from_pT(f.fs.gas);
                break;
            case InterpolateOption.rhou:
                gmodel.update_thermo_from_rhou(f.fs.gas);
                break;
            case InterpolateOption.rhop:
                gmodel.update_thermo_from_rhop(f.fs.gas);
                break;
            case InterpolateOption.rhot:
                gmodel.update_thermo_from_rhoT(f.fs.gas);
                break;
            }
            gmodel.update_trans_coeffs(f.fs.gas);
        }
    } // end apply_unstructured_grid()

    override void apply_for_interface_structured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        BoundaryCondition bc = blk.bc[which_boundary];
        auto gmodel = blk.myConfig.gmodel;
        final switch (mInterpOpt) {
        case InterpolateOption.pt:
            gmodel.update_thermo_from_pT(f.fs.gas);
            break;
        case InterpolateOption.rhou:
            gmodel.update_thermo_from_rhou(f.fs.gas);
            break;
        case InterpolateOption.rhop:
            gmodel.update_thermo_from_rhop(f.fs.gas);
            break;
        case InterpolateOption.rhot:
            gmodel.update_thermo_from_rhoT(f.fs.gas);
            break;
        }
        gmodel.update_trans_coeffs(f.fs.gas);
    }

    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");
        BoundaryCondition bc = blk.bc[which_boundary];
        auto gmodel = blk.myConfig.gmodel;
        foreach (i, f; bc.faces) {
            final switch (mInterpOpt) {
            case InterpolateOption.pt:
                gmodel.update_thermo_from_pT(f.fs.gas);
                break;
            case InterpolateOption.rhou:
                gmodel.update_thermo_from_rhou(f.fs.gas);
                break;
            case InterpolateOption.rhop:
                gmodel.update_thermo_from_rhop(f.fs.gas);
                break;
            case InterpolateOption.rhot:
                gmodel.update_thermo_from_rhoT(f.fs.gas);
                break;
            }
            gmodel.update_trans_coeffs(f.fs.gas);
        }
    } // end apply_structured_grid()

    private:
    InterpolateOption mInterpOpt;
} // end class BIE_UpdateThermoTransCoeffs

class BIE_WallTurbulent : BoundaryInterfaceEffect {
    this(int id, int boundary)
    {
        super(id, boundary, "WallTurbulent");
    }

    override string toString() const
    {
        return "WallTurbulent()";
    }

    override void apply_for_interface_unstructured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        BoundaryCondition bc = blk.bc[which_boundary];
        version(turbulence) {
	    if (bc.outsigns[f.i_bndry] == 1) {
                blk.myConfig.turb_model.set_flowstate_at_wall(gtl, f, f.left_cell, *(f.fs));
	    } else {
                blk.myConfig.turb_model.set_flowstate_at_wall(gtl, f, f.right_cell, *(f.fs));
            }
        }
    } // end apply_unstructured_grid()

    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
        BoundaryCondition bc = blk.bc[which_boundary];
        version(turbulence) {
            foreach (i, f; bc.faces) {
                if (bc.outsigns[i] == 1) {
                    blk.myConfig.turb_model.set_flowstate_at_wall(gtl, f, f.left_cell, *(f.fs));
                } else {
                    blk.myConfig.turb_model.set_flowstate_at_wall(gtl, f, f.right_cell, *(f.fs));
                }
            }
        }
    } // end apply_unstructured_grid()

    override void apply_for_interface_structured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        BoundaryCondition bc = blk.bc[which_boundary];
        version(turbulence) {
	    if (bc.outsigns[f.i_bndry] == 1) {
                blk.myConfig.turb_model.set_flowstate_at_wall(gtl, f, f.left_cell, *(f.fs));
	    } else {
                blk.myConfig.turb_model.set_flowstate_at_wall(gtl, f, f.right_cell, *(f.fs));
            }
        }
    } // end apply_unstructured_grid()

    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");
        BoundaryCondition bc = blk.bc[which_boundary];
        version(turbulence) {
            foreach (i, f; bc.faces) {
                if (bc.outsigns[i] == 1) {
                    blk.myConfig.turb_model.set_flowstate_at_wall(gtl, f, f.left_cell, *(f.fs));
                } else {
                    blk.myConfig.turb_model.set_flowstate_at_wall(gtl, f, f.right_cell, *(f.fs));
                }
            }
        }
    } // end apply_structured_grid()


    //@nogc
    //number ideal_omega_at_wall(in FluidFVCell cell, number d0)
    //// As recommended by Wilson Chan, we use Menter's correction
    //// for omega values at the wall. This appears as Eqn A12 in
    //// Menter's paper.
    //// Reference:
    //// Menter (1994)
    //// Two-Equation Eddy-Viscosity Turbulence Models for
    //// Engineering Applications.
    //// AIAA Journal, 32:8, pp. 1598--1605
    //{
    //    auto wall_gas = cell.fs.gas;
    //    // Note: d0 is half_cell_width_at_wall.
    //    number nu = wall_gas.mu / wall_gas.rho;
    //    double beta1 = 0.075;
    //    return 10 * (6 * nu) / (beta1 * d0 * d0);
    //}
} // end class BIE_Turbulent

class BIE_WallFunction : BoundaryInterfaceEffect {
    this(int id, int boundary)
    {
        super(id, boundary, "WallFunction_InterfaceEffect");
        _faces_need_to_be_flagged = true;
    }

    override string toString() const
    {
        return "WallFunction_InterfaceEffect()";
    }

    override void apply_for_interface_unstructured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        throw new FlowSolverException("WallFunction_InterfaceEffect bc not implemented for unstructured grids.");
    }

    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
        throw new FlowSolverException("WallFunction_InterfaceEffect bc not implemented for unstructured grids.");
    }

    override void apply_for_interface_structured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        f.use_wall_function_shear_and_heat_flux = true; // Do this every time; Should ask Wilson about it. PJ 2021-06-19
        BoundaryCondition bc = blk.bc[which_boundary];
        auto c = (bc.outsigns[f.i_bndry] == 1) ? f.left_cells[0] : f.right_cells[0];
        wall_function(c, f);
    }

    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");
        BoundaryCondition bc = blk.bc[which_boundary];
        if (_faces_need_to_be_flagged) {
            // Flag faces, just once.
            foreach (i, f; bc.faces) {
                f.use_wall_function_shear_and_heat_flux = true;
            }
            _faces_need_to_be_flagged = false;
        }
        // Do some real work.
        foreach (i, f; bc.faces) {
            auto c = (bc.outsigns[i] == 1) ? f.left_cells[0] : f.right_cells[0];
            wall_function(c, f);
        }
    } // end apply_structured_grid()

    void wall_function(const FluidFVCell cell, FVInterface IFace)
    // Implement Nichols' and Nelson's wall function boundary condition
    // Reference:
    //  Nichols RH & Nelson CC (2004)
    //  Wall Function Boundary Conditions Inclding Heat Transfer
    //  and Compressibility.
    //  AIAA Journal, 42:6, pp. 1107--1114
    // Authors N. Gibbons and W. Chan
    // NOTE: IFace.fs will receive updated values of tke and omega for later
    //       copying to boundary cells.
    {
        auto gmodel = blk.myConfig.gmodel;
        // Compute tangent velocity at nearest interior point and wall interface
        number du, vt1_2_angle;
        number cell_tangent0, cell_tangent1, face_tangent0, face_tangent1;
        Vector3 cellVel = cell.fs.vel;
        Vector3 faceVel = IFace.fs.vel;
        cellVel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
        faceVel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
        if ( blk.myConfig.dimensions == 2 ) {
            cell_tangent0 = sqrt( pow(cellVel.y, 2.0) + pow(cellVel.z, 2.0) );
            face_tangent0 = sqrt( pow(faceVel.y, 2.0) + pow(faceVel.z, 2.0) );
            du = fabs(cell_tangent0 - face_tangent0);
        } else {
            cell_tangent0 = cellVel.y;
            cell_tangent1 = cellVel.z;
            face_tangent0 = faceVel.y;
            face_tangent1 = faceVel.z;
            vt1_2_angle = atan2(fabs(cellVel.z - faceVel.z), fabs(cellVel.y - faceVel.y));
            du = sqrt( pow((cell_tangent0-face_tangent0),2.0) + pow((cell_tangent1-face_tangent1),2.0) );
        }

        number tau_wall, q_wall;
        SolveShearStressAndHeatTransfer(cell, IFace, gmodel, du, tau_wall, q_wall);

        // Store wall shear stress and heat flux to be used later to replace viscous
        // stress in flux calculations. Also, for wall shear stress, transform value
        // back to the global frame of reference.
        double reverse_flag0 = 1.0; double reverse_flag1 = 1.0;
        Vector3 local_tau_wall;
        if ( blk.myConfig.dimensions == 2 ) {
            if ( face_tangent0 > cell_tangent0 ) reverse_flag0 = -1.0;
        } else {
            if ( face_tangent0 > cell_tangent0 ) reverse_flag0 = -1.0;
            if ( face_tangent1 > cell_tangent1 ) reverse_flag1 = -1.0;
        }
        if ( IFace.bc_id == Face.north || IFace.bc_id == Face.east || IFace.bc_id == Face.top ) {
            if ( blk.myConfig.dimensions == 2 ) {
                IFace.tau_wall_x = -1.0 * reverse_flag0 * tau_wall * IFace.n.y;
                IFace.tau_wall_y = -1.0 * reverse_flag0 * tau_wall * IFace.n.x;
                IFace.tau_wall_z = 0.0;
            } else {
                local_tau_wall = Vector3(to!number(0.0),
                                         -1.0 * reverse_flag0 * tau_wall * cos(vt1_2_angle),
                                         -1.0 * reverse_flag1 * tau_wall * sin(vt1_2_angle));
            }
            IFace.q = -1.0 * q_wall;
        } else { // South, West and Bottom
            if ( blk.myConfig.dimensions == 2 ) {
                IFace.tau_wall_x = reverse_flag0 * tau_wall * IFace.n.y;
                IFace.tau_wall_y = reverse_flag0 * tau_wall * IFace.n.x;
                IFace.tau_wall_z = 0.0;
            } else {
                local_tau_wall = Vector3(to!number(0.0),
                                         reverse_flag0 * tau_wall * cos(vt1_2_angle),
                                         reverse_flag1 * tau_wall * sin(vt1_2_angle));
            }
            IFace.q = q_wall;
        }
        if ( blk.myConfig.dimensions == 3 ) {
            local_tau_wall.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
            IFace.tau_wall_x = local_tau_wall.x;
            IFace.tau_wall_y = local_tau_wall.y;
            IFace.tau_wall_z = local_tau_wall.z;
        }
        number tke, omega, mu_t;
        mu_t = TurbulentViscosity(cell, IFace, gmodel, du, tau_wall, q_wall);
        komegaVariables(cell, IFace, du, tau_wall, mu_t, tke, omega); // TODO: Generalise

        version(turbulence) {
            // Assign updated values of tke and omega to IFace.fs for
            // later copying to boundary cells.
            IFace.fs.turb[0] = tke;
            IFace.fs.turb[1] = omega;
        }
        return;
    } // end wall_function()

protected:
    void SolveShearStressAndHeatTransfer(const FluidFVCell cell, FVInterface IFace,
                                         GasModel gmodel, const number du,
                                         ref number tau_wall, ref number q_wall){
        // Compute wall gas properties from either ...
        number cp = gmodel.Cp(cell.fs.gas);
        number Pr = cell.fs.gas.mu * cp / cell.fs.gas.k;
        number recovery = pow(Pr, (1.0/3.0));

        number T = cell.fs.gas.T;
        number T_wall = IFace.fs.gas.T;
        IFace.fs.gas.p = cell.fs.gas.p;
        gmodel.update_thermo_from_pT(IFace.fs.gas);
        gmodel.update_trans_coeffs(IFace.fs.gas);
        number rho_wall = IFace.fs.gas.rho;
        number k_lam_wall = IFace.fs.gas.k;
        number mu_lam_wall = IFace.fs.gas.mu;

        number wall_dist = distance_between(cell.pos[0], IFace.pos);
        number dudy = du / wall_dist;
        tau_wall = mu_lam_wall * dudy;

        number dTdy = fabs(T - T_wall)/wall_dist;
        q_wall = k_lam_wall * dTdy;

        // Two variable Newton's method to solve equations (8) and (10) simulteneously
        const number tolerance = 1.0e-10;
        number error = 1e32;
        ulong iterations = 0;
        number f1, f2, df1_dtau, df1_dq, df2_dtau, df2_dq;
        number diff_tau; number diff_q;
        while (error>tolerance) {
            f1 = WhiteChristophError(tau_wall, q_wall, rho_wall, T_wall, mu_lam_wall,
                                     k_lam_wall, recovery, du, wall_dist, cp);
            f2 = CroccoBusemanError(tau_wall, q_wall, T_wall, rho_wall, mu_lam_wall,
                                    k_lam_wall, du, T, recovery, cp);
            error = sqrt(f1*f1 + f2*f2);

            WhiteChristophJacobian(tau_wall, q_wall, rho_wall, T_wall, mu_lam_wall,
                                   k_lam_wall, recovery, du, wall_dist, cp, df1_dtau, df1_dq);
            CroccoBusemanJacobian(tau_wall, q_wall, mu_lam_wall, k_lam_wall, du, df2_dtau, df2_dq);

            diff_tau = (f1*df2_dq/df1_dq - f2)/(df2_dtau - df2_dq*df1_dtau/df1_dq);
            diff_q = (f1*df2_dtau/df1_dtau - f2)/(df2_dq - df2_dtau*df1_dq/df1_dtau);

            tau_wall += diff_tau;
            q_wall += diff_q;
            iterations += 1;
            if (iterations>1000) throw new Exception("Convergence failure in turbulent wall function");
        }
        return;
    }

    void komegaVariables(const FluidFVCell cell, const FVInterface IFace, const number du,
                         const number tau_wall, const number mu_t, ref number tke, ref number omega){
        // Compute omega (Eq 19 - 21)
        const number C_mu = 0.09;
        const number kappa = 0.4;

        number wall_dist = distance_between(cell.pos[0], IFace.pos);
        number rho_wall = IFace.fs.gas.rho;
        number mu_lam_wall = IFace.fs.gas.mu;
        number rho = cell.fs.gas.rho;

        number u_tau = sqrt( tau_wall / rho_wall );
        number omega_i = 6.0*mu_lam_wall / (0.075*rho_wall*wall_dist*wall_dist);
        number omega_o = u_tau / (sqrt(C_mu)*kappa*wall_dist);

        omega = sqrt(omega_i*omega_i + omega_o*omega_o);  // Compute omega (Eq 21)
        tke = omega * mu_t / rho;                         // Compute tke (Eq 22)
    }

    number TurbulentViscosity(const FluidFVCell cell, const FVInterface IFace, GasModel gmodel,
                              const number u, const number tau_wall, const number q_wall){
        // Turbulence model boundary conditions (Eq 15 & 14)
        // Note that the formulation of y_white_y_plus (Eq 15) is now directly
        // affected by the limiter which was set earlier to help get past large
        // velocity gradients at the wall for cases with initialised with high
        // velocity in flow domain.
        number cp = gmodel.Cp(cell.fs.gas);
        number Pr = cell.fs.gas.mu * cp / cell.fs.gas.k;
        number recovery = pow(Pr, (1.0/3.0));
        number T_wall = IFace.fs.gas.T;
        number rho_wall = IFace.fs.gas.rho;
        number k_lam_wall = IFace.fs.gas.k;
        number mu_lam_wall = IFace.fs.gas.mu;
        number mu_lam = cell.fs.gas.mu;
        number y = distance_between(cell.pos[0], IFace.pos);

        const number kappa = 0.4;
        const number B = 5.5;

        number u_tau = sqrt( tau_wall / rho_wall );
        number u_plus = u / u_tau;
        number y_plus = rho_wall*u_tau*y/mu_lam_wall;
        number Gam = recovery * u_tau * u_tau / (2.0 * cp * T_wall);
        number Beta = q_wall * mu_lam_wall / (rho_wall*T_wall*k_lam_wall*u_tau);
        number Q = sqrt(Beta*Beta + 4.0*Gam);
        number Phi = asin(-1.0 * Beta / Q);
        number alpha = (2.0*Gam*u_plus - Beta)/Q;
        if (alpha > 0.0) alpha = fmin(alpha, 1-1e-12);
        else alpha = fmax(alpha, -1+1e-12);

        number y_plus_white = exp((kappa/sqrt(Gam))*(asin(alpha) - Phi))*exp(-1.0*kappa*B);

        number y_white_y_plus = 2.0 * y_plus_white * kappa*sqrt(Gam)/Q
                * pow((1.0 - pow(alpha,2.0)), 0.5);
        number mu_coeff = 1.0 + y_white_y_plus
                - kappa*exp(-1.0*kappa*B) * (1.0 + kappa*u_plus + kappa*u_plus*kappa*u_plus/2.0)
                - mu_lam/mu_lam_wall;
        // Limit turbulent-to-laminar viscosity ratio between zero and the global limit.
        if ( mu_coeff < 0.0 ) mu_coeff = 0.0;
        mu_coeff = fmin(mu_coeff, blk.myConfig.max_mu_t_factor);
        // Compute turbulent viscosity; forcing mu_t to zero, if result is negative.
        return mu_lam_wall * mu_coeff;
    }

    number WhiteChristophError(number tau_wall, number q_wall, number rho_wall, number T_wall,
                               number mu_lam_wall, number k_lam_wall, number recovery, number u,
                               number y, number cp){
        /*
           Equation (8) + (9) from Nichols paper

           Note:
               In the calculation of y+ defined by White and Christoph
               (Eq 9), the equation breaks down when the value of
               asin((2.0*Gam*u_plus - Beta)/Q) goes larger than 1.0 or
               smaller than -1.0. For cases where we initialise the flow
               solution with high flow velocity, du (and hence u_plus)
               becomes large enough to exceed this limit. A limiter is
               therefore implemented here to help get past this initially
               large velocity gradient at the wall. Note that this limiter
               is not in Nichols and Nelson's paper. We set the limit to
               either a value just below 1.0, or just above -1.0, to avoid
               the calculation of y_white_y_plus (Eq. 15) from blowing up.
        */
        const number kappa = 0.4;
        const number B = 5.5;

        number u_tau = sqrt( tau_wall / rho_wall );
        number u_plus = u / u_tau;
        number y_plus = rho_wall*u_tau*y/mu_lam_wall;
        number Gam = recovery * u_tau * u_tau / (2.0 * cp * T_wall);
        number Beta = q_wall * mu_lam_wall / (rho_wall*T_wall*k_lam_wall*u_tau);
        number Q = sqrt(Beta*Beta + 4.0*Gam);
        number Phi = asin(-1.0 * Beta / Q);
        number alpha = (2.0*Gam*u_plus - Beta)/Q;
        if (alpha > 0.0) alpha = fmin(alpha, 1-1e-12);
        else alpha = fmax(alpha, -1+1e-12);

        number y_plus_white = exp((kappa/sqrt(Gam))*(asin(alpha) - Phi))*exp(-1.0*kappa*B);
        number thing = exp(-1.0*kappa*B) * ( 1.0 + kappa*u_plus + pow((kappa*u_plus), 2.0) / 2.0
                                                                + pow((kappa*u_plus), 3.0) / 6.0 );
        number error = u_plus + y_plus_white - thing - y_plus;
        return error;
    }

    void WhiteChristophJacobian(number tau_wall, number q_wall, number rho_wall, number T_wall,
                                number mu_lam_wall, number k_lam_wall, number recovery, number u,
                                number y, number cp, ref number dWC_dtauw, ref number dWC_dqw){

        const number kappa = 0.4;
        const number B = 5.5;
        number u_tau = sqrt( tau_wall / rho_wall );
        number u_plus = u / u_tau;
        number y_plus = rho_wall*u_tau*y/mu_lam_wall;
        number Gam = recovery * u_tau * u_tau / (2.0 * cp * T_wall);
        number Beta = q_wall * mu_lam_wall / (rho_wall*T_wall*k_lam_wall*u_tau);
        number Q = sqrt(Beta*Beta + 4.0*Gam);
        number Phi = asin(-1.0 * Beta / Q);
        number alpha = (2.0*Gam*u_plus - Beta)/Q;
        if (alpha > 0.0) alpha = fmin(alpha, 1-1e-12);
        else alpha = fmax(alpha, -1+1e-12);
        number P = asin(alpha);

        number y_plus_white = exp((kappa/sqrt(Gam))*(P - Phi))*exp(-1.0*kappa*B);

        number duplus_dtauw = -u/2.0*sqrt(rho_wall/tau_wall/tau_wall/tau_wall);
        number dyplus_dtauw = y/2.0/mu_lam_wall*sqrt(rho_wall/tau_wall);
        number dthing_duplus= exp(-1.0*kappa*B)*(kappa + kappa*kappa*u_plus +
                                                 0.5*kappa*kappa*kappa*u_plus*u_plus);
        number dthing_dtauw = dthing_duplus*duplus_dtauw;

        number dBeta_dtauw  = -q_wall*mu_lam_wall/2.0/T_wall/k_lam_wall/sqrt(rho_wall*tau_wall*tau_wall*tau_wall);
        number dGam_dtauw   = recovery/2.0/cp/T_wall/rho_wall;
        number dQ_dtauw     = 0.5/Q*(2.0*Beta*dBeta_dtauw + 4.0*dGam_dtauw);
        number dalpha_dtauw = (2*Gam*duplus_dtauw + 2*u_plus*dGam_dtauw - dBeta_dtauw)/Q - alpha/Q*dQ_dtauw;
        number dP_dtauw     = dalpha_dtauw/sqrt(1-alpha*alpha);
        number dPhi_dtauw   = (Beta*dQ_dtauw - Q*dBeta_dtauw)/Q/Q/sqrt(1.0 - Beta*Beta/Q/Q);

        number dyplusWhite_dtauw = y_plus_white*kappa/Gam*(dP_dtauw*sqrt(Gam) - dPhi_dtauw*sqrt(Gam) -
                                                           (P-Phi)*dGam_dtauw/2.0/sqrt(Gam));

        dWC_dtauw = duplus_dtauw + dyplusWhite_dtauw - dthing_dtauw - dyplus_dtauw;

        number duplus_dqw = 0.0;
        number dyplus_dqw = 0.0;
        number dthing_dqw = 0.0;

        number dBeta_dqw = mu_lam_wall/T_wall/k_lam_wall/sqrt(rho_wall*tau_wall);
        number dGam_dqw = 0.0;
        number dQ_dqw     = 0.5/Q*(2.0*Beta*dBeta_dqw + 4.0*dGam_dqw);
        number dalpha_dqw = (2*Gam*duplus_dqw + 2*u_plus*dGam_dqw - dBeta_dqw)/Q - alpha/Q*dQ_dqw;
        number dP_dqw     = dalpha_dqw/sqrt(1-alpha*alpha);
        number dPhi_dqw   = (Beta*dQ_dqw - Q*dBeta_dqw)/Q/Q/sqrt(1.0 - Beta*Beta/Q/Q);

        number dyplusWhite_dqw = y_plus_white*kappa/Gam*(dP_dqw*sqrt(Gam) - dPhi_dqw*sqrt(Gam) -
                                                         (P-Phi)*dGam_dqw/2.0/sqrt(Gam));

        dWC_dqw = dyplusWhite_dqw;
    }

    number CroccoBusemanError(number tau_wall, number q_wall, number T_wall, number rho_wall,
                              number mu_lam_wall, number k_lam_wall, number u, number T,
                              number recovery, number cp){
        /*
           Equation (10) from Nichols Paper

       */
        number u_tau = sqrt( tau_wall / rho_wall );
        number u_plus = u / u_tau;
        number Gam = recovery * u_tau * u_tau / (2.0 * cp * T_wall);
        number Beta = q_wall * mu_lam_wall / (rho_wall*T_wall*k_lam_wall*u_tau);
        number thing = 1.0 + Beta*u_plus - Gam*u_plus*u_plus;
        number error = T_wall*thing - T;
        return error;
    }

    void CroccoBusemanJacobian(number tau_wall, number q_wall, number mu_lam_wall, number k_lam_wall,
                               number u, ref number dCB_dtauw, ref number dCB_dqw){

        dCB_dtauw = -q_wall*mu_lam_wall*u/(k_lam_wall*tau_wall*tau_wall);
        dCB_dqw = mu_lam_wall*u/k_lam_wall/tau_wall;
        return;
    }
private:
    bool _faces_need_to_be_flagged = true;

} // end class BIE_WallFunction

class BIE_AdiabaticWallFunction : BIE_WallFunction {
    this(int id, int boundary)
    {
        super(id, boundary);
        desc = "AdiabaticWallFunction_InterfaceEffect";
    }

    override string toString() const
    {
        return "AdiabaticWallFunction_InterfaceEffect()";
    }

protected:
    override void SolveShearStressAndHeatTransfer(const FluidFVCell cell, FVInterface IFace, GasModel gmodel,
                                                  const number du, ref number tau_wall, ref number q_wall){
        number cp = gmodel.Cp(cell.fs.gas);
        number Pr = cell.fs.gas.mu * cp / cell.fs.gas.k;
        number recovery = pow(Pr, (1.0/3.0));

        // set wall properties using the Crocco-Busemann relation (Eq 11)
        number T = cell.fs.gas.T;
        number T_wall = T + recovery * du * du / (2.0 * cp);
        IFace.fs.gas.p = cell.fs.gas.p;
        IFace.fs.gas.T = T_wall;
        gmodel.update_thermo_from_pT(IFace.fs.gas);
        gmodel.update_trans_coeffs(IFace.fs.gas);
        number rho_wall = IFace.fs.gas.rho;
        number k_lam_wall = IFace.fs.gas.k;
        number mu_lam_wall = IFace.fs.gas.mu;

        number wall_dist = distance_between(cell.pos[0], IFace.pos);
        number dudy = du / wall_dist;
        tau_wall = mu_lam_wall * dudy;

        q_wall = 0.0;

        // One variable Newton's method to solve equation (8) with q_wall = 0.0
        const number tolerance = 1.0e-10;
        number error = 1e32;
        ulong iterations = 0;
        number f1, df1_dtau, df1_dq;
        number diff_tau;
        while (error>tolerance) {
            f1 = WhiteChristophError(tau_wall, q_wall, rho_wall, T_wall, mu_lam_wall, k_lam_wall,
                                     recovery, du, wall_dist, cp);
            error = sqrt(f1*f1);

            WhiteChristophJacobian(tau_wall, q_wall, rho_wall, T_wall, mu_lam_wall, k_lam_wall,
                                   recovery, du, wall_dist, cp, df1_dtau, df1_dq);

            diff_tau = -f1/df1_dtau;
            tau_wall += diff_tau;
            iterations += 1;
            if (iterations>1000) throw new Exception("Convergence failure in turbulent wall function");
        }
        return;
    }

}

// NOTE: This GAS DOMAIN boundary effect has a large
//       and important side-effect:
//       IT ALSO SETS THE FLUX IN THE ADJACENT SOLID DOMAIN
//       AT THE TIME IT IS CALLED.
// TODO: We need to work out a way to coordinate this
//       interface effect (ie. the setting of temperature)
//       with the flux effect. Ideally, we only want to compute
//       the temperature/flux once per update. This will require
//       some storage at the block level, or in the in the
//       gas/solid interface module since we can't share information
//       (easily) between different types of boundary condition
//       objects. We need to store the energy flux somewhere where it
//       so that we can use it again in the boundary flux effect.
//       It's no good storing the flux in the interface object since
//       that will be changed during the diffusive flux calculation.

class BIE_TemperatureFromGasSolidInterface : BoundaryInterfaceEffect {
public:
    int neighbourSolidBlk;
    int neighbourSolidFace;
    int neighbourOrientation;

    this(int id, int boundary,
         int otherBlock, int otherFace, int orient)
    {
        super(id, boundary, "TemperatureFromGasSolidInterface");
        neighbourSolidBlk = otherBlock;
        neighbourSolidFace = otherFace;
        neighbourOrientation = orient;
    }

    override string toString() const
    {
        return "TemperatureFromGasSolidInterface()";
    }

    @nogc
    override void apply_for_interface_unstructured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        throw new Error("BIE_TemperatureFromGasSolidInterface.apply_unstructured_grid() not implemented yet");
    }

    @nogc
    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
        throw new Error("BIE_TemperatureFromGasSolidInterface.apply_unstructured_grid() not implemented yet");
    }

    @nogc
    override void apply_for_interface_structured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        auto myBC = blk.bc[which_boundary];
        number dxG, dyG, dzG, dnG, dxS, dyS, dzS, dnS;
        number kG_dnG, kS_dnS, cosA, cosB, cosC;
        number T, q;

        cosA = myBC.ifaces[f.i_bndry].n.x;
        cosB = myBC.ifaces[f.i_bndry].n.y;
        cosC = myBC.ifaces[f.i_bndry].n.z;

        dxG = myBC.ifaces[f.i_bndry].pos.x - myBC.gasCells[f.i_bndry].pos[0].x;
        dyG = myBC.ifaces[f.i_bndry].pos.y - myBC.gasCells[f.i_bndry].pos[0].y;
        dzG = myBC.ifaces[f.i_bndry].pos.z - myBC.gasCells[f.i_bndry].pos[0].z;
        dnG = fabs(cosA*dxG + cosB*dyG + cosC*dzG);

        dxS = myBC.ifaces[f.i_bndry].pos.x - myBC.solidCells[f.i_bndry].pos.x;
        dyS = myBC.ifaces[f.i_bndry].pos.y - myBC.solidCells[f.i_bndry].pos.y;
        dzS = myBC.ifaces[f.i_bndry].pos.z - myBC.solidCells[f.i_bndry].pos.z;
        dnS = fabs(cosA*dxS + cosB*dyS + cosC*dzS);

        kG_dnG = myBC.gasCells[f.i_bndry].fs.gas.k / dnG;
        kS_dnS = myBC.solidCells[f.i_bndry].ss.k / dnS;

        T = (myBC.gasCells[f.i_bndry].fs.gas.T*kG_dnG + myBC.solidCells[f.i_bndry].T*kS_dnS) / (kG_dnG + kS_dnS);

        // Finally update properties in interfaces
        myBC.ifaces[f.i_bndry].fs.gas.T = T;
    }

    @nogc
    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        auto myBC = blk.bc[which_boundary];
        number dxG, dyG, dzG, dnG, dxS, dyS, dzS, dnS;
        number kG_dnG, kS_dnS, cosA, cosB, cosC;
        number T, q;

        foreach ( i; 0 .. myBC.ifaces.length ) {
            cosA = myBC.ifaces[i].n.x;
            cosB = myBC.ifaces[i].n.y;
            cosC = myBC.ifaces[i].n.z;

            dxG = myBC.ifaces[i].pos.x - myBC.gasCells[i].pos[0].x;
            dyG = myBC.ifaces[i].pos.y - myBC.gasCells[i].pos[0].y;
            dzG = myBC.ifaces[i].pos.z - myBC.gasCells[i].pos[0].z;
            dnG = fabs(cosA*dxG + cosB*dyG + cosC*dzG);

            dxS = myBC.ifaces[i].pos.x - myBC.solidCells[i].pos.x;
            dyS = myBC.ifaces[i].pos.y - myBC.solidCells[i].pos.y;
            dzS = myBC.ifaces[i].pos.z - myBC.solidCells[i].pos.z;
            dnS = fabs(cosA*dxS + cosB*dyS + cosC*dzS);

            kG_dnG = myBC.gasCells[i].fs.gas.k / dnG;
            kS_dnS = myBC.solidCells[i].ss.k / dnS;

            T = (myBC.gasCells[i].fs.gas.T*kG_dnG + myBC.solidCells[i].T*kS_dnS) / (kG_dnG + kS_dnS);

            // Finally update properties in interfaces
            myBC.ifaces[i].fs.gas.T = T;
        }
    }

} // end class BIE_TemperatureFromGasSolidInterface

class BIE_TemperatureFromGasSolidInterface2 : BoundaryInterfaceEffect {
public:
    int neighbourSolidBlk;
    int neighbourSolidFace;
    int neighbourOrientation;

    this(int id, int boundary,
         int otherBlock, int otherFace, int orient)
    {
        super(id, boundary, "TemperatureFromGasSolidInterface2");
        neighbourSolidBlk = otherBlock;
        neighbourSolidFace = otherFace;
        neighbourOrientation = orient;
    }

    override string toString() const
    {
        return "TemperatureFromGasSolidInterface2()";
    }

    @nogc
    override void apply_for_interface_unstructured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        throw new Error("BIE_TemperatureFromGasSolidInterface2.apply_unstructured_grid() not implemented yet");
    }

    @nogc
    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
        throw new Error("BIE_TemperatureFromGasSolidInterface2.apply_unstructured_grid() not implemented yet");
    }

    @nogc
    override void apply_for_interface_structured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        auto myBC = blk.bc[which_boundary];
        number Twall = myBC.solidCells[f.i_bndry].T;
        f.fs.gas.T = Twall;
        version(multi_T_gas) {
            foreach(ref elem; f.fs.gas.T_modes) { elem = Twall; }
        }
    }

    @nogc
    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        auto myBC = blk.bc[which_boundary];
        number Twall;
        foreach ( i; 0 .. myBC.ifaces.length ) {
            Twall = myBC.solidCells[i].T;
            myBC.ifaces[i].fs.gas.T = Twall;
            version(multi_T_gas) {
                foreach(ref elem; myBC.ifaces[i].fs.gas.T_modes) { elem = Twall; }
            }
        }
    }

} // end class BIE_TemperatureFromGasSolidInterface2

class BIE_ThermionicRadiativeEquilibrium : BoundaryInterfaceEffect {
    this(int id, int boundary, double emissivity, double Ar, double phi, int ThermionicEmissionActive, size_t nsp,
         const string gas_file_name, const string catalytic_type, double[] massfAtWall)
    {
        super(id, boundary, "ThermionicRadiativeEquilibrium");
        this.emissivity = emissivity;
        this.Ar = Ar;
        this.phi = phi*Qe;  // Convert phi from input 'eV' to 'J'
        this.ThermionicEmissionActive = ThermionicEmissionActive;
        this.nsp = nsp;
        switch(catalytic_type){
            case "none":
                break;
            case "equilibrium":
                eqcalc = new EquilibriumCalculator(gas_file_name);
                catalytic = true;
                break;
            case "fixed_composition":
                if (massfAtWall.length!=nsp)
                    throw new Error("massfAtWall.length does not match gmodel.n_species");
                fmassf = massfAtWall;
                setmassf = true;
                catalytic = true;
                break;
            default:
                throw new Error("TRE Boundary Condition does not support catalytic type: " ~ catalytic_type);
        }

        if (catalytic) {
            jx.length = nsp;
            jy.length = nsp;
            jz.length = nsp;
        }

        if ((catalytic) && (GlobalConfig.mass_diffusion_model == MassDiffusionModel.none))
            throw new Error("Catalytic Wall requires a mass diffusion model.");
    }

    override string toString() const
    {
        return "BIE_ThermionicRadiativeEquilibrium(ThermionicEmissionActive=" ~
            to!string(ThermionicEmissionActive) ~
            ", Work Function =" ~ to!string(phi/Qe) ~
            "eV , emissivity=" ~ to!string(emissivity) ~
            ", Richardson Constant=" ~ to!string(Ar) ~
            ")";
    }

    @nogc
    override void apply_for_interface_unstructured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        BoundaryCondition bc = blk.bc[which_boundary];
        int outsign = bc.outsigns[f.i_bndry];
        auto c = (outsign == 1) ? f.left_cell : f.right_cell;
        solve_for_wall_temperature_and_energy_flux(c, f, outsign);
    }

    @nogc
    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
        BoundaryCondition bc = blk.bc[which_boundary];
        foreach (i, f; bc.faces) {
            int outsign = bc.outsigns[i];
            auto c = (outsign == 1) ? f.left_cell : f.right_cell;
            solve_for_wall_temperature_and_energy_flux(c, f, outsign);
        }
    } // end apply_unstructured_grid()

    @nogc
    override void apply_for_interface_structured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        BoundaryCondition bc = blk.bc[which_boundary];
        int outsign = bc.outsigns[f.i_bndry];
        auto c = (outsign == 1) ? f.left_cells[0] : f.right_cells[0];
        solve_for_wall_temperature_and_energy_flux(c, f, outsign);
    }

    @nogc
    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");
        BoundaryCondition bc = blk.bc[which_boundary];
        foreach (i, f; bc.faces) {
            int outsign = bc.outsigns[i];
            auto c = (outsign == 1) ? f.left_cells[0] : f.right_cells[0];
            solve_for_wall_temperature_and_energy_flux(c, f, outsign);
        }
    } // end apply_structured_grid()

protected:
    // Function inputs from Eilmer4 .lua simulation input
    double emissivity;  // Input emissivity, 0<e<=1.0. Assumed black body radiation out from wall
    double Ar;          // Richardson constant, material-dependent
    double phi;         // Work function, material dependent. Input units in eV,
                        // this gets converted to Joules by multiplying by Elementary charge, Qe
    int ThermionicEmissionActive;  // Whether or not Thermionic Emission is active. Default is 'on'

    // Pieces needed for catalytic heat transfer
    number[] jx, jy, jz;
    double[] fmassf;
    bool setmassf, catalytic;
    size_t nsp;
    EquilibriumCalculator eqcalc;

    // Constants used in analysis
    immutable double SB_sigma = 5.670373e-8;  // Stefan-Boltzmann constant.   Units: W/(m^2 K^4)
    immutable double kb = 1.38064852e-23;     // Boltzmann constant.          Units: (m^2 kg)/(s^2 K^1)
    immutable double Qe = 1.60217662e-19;     // Elementary charge.           Units: C

    @nogc
    void solve_for_wall_temperature_and_energy_flux(FluidFVCell cell, FVInterface IFace, int outsign)
    {
    /*
        Set the temperature of the wall interface to satisfy the thermionic+radiative energy
        balance equations. Note the the outsign parameter is used to flip the faces normal
        vector to be pointing into the domain, if needed.
    */
        auto gmodel = blk.myConfig.gmodel;
        uint n_modes = blk.myConfig.n_modes;
        double viscous_factor = blk.myConfig.viscous_factor;

        IFace.fs.gas.p = cell.fs.gas.p;
        if (setmassf) foreach (isp; 0 .. nsp) IFace.fs.gas.massf[isp] = fmassf[isp];
        number Twall;
        number T0 = IFace.fs.gas.T;

        Twall = T0;
        // One variable Newton's method to solve heat balance equations
        immutable number eta = 1e-9;
        immutable number tolerance = 1.0e-10;
        number error = 1e32;
        ulong iterations = 0;
        number f, dfdT,df;
        number diff_T;
        while (error>tolerance) {
            f = ThermionicRadiativeEnergyBalance(gmodel, cell, IFace, n_modes, outsign, Twall);
            df= ThermionicRadiativeEnergyBalance(gmodel, cell, IFace, n_modes, outsign, Twall+eta);
            dfdT = (df-f)/eta;

            error = sqrt(f*f);
            diff_T = -f/dfdT;
            Twall += diff_T;
            iterations += 1;
            if (iterations>100) {
                string msg = "Convergence failure in ThermionicRadiative boundary condition.";
                debug{
                    msg ~= format("\ndiff T %f Twall %f error %f f %f\n", diff_T, Twall, error, f);
                    msg ~= format("Face.fs.gas: %s\n", IFace.fs.gas);
                    msg ~= format("Face.id: %s Face.pos %s", IFace.id, IFace.pos);
                }
                throw new Exception(msg);
            }
        }

        IFace.fs.gas.T = Twall;
        version(multi_T_gas) { foreach (imode; 0 .. n_modes) IFace.fs.gas.T_modes[imode] = Twall; }
        if (eqcalc) eqcalc.set_massf_from_pT(IFace.fs.gas);
        gmodel.update_thermo_from_pT(IFace.fs.gas);
        gmodel.update_trans_coeffs(IFace.fs.gas);
        return;

    } // end solve_for_wall_temperature_and_energy_flux()

    @nogc
    number ThermionicRadiativeEnergyBalance(GasModel gmodel, FluidFVCell cell, FVInterface IFace,
                                            uint n_modes, int outsign, number Twall){
    /*
        Energy flux balance at the wall, from Alkandry, 2014 equations 6 and 10.
        This version uses a hardcoded gradients_leastsq call, because the viscous flux
        routine was causing weird side effects. You will have to FIXME for catalytic walls.
        @author: Nick Gibbons
    */
        IFace.fs.gas.T = Twall;
        version(multi_T_gas) { foreach (imode; 0 .. n_modes) IFace.fs.gas.T_modes[imode] = Twall; }
        if (eqcalc) eqcalc.set_massf_from_pT(IFace.fs.gas);
        gmodel.update_thermo_from_pT(IFace.fs.gas);
        gmodel.update_trans_coeffs(IFace.fs.gas);

        cell.grad.gradients_leastsq(blk.myConfig, cell.cloud_fs, cell.cloud_pos, *(cell.ws_grad));
        IFace.grad.copy_values_from(*(cell.grad));
        number qx = IFace.fs.gas.k*cell.grad.T[0];
        number qy = IFace.fs.gas.k*cell.grad.T[1];
        number qz = IFace.fs.gas.k*cell.grad.T[2];
        version(multi_T_gas) {
            qx += IFace.fs.gas.k_modes[0]*cell.grad.T_modes[0][0];
            qy += IFace.fs.gas.k_modes[0]*cell.grad.T_modes[0][1];
            qz += IFace.fs.gas.k_modes[0]*cell.grad.T_modes[0][2];
        }
        // Negative sign is because heat flows along temperature gradient from hot to cold
        number q_conduction = -1.0*(IFace.n.x*qx + IFace.n.y*qy + IFace.n.z*qz);
        // We then correct for the direction of IFace.n using outsign
        q_conduction *= outsign;

        // Species diffusion has the opposite sign, for some reason.
        number q_diffusion = to!number(0.0);
        if (catalytic){
            blk.myConfig.massDiffusion.update_mass_fluxes(*(IFace.fs), *(cell.grad), jx, jy, jz);
            foreach (isp; 0 .. nsp) {
                number h = gmodel.enthalpy(IFace.fs.gas, cast(int)isp);
                q_diffusion += (jx[isp]*h*IFace.n.x + jy[isp]*h*IFace.n.y + jz[isp]*h*IFace.n.z);
            }
        }
        q_diffusion *= outsign;

        number q_rad = emissivity*SB_sigma*Twall*Twall;
        number q_thermionic = to!number(0.0);
        if (ThermionicEmissionActive == 1) {
            q_thermionic = Ar/Qe*exp(-phi/kb/Twall)*(phi + 2*kb*Twall);
        }

        // We actually solve for the energy divided by T^2 to reduce the size of the numbers
        // and reduce floating point round off error.
        return q_rad + q_thermionic - q_conduction/Twall/Twall - q_diffusion/Twall/Twall;
    }
} // end class BIE_ThermionicRadiativeEquilibrium

class BIE_EquilibriumComposition : BoundaryInterfaceEffect {
/*
    Equilibrium Boundary Interface effect for use in catalytic
    wall studies. Note that this effect assumes that the pressure
    and temperature at the interface face been set correctly,
    and that the transprops will be recomputed after it is called.

    @author: Nick Gibbons
*/
    this(int id, int boundary, string gas_file_name)
    {
        super(id, boundary, "EquilibriumComposition");
        eqcalc = new EquilibriumCalculator(gas_file_name);
    }

    override string toString() const
    {
        return "BIE_EquilibriumComposition()";
    }

    @nogc
    override void apply_for_interface_unstructured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        set_equilibrium_composition(f);
    }

    @nogc
    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
        BoundaryCondition bc = blk.bc[which_boundary];
        foreach (f; bc.faces) {
            set_equilibrium_composition(f);
        }
    }

    @nogc
    override void apply_for_interface_structured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        set_equilibrium_composition(f);
    }

    @nogc
    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");
        BoundaryCondition bc = blk.bc[which_boundary];
        foreach (f; bc.faces) {
            set_equilibrium_composition(f);
        }
    }

protected:
    EquilibriumCalculator eqcalc;

    @nogc
    void set_equilibrium_composition(FVInterface IFace) {
    /*
        Set IFace.gs.gas.massf to the equilibrium results, assuming a later interface effect
        is going to tidy up the thermodynamic state and transport coefficients.

        See file src/kinetics/equilibrium_update for the
        EquilibriumCalculator source code.
    */
        eqcalc.set_massf_from_pT(IFace.fs.gas);
    }


} // end class BIE_EquilibriumComposition
