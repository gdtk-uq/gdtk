
module bc.ghost_cell_effect.porous_wall;

import std.json;
import std.string;
import std.conv;
import std.stdio;
import std.math;
import nm.complex;
import nm.number;
import util.lua;
import util.lua_service;

import geom;
import globalconfig;
import globaldata;
import flowstate;
import fvinterface;
import fvcell;
import fluidblock;
import sfluidblock;
import gas;
import bc;


class PorousWallGhostCellEffect : GhostCellEffect {
/*
    Steady-state Darcy-Forchheimer porous wall equations.
    Based on:
        "Combustion enhancement in a scramjet engine using oxygen
        enrichment and porous fuel injection"
        Bianca R. Capra, R. R. Boyce, M. Kuhn and H. Hald
        J. Fluid Mech. (2015), vol. 767, pp. 173–198.
        doi:10.1017/jfm.2015.43

    We also reference:
        "Permeability of ceramic foams to compressible and in-
        compressible flow"
        E. A. Moreira amd M. D. M. Innocentini and J. R. Coury
        Journal of the European Ceramic Society 24 (2004) 3209–3218

    @author: Nick Gibbons (July 2024)
*/
public:
    this(int id, int boundary, double[] injectant_massf, double kF, double kD, double wall_temperature, double plenum_pressure, double porous_plate_thickness)
    {
        super(id, boundary, "PorousWallInGhostCell");
        this.kF = kF;
        this.kD = kD;
        this.wall_temperature = wall_temperature;
        this.plenum_pressure = plenum_pressure;
        this.porous_plate_thickness = porous_plate_thickness;

        auto gmodel = GlobalConfig.gmodel_master;
        auto nturb = GlobalConfig.turb_model.nturb;
        this.injectant = FlowState(gmodel, nturb);
        this.plenum_state = GasState(gmodel);
        set_plenum_state(plenum_pressure, wall_temperature, injectant_massf);
    }

    override string toString() const
    {
        return "PorousFlowInGhostCell(kF=" ~ to!string(kF) ~
            ", kD=" ~ to!string(kD) ~
            ", Twall=" ~ to!string(wall_temperature) ~
            ", plenum_pressure=" ~ to!string(plenum_pressure) ~ ")";
    }


    override void apply_for_interface_unstructured_grid(double t, int gtl, int ftl, FVInterface face)
    {
        BoundaryCondition bc = blk.bc[which_boundary];
        int outsign = bc.outsigns[face.i_bndry];
        FVCell cell = (outsign == 1) ? face.left_cell : face.right_cell;
        FVCell ghost0 = (outsign == 1) ? face.right_cell : face.left_cell;

        set_injectant_from_flow_pressure(cell.fs.gas.p, face.n, outsign);
        ghost0.fs.copy_values_from(injectant);
    }

    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
        BoundaryCondition bc = blk.bc[which_boundary];
        // Now, apply the ghost-cell conditions
        foreach (i, face; bc.faces) {
            int outsign = bc.outsigns[i];
            FVCell cell = (outsign == 1) ? face.left_cell : face.right_cell;
            FVCell ghost0 = (outsign == 1) ? face.right_cell : face.left_cell;
            set_injectant_from_flow_pressure(cell.fs.gas.p, face.n, outsign);
            ghost0.fs.copy_values_from(injectant);
        }
    } // end apply_unstructured_grid()

    override void apply_for_interface_structured_grid(double t, int gtl, int ftl, FVInterface face)
    {
        auto bc = blk.bc[which_boundary];
        int outsign = bc.outsigns[face.i_bndry];
        FVCell cell = (outsign == 1) ? face.left_cells[0] : face.right_cells[0];
        set_injectant_from_flow_pressure(cell.fs.gas.p, face.n, outsign);

        foreach (n; 0 .. blk.n_ghost_cell_layers) {
            FVCell ghost = (outsign == 1) ? face.right_cells[n] : face.left_cells[n];
            ghost.fs.copy_values_from(injectant);
        }
    }

    // not @nogc
    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        BoundaryCondition bc = blk.bc[which_boundary];

        foreach (i, face; bc.faces) {
            int outsign = bc.outsigns[i];
            FVCell cell = (outsign == 1) ? face.left_cells[0] : face.right_cells[0];
            set_injectant_from_flow_pressure(cell.fs.gas.p, face.n, outsign);

            foreach (n; 0 .. blk.n_ghost_cell_layers) {
                FVCell ghost = (outsign == 1) ? face.right_cells[n] : face.left_cells[n];
                ghost.fs.copy_values_from(injectant);
            }
        }
    } // end apply_structured_grid()
private:
    double kF;
    double kD;
    double wall_temperature;
    double plenum_pressure;
    double porous_plate_thickness;
    GasState plenum_state;
    FlowState injectant;

    @nogc
    void set_plenum_state(double pressure, double temperature, double[] mass_fractions)
    {
        auto gmodel = GlobalConfig.gmodel_master;
        plenum_state.p = pressure;
        plenum_state.T = temperature;
        foreach(i, Y; mass_fractions) plenum_state.massf[i] = Y;
        gmodel.update_thermo_from_pT(plenum_state);
        gmodel.update_trans_coeffs(plenum_state);
    }

    @nogc
    void set_injectant_from_flow_pressure(number pout, Vector3 n, int outsign)
    {
    /*
        Equation 3.1 from Capra et al. 2015. We also use the weird expression
        for delta p/L from Moreira et al. 2004, equation 1b.
    */

        number rhoin = plenum_state.rho;
        number pin = plenum_state.p;
        number muin = plenum_state.mu;

        number A = rhoin/kF;
        number B = muin/kD;
        number C = -(pin*pin - pout*pout)/2.0/pin/porous_plate_thickness;

        number v = (-B + sqrt(B*B - 4.0*A*C))/(2.0*A);
        assert(v>0.0);

        // Note that this is the velocity on the plenum side, the velocity at the
        // outflow can be computed using mass conservation.
        auto gmodel = blk.myConfig.gmodel;
        injectant.gas.p = pout;
        injectant.gas.T = wall_temperature;
        injectant.gas.massf[] = plenum_state.massf[];
        gmodel.update_thermo_from_pT(injectant.gas);

        number mdot_per_m2 = v*rhoin;
        number vout = mdot_per_m2/injectant.gas.rho;

        injectant.vel = -1.0*outsign*n*vout;

        debug{
            writefln("GC pressureout %e. Computed mass flux %e vout %e", pout, mdot_per_m2, vout);
        }
        return;
    }
} // end class PorousWallInGhostCell
