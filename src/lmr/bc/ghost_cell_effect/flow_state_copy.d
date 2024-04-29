// flow_state_copy.d

module bc.ghost_cell_effect.flow_state_copy;

import std.json;
import std.string;
import std.conv;
import std.stdio;
import std.math;

import nm.number;
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


class GhostCellFlowStateCopy : GhostCellEffect {
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

    @nogc
    override void apply_for_interface_unstructured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        BoundaryCondition bc = blk.bc[which_boundary];
        auto ghost = (bc.outsigns[f.i_bndry] == 1) ? f.right_cell : f.left_cell;
        ghost.fs.copy_values_from(fstate);
        if (r > 0.0) { compute_source_flow(ghost); }
        if (blk.omegaz != 0.0) {
            into_rotating_frame(ghost.fs.vel, ghost.pos[gtl], blk.omegaz);
        }
    } // end apply_for_interface_unstructured_grid()

    @nogc
    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
        BoundaryCondition bc = blk.bc[which_boundary];
        foreach (i, f; bc.faces) {
            auto ghost = (bc.outsigns[i] == 1) ? f.right_cell : f.left_cell;
            ghost.fs.copy_values_from(fstate);
            if (r > 0.0) { compute_source_flow(ghost); }
            if (blk.omegaz != 0.0) {
                into_rotating_frame(ghost.fs.vel, ghost.pos[gtl], blk.omegaz);
            }
        }
    } // end apply_unstructured_grid()

    @nogc
    override void apply_for_interface_structured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");
        BoundaryCondition bc = blk.bc[which_boundary];
        foreach (n; 0 .. blk.n_ghost_cell_layers) {
            auto ghost = (bc.outsigns[f.i_bndry] == 1) ? f.right_cells[n] : f.left_cells[n];
            ghost.fs.copy_values_from(fstate);
            if (r > 0.0) { compute_source_flow(ghost); }
            if (blk.omegaz != 0.0) {
                into_rotating_frame(ghost.fs.vel, ghost.pos[gtl], blk.omegaz);
            }
        }
    } // end apply_for_interface_structured_grid()

    @nogc
    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");
        BoundaryCondition bc = blk.bc[which_boundary];
        foreach (i, f; bc.faces) {
            foreach (n; 0 .. blk.n_ghost_cell_layers) {
                auto ghost = (bc.outsigns[i] == 1) ? f.right_cells[n] : f.left_cells[n];
                ghost.fs.copy_values_from(fstate);
                if (r > 0.0) { compute_source_flow(ghost); }
                if (blk.omegaz != 0.0) {
                    into_rotating_frame(ghost.fs.vel, ghost.pos[gtl], blk.omegaz);
                }
            }
        }
    } // end apply_structured_grid()

private:
    @nogc
    void compute_source_flow(FluidFVCell c)
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
        double dx = c.pos[0].x.re - x0;
        double dy = c.pos[0].y.re - y0;
        double dz = (cqi.threeD) ? c.pos[0].z.re - z0 : 0.0;
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
        FlowState* fs = c.fs; // local pointer
        fs.gas.p = p; // not really needed because we use rhou
        fs.gas.rho = rho;
        fs.gas.u = u;
        gmodel.update_thermo_from_rhou(fs.gas);
        fs.vel.set(velx, vely, velz);
    } // end apply_to_single_face()

} // end class GhostCellFlowStateCopy
