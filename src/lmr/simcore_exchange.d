// simcore_exchange.d
// 2021-04-15: extracted from simcore.d
//

module simcore_exchange;

import std.math;
import std.stdio;
import std.file;
import std.conv;
import std.array;
import std.format;
import std.string;
import std.algorithm;
import std.parallelism;
import ntypes.complex;
import nm.number;

import geom;
import geom.misc.kdtree;
import gas;
import globalconfig;
import globaldata;
import flowstate;
import fluidblock;
import sfluidblock;
import ufluidblock;
import ssolidblock;
import solidfvinterface;
import solid_full_face_copy;
import solid_gas_full_face_copy;
import bc.ghost_cell_effect.gas_solid_full_face_copy;
import bc;
import user_defined_source_terms;
import solid_udf_source_terms;
import grid_motion;
import grid_motion_udf;
import grid_motion_shock_fitting;
version(mpi_parallel) {
    import mpi;
}

void exchange_ghost_cell_geometry_data()
/*
    This routine copies the following cell geometric data from one block
    into the neighbouring block's ghost cells:
     - cell.pos.x
     - cell.volume[j]
     - cell.areaxy[j]
     - cell.iLength
     - cell.jLength
     - cell.kLength
     - cell.L_min
     - cell.L_max

     Notes: - This code was originally present in exchange_ghost_cell_boundary_data, but
              was disintered into a separate routine for the geometry (this one) and the
              flowstate (handled by the original routine).
     @author: Nick Gibbons (19/08/21)
*/
{
    foreach (blk; localFluidBlocks) {
        foreach(bc; blk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce1 = cast(GhostCellMappedCellCopy) gce;
                if (mygce1) { mygce1.exchange_geometry_phase0(); }
                auto mygce2 = cast(GhostCellFullFaceCopy) gce;
                if (mygce2) { mygce2.exchange_geometry_phase0(); }
            }
        }
    }
    foreach (blk; localFluidBlocks) {
        foreach(bc; blk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce1 = cast(GhostCellMappedCellCopy) gce;
                if (mygce1) { mygce1.exchange_geometry_phase1(); }
                auto mygce2 = cast(GhostCellFullFaceCopy) gce;
                if (mygce2) { mygce2.exchange_geometry_phase1(); }
            }
        }
    }
    foreach (blk; localFluidBlocks) {
        foreach(bc; blk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce1 = cast(GhostCellMappedCellCopy) gce;
                if (mygce1) { mygce1.exchange_geometry_phase2(); }
                auto mygce2 = cast(GhostCellFullFaceCopy) gce;
                if (mygce2) { mygce2.exchange_geometry_phase2(); }
            }
        }
    }
}


void exchange_ghost_cell_boundary_data(double t, int gtl, int ftl)
// We have hoisted the exchange of ghost-cell data out of the GhostCellEffect class
// that used to live only inside the boundary condition attached to a block.
// The motivation for allowing this leakage of abstraction is that the MPI
// exchange of messages requires a coordination of actions that spans blocks.
// Sigh...  2017-01-24 PJ
// p.s. The data for that coordination is still buried in the FullFaceCopy class.
// No need to have all its guts hanging out.
{
    foreach (blk; localFluidBlocks) {
        foreach(bc; blk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce1 = cast(GhostCellMappedCellCopy) gce;
                if (mygce1) { mygce1.exchange_flowstate_phase0(t, gtl, ftl); }
                auto mygce2 = cast(GhostCellFullFaceCopy) gce;
                if (mygce2) { mygce2.exchange_flowstate_phase0(t, gtl, ftl); }
            }
        }
    }
    foreach (blk; localFluidBlocks) {
        foreach(bc; blk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce1 = cast(GhostCellMappedCellCopy) gce;
                if (mygce1) { mygce1.exchange_flowstate_phase1(t, gtl, ftl); }
                auto mygce2 = cast(GhostCellFullFaceCopy) gce;
                if (mygce2) { mygce2.exchange_flowstate_phase1(t, gtl, ftl); }
            }
        }
    }
    foreach (blk; localFluidBlocks) {
        foreach(bc; blk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce1 = cast(GhostCellMappedCellCopy) gce;
                if (mygce1) { mygce1.exchange_flowstate_phase2(t, gtl, ftl); }
                auto mygce2 = cast(GhostCellFullFaceCopy) gce;
                if (mygce2) { mygce2.exchange_flowstate_phase2(t, gtl, ftl); }
            }
        }
    }
} // end exchange_ghost_cell_boundary_data()

void exchange_ghost_cell_turbulent_viscosity()
{
    /*
        Some turbulence models (such as k-w) require the viscous gradients to
        compute mu_t and k_t. This means that we need to do a second exchange
        in the middle of calculating the viscous terms, to make sure that the
        ghost cell mu_t and k_t are up to date, or else we lose second order
        convergence at the boundaries. Rather than exchange the whole
        flowstate, this routine just does mu_t and k_t. It should be invoked
        only when we are sure that it is needed, by asking the turbulence model.

        @author: Nick Gibbons (20/08/21)
    */
    if (!GlobalConfig.turb_model.isTurbulent) return;

    foreach (blk; localFluidBlocks) {
        foreach(bc; blk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce1 = cast(GhostCellMappedCellCopy) gce;
                if (mygce1) { mygce1.exchange_turbulent_transprops_phase0(); }
                auto mygce2 = cast(GhostCellFullFaceCopy) gce;
                if (mygce2) { mygce2.exchange_turbulent_transprops_phase0(); }
            }
        }
    }
    foreach (blk; localFluidBlocks) {
        foreach(bc; blk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce1 = cast(GhostCellMappedCellCopy) gce;
                if (mygce1) { mygce1.exchange_turbulent_transprops_phase1(); }
                auto mygce2 = cast(GhostCellFullFaceCopy) gce;
                if (mygce2) { mygce2.exchange_turbulent_transprops_phase1(); }
            }
        }
    }
    foreach (blk; localFluidBlocks) {
        foreach(bc; blk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce1 = cast(GhostCellMappedCellCopy) gce;
                if (mygce1) { mygce1.exchange_turbulent_transprops_phase2(); }
                auto mygce2 = cast(GhostCellFullFaceCopy) gce;
                if (mygce2) { mygce2.exchange_turbulent_transprops_phase2(); }
            }
        }
    }
} // end exchange_ghost_cell_boundary_data()


void exchange_ghost_cell_shock_data(double t, int gtl, int ftl)
// exchange shock detector data
{
    foreach (blk; localFluidBlocks) {
        foreach(bc; blk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce = cast(GhostCellFullFaceCopy) gce;
                if (mygce) { mygce.exchange_shock_phase0(t, gtl, ftl); }
            }
        }
    }
    foreach (blk; localFluidBlocks) {
        foreach(bc; blk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce = cast(GhostCellFullFaceCopy) gce;
                if (mygce) { mygce.exchange_shock_phase1(t, gtl, ftl); }
            }
        }
    }
    foreach (blk; localFluidBlocks) {
        foreach(bc; blk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce = cast(GhostCellFullFaceCopy) gce;
                if (mygce) { mygce.exchange_shock_phase2(t, gtl, ftl); }
            }
        }
    }
} // end exchange_ghost_cell_boundary_data()

void exchange_ghost_cell_boundary_convective_gradient_data(double t, int gtl, int ftl)
// We have hoisted the exchange of ghost-cell data out of the GhostCellEffect class
// that used to live only inside the boundary condition attached to a block.
// The motivation for allowing this leakage of abstraction is that the MPI
// exchange of messages requires a coordination of actions that spans blocks.
// Sigh...  2017-01-24 PJ
// p.s. The data for that coordination is still buried in the FullFaceCopy class.
// No need to have all its guts hanging out.
{
    foreach (blk; localFluidBlocks) {
        foreach(bc; blk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce1 = cast(GhostCellMappedCellCopy) gce;
                if (mygce1) { mygce1.exchange_convective_gradient_phase0(t, gtl, ftl); }
            }
        }
    }
    foreach (blk; localFluidBlocks) {
        foreach(bc; blk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce1 = cast(GhostCellMappedCellCopy) gce;
                if (mygce1) { mygce1.exchange_convective_gradient_phase1(t, gtl, ftl); }
            }
        }
    }
    foreach (blk; localFluidBlocks) {
        foreach(bc; blk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce1 = cast(GhostCellMappedCellCopy) gce;
                if (mygce1) { mygce1.exchange_convective_gradient_phase2(t, gtl, ftl); }
            }
        }
    }
} // end exchange_ghost_cell_boundary_convective_gradient_data()

void exchange_ghost_cell_boundary_viscous_gradient_data(double t, int gtl, int ftl)
// We have hoisted the exchange of ghost-cell data out of the GhostCellEffect class
// that used to live only inside the boundary condition attached to a block.
// The motivation for allowing this leakage of abstraction is that the MPI
// exchange of messages requires a coordination of actions that spans blocks.
// Sigh...  2017-01-24 PJ
// p.s. The data for that coordination is still buried in the FullFaceCopy class.
// No need to have all its guts hanging out.
{
    foreach (blk; localFluidBlocks) {
        foreach(bc; blk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce1 = cast(GhostCellMappedCellCopy) gce;
                if (mygce1) { mygce1.exchange_viscous_gradient_phase0(t, gtl, ftl); }
                auto mygce2 = cast(GhostCellFullFaceCopy) gce;
                if (mygce2) { mygce2.exchange_viscous_gradient_phase0(t, gtl, ftl); }
            }
        }
    }
    foreach (blk; localFluidBlocks) {
        foreach(bc; blk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce1 = cast(GhostCellMappedCellCopy) gce;
                if (mygce1) { mygce1.exchange_viscous_gradient_phase1(t, gtl, ftl); }
                auto mygce2 = cast(GhostCellFullFaceCopy) gce;
                if (mygce2) { mygce2.exchange_viscous_gradient_phase1(t, gtl, ftl); }
            }
        }
    }
    foreach (blk; localFluidBlocks) {
        foreach(bc; blk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce1 = cast(GhostCellMappedCellCopy) gce;
                if (mygce1) { mygce1.exchange_viscous_gradient_phase2(t, gtl, ftl); }
                auto mygce2 = cast(GhostCellFullFaceCopy) gce;
                if (mygce2) { mygce2.exchange_viscous_gradient_phase2(t, gtl, ftl); }
            }
        }
    }
} // end exchange_ghost_cell_boundary_viscous_gradient_data()

void exchange_ghost_cell_solid_boundary_data()
{
    foreach (myblk; localSolidBlocks) {
        foreach (bc; myblk.bc) {
            foreach (gce; bc.preSpatialDerivActionAtBndryCells) {
                auto mygce = cast(SolidGCE_SolidGhostCellFullFaceCopy)gce;
                if (mygce) { mygce.exchange_solid_data_phase0(); }
            }
        }
    }
    foreach (myblk; localSolidBlocks) {
        foreach (bc; myblk.bc) {
            foreach (gce; bc.preSpatialDerivActionAtBndryCells) {
                auto mygce = cast(SolidGCE_SolidGhostCellFullFaceCopy)gce;
                if (mygce) { mygce.exchange_solid_data_phase1(); }
            }
        }
    }
    foreach (myblk; localSolidBlocks) {
        foreach (bc; myblk.bc) {
            foreach (gce; bc.preSpatialDerivActionAtBndryCells) {
                auto mygce = cast(SolidGCE_SolidGhostCellFullFaceCopy)gce;
                if (mygce) { mygce.exchange_solid_data_phase2(); }
            }
        }
    }
} // end exchange_ghost_cell_solid_boundary_data()


void exchange_ghost_cell_gas_solid_boundary_data()
{
    foreach (myblk; localFluidBlocks) {
        foreach (bc; myblk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce = cast(GhostCellGasSolidFullFaceCopy)gce;
                if (mygce) { mygce.exchange_solidstate_phase0(); }
            }
        }
    }
    foreach (myblk; localSolidBlocks) {
        foreach (bc; myblk.bc) {
            foreach (gce; bc.preSpatialDerivActionAtBndryCells) {
                auto mygce = cast(GhostCellSolidGasFullFaceCopy)gce;
                if (mygce) { mygce.exchange_fluidstate_phase0(); }
            }
        }
    }
    foreach (myblk; localFluidBlocks) {
        foreach (bc; myblk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce = cast(GhostCellGasSolidFullFaceCopy)gce;
                if (mygce) { mygce.exchange_fluidstate_phase1(); }
            }
        }
    }
    foreach (myblk; localSolidBlocks) {
        foreach (bc; myblk.bc) {
            foreach (gce; bc.preSpatialDerivActionAtBndryCells) {
                auto mygce = cast(GhostCellSolidGasFullFaceCopy)gce;
                if (mygce) { mygce.exchange_solidstate_phase1(); }
            }
        }
    }
    foreach (myblk; localFluidBlocks) {
        foreach (bc; myblk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce = cast(GhostCellGasSolidFullFaceCopy)gce;
                if (mygce) { mygce.exchange_solidstate_phase2(); }
            }
        }
    }
    foreach (myblk; localSolidBlocks) {
        foreach (bc; myblk.bc) {
            foreach (gce; bc.preSpatialDerivActionAtBndryCells) {
                auto mygce = cast(GhostCellSolidGasFullFaceCopy)gce;
                if (mygce) { mygce.exchange_fluidstate_phase2(); }
            }
        }
    }
}

void send_gas_domain_boundary_heat_flux_data_to_solid_domain()
{
    foreach (myblk; localSolidBlocks) {
        foreach (bc; myblk.bc) {
            foreach (gce; bc.preSpatialDerivActionAtBndryCells) {
                auto mygce = cast(GhostCellSolidGasFullFaceCopy)gce;
                if (mygce) { mygce.exchange_fluidstate_heat_flux_phase0(); }
            }
        }
    }
    foreach (myblk; localFluidBlocks) {
        foreach (bc; myblk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce = cast(GhostCellGasSolidFullFaceCopy)gce;
                if (mygce) { mygce.exchange_fluidstate_heat_flux_phase1(); }
            }
        }
    }
    foreach (myblk; localSolidBlocks) {
        foreach (bc; myblk.bc) {
            foreach (gce; bc.preSpatialDerivActionAtBndryCells) {
                auto mygce = cast(GhostCellSolidGasFullFaceCopy)gce;
                if (mygce) { mygce.exchange_fluidstate_heat_flux_phase2(); }
            }
        }
    }
}

void send_solid_domain_boundary_temperature_data_to_gas_domain()
{
    foreach (myblk; localFluidBlocks) {
        foreach (bc; myblk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce = cast(GhostCellGasSolidFullFaceCopy)gce;
                if (mygce) { mygce.exchange_solidstate_temperature_phase0(); }
            }
        }
    }
    foreach (myblk; localSolidBlocks) {
        foreach (bc; myblk.bc) {
            foreach (gce; bc.preSpatialDerivActionAtBndryCells) {
                auto mygce = cast(GhostCellSolidGasFullFaceCopy)gce;
                if (mygce) { mygce.exchange_solidstate_temperature_phase1(); }
            }
        }
    }
    foreach (myblk; localFluidBlocks) {
        foreach (bc; myblk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce = cast(GhostCellGasSolidFullFaceCopy)gce;
                if (mygce) { mygce.exchange_solidstate_temperature_phase2(); }
            }
        }
    }
}
