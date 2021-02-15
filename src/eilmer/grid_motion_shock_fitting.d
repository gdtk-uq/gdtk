// grid_motion_shock_fitting.d
// Module for implementing shock-fitting with a moving grid in Eilmer4.
//
// 2015-Nov Kyle D. original implmentation (moving grid, shock fitting).
// 2019 Lachlan Whyborn multiple blocks and MPI
// 2019-Nov PJ tidy-up to make the code more readable (by PJ, at least).
// 2021-02-15 PJ complete rework.

module grid_motion_shock_fitting;

import std.math;
import nm.complex;
import nm.number;
version(mpi_parallel) {
    import mpi;
    import mpi.util;
}

import globalconfig;
import globaldata;
import fvcore;
import fvvertex;
import fvinterface;
import fvcell;
import bc;
import fluidblock;
import sfluidblock;
import geom;
import grid_motion;
import bc;
version(mpi_parallel) {
    import bc.ghost_cell_effect.full_face_copy : MPI_Wait_a_while, make_mpi_tag;
}

version(mpi_parallel) {
    // MPI-parallel flavour of shock-fitting functions.
    // [TODO] 2021-02-15 PJ
} else {
    // Shared-memory-parallel.
    // [TODO] 2021-02-15 PJ
}
