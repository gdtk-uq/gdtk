// bcs.cu
// This file is included into simulate.cu, somewhere around the middle of that file.
//
// For a given block and boundary, we work across all faces in the boundary
// and set the ghost-cell flow states appropriate for the boundary condition.
// These functions are intended to be used just for the CPU flavour of the code.
//
// PJ 2022-10-03

#ifndef BCS_INCLUDED
#define BCS_INCLUDED

// Part A.
// Set the ghost-cell quantities to effect the boundary-conditions for convective fluxes.

__host__
void bc_wall_reflect_normal_velocity(int iblk, int ibc)
// Copy data, reflecting velocity.
{
    BConfig& cfg = blk_configs[iblk];
    Block& blk = fluidBlocks[iblk];
    switch (ibc) {
    case Face::iminus: // jk across face
        for (int k=0; k < cfg.nkc; k++) {
            for (int j=0; j < cfg.njc; j++) {
                FVFace& f = blk.faces[cfg.iFaceIndex(0, j, k)];
                FVCell& c = blk.cells[f.right_cells[0]];
                FlowState& fs0 = blk.cells[f.left_cells[0]].fs;
                fs0 = c.fs;
                fs0.vel.transform_to_local_frame(f.n, f.t1, f.t2);
                fs0.vel.x = -(fs0.vel.x);
                fs0.vel.transform_to_global_frame(f.n, f.t1, f.t2);
                FlowState& fs1 = blk.cells[f.left_cells[1]].fs;
                fs1 = c.fs;
                fs1.vel.transform_to_local_frame(f.n, f.t1, f.t2);
                fs1.vel.x = -(fs1.vel.x);
                fs1.vel.transform_to_global_frame(f.n, f.t1, f.t2);
            } // end for j
        } // end for k
        break;
    case Face::iplus: // jk across face
        for (int k=0; k < cfg.nkc; k++) {
            for (int j=0; j < cfg.njc; j++) {
                FVFace& f = blk.faces[cfg.iFaceIndex(cfg.nic, j, k)];
                FVCell& c = blk.cells[f.left_cells[0]];
                FlowState& fs0 = blk.cells[f.right_cells[0]].fs;
                fs0 = c.fs;
                fs0.vel.transform_to_local_frame(f.n, f.t1, f.t2);
                fs0.vel.x = -(fs0.vel.x);
                fs0.vel.transform_to_global_frame(f.n, f.t1, f.t2);
                FlowState& fs1 = blk.cells[f.right_cells[1]].fs;
                fs1 = c.fs;
                fs1.vel.transform_to_local_frame(f.n, f.t1, f.t2);
                fs1.vel.x = -(fs1.vel.x);
                fs1.vel.transform_to_global_frame(f.n, f.t1, f.t2);
            } // end for j
        } // end for k
        break;
    case Face::jminus: // ik across face
        for (int k=0; k < cfg.nkc; k++) {
            for (int i=0; i < cfg.nic; i++) {
                FVFace& f = blk.faces[cfg.jFaceIndex(i, 0, k)];
                FVCell& c = blk.cells[f.right_cells[0]];
                FlowState& fs0 = blk.cells[f.left_cells[0]].fs;
                fs0 = c.fs;
                fs0.vel.transform_to_local_frame(f.n, f.t1, f.t2);
                fs0.vel.x = -(fs0.vel.x);
                fs0.vel.transform_to_global_frame(f.n, f.t1, f.t2);
                FlowState& fs1 = blk.cells[f.left_cells[1]].fs;
                fs1 = c.fs;
                fs1.vel.transform_to_local_frame(f.n, f.t1, f.t2);
                fs1.vel.x = -(fs1.vel.x);
                fs1.vel.transform_to_global_frame(f.n, f.t1, f.t2);
            } // end for i
        } // end for k
        break;
    case Face::jplus: // ik across face
        for (int k=0; k < cfg.nkc; k++) {
            for (int i=0; i < cfg.nic; i++) {
                FVFace& f = blk.faces[cfg.jFaceIndex(i, cfg.njc, k)];
                FVCell& c = blk.cells[f.left_cells[0]];
                FlowState& fs0 = blk.cells[f.right_cells[0]].fs;
                fs0 = c.fs;
                fs0.vel.transform_to_local_frame(f.n, f.t1, f.t2);
                fs0.vel.x = -(fs0.vel.x);
                fs0.vel.transform_to_global_frame(f.n, f.t1, f.t2);
                FlowState& fs1 = blk.cells[f.right_cells[1]].fs;
                fs1 = c.fs;
                fs1.vel.transform_to_local_frame(f.n, f.t1, f.t2);
                fs1.vel.x = -(fs1.vel.x);
                fs1.vel.transform_to_global_frame(f.n, f.t1, f.t2);
            } // end for i
        } // end for k
        break;
    case Face::kminus: // ij across face
        for (int j=0; j < cfg.njc; j++) {
            for (int i=0; i < cfg.nic; i++) {
                FVFace& f = blk.faces[cfg.kFaceIndex(i, j, 0)];
                FVCell& c = blk.cells[f.right_cells[0]];
                FlowState& fs0 = blk.cells[f.left_cells[0]].fs;
                fs0 = c.fs;
                fs0.vel.transform_to_local_frame(f.n, f.t1, f.t2);
                fs0.vel.x = -(fs0.vel.x);
                fs0.vel.transform_to_global_frame(f.n, f.t1, f.t2);
                FlowState& fs1 = blk.cells[f.left_cells[1]].fs;
                fs1 = c.fs;
                fs1.vel.transform_to_local_frame(f.n, f.t1, f.t2);
                fs1.vel.x = -(fs1.vel.x);
                fs1.vel.transform_to_global_frame(f.n, f.t1, f.t2);
            } // end for i
        } // end for j
        break;
    case Face::kplus: // ij across face
        for (int j=0; j < cfg.njc; j++) {
            for (int i=0; i < cfg.nic; i++) {
                FVFace& f = blk.faces[cfg.kFaceIndex(i, j, cfg.nkc)];
                FVCell& c = blk.cells[f.left_cells[0]];
                FlowState& fs0 = blk.cells[f.right_cells[0]].fs;
                fs0 = c.fs;
                fs0.vel.transform_to_local_frame(f.n, f.t1, f.t2);
                fs0.vel.x = -(fs0.vel.x);
                fs0.vel.transform_to_global_frame(f.n, f.t1, f.t2);
                FlowState& fs1 = blk.cells[f.right_cells[1]].fs;
                fs1 = c.fs;
                fs1.vel.transform_to_local_frame(f.n, f.t1, f.t2);
                fs1.vel.x = -(fs1.vel.x);
                fs1.vel.transform_to_global_frame(f.n, f.t1, f.t2);
            } // end for j
        } // end for k
        break;
    }
} // end bc_wall_reflect_normal_velocity()


__host__
void bc_exchange(int iblk, int ibc)
// Depending on which face, look up the adjacent block in ijk indices and copy the flow data.
// The other block will have a corresponding boundary with the same type of boundary condition,
// so the "exchange" occurs is two phases.
{
    BConfig& cfg = blk_configs[iblk];
    Block& blk = fluidBlocks[iblk];
    //
    switch (ibc) {
    case Face::iminus: { // jk across face
        int other_i = cfg.i - 1;
        if (other_i < 0) { other_i = Config::nib-1; } // Wrap around.
        int other_j = cfg.j;
        int other_k = cfg.k;
        int other_id = Config::blk_ids[other_i][other_j][other_k];
        Block& other_blk = fluidBlocks[other_id];
        BConfig& other_cfg = blk_configs[other_id];
        for (int k=0; k < cfg.nkc; k++) {
            for (int j=0; j < cfg.njc; j++) {
                FVFace& f = blk.faces[cfg.iFaceIndex(0, j, k)];
                FVFace& other_f = other_blk.faces[other_cfg.iFaceIndex(other_cfg.nic, j, k)];
                blk.cells[f.left_cells[0]].fs = other_blk.cells[other_f.left_cells[0]].fs;
                blk.cells[f.left_cells[1]].fs = other_blk.cells[other_f.left_cells[1]].fs;
            } // end for j
        } // end for k
        break;
    }
    case Face::iplus: { // jk across face
        int other_i = cfg.i + 1;
        if (other_i >= Config::nib) { other_i = 0; } // Wrap around.
        int other_j = cfg.j;
        int other_k = cfg.k;
        int other_id = Config::blk_ids[other_i][other_j][other_k];
        Block& other_blk = fluidBlocks[other_id];
        BConfig& other_cfg = blk_configs[other_id];
        for (int k=0; k < cfg.nkc; k++) {
            for (int j=0; j < cfg.njc; j++) {
                FVFace& f = blk.faces[cfg.iFaceIndex(cfg.nic, j, k)];
                FVFace& other_f = other_blk.faces[other_cfg.iFaceIndex(0, j, k)];
                blk.cells[f.right_cells[0]].fs = other_blk.cells[other_f.right_cells[0]].fs;
                blk.cells[f.right_cells[1]].fs = other_blk.cells[other_f.right_cells[1]].fs;
            } // end for j
        } // end for k
        break;
    }
    case Face::jminus: { // ik across face
        int other_i = cfg.i;
        int other_j = cfg.j - 1;
        if (other_j < 0) { other_j = Config::njb-1; } // Wrap around.
        int other_k = cfg.k;
        int other_id = Config::blk_ids[other_i][other_j][other_k];
        Block& other_blk = fluidBlocks[other_id];
        BConfig& other_cfg = blk_configs[other_id];
        for (int k=0; k < cfg.nkc; k++) {
            for (int i=0; i < cfg.nic; i++) {
                FVFace& f = blk.faces[cfg.jFaceIndex(i, 0, k)];
                FVFace& other_f = other_blk.faces[other_cfg.jFaceIndex(i, other_cfg.njc, k)];
                blk.cells[f.left_cells[0]].fs = other_blk.cells[other_f.left_cells[0]].fs;
                blk.cells[f.left_cells[1]].fs = other_blk.cells[other_f.left_cells[1]].fs;
            } // end for i
        } // end for k
        break;
    }
    case Face::jplus: { // ik across face
        int other_i = cfg.i;
        int other_j = cfg.j + 1;
        if (other_j >= Config::njb) { other_j = 0; } // Wrap around.
        int other_k = cfg.k;
        int other_id = Config::blk_ids[other_i][other_j][other_k];
        Block& other_blk = fluidBlocks[other_id];
        BConfig& other_cfg = blk_configs[other_id];
        for (int k=0; k < cfg.nkc; k++) {
            for (int i=0; i < cfg.nic; i++) {
                FVFace& f = blk.faces[cfg.jFaceIndex(i, cfg.njc, k)];
                FVFace& other_f = other_blk.faces[other_cfg.jFaceIndex(i, 0, k)];
                blk.cells[f.right_cells[0]].fs = other_blk.cells[other_f.right_cells[0]].fs;
                blk.cells[f.right_cells[1]].fs = other_blk.cells[other_f.right_cells[1]].fs;
            } // end for i
        } // end for k
        break;
    }
    case Face::kminus: { // ij across face
        int other_i = cfg.i;
        int other_j = cfg.j;
        int other_k = cfg.k - 1;
        if (other_k < 0) { other_k = Config::nkb-1; } // Wrap around.
        int other_id = Config::blk_ids[other_i][other_j][other_k];
        Block& other_blk = fluidBlocks[other_id];
        BConfig& other_cfg = blk_configs[other_id];
        for (int j=0; j < cfg.njc; j++) {
            for (int i=0; i < cfg.nic; i++) {
                FVFace& f = blk.faces[cfg.kFaceIndex(i, j, 0)];
                FVFace& other_f = other_blk.faces[other_cfg.kFaceIndex(i, j, other_cfg.nkc)];
                blk.cells[f.left_cells[0]].fs = other_blk.cells[other_f.left_cells[0]].fs;
                blk.cells[f.left_cells[1]].fs = other_blk.cells[other_f.left_cells[1]].fs;
            } // end for i
        } // end for j
        break;
    }
    case Face::kplus: { // ij across face
        int other_i = cfg.i;
        int other_j = cfg.j;
        int other_k = cfg.k + 1;
        if (other_k >= Config::nkb) { other_k = 0; } // Wrap around.
        int other_id = Config::blk_ids[other_i][other_j][other_k];
        Block& other_blk = fluidBlocks[other_id];
        BConfig& other_cfg = blk_configs[other_id];
        for (int j=0; j < cfg.njc; j++) {
            for (int i=0; i < cfg.nic; i++) {
                FVFace& f = blk.faces[cfg.kFaceIndex(i, j, cfg.nkc)];
                FVFace& other_f = other_blk.faces[other_cfg.kFaceIndex(i, j, 0)];
                blk.cells[f.right_cells[0]].fs = other_blk.cells[other_f.right_cells[0]].fs;
                blk.cells[f.right_cells[1]].fs = other_blk.cells[other_f.right_cells[1]].fs;
            } // end for j
        } // end for k
        break;
    }
    } // end switch()
} // end bc_exchange()


__host__
void bc_inflow(int iblk, int ibc, FlowState& inflow)
// Copy the associated flow state data into the ghost cells.
{
    BConfig& cfg = blk_configs[iblk];
    Block& blk = fluidBlocks[iblk];
    switch (ibc) {
    case Face::iminus: // jk across face
        for (int k=0; k < cfg.nkc; k++) {
            for (int j=0; j < cfg.njc; j++) {
                FVFace& f = blk.faces[cfg.iFaceIndex(0, j, k)];
                blk.cells[f.left_cells[0]].fs = inflow;
                blk.cells[f.left_cells[1]].fs = inflow;
            } // end for j
        } // end for k
        break;
    case Face::iplus: // jk across face
        for (int k=0; k < cfg.nkc; k++) {
            for (int j=0; j < cfg.njc; j++) {
                FVFace& f = blk.faces[cfg.iFaceIndex(cfg.nic, j, k)];
                blk.cells[f.right_cells[0]].fs = inflow;
                blk.cells[f.right_cells[1]].fs = inflow;
            } // end for j
        } // end for k
        break;
    case Face::jminus: // ik across face
        for (int k=0; k < cfg.nkc; k++) {
            for (int i=0; i < cfg.nic; i++) {
                FVFace& f = blk.faces[cfg.jFaceIndex(i, 0, k)];
                blk.cells[f.left_cells[0]].fs = inflow;
                blk.cells[f.left_cells[1]].fs = inflow;
            } // end for i
        } // end for k
        break;
    case Face::jplus: // ik across face
        for (int k=0; k < cfg.nkc; k++) {
            for (int i=0; i < cfg.nic; i++) {
                FVFace& f = blk.faces[cfg.jFaceIndex(i, cfg.njc, k)];
                blk.cells[f.right_cells[0]].fs = inflow;
                blk.cells[f.right_cells[1]].fs = inflow;
            } // end for i
        } // end for k
        break;
    case Face::kminus: // ij across face
        for (int j=0; j < cfg.njc; j++) {
            for (int i=0; i < cfg.nic; i++) {
                FVFace& f = blk.faces[cfg.kFaceIndex(i, j, 0)];
                blk.cells[f.left_cells[0]].fs = inflow;
                blk.cells[f.left_cells[1]].fs = inflow;
            } // end for i
        } // end for j
        break;
    case Face::kplus: // ij across face
        for (int j=0; j < cfg.njc; j++) {
            for (int i=0; i < cfg.nic; i++) {
                FVFace& f = blk.faces[cfg.kFaceIndex(i, j, cfg.nkc)];
                blk.cells[f.right_cells[0]].fs = inflow;
                blk.cells[f.right_cells[1]].fs = inflow;
            } // end for j
        } // end for k
        break;
    }
} // end bc_inflow()


__host__ __device__
FlowState compute_inflow_state_from_function(int ifn, Vector3 pos)
{
    FlowState inflow;
    // Dummy value.
    inflow.gas.p = 100.0e3;
    inflow.gas.T = 300.0;
    inflow.gas.update_from_pT();
    inflow.vel.set(0.0, 0.0, 0.0);
    //
    switch (ifn) {
    case BCFunction::none:
        break;
    case BCFunction::supersonic_vortex:
        // [TODO]
        break;
    case BCFunction::laminar_boundary_layer:
        // [TODO]
        break;
    case BCFunction::manufactured_solution:
        // [TODO]
        break;
    default:
        break;
    }
    return inflow;
}


__host__
void bc_inflow_function(int iblk, int ibc, int ifn)
// Copy the associated flow state data into the ghost cells.
{
    BConfig& cfg = blk_configs[iblk];
    Block& blk = fluidBlocks[iblk];
    //
    switch (ibc) {
    case Face::iminus: // jk across face
        for (int k=0; k < cfg.nkc; k++) {
            for (int j=0; j < cfg.njc; j++) {
                FVFace& f = blk.faces[cfg.iFaceIndex(0, j, k)];
                FVCell& g0 = blk.cells[f.left_cells[0]];
                g0.fs = compute_inflow_state_from_function(ifn, g0.pos);
                FVCell& g1 = blk.cells[f.left_cells[1]];
                g1.fs = compute_inflow_state_from_function(ifn, g1.pos);
            } // end for j
        } // end for k
        break;
    case Face::iplus: // jk across face
        for (int k=0; k < cfg.nkc; k++) {
            for (int j=0; j < cfg.njc; j++) {
                FVFace& f = blk.faces[cfg.iFaceIndex(cfg.nic, j, k)];
                FVCell& g0 = blk.cells[f.right_cells[0]];
                g0.fs = compute_inflow_state_from_function(ifn, g0.pos);
                FVCell& g1 = blk.cells[f.right_cells[1]];
                g1.fs = compute_inflow_state_from_function(ifn, g1.pos);
            } // end for j
        } // end for k
        break;
    case Face::jminus: // ik across face
        for (int k=0; k < cfg.nkc; k++) {
            for (int i=0; i < cfg.nic; i++) {
                FVFace& f = blk.faces[cfg.jFaceIndex(i, 0, k)];
                FVCell& g0 = blk.cells[f.left_cells[0]];
                g0.fs = compute_inflow_state_from_function(ifn, g0.pos);
                FVCell& g1 = blk.cells[f.left_cells[1]];
                g1.fs = compute_inflow_state_from_function(ifn, g1.pos);
            } // end for i
        } // end for k
        break;
    case Face::jplus: // ik across face
        for (int k=0; k < cfg.nkc; k++) {
            for (int i=0; i < cfg.nic; i++) {
                FVFace& f = blk.faces[cfg.jFaceIndex(i, cfg.njc, k)];
                FVCell& g0 = blk.cells[f.right_cells[0]];
                g0.fs = compute_inflow_state_from_function(ifn, g0.pos);
                FVCell& g1 = blk.cells[f.right_cells[1]];
                g1.fs = compute_inflow_state_from_function(ifn, g1.pos);
            } // end for i
        } // end for k
        break;
    case Face::kminus: // ij across face
        for (int j=0; j < cfg.njc; j++) {
            for (int i=0; i < cfg.nic; i++) {
                FVFace& f = blk.faces[cfg.kFaceIndex(i, j, 0)];
                FVCell& g0 = blk.cells[f.left_cells[0]];
                g0.fs = compute_inflow_state_from_function(ifn, g0.pos);
                FVCell& g1 = blk.cells[f.left_cells[1]];
                g1.fs = compute_inflow_state_from_function(ifn, g1.pos);
            } // end for i
        } // end for j
        break;
    case Face::kplus: // ij across face
        for (int j=0; j < cfg.njc; j++) {
            for (int i=0; i < cfg.nic; i++) {
                FVFace& f = blk.faces[cfg.kFaceIndex(i, j, cfg.nkc)];
                FVCell& g0 = blk.cells[f.right_cells[0]];
                g0.fs = compute_inflow_state_from_function(ifn, g0.pos);
                FVCell& g1 = blk.cells[f.right_cells[1]];
                g1.fs = compute_inflow_state_from_function(ifn, g1.pos);
            } // end for j
        } // end for k
        break;
    }
} // end bc_inflow_function()


__host__
void bc_outflow(int iblk, int ibc)
// Copy the interior flow states to the ghost cells.
{
    BConfig& cfg = blk_configs[iblk];
    Block& blk = fluidBlocks[iblk];
    switch (ibc) {
    case Face::iminus: // jk across face
        for (int k=0; k < cfg.nkc; k++) {
            for (int j=0; j < cfg.njc; j++) {
                FVFace& f = blk.faces[cfg.iFaceIndex(0, j, k)];
                FVCell& c = blk.cells[f.right_cells[0]];
                blk.cells[f.left_cells[0]].fs = c.fs;
                blk.cells[f.left_cells[1]].fs = c.fs;
            } // end for j
        } // end for k
        break;
    case Face::iplus: // jk across face
        for (int k=0; k < cfg.nkc; k++) {
            for (int j=0; j < cfg.njc; j++) {
                FVFace& f = blk.faces[cfg.iFaceIndex(cfg.nic, j, k)];
                FVCell& c = blk.cells[f.left_cells[0]];
                blk.cells[f.right_cells[0]].fs = c.fs;
                blk.cells[f.right_cells[1]].fs = c.fs;
            } // end for j
        } // end for k
        break;
    case Face::jminus: // ik across face
        for (int k=0; k < cfg.nkc; k++) {
            for (int i=0; i < cfg.nic; i++) {
                FVFace& f = blk.faces[cfg.jFaceIndex(i, 0, k)];
                FVCell& c = blk.cells[f.right_cells[0]];
                blk.cells[f.left_cells[0]].fs = c.fs;
                blk.cells[f.left_cells[1]].fs = c.fs;
            } // end for i
        } // end for k
        break;
    case Face::jplus: // ik across face
        for (int k=0; k < cfg.nkc; k++) {
            for (int i=0; i < cfg.nic; i++) {
                FVFace& f = blk.faces[cfg.jFaceIndex(i, cfg.njc, k)];
                FVCell& c = blk.cells[f.left_cells[0]];
                blk.cells[f.right_cells[0]].fs = c.fs;
                blk.cells[f.right_cells[1]].fs = c.fs;
            } // end for i
        } // end for k
        break;
    case Face::kminus: // ij across face
        for (int j=0; j < cfg.njc; j++) {
            for (int i=0; i < cfg.nic; i++) {
                FVFace& f = blk.faces[cfg.kFaceIndex(i, j, 0)];
                FVCell& c = blk.cells[f.right_cells[0]];
                blk.cells[f.left_cells[0]].fs = c.fs;
                blk.cells[f.left_cells[1]].fs = c.fs;
            } // end for i
        } // end for j
        break;
    case Face::kplus: // ij across face
        for (int j=0; j < cfg.njc; j++) {
            for (int i=0; i < cfg.nic; i++) {
                FVFace& f = blk.faces[cfg.kFaceIndex(i, j, cfg.nkc)];
                FVCell& c = blk.cells[f.left_cells[0]];
                blk.cells[f.right_cells[0]].fs = c.fs;
                blk.cells[f.right_cells[1]].fs = c.fs;
            } // end for j
        } // end for k
        break;
    }
} // end bc_outflow()


__host__
void apply_boundary_conditions_for_convective_fluxes()
// Work through all fluid blocks and fill in the flow properties for ghost cells
// that sit around the block boundaries.
//
// Since the boundary-condition code needs a view of all blocks and
// most of the coperations are switching between code to copy specific data,
// we expect the CPU to apply the boundary conditions more effectively than the GPU.
// Measurements might tell us otherwise -- and they did (PJ 2022-10-27).
{
    #pragma omp parallel for
    for (int iblk=0; iblk < Config::nFluidBlocks; iblk++) {
        BConfig& blk_config = blk_configs[iblk];
        if (!blk_config.active) continue;
        for (int ibc=0; ibc < 6; ibc++) {
            switch (blk_config.bcCodes[ibc]) {
            case BCCode::wall_with_slip:
            case BCCode::wall_no_slip_adiabatic:
            case BCCode::wall_no_slip_fixed_T:
                bc_wall_reflect_normal_velocity(iblk, ibc);
                break;
            case BCCode::exchange:
                bc_exchange(iblk, ibc);
                break;
            case BCCode::inflow:
                bc_inflow(iblk, ibc, Config::flow_states[blk_config.bc_fs[ibc]]);
                break;
            case BCCode::inflow_function:
                bc_inflow_function(iblk, ibc, blk_config.bc_fun[ibc]);
                break;
            case BCCode::outflow:
                bc_outflow(iblk, ibc);
                break;
            default:
                throw runtime_error("Invalid bcCode: "+to_string(blk_config.bcCodes[ibc]));
            }
        } // end for ibc
    } // end for iblk
} // end apply_boundary_conditions_for_convective_fluxes()


//-------------------------------------------------------------------------------------------------
// Part B.
// Set the ghost-cell quantities to effect the boundary-conditions for convective fluxes on the GPU.

__host__
void configure_exchange_info(vector<Block>& blks, vector<BConfig>& cfgs)
// Set up the per-face information for the exchange boundary condition.
// Do this after all blocks have been configured wecause we need to dip into
// the other block to get the cell indices on the corresponding boundary face.
{
    for (int iblk=0; iblk < Config::nFluidBlocks; ++iblk) {
        BConfig& cfg = cfgs[iblk];
        if (!cfg.active) continue;
        Block& blk = blks[iblk];
        if (cfg.bcCodes[Face::iminus] == BCCode::exchange) { // jk across face
            int other_i = cfg.i - 1;
            if (other_i < 0) { other_i = Config::nib-1; } // Wrap around.
            int other_j = cfg.j;
            int other_k = cfg.k;
            int other_id = Config::blk_ids[other_i][other_j][other_k];
            Block& other_blk = blks[other_id];
            BConfig& other_cfg = cfgs[other_id];
            for (int k=0; k < cfg.nkc; k++) {
                for (int j=0; j < cfg.njc; j++) {
                    FVFace& f = blk.faces[cfg.iFaceIndex(0, j, k)];
                    FVFace& other_f = other_blk.faces[other_cfg.iFaceIndex(other_cfg.nic, j, k)];
                    f.other_blkId = other_id;
                    f.other_cells[0] = other_f.left_cells[0];
                    f.other_cells[1] = other_f.left_cells[1];
                } // end for j
            } // end for k
        }
        if (cfg.bcCodes[Face::iplus] == BCCode::exchange) { // jk across face
            int other_i = cfg.i + 1;
            if (other_i >= Config::nib) { other_i = 0; } // Wrap around.
            int other_j = cfg.j;
            int other_k = cfg.k;
            int other_id = Config::blk_ids[other_i][other_j][other_k];
            Block& other_blk = fluidBlocks[other_id];
            BConfig& other_cfg = blk_configs[other_id];
            for (int k=0; k < cfg.nkc; k++) {
                for (int j=0; j < cfg.njc; j++) {
                    FVFace& f = blk.faces[cfg.iFaceIndex(cfg.nic, j, k)];
                    FVFace& other_f = other_blk.faces[other_cfg.iFaceIndex(0, j, k)];
                    f.other_blkId = other_id;
                    f.other_cells[0] = other_f.right_cells[0];
                    f.other_cells[1] = other_f.right_cells[1];
                } // end for j
            } // end for k
        }
        if (cfg.bcCodes[Face::jminus] == BCCode::exchange) { // ik across face
            int other_i = cfg.i;
            int other_j = cfg.j - 1;
            if (other_j < 0) { other_j = Config::njb-1; } // Wrap around.
            int other_k = cfg.k;
            int other_id = Config::blk_ids[other_i][other_j][other_k];
            Block& other_blk = fluidBlocks[other_id];
            BConfig& other_cfg = blk_configs[other_id];
            for (int k=0; k < cfg.nkc; k++) {
                for (int i=0; i < cfg.nic; i++) {
                    FVFace& f = blk.faces[cfg.jFaceIndex(i, 0, k)];
                    FVFace& other_f = other_blk.faces[other_cfg.jFaceIndex(i, other_cfg.njc, k)];
                    f.other_blkId = other_id;
                    f.other_cells[0] = other_f.left_cells[0];
                    f.other_cells[1] = other_f.left_cells[1];
                } // end for i
            } // end for k
        }
        if (cfg.bcCodes[Face::jplus] == BCCode::exchange) { // ik across face
            int other_i = cfg.i;
            int other_j = cfg.j + 1;
            if (other_j >= Config::njb) { other_j = 0; } // Wrap around.
            int other_k = cfg.k;
            int other_id = Config::blk_ids[other_i][other_j][other_k];
            Block& other_blk = fluidBlocks[other_id];
            BConfig& other_cfg = blk_configs[other_id];
            for (int k=0; k < cfg.nkc; k++) {
                for (int i=0; i < cfg.nic; i++) {
                    FVFace& f = blk.faces[cfg.jFaceIndex(i, cfg.njc, k)];
                    FVFace& other_f = other_blk.faces[other_cfg.jFaceIndex(i, 0, k)];
                    f.other_blkId = other_id;
                    f.other_cells[0] = other_f.right_cells[0];
                    f.other_cells[1] = other_f.right_cells[1];
                } // end for i
            } // end for k
        }
        if (cfg.bcCodes[Face::kminus] == BCCode::exchange) { // ij across face
            int other_i = cfg.i;
            int other_j = cfg.j;
            int other_k = cfg.k - 1;
            if (other_k < 0) { other_k = Config::nkb-1; } // Wrap around.
            int other_id = Config::blk_ids[other_i][other_j][other_k];
            Block& other_blk = fluidBlocks[other_id];
            BConfig& other_cfg = blk_configs[other_id];
            for (int j=0; j < cfg.njc; j++) {
                for (int i=0; i < cfg.nic; i++) {
                    FVFace& f = blk.faces[cfg.kFaceIndex(i, j, 0)];
                    FVFace& other_f = other_blk.faces[other_cfg.kFaceIndex(i, j, other_cfg.nkc)];
                    f.other_blkId = other_id;
                    f.other_cells[0] = other_f.left_cells[0];
                    f.other_cells[1] = other_f.left_cells[1];
                } // end for i
            } // end for j
        }
        if (cfg.bcCodes[Face::kplus] == BCCode::exchange) { // ij across face
            int other_i = cfg.i;
            int other_j = cfg.j;
            int other_k = cfg.k + 1;
            if (other_k >= Config::nkb) { other_k = 0; } // Wrap around.
            int other_id = Config::blk_ids[other_i][other_j][other_k];
            Block& other_blk = fluidBlocks[other_id];
            BConfig& other_cfg = blk_configs[other_id];
            for (int j=0; j < cfg.njc; j++) {
                for (int i=0; i < cfg.nic; i++) {
                    FVFace& f = blk.faces[cfg.kFaceIndex(i, j, cfg.nkc)];
                    FVFace& other_f = other_blk.faces[other_cfg.kFaceIndex(i, j, 0)];
                    f.other_blkId = other_id;
                    f.other_cells[0] = other_f.right_cells[0];
                    f.other_cells[1] = other_f.right_cells[1];
                } // end for j
            } // end for k
        }
    } // end for iblk
} // end configure_exchange_info()


__device__
void apply_convective_boundary_condition(FVFace& f, FVCell cells[],
                                         FlowState flowStates[], Block blks[])
// For a given FVFace, set the ghost cell FlowStates according to the type of boundary condition.
// Input:
// f:          reference to the FVFace
// cells:      array of cells in the current block
// flowStates: array of FlowStates needed by inflow boundary condition
// blks:       array of Block objects needed by the exchange boundary condition
//
{
    if (f.bcCode < 0) return; // Interior face, leave now.
    //
    switch (f.bcCode) {
    case BCCode::wall_with_slip:
    case BCCode::wall_no_slip_adiabatic:
    case BCCode::wall_no_slip_fixed_T: {
        // Copy data, reflecting velocity.
        if (f.bcId == Face::iplus || f.bcId == Face::jplus || f.bcId == Face::kplus) {
            // [TODO] PJ 2022-11-11
            // Note that the following fill of the ghost cell FlowStates reuses the same interior data.
            // It may be better to also use the data from the second interior cell.
            FlowState fs0 = cells[f.left_cells[0]].fs;
            fs0.vel.transform_to_local_frame(f.n, f.t1, f.t2);
            fs0.vel.x = -(fs0.vel.x);
            fs0.vel.transform_to_global_frame(f.n, f.t1, f.t2);
            cells[f.right_cells[0]].fs = fs0;
            cells[f.right_cells[1]].fs = fs0;
        } else {
            FlowState fs0 = cells[f.right_cells[0]].fs;
            fs0.vel.transform_to_local_frame(f.n, f.t1, f.t2);
            fs0.vel.x = -(fs0.vel.x);
            fs0.vel.transform_to_global_frame(f.n, f.t1, f.t2);
            cells[f.left_cells[0]].fs = fs0;
            cells[f.left_cells[1]].fs = fs0;
        }
        break;
    }
    case BCCode::exchange: {
        Block& other_blk = blks[f.other_blkId];
        FVCell* other_cells = other_blk.cells_on_gpu;
        // Note that this function only works from within a kernel.
        // On the CPU, we ould like to do the following but it is not allowed from a __device__ function.
        // FVCell* other_cells = other_blk.cells.data();
        if (f.bcId == Face::iplus || f.bcId == Face::jplus || f.bcId == Face::kplus) {
            cells[f.right_cells[0]].fs = other_cells[f.other_cells[0]].fs;
            cells[f.right_cells[1]].fs = other_cells[f.other_cells[1]].fs;
        } else {
            cells[f.left_cells[0]].fs = other_cells[f.other_cells[0]].fs;
            cells[f.left_cells[1]].fs = other_cells[f.other_cells[1]].fs;
        }
        break;
    }
    case BCCode::inflow: {
        FlowState inflow = flowStates[f.inflowId];
        if (f.bcId == Face::iplus || f.bcId == Face::jplus || f.bcId == Face::kplus) {
            cells[f.right_cells[0]].fs = inflow;
            cells[f.right_cells[1]].fs = inflow;
        } else {
            cells[f.left_cells[0]].fs = inflow;
            cells[f.left_cells[1]].fs = inflow;
        }
        break;
    }
    case BCCode::inflow_function: {
        if (f.bcId == Face::iplus || f.bcId == Face::jplus || f.bcId == Face::kplus) {
            FVCell& g0 = cells[f.right_cells[0]];
            g0.fs = compute_inflow_state_from_function(f.bcFun, g0.pos);
            FVCell& g1 = cells[f.right_cells[1]];
            g1.fs = compute_inflow_state_from_function(f.bcFun, g1.pos);
        } else {
            FVCell& g0 = cells[f.left_cells[0]];
            g0.fs = compute_inflow_state_from_function(f.bcFun, g0.pos);
            FVCell& g1 = cells[f.left_cells[1]];
            g1.fs = compute_inflow_state_from_function(f.bcFun, g1.pos);
        }
        break;
    }
    case BCCode::outflow: {
        if (f.bcId == Face::iplus || f.bcId == Face::jplus || f.bcId == Face::kplus) {
            FlowState fs0 = cells[f.left_cells[0]].fs;
            cells[f.right_cells[0]].fs = fs0;
            cells[f.right_cells[1]].fs = fs0;
        } else {
            FlowState fs0 = cells[f.right_cells[0]].fs;
            cells[f.left_cells[0]].fs = fs0;
            cells[f.left_cells[1]].fs = fs0;
        }
        break;
    }
    default:
        // Do nothing.
        break;
    }
} // end apply_convective_boundary_condition()


//-------------------------------------------------------------------------------------------------
// Part C.
// Set the FlowStates at the actual faces for effecting the boundary conditions for viscous fluxes.


__host__ __device__
void apply_viscous_boundary_condition(FVFace& f)
// Set the FlowState according to the type of boundary condition.
// Will overwrite some of the FlowState properties computed earlier
// in the convective-flux calculation.
{
    switch (f.bcCode) {
    case BCCode::wall_no_slip_adiabatic:
        f.fs.vel.set(0.0, 0.0, 0.0);
        break;
    case BCCode::wall_no_slip_fixed_T:
        f.fs.vel.set(0.0, 0.0, 0.0);
        f.fs.gas.T = f.TWall;
        break;
    case BCCode::inflow_function:
        f.fs = compute_inflow_state_from_function(f.bcFun, f.pos);
    default:
        // Do nothing.
        break;
    }
} // end apply_viscous_boundary_condition()


__host__
void apply_viscous_boundary_conditions(Block& blk)
{
    for (auto& face : blk.faces) {
        apply_viscous_boundary_condition(face);
    }
}

#endif
