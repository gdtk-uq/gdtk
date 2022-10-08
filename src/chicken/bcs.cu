// bcs.cu
// This file is included into simulate.cu, somewhere around the middle of that file.
//
// For a given block and boundary, we work across all faces in the boundary
// and set the ghost-cell flow states appropriate for the boundary condition.
//
// PJ 2022-10-03

#ifndef BCS_INCLUDED
#define BCS_INCLUDED

__host__
void bc_wall_with_slip(int iblk, int ibc)
// Copy data, reflecting velocity.
{
    Block& blk = fluidBlocks[iblk];
    switch (ibc) {
    case Face::iminus: // jk across face
        for (int k=0; k < blk.nkc; k++) {
            for (int j=0; j < blk.njc; j++) {
                FVFace& f = blk.iFaces[blk.iFaceIndex(0, j, k)];
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
        for (int k=0; k < blk.nkc; k++) {
            for (int j=0; j < blk.njc; j++) {
                FVFace& f = blk.iFaces[blk.iFaceIndex(blk.nic, j, k)];
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
        for (int k=0; k < blk.nkc; k++) {
            for (int i=0; i < blk.nic; i++) {
                FVFace& f = blk.jFaces[blk.jFaceIndex(i, 0, k)];
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
        for (int k=0; k < blk.nkc; k++) {
            for (int i=0; i < blk.nic; i++) {
                FVFace& f = blk.jFaces[blk.jFaceIndex(i, blk.njc, k)];
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
        for (int j=0; j < blk.njc; j++) {
            for (int i=0; i < blk.nic; i++) {
                FVFace& f = blk.kFaces[blk.kFaceIndex(i, j, 0)];
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
        for (int j=0; j < blk.njc; j++) {
            for (int i=0; i < blk.nic; i++) {
                FVFace& f = blk.kFaces[blk.kFaceIndex(i, j, blk.nkc)];
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
} // end bc_wall_with_slip()


__host__
void bc_wall_no_slip(int iblk, int ibc)
// Copy data, reflecting velocity, [TODO] then set the face velocity to zero.
{
    Block& blk = fluidBlocks[iblk];
    switch (ibc) {
    case Face::iminus: // jk across face
        for (int k=0; k < blk.nkc; k++) {
            for (int j=0; j < blk.njc; j++) {
                FVFace& f = blk.iFaces[blk.iFaceIndex(0, j, k)];
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
        for (int k=0; k < blk.nkc; k++) {
            for (int j=0; j < blk.njc; j++) {
                FVFace& f = blk.iFaces[blk.iFaceIndex(blk.nic, j, k)];
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
        for (int k=0; k < blk.nkc; k++) {
            for (int i=0; i < blk.nic; i++) {
                FVFace& f = blk.jFaces[blk.jFaceIndex(i, 0, k)];
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
        for (int k=0; k < blk.nkc; k++) {
            for (int i=0; i < blk.nic; i++) {
                FVFace& f = blk.jFaces[blk.jFaceIndex(i, blk.njc, k)];
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
        for (int j=0; j < blk.njc; j++) {
            for (int i=0; i < blk.nic; i++) {
                FVFace& f = blk.kFaces[blk.kFaceIndex(i, j, 0)];
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
        for (int j=0; j < blk.njc; j++) {
            for (int i=0; i < blk.nic; i++) {
                FVFace& f = blk.kFaces[blk.kFaceIndex(i, j, blk.nkc)];
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
} // end bc_wall_no_slip()


__host__
void bc_exchange(int iblk, int ibc)
// Depending on which face, look up the adjacent block in ijk indices and copy the flow data.
// The other block will have a corresponding boundary with the same type of boundary condition,
// so the "exchange" occurs is two phases.
{
    BConfig& blk_config = Config::blk_configs[iblk];
    Block& blk = fluidBlocks[iblk];
    //
    switch (ibc) {
    case Face::iminus: { // jk across face
        int other_i = blk_config.i - 1;
        if (other_i < 0) { other_i = Config::nib-1; } // Wrap around.
        int other_j = blk_config.j;
        int other_k = blk_config.k;
        int other_id = Config::blk_ids[other_i][other_j][other_k];
        Block& other_blk = fluidBlocks[other_id];
        for (int k=0; k < blk.nkc; k++) {
            for (int j=0; j < blk.njc; j++) {
                FVFace& f = blk.iFaces[blk.iFaceIndex(0, j, k)];
                FVFace& other_f = other_blk.iFaces[other_blk.iFaceIndex(other_blk.nic, j, k)];
                blk.cells[f.left_cells[0]].fs = other_blk.cells[other_f.left_cells[0]].fs;
                blk.cells[f.left_cells[1]].fs = other_blk.cells[other_f.left_cells[1]].fs;
            } // end for j
        } // end for k
        break;
    }
    case Face::iplus: { // jk across face
        int other_i = blk_config.i + 1;
        if (other_i >= Config::nib) { other_i = 0; } // Wrap around.
        int other_j = blk_config.j;
        int other_k = blk_config.k;
        int other_id = Config::blk_ids[other_i][other_j][other_k];
        Block& other_blk = fluidBlocks[other_id];
        for (int k=0; k < blk.nkc; k++) {
            for (int j=0; j < blk.njc; j++) {
                FVFace& f = blk.iFaces[blk.iFaceIndex(blk.nic, j, k)];
                FVFace& other_f = other_blk.iFaces[other_blk.iFaceIndex(0, j, k)];
                blk.cells[f.right_cells[0]].fs = other_blk.cells[other_f.right_cells[0]].fs;
                blk.cells[f.right_cells[1]].fs = other_blk.cells[other_f.right_cells[1]].fs;
            } // end for j
        } // end for k
        break;
    }
    case Face::jminus: { // ik across face
        int other_i = blk_config.i;
        int other_j = blk_config.j - 1;
        if (other_j < 0) { other_j = Config::njb-1; } // Wrap around.
        int other_k = blk_config.k;
        int other_id = Config::blk_ids[other_i][other_j][other_k];
        Block& other_blk = fluidBlocks[other_id];
        for (int k=0; k < blk.nkc; k++) {
            for (int i=0; i < blk.nic; i++) {
                FVFace& f = blk.jFaces[blk.jFaceIndex(i, 0, k)];
                FVFace& other_f = other_blk.jFaces[other_blk.jFaceIndex(i, other_blk.njc, k)];
                blk.cells[f.left_cells[0]].fs = other_blk.cells[other_f.left_cells[0]].fs;
                blk.cells[f.left_cells[1]].fs = other_blk.cells[other_f.left_cells[1]].fs;
            } // end for i
        } // end for k
        break;
    }
    case Face::jplus: { // ik across face
        int other_i = blk_config.i;
        int other_j = blk_config.j + 1;
        if (other_j >= Config::njb) { other_j = 0; } // Wrap around.
        int other_k = blk_config.k;
        int other_id = Config::blk_ids[other_i][other_j][other_k];
        Block& other_blk = fluidBlocks[other_id];
        for (int k=0; k < blk.nkc; k++) {
            for (int i=0; i < blk.nic; i++) {
                FVFace& f = blk.jFaces[blk.jFaceIndex(i, blk.njc, k)];
                FVFace& other_f = other_blk.jFaces[other_blk.jFaceIndex(i, 0, k)];
                blk.cells[f.right_cells[0]].fs = other_blk.cells[other_f.right_cells[0]].fs;
                blk.cells[f.right_cells[1]].fs = other_blk.cells[other_f.right_cells[1]].fs;
            } // end for i
        } // end for k
        break;
    }
    case Face::kminus: { // ij across face
        int other_i = blk_config.i;
        int other_j = blk_config.j;
        int other_k = blk_config.k - 1;
        if (other_k < 0) { other_k = Config::nkb-1; } // Wrap around.
        int other_id = Config::blk_ids[other_i][other_j][other_k];
        Block& other_blk = fluidBlocks[other_id];
        for (int j=0; j < blk.njc; j++) {
            for (int i=0; i < blk.nic; i++) {
                FVFace& f = blk.kFaces[blk.kFaceIndex(i, j, 0)];
                FVFace& other_f = other_blk.kFaces[other_blk.kFaceIndex(i, j, other_blk.nkc)];
                blk.cells[f.left_cells[0]].fs = other_blk.cells[other_f.left_cells[0]].fs;
                blk.cells[f.left_cells[1]].fs = other_blk.cells[other_f.left_cells[1]].fs;
            } // end for i
        } // end for j
        break;
    }
    case Face::kplus: { // ij across face
        int other_i = blk_config.i;
        int other_j = blk_config.j;
        int other_k = blk_config.k + 1;
        if (other_k >= Config::nkb) { other_k = 0; } // Wrap around.
        int other_id = Config::blk_ids[other_i][other_j][other_k];
        Block& other_blk = fluidBlocks[other_id];
        for (int j=0; j < blk.njc; j++) {
            for (int i=0; i < blk.nic; i++) {
                FVFace& f = blk.kFaces[blk.kFaceIndex(i, j, blk.nkc)];
                FVFace& other_f = other_blk.kFaces[other_blk.kFaceIndex(i, j, 0)];
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
    Block& blk = fluidBlocks[iblk];
    switch (ibc) {
    case Face::iminus: // jk across face
        for (int k=0; k < blk.nkc; k++) {
            for (int j=0; j < blk.njc; j++) {
                FVFace& f = blk.iFaces[blk.iFaceIndex(0, j, k)];
                FlowState& fs0 = blk.cells[f.left_cells[0]].fs;
                fs0 = inflow;
                FlowState& fs1 = blk.cells[f.left_cells[1]].fs;
                fs1 = inflow;
            } // end for j
        } // end for k
        break;
    case Face::iplus: // jk across face
        for (int k=0; k < blk.nkc; k++) {
            for (int j=0; j < blk.njc; j++) {
                FVFace& f = blk.iFaces[blk.iFaceIndex(blk.nic, j, k)];
                FlowState& fs0 = blk.cells[f.right_cells[0]].fs;
                fs0 = inflow;
                FlowState& fs1 = blk.cells[f.right_cells[1]].fs;
                fs1 = inflow;
            } // end for j
        } // end for k
        break;
    case Face::jminus: // ik across face
        for (int k=0; k < blk.nkc; k++) {
            for (int i=0; i < blk.nic; i++) {
                FVFace& f = blk.jFaces[blk.jFaceIndex(i, 0, k)];
                FlowState& fs0 = blk.cells[f.left_cells[0]].fs;
                fs0 = inflow;
                FlowState& fs1 = blk.cells[f.left_cells[1]].fs;
                fs1 = inflow;
            } // end for i
        } // end for k
        break;
    case Face::jplus: // ik across face
        for (int k=0; k < blk.nkc; k++) {
            for (int i=0; i < blk.nic; i++) {
                FVFace& f = blk.jFaces[blk.jFaceIndex(i, blk.njc, k)];
                FlowState& fs0 = blk.cells[f.right_cells[0]].fs;
                fs0 = inflow;
                FlowState& fs1 = blk.cells[f.right_cells[1]].fs;
                fs1 = inflow;
            } // end for i
        } // end for k
        break;
    case Face::kminus: // ij across face
        for (int j=0; j < blk.njc; j++) {
            for (int i=0; i < blk.nic; i++) {
                FVFace& f = blk.kFaces[blk.kFaceIndex(i, j, 0)];
                FlowState& fs0 = blk.cells[f.left_cells[0]].fs;
                fs0 = inflow;
                FlowState& fs1 = blk.cells[f.left_cells[1]].fs;
                fs1 = inflow;
            } // end for i
        } // end for j
        break;
    case Face::kplus: // ij across face
        for (int j=0; j < blk.njc; j++) {
            for (int i=0; i < blk.nic; i++) {
                FVFace& f = blk.kFaces[blk.kFaceIndex(i, j, blk.nkc)];
                FlowState& fs0 = blk.cells[f.right_cells[0]].fs;
                fs0 = inflow;
                FlowState& fs1 = blk.cells[f.right_cells[1]].fs;
                fs1 = inflow;
            } // end for j
        } // end for k
        break;
    }
} // end bc_inflow()


__host__
void bc_outflow(int iblk, int ibc)
// Copy the interior flow states to the ghost cells.
{
    Block& blk = fluidBlocks[iblk];
    switch (ibc) {
    case Face::iminus: // jk across face
        for (int k=0; k < blk.nkc; k++) {
            for (int j=0; j < blk.njc; j++) {
                FVFace& f = blk.iFaces[blk.iFaceIndex(0, j, k)];
                FVCell& c = blk.cells[f.right_cells[0]];
                FlowState& fs0 = blk.cells[f.left_cells[0]].fs;
                fs0 = c.fs;
                FlowState& fs1 = blk.cells[f.left_cells[1]].fs;
                fs1 = c.fs;
            } // end for j
        } // end for k
        break;
    case Face::iplus: // jk across face
        for (int k=0; k < blk.nkc; k++) {
            for (int j=0; j < blk.njc; j++) {
                FVFace& f = blk.iFaces[blk.iFaceIndex(blk.nic, j, k)];
                FVCell& c = blk.cells[f.left_cells[0]];
                FlowState& fs0 = blk.cells[f.right_cells[0]].fs;
                fs0 = c.fs;
                FlowState& fs1 = blk.cells[f.right_cells[1]].fs;
                fs1 = c.fs;
            } // end for j
        } // end for k
        break;
    case Face::jminus: // ik across face
        for (int k=0; k < blk.nkc; k++) {
            for (int i=0; i < blk.nic; i++) {
                FVFace& f = blk.jFaces[blk.jFaceIndex(i, 0, k)];
                FVCell& c = blk.cells[f.right_cells[0]];
                FlowState& fs0 = blk.cells[f.left_cells[0]].fs;
                fs0 = c.fs;
                FlowState& fs1 = blk.cells[f.left_cells[1]].fs;
                fs1 = c.fs;
            } // end for i
        } // end for k
        break;
    case Face::jplus: // ik across face
        for (int k=0; k < blk.nkc; k++) {
            for (int i=0; i < blk.nic; i++) {
                FVFace& f = blk.jFaces[blk.jFaceIndex(i, blk.njc, k)];
                FVCell& c = blk.cells[f.left_cells[0]];
                FlowState& fs0 = blk.cells[f.right_cells[0]].fs;
                fs0 = c.fs;
                FlowState& fs1 = blk.cells[f.right_cells[1]].fs;
                fs1 = c.fs;
            } // end for i
        } // end for k
        break;
    case Face::kminus: // ij across face
        for (int j=0; j < blk.njc; j++) {
            for (int i=0; i < blk.nic; i++) {
                FVFace& f = blk.kFaces[blk.kFaceIndex(i, j, 0)];
                FVCell& c = blk.cells[f.right_cells[0]];
                FlowState& fs0 = blk.cells[f.left_cells[0]].fs;
                fs0 = c.fs;
                FlowState& fs1 = blk.cells[f.left_cells[1]].fs;
                fs1 = c.fs;
            } // end for i
        } // end for j
        break;
    case Face::kplus: // ij across face
        for (int j=0; j < blk.njc; j++) {
            for (int i=0; i < blk.nic; i++) {
                FVFace& f = blk.kFaces[blk.kFaceIndex(i, j, blk.nkc)];
                FVCell& c = blk.cells[f.left_cells[0]];
                FlowState& fs0 = blk.cells[f.right_cells[0]].fs;
                fs0 = c.fs;
                FlowState& fs1 = blk.cells[f.right_cells[1]].fs;
                fs1 = c.fs;
            } // end for j
        } // end for k
        break;
    }
} // end bc_outflow()

#endif
