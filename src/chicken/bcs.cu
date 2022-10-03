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
void bc_wall_with_slip(Block* blk_ptr, int ibc)
// Copy data, reflecting velocity.
{
    switch (ibc) {
    case Face::iminus: // jk across face
        for (int k=0; k < blk_ptr->nkc; k++) {
            for (int j=0; j < blk_ptr->njc; j++) {
                FVFace& f = blk_ptr->iFaces[blk_ptr->iFaceIndex(0, j, k)];
                FVCell& c = blk_ptr->cells[f.right_cells[0]];
                FlowState& fs0 = blk_ptr->cells[f.left_cells[0]].fs;
                fs0 = c.fs;
                fs0.vel.transform_to_local_frame(f.n, f.t1, f.t2);
                fs0.vel.x = -(fs0.vel.x);
                fs0.vel.transform_to_global_frame(f.n, f.t1, f.t2);
                FlowState& fs1 = blk_ptr->cells[f.left_cells[1]].fs;
                fs1 = c.fs;
                fs1.vel.transform_to_local_frame(f.n, f.t1, f.t2);
                fs1.vel.x = -(fs1.vel.x);
                fs1.vel.transform_to_global_frame(f.n, f.t1, f.t2);
            } // end for j
        } // end for k
        break;
    case Face::iplus: // jk across face
        for (int k=0; k < blk_ptr->nkc; k++) {
            for (int j=0; j < blk_ptr->njc; j++) {
                FVFace& f = blk_ptr->iFaces[blk_ptr->iFaceIndex(blk_ptr->nic, j, k)];
                FVCell& c = blk_ptr->cells[f.left_cells[0]];
                FlowState& fs0 = blk_ptr->cells[f.right_cells[0]].fs;
                fs0 = c.fs;
                fs0.vel.transform_to_local_frame(f.n, f.t1, f.t2);
                fs0.vel.x = -(fs0.vel.x);
                fs0.vel.transform_to_global_frame(f.n, f.t1, f.t2);
                FlowState& fs1 = blk_ptr->cells[f.right_cells[1]].fs;
                fs1 = c.fs;
                fs1.vel.transform_to_local_frame(f.n, f.t1, f.t2);
                fs1.vel.x = -(fs1.vel.x);
                fs1.vel.transform_to_global_frame(f.n, f.t1, f.t2);
            } // end for j
        } // end for k
        break;
    case Face::jminus: // ik across face
        for (int k=0; k < blk_ptr->nkc; k++) {
            for (int i=0; i < blk_ptr->nic; i++) {
                FVFace& f = blk_ptr->iFaces[blk_ptr->iFaceIndex(i, 0, k)];
                FVCell& c = blk_ptr->cells[f.right_cells[0]];
                FlowState& fs0 = blk_ptr->cells[f.left_cells[0]].fs;
                fs0 = c.fs;
                fs0.vel.transform_to_local_frame(f.n, f.t1, f.t2);
                fs0.vel.x = -(fs0.vel.x);
                fs0.vel.transform_to_global_frame(f.n, f.t1, f.t2);
                FlowState& fs1 = blk_ptr->cells[f.left_cells[1]].fs;
                fs1 = c.fs;
                fs1.vel.transform_to_local_frame(f.n, f.t1, f.t2);
                fs1.vel.x = -(fs1.vel.x);
                fs1.vel.transform_to_global_frame(f.n, f.t1, f.t2);
            } // end for i
        } // end for k
        break;
    case Face::jplus: // ik across face
        for (int k=0; k < blk_ptr->nkc; k++) {
            for (int i=0; i < blk_ptr->nic; i++) {
                FVFace& f = blk_ptr->iFaces[blk_ptr->iFaceIndex(i, blk_ptr->njc, k)];
                FVCell& c = blk_ptr->cells[f.left_cells[0]];
                FlowState& fs0 = blk_ptr->cells[f.right_cells[0]].fs;
                fs0 = c.fs;
                fs0.vel.transform_to_local_frame(f.n, f.t1, f.t2);
                fs0.vel.x = -(fs0.vel.x);
                fs0.vel.transform_to_global_frame(f.n, f.t1, f.t2);
                FlowState& fs1 = blk_ptr->cells[f.right_cells[1]].fs;
                fs1 = c.fs;
                fs1.vel.transform_to_local_frame(f.n, f.t1, f.t2);
                fs1.vel.x = -(fs1.vel.x);
                fs1.vel.transform_to_global_frame(f.n, f.t1, f.t2);
            } // end for i
        } // end for k
        break;
    case Face::kminus: // ij across face
        for (int j=0; j < blk_ptr->njc; j++) {
            for (int i=0; i < blk_ptr->nic; i++) {
                FVFace& f = blk_ptr->iFaces[blk_ptr->iFaceIndex(i, j, 0)];
                FVCell& c = blk_ptr->cells[f.right_cells[0]];
                FlowState& fs0 = blk_ptr->cells[f.left_cells[0]].fs;
                fs0 = c.fs;
                fs0.vel.transform_to_local_frame(f.n, f.t1, f.t2);
                fs0.vel.x = -(fs0.vel.x);
                fs0.vel.transform_to_global_frame(f.n, f.t1, f.t2);
                FlowState& fs1 = blk_ptr->cells[f.left_cells[1]].fs;
                fs1 = c.fs;
                fs1.vel.transform_to_local_frame(f.n, f.t1, f.t2);
                fs1.vel.x = -(fs1.vel.x);
                fs1.vel.transform_to_global_frame(f.n, f.t1, f.t2);
            } // end for i
        } // end for j
        break;
    case Face::kplus: // ij across face
        for (int j=0; j < blk_ptr->njc; j++) {
            for (int i=0; i < blk_ptr->nic; i++) {
                FVFace& f = blk_ptr->iFaces[blk_ptr->iFaceIndex(i, j, blk_ptr->nkc)];
                FVCell& c = blk_ptr->cells[f.left_cells[0]];
                FlowState& fs0 = blk_ptr->cells[f.right_cells[0]].fs;
                fs0 = c.fs;
                fs0.vel.transform_to_local_frame(f.n, f.t1, f.t2);
                fs0.vel.x = -(fs0.vel.x);
                fs0.vel.transform_to_global_frame(f.n, f.t1, f.t2);
                FlowState& fs1 = blk_ptr->cells[f.right_cells[1]].fs;
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
void bc_wall_no_slip(Block* blk_ptr, int ibc)
// Copy data, reflecting velocity, [TODO] then set the face velocity to zero.
{
    switch (ibc) {
    case Face::iminus: // jk across face
        for (int k=0; k < blk_ptr->nkc; k++) {
            for (int j=0; j < blk_ptr->njc; j++) {
                FVFace& f = blk_ptr->iFaces[blk_ptr->iFaceIndex(0, j, k)];
                FVCell& c = blk_ptr->cells[f.right_cells[0]];
                FlowState& fs0 = blk_ptr->cells[f.left_cells[0]].fs;
                fs0 = c.fs;
                fs0.vel.transform_to_local_frame(f.n, f.t1, f.t2);
                fs0.vel.x = -(fs0.vel.x);
                fs0.vel.transform_to_global_frame(f.n, f.t1, f.t2);
                FlowState& fs1 = blk_ptr->cells[f.left_cells[1]].fs;
                fs1 = c.fs;
                fs1.vel.transform_to_local_frame(f.n, f.t1, f.t2);
                fs1.vel.x = -(fs1.vel.x);
                fs1.vel.transform_to_global_frame(f.n, f.t1, f.t2);
            } // end for j
        } // end for k
        break;
    case Face::iplus: // jk across face
        for (int k=0; k < blk_ptr->nkc; k++) {
            for (int j=0; j < blk_ptr->njc; j++) {
                FVFace& f = blk_ptr->iFaces[blk_ptr->iFaceIndex(blk_ptr->nic, j, k)];
                FVCell& c = blk_ptr->cells[f.left_cells[0]];
                FlowState& fs0 = blk_ptr->cells[f.right_cells[0]].fs;
                fs0 = c.fs;
                fs0.vel.transform_to_local_frame(f.n, f.t1, f.t2);
                fs0.vel.x = -(fs0.vel.x);
                fs0.vel.transform_to_global_frame(f.n, f.t1, f.t2);
                FlowState& fs1 = blk_ptr->cells[f.right_cells[1]].fs;
                fs1 = c.fs;
                fs1.vel.transform_to_local_frame(f.n, f.t1, f.t2);
                fs1.vel.x = -(fs1.vel.x);
                fs1.vel.transform_to_global_frame(f.n, f.t1, f.t2);
            } // end for j
        } // end for k
        break;
    case Face::jminus: // ik across face
        for (int k=0; k < blk_ptr->nkc; k++) {
            for (int i=0; i < blk_ptr->nic; i++) {
                FVFace& f = blk_ptr->iFaces[blk_ptr->iFaceIndex(i, 0, k)];
                FVCell& c = blk_ptr->cells[f.right_cells[0]];
                FlowState& fs0 = blk_ptr->cells[f.left_cells[0]].fs;
                fs0 = c.fs;
                fs0.vel.transform_to_local_frame(f.n, f.t1, f.t2);
                fs0.vel.x = -(fs0.vel.x);
                fs0.vel.transform_to_global_frame(f.n, f.t1, f.t2);
                FlowState& fs1 = blk_ptr->cells[f.left_cells[1]].fs;
                fs1 = c.fs;
                fs1.vel.transform_to_local_frame(f.n, f.t1, f.t2);
                fs1.vel.x = -(fs1.vel.x);
                fs1.vel.transform_to_global_frame(f.n, f.t1, f.t2);
            } // end for i
        } // end for k
        break;
    case Face::jplus: // ik across face
        for (int k=0; k < blk_ptr->nkc; k++) {
            for (int i=0; i < blk_ptr->nic; i++) {
                FVFace& f = blk_ptr->iFaces[blk_ptr->iFaceIndex(i, blk_ptr->njc, k)];
                FVCell& c = blk_ptr->cells[f.left_cells[0]];
                FlowState& fs0 = blk_ptr->cells[f.right_cells[0]].fs;
                fs0 = c.fs;
                fs0.vel.transform_to_local_frame(f.n, f.t1, f.t2);
                fs0.vel.x = -(fs0.vel.x);
                fs0.vel.transform_to_global_frame(f.n, f.t1, f.t2);
                FlowState& fs1 = blk_ptr->cells[f.right_cells[1]].fs;
                fs1 = c.fs;
                fs1.vel.transform_to_local_frame(f.n, f.t1, f.t2);
                fs1.vel.x = -(fs1.vel.x);
                fs1.vel.transform_to_global_frame(f.n, f.t1, f.t2);
            } // end for i
        } // end for k
        break;
    case Face::kminus: // ij across face
        for (int j=0; j < blk_ptr->njc; j++) {
            for (int i=0; i < blk_ptr->nic; i++) {
                FVFace& f = blk_ptr->iFaces[blk_ptr->iFaceIndex(i, j, 0)];
                FVCell& c = blk_ptr->cells[f.right_cells[0]];
                FlowState& fs0 = blk_ptr->cells[f.left_cells[0]].fs;
                fs0 = c.fs;
                fs0.vel.transform_to_local_frame(f.n, f.t1, f.t2);
                fs0.vel.x = -(fs0.vel.x);
                fs0.vel.transform_to_global_frame(f.n, f.t1, f.t2);
                FlowState& fs1 = blk_ptr->cells[f.left_cells[1]].fs;
                fs1 = c.fs;
                fs1.vel.transform_to_local_frame(f.n, f.t1, f.t2);
                fs1.vel.x = -(fs1.vel.x);
                fs1.vel.transform_to_global_frame(f.n, f.t1, f.t2);
            } // end for i
        } // end for j
        break;
    case Face::kplus: // ij across face
        for (int j=0; j < blk_ptr->njc; j++) {
            for (int i=0; i < blk_ptr->nic; i++) {
                FVFace& f = blk_ptr->iFaces[blk_ptr->iFaceIndex(i, j, blk_ptr->nkc)];
                FVCell& c = blk_ptr->cells[f.left_cells[0]];
                FlowState& fs0 = blk_ptr->cells[f.right_cells[0]].fs;
                fs0 = c.fs;
                fs0.vel.transform_to_local_frame(f.n, f.t1, f.t2);
                fs0.vel.x = -(fs0.vel.x);
                fs0.vel.transform_to_global_frame(f.n, f.t1, f.t2);
                FlowState& fs1 = blk_ptr->cells[f.right_cells[1]].fs;
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
    auto* blk_config = &(Config::blk_configs[iblk]);
    auto* blk_ptr = fluidBlocks[iblk];
    //
    switch (ibc) {
    case Face::iminus: { // jk across face
        int other_i = blk_config->i - 1;
        if (other_i < 0) { other_i = Config::nib; } // Wrap around.
        int other_j = blk_config->j;
        int other_k = blk_config->k;
        int other_id = Config::blk_ids[other_i][other_j][other_k];
        Block* other_ptr = fluidBlocks[other_id];
        for (int k=0; k < blk_ptr->nkc; k++) {
            for (int j=0; j < blk_ptr->njc; j++) {
                FVFace& f = blk_ptr->iFaces[blk_ptr->iFaceIndex(0, j, k)];
                FVFace& other_f = blk_ptr->iFaces[blk_ptr->iFaceIndex(blk_ptr->nic, j, k)];
                blk_ptr->cells[f.left_cells[0]].fs = other_ptr->cells[f.left_cells[0]].fs;
                blk_ptr->cells[f.left_cells[1]].fs = other_ptr->cells[f.left_cells[1]].fs;
            } // end for j
        } // end for k
        break;
    }
    case Face::iplus: { // jk across face
        int other_i = blk_config->i + 1;
        if (other_i >= Config::nib) { other_i = 0; } // Wrap around.
        int other_j = blk_config->j;
        int other_k = blk_config->k;
        int other_id = Config::blk_ids[other_i][other_j][other_k];
        Block* other_ptr = fluidBlocks[other_id];
        for (int k=0; k < blk_ptr->nkc; k++) {
            for (int j=0; j < blk_ptr->njc; j++) {
                FVFace& f = blk_ptr->iFaces[blk_ptr->iFaceIndex(blk_ptr->nic, j, k)];
                FVFace& other_f = blk_ptr->iFaces[blk_ptr->iFaceIndex(0, j, k)];
                blk_ptr->cells[f.right_cells[0]].fs = other_ptr->cells[f.right_cells[0]].fs;
                blk_ptr->cells[f.right_cells[1]].fs = other_ptr->cells[f.right_cells[1]].fs;
            } // end for j
        } // end for k
        break;
    }
    case Face::jminus: { // ik across face
        int other_i = blk_config->i;
        int other_j = blk_config->j - 1;
        if (other_j < 0) { other_j = Config::njb; } // Wrap around.
        int other_k = blk_config->k;
        int other_id = Config::blk_ids[other_i][other_j][other_k];
        Block* other_ptr = fluidBlocks[other_id];
        for (int k=0; k < blk_ptr->nkc; k++) {
            for (int i=0; i < blk_ptr->nic; i++) {
                FVFace& f = blk_ptr->iFaces[blk_ptr->iFaceIndex(i, 0, k)];
                FVFace& other_f = blk_ptr->iFaces[blk_ptr->iFaceIndex(i, blk_ptr->njc, k)];
                blk_ptr->cells[f.left_cells[0]].fs = other_ptr->cells[f.left_cells[0]].fs;
                blk_ptr->cells[f.left_cells[1]].fs = other_ptr->cells[f.left_cells[1]].fs;
            } // end for i
        } // end for k
        break;
    }
    case Face::jplus: { // ik across face
        int other_i = blk_config->i;
        int other_j = blk_config->j + 1;
        if (other_j >= Config::njb) { other_j = 0; } // Wrap around.
        int other_k = blk_config->k;
        int other_id = Config::blk_ids[other_i][other_j][other_k];
        Block* other_ptr = fluidBlocks[other_id];
        for (int k=0; k < blk_ptr->nkc; k++) {
            for (int i=0; i < blk_ptr->nic; i++) {
                FVFace& f = blk_ptr->iFaces[blk_ptr->iFaceIndex(i, blk_ptr->njc, k)];
                FVFace& other_f = blk_ptr->iFaces[blk_ptr->iFaceIndex(i, 0, k)];
                blk_ptr->cells[f.right_cells[0]].fs = other_ptr->cells[f.right_cells[0]].fs;
                blk_ptr->cells[f.right_cells[1]].fs = other_ptr->cells[f.right_cells[1]].fs;
            } // end for i
        } // end for k
        break;
    }
    case Face::kminus: { // ij across face
        int other_i = blk_config->i;
        int other_j = blk_config->j;
        int other_k = blk_config->k - 1;
        if (other_k < 0) { other_k = Config::nkb; } // Wrap around.
        int other_id = Config::blk_ids[other_i][other_j][other_k];
        Block* other_ptr = fluidBlocks[other_id];
        for (int j=0; j < blk_ptr->njc; j++) {
            for (int i=0; i < blk_ptr->nic; i++) {
                FVFace& f = blk_ptr->iFaces[blk_ptr->iFaceIndex(i, j, 0)];
                FVFace& other_f = blk_ptr->iFaces[blk_ptr->iFaceIndex(i, j, blk_ptr->nkc)];
                blk_ptr->cells[f.left_cells[0]].fs = other_ptr->cells[f.left_cells[0]].fs;
                blk_ptr->cells[f.left_cells[1]].fs = other_ptr->cells[f.left_cells[1]].fs;
            } // end for i
        } // end for j
        break;
    }
    case Face::kplus: { // ij across face
        int other_i = blk_config->i;
        int other_j = blk_config->j;
        int other_k = blk_config->k + 1;
        if (other_k >= Config::nkb) { other_k = 0; } // Wrap around.
        int other_id = Config::blk_ids[other_i][other_j][other_k];
        Block* other_ptr = fluidBlocks[other_id];
        for (int j=0; j < blk_ptr->njc; j++) {
            for (int i=0; i < blk_ptr->nic; i++) {
                FVFace& f = blk_ptr->iFaces[blk_ptr->iFaceIndex(i, j, blk_ptr->nkc)];
                FVFace& other_f = blk_ptr->iFaces[blk_ptr->iFaceIndex(i, j, 0)];
                blk_ptr->cells[f.right_cells[0]].fs = other_ptr->cells[f.right_cells[0]].fs;
                blk_ptr->cells[f.right_cells[1]].fs = other_ptr->cells[f.right_cells[1]].fs;
            } // end for j
        } // end for k
        break;
    }
    } // end switch()
} // end bc_exchange()


__host__
void bc_inflow(Block* blk_ptr, int ibc, FlowState& inflow)
// Copy the associated flow state data into the ghost cells.
{
    switch (ibc) {
    case Face::iminus: // jk across face
        for (int k=0; k < blk_ptr->nkc; k++) {
            for (int j=0; j < blk_ptr->njc; j++) {
                FVFace& f = blk_ptr->iFaces[blk_ptr->iFaceIndex(0, j, k)];
                FlowState& fs0 = blk_ptr->cells[f.left_cells[0]].fs;
                fs0 = inflow;
                FlowState& fs1 = blk_ptr->cells[f.left_cells[1]].fs;
                fs1 = inflow;
            } // end for j
        } // end for k
        break;
    case Face::iplus: // jk across face
        for (int k=0; k < blk_ptr->nkc; k++) {
            for (int j=0; j < blk_ptr->njc; j++) {
                FVFace& f = blk_ptr->iFaces[blk_ptr->iFaceIndex(blk_ptr->nic, j, k)];
                FlowState& fs0 = blk_ptr->cells[f.right_cells[0]].fs;
                fs0 = inflow;
                FlowState& fs1 = blk_ptr->cells[f.right_cells[1]].fs;
                fs1 = inflow;
            } // end for j
        } // end for k
        break;
    case Face::jminus: // ik across face
        for (int k=0; k < blk_ptr->nkc; k++) {
            for (int i=0; i < blk_ptr->nic; i++) {
                FVFace& f = blk_ptr->iFaces[blk_ptr->iFaceIndex(i, 0, k)];
                FlowState& fs0 = blk_ptr->cells[f.left_cells[0]].fs;
                fs0 = inflow;
                FlowState& fs1 = blk_ptr->cells[f.left_cells[1]].fs;
                fs1 = inflow;
            } // end for i
        } // end for k
        break;
    case Face::jplus: // ik across face
        for (int k=0; k < blk_ptr->nkc; k++) {
            for (int i=0; i < blk_ptr->nic; i++) {
                FVFace& f = blk_ptr->iFaces[blk_ptr->iFaceIndex(i, blk_ptr->njc, k)];
                FlowState& fs0 = blk_ptr->cells[f.right_cells[0]].fs;
                fs0 = inflow;
                FlowState& fs1 = blk_ptr->cells[f.right_cells[1]].fs;
                fs1 = inflow;
            } // end for i
        } // end for k
        break;
    case Face::kminus: // ij across face
        for (int j=0; j < blk_ptr->njc; j++) {
            for (int i=0; i < blk_ptr->nic; i++) {
                FVFace& f = blk_ptr->iFaces[blk_ptr->iFaceIndex(i, j, 0)];
                FlowState& fs0 = blk_ptr->cells[f.left_cells[0]].fs;
                fs0 = inflow;
                FlowState& fs1 = blk_ptr->cells[f.left_cells[1]].fs;
                fs1 = inflow;
            } // end for i
        } // end for j
        break;
    case Face::kplus: // ij across face
        for (int j=0; j < blk_ptr->njc; j++) {
            for (int i=0; i < blk_ptr->nic; i++) {
                FVFace& f = blk_ptr->iFaces[blk_ptr->iFaceIndex(i, j, blk_ptr->nkc)];
                FlowState& fs0 = blk_ptr->cells[f.right_cells[0]].fs;
                fs0 = inflow;
                FlowState& fs1 = blk_ptr->cells[f.right_cells[1]].fs;
                fs1 = inflow;
            } // end for j
        } // end for k
        break;
    }
} // end bc_inflow()


__host__
void bc_outflow(Block* blk_ptr, int ibc)
// Copy the interior flow states to the ghost cells.
{
    switch (ibc) {
    case Face::iminus: // jk across face
        for (int k=0; k < blk_ptr->nkc; k++) {
            for (int j=0; j < blk_ptr->njc; j++) {
                FVFace& f = blk_ptr->iFaces[blk_ptr->iFaceIndex(0, j, k)];
                FVCell& c = blk_ptr->cells[f.right_cells[0]];
                FlowState& fs0 = blk_ptr->cells[f.left_cells[0]].fs;
                fs0 = c.fs;
                FlowState& fs1 = blk_ptr->cells[f.left_cells[1]].fs;
                fs1 = c.fs;
            } // end for j
        } // end for k
        break;
    case Face::iplus: // jk across face
        for (int k=0; k < blk_ptr->nkc; k++) {
            for (int j=0; j < blk_ptr->njc; j++) {
                FVFace& f = blk_ptr->iFaces[blk_ptr->iFaceIndex(blk_ptr->nic, j, k)];
                FVCell& c = blk_ptr->cells[f.left_cells[0]];
                FlowState& fs0 = blk_ptr->cells[f.right_cells[0]].fs;
                fs0 = c.fs;
                FlowState& fs1 = blk_ptr->cells[f.right_cells[1]].fs;
                fs1 = c.fs;
            } // end for j
        } // end for k
        break;
    case Face::jminus: // ik across face
        for (int k=0; k < blk_ptr->nkc; k++) {
            for (int i=0; i < blk_ptr->nic; i++) {
                FVFace& f = blk_ptr->iFaces[blk_ptr->iFaceIndex(i, 0, k)];
                FVCell& c = blk_ptr->cells[f.right_cells[0]];
                FlowState& fs0 = blk_ptr->cells[f.left_cells[0]].fs;
                fs0 = c.fs;
                FlowState& fs1 = blk_ptr->cells[f.left_cells[1]].fs;
                fs1 = c.fs;
            } // end for i
        } // end for k
        break;
    case Face::jplus: // ik across face
        for (int k=0; k < blk_ptr->nkc; k++) {
            for (int i=0; i < blk_ptr->nic; i++) {
                FVFace& f = blk_ptr->iFaces[blk_ptr->iFaceIndex(i, blk_ptr->njc, k)];
                FVCell& c = blk_ptr->cells[f.left_cells[0]];
                FlowState& fs0 = blk_ptr->cells[f.right_cells[0]].fs;
                fs0 = c.fs;
                FlowState& fs1 = blk_ptr->cells[f.right_cells[1]].fs;
                fs1 = c.fs;
            } // end for i
        } // end for k
        break;
    case Face::kminus: // ij across face
        for (int j=0; j < blk_ptr->njc; j++) {
            for (int i=0; i < blk_ptr->nic; i++) {
                FVFace& f = blk_ptr->iFaces[blk_ptr->iFaceIndex(i, j, 0)];
                FVCell& c = blk_ptr->cells[f.right_cells[0]];
                FlowState& fs0 = blk_ptr->cells[f.left_cells[0]].fs;
                fs0 = c.fs;
                FlowState& fs1 = blk_ptr->cells[f.left_cells[1]].fs;
                fs1 = c.fs;
            } // end for i
        } // end for j
        break;
    case Face::kplus: // ij across face
        for (int j=0; j < blk_ptr->njc; j++) {
            for (int i=0; i < blk_ptr->nic; i++) {
                FVFace& f = blk_ptr->iFaces[blk_ptr->iFaceIndex(i, j, blk_ptr->nkc)];
                FVCell& c = blk_ptr->cells[f.left_cells[0]];
                FlowState& fs0 = blk_ptr->cells[f.right_cells[0]].fs;
                fs0 = c.fs;
                FlowState& fs1 = blk_ptr->cells[f.right_cells[1]].fs;
                fs1 = c.fs;
            } // end for j
        } // end for k
        break;
    }
} // end bc_outflow()

#endif
