// block.cu
// Include file for chicken.
// PJ 2022-09-11

#ifndef BLOCK_INCLUDED
#define BLOCK_INCLUDED

#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <limits>
#include "include/bxzstr/bxzstr.hpp"

#include "number.cu"
#include "rsla.cu"
#include "vector3.cu"
#include "config.cu"
#include "gas.cu"
#include "vertex.cu"
#include "flow.cu"
#include "face.cu"
#include "cell.cu"

using namespace std;


__host__ __device__
int setup_LSQ_arrays_at_face(FVFace& f, FVCell cells[], FVFace faces[])
// Prepare the inverse of the least-squares design matrix and use it to
// prepare the final weights of the points in the cloud at each face.
// The evaluation of each flow gradient becomes a matrix*vector product.
//
// Returns 0 if all successful, else 1 if a singular matrix is encountered.
//
// Adapted from the Eilmer4 code PJ, 2022-10-22.
// Placed in block.cu so that FVCell and FVFace struct definitions can be seen.
{
    // Get pointers to all of the cloud positions.
    Vector3* cloud_pos[cloud_nmax];
    for (int i=0; i < f.cloud_nc; i++) { cloud_pos[i] = &(cells[f.cells_in_cloud[i]].pos); }
    for (int i=0; i < f.cloud_nf; i++) { cloud_pos[f.cloud_nc+i] = &(faces[f.faces_in_cloud[i]].pos); }
    int cloud_n = f.cloud_nc + f.cloud_nf;
    //
    // Calculate the weights used in the least-squares gradient calculation.
    // These are the square of the weights on the original linear constraint eqns
    // and are calculated with the face centre as the reference point.
    number weights2[cloud_nmax];
    number dx[cloud_nmax], dy[cloud_nmax], dz[cloud_nmax];
    number x0 = f.pos.x; number y0 = f.pos.y; number z0 = f.pos.z;
    for (int i=0; i < cloud_n; i++) {
        dx[i] = cloud_pos[i]->x - x0;
        dy[i] = cloud_pos[i]->y - y0;
        dz[i] = cloud_pos[i]->z - z0;
        weights2[i] = 1.0/(dx[i]*dx[i] + dy[i]*dy[i] + dz[i]*dz[i]);
    }
    //
    // Set up the matrix for the normal equations.
    //
    number xx = 0.0; number xy = 0.0; number xz = 0.0;
    number yy = 0.0; number yz = 0.0; number zz = 0.0;
    for (int i=0; i < cloud_n; i++) {
        xx += weights2[i]*dx[i]*dx[i];
        xy += weights2[i]*dx[i]*dy[i];
        xz += weights2[i]*dx[i]*dz[i];
        yy += weights2[i]*dy[i]*dy[i];
        yz += weights2[i]*dy[i]*dz[i];
        zz += weights2[i]*dz[i]*dz[i];
    }
    number xTx[3][3]; // normal matrix
    xTx[0][0] = xx; xTx[0][1] = xy; xTx[0][2] = xz;
    xTx[1][0] = xy; xTx[1][1] = yy; xTx[1][2] = yz;
    xTx[2][0] = xz; xTx[2][1] = yz; xTx[2][2] = zz;
    //
    number xTxInv[3][3]; // Inverse of normal matrix.
    number very_small_value = 1.0e-16; // Should be 1.0e-32 (normInf(xTx))^^3;
    if (0 != rsla::MInverse(xTx, xTxInv, very_small_value)) {
        return 1;
    }
    // Prepare final weights for later use in the reconstruction phase.
    for (int i=0; i < cloud_n; i++) {
        f.wx[i] = weights2[i]*(xTxInv[0][0]*dx[i] + xTxInv[0][1]*dy[i] + xTxInv[0][2]*dz[i]);
        f.wy[i] = weights2[i]*(xTxInv[1][0]*dx[i] + xTxInv[1][1]*dy[i] + xTxInv[1][2]*dz[i]);
        f.wz[i] = weights2[i]*(xTxInv[2][0]*dx[i] + xTxInv[2][1]*dy[i] + xTxInv[2][2]*dz[i]);
    }
    return 0; // All weights successfully computed.
} // end setup_LSQ_arrays()


__host__ __device__
void add_viscous_fluxes_at_face(FVFace& f, FVCell cells[], FVFace faces[])
// Compute the flow quantity gradients at the face centre,
// making use of the least-squares coefficients prepared at the start of stepping.
// Then add the viscous component of the fluxes of mass, momentum and energy
// to the convective flux values that were computed eariler.
{
    // Get local copies of the cloud FlowStates.
    FlowState cloud_fs[cloud_nmax];
    for (int i=0; i < f.cloud_nc; i++) { cloud_fs[i] = cells[f.cells_in_cloud[i]].fs; }
    for (int i=0; i < f.cloud_nf; i++) { cloud_fs[f.cloud_nc+i] = faces[f.faces_in_cloud[i]].fs; }
    int cloud_n = f.cloud_nc + f.cloud_nf;
    FlowState fs0 = f.fs;
    // Now, compute the gradients, one flow quantity at a time.
    // On device, local memory will be faster than accumulating results in global memory.
    Vector3 grad_T{0.0, 0.0, 0.0};
    Vector3 grad_vx{0.0, 0.0, 0.0};
    Vector3 grad_vy{0.0, 0.0, 0.0};
    Vector3 grad_vz{0.0, 0.0, 0.0};
    number dq = 0.0;
    for (int i=0; i < cloud_n; i++) {
	number wx = f.wx[i];
	number wy = f.wy[i];
	number wz = f.wz[i];
	// temperature
        dq = cloud_fs[i].gas.T - fs0.gas.T;
        grad_T.x += wx * dq;
        grad_T.y += wy * dq;
        grad_T.z += wz * dq;
	// x velocity
        dq = cloud_fs[i].vel.x - fs0.vel.x;
        grad_vx.x += wx * dq;
        grad_vx.y += wy * dq;
        grad_vx.z += wz * dq;
	// y velocity
        dq = cloud_fs[i].vel.y - fs0.vel.y;
        grad_vy.x += wx * dq;
        grad_vy.y += wy * dq;
        grad_vy.z += wz * dq;
	// z velocity
        dq = cloud_fs[i].vel.z - fs0.vel.z;
        grad_vz.x += wx * dq;
        grad_vz.y += wy * dq;
        grad_vz.z += wz * dq;
    }
    // Calculate the viscous fluxes of mass, momentum and energy by
    // Combining the flow-quantity gradients with the transport coefficients.
    number mu, k;
    fs0.gas.trans_coeffs(mu, k);
    number lmbda = -2.0/3.0 * mu;
    // Shear stresses.
    number tau_xx = 2.0*mu*grad_vx.x + lmbda*(grad_vx.x + grad_vy.y + grad_vz.z);
    number tau_yy = 2.0*mu*grad_vy.y + lmbda*(grad_vx.x + grad_vy.y + grad_vz.z);
    number tau_zz = 2.0*mu*grad_vz.z + lmbda*(grad_vx.x + grad_vy.y + grad_vz.z);
    number tau_xy = mu * (grad_vx.y + grad_vy.x);
    number tau_xz = mu * (grad_vx.z + grad_vz.x);
    number tau_yz = mu * (grad_vy.z + grad_vz.y);
    // Thermal conduction.
    number qx = k * grad_T.x;
    number qy = k * grad_T.y;
    number qz = k * grad_T.z;
    // Combine into fluxes: store as the dot product (F.n).
    Vector3 n = f.n;
    ConservedQuantities F = f.F; // face holds convective fluxes already.
    // Mass flux -- NO CONTRIBUTION
    F[CQI::xMom] -= tau_xx*n.x + tau_xy*n.y + tau_xz*n.z;
    F[CQI::yMom] -= tau_xy*n.x + tau_yy*n.y + tau_yz*n.z;
    F[CQI::zMom] -= tau_xz*n.x + tau_yz*n.y + tau_zz*n.z;
    F[CQI::totEnergy] -=
        (tau_xx*fs0.vel.x + tau_xy*fs0.vel.y + tau_xz*fs0.vel.z + qx)*n.x +
        (tau_xy*fs0.vel.x + tau_yy*fs0.vel.y + tau_yz*fs0.vel.z + qy)*n.y +
        (tau_xz*fs0.vel.x + tau_yz*fs0.vel.y + tau_zz*fs0.vel.z + qz)*n.z;
    f.F = F;
} // end add_viscous_fluxes_at_face()

//-----------------------------------------------------------------------------------

struct Block {
    // Storage for active cells and ghost cells.
    vector<FVCell> cells;
    FVCell* cells_on_gpu;
    //
    // Active cells have conserved quantities data, along with the time derivatives.
    vector<ConservedQuantities> Q;
    vector<ConservedQuantities> dQdt;
    ConservedQuantities* Q_on_gpu;
    ConservedQuantities* dQdt_on_gpu;
    //
    // Collection of faces which bound the active cells.
    // We compute fluxes of conserved flow properties across these faces.
    vector<FVFace> faces;
    FVFace* faces_on_gpu;
    //
    // The vertices are used to define the locations and geometric properties
    // of faces and cells.
    vector<Vector3> vertices;
    Vector3* vertices_on_gpu;


    __host__
    string toString() {
        string repr = "Block()";
        return repr;
    }

    __host__
    size_t configure(const BConfig& cfg)
    // Set up the block to hold the grid and flow data.
    // Do this before reading a grid or flow file.
    {
        size_t bytes_allocated = 0;
        // Now that we know the numbers of cells, resize the data store to fit them all.
        cells.resize(cfg.nActiveCells + cfg.nTotalGhostCells);
        bytes_allocated += cells.size()*sizeof(FVCell);
        if (cfg.active) {
            Q.resize(cfg.nActiveCells*2);
            dQdt.resize(cfg.nActiveCells*3);
        }
        bytes_allocated += (Q.size()+dQdt.size())*sizeof(ConservedQuantities);
        //
        // Each set of finite-volume faces is in the index-plane of the corresponding vertices
        // but we pack them all into the one vector.
        faces.resize(cfg.nFaces);
        bytes_allocated += faces.size()*sizeof(FVFace);
        //
        // And the vertices.
        vertices.resize((cfg.nic+1)*(cfg.njc+1)*(cfg.nkc+1));
        bytes_allocated += vertices.size()*sizeof(Vector3);
        //
#ifdef CUDA
        // We need to allocate corresponding memory space on the GPU.
        auto status = cudaMalloc(&cells_on_gpu, cells.size()*sizeof(FVCell));
        if (status) throw runtime_error("Could not allocate cells on gpu.");
        status = cudaMalloc(&Q_on_gpu, Q.size()*sizeof(ConservedQuantities));
        if (status) throw runtime_error("Could not allocate Q on gpu.");
        status = cudaMalloc(&dQdt_on_gpu, dQdt.size()*sizeof(ConservedQuantities));
        if (status) throw runtime_error("Could not allocate dQdt on gpu.");
        status = cudaMalloc(&faces_on_gpu, faces.size()*sizeof(FVFace));
        if (status) throw runtime_error("Could not allocate faces on gpu.");
        status = cudaMalloc(&vertices_on_gpu, vertices.size()*sizeof(Vector3));
        if (status) throw runtime_error("Could not allocate vertices on gpu.");
#endif
        //
        // Make connections from cells to faces and vertices.
        for (int k=0; k < cfg.nkc; k++) {
            for (int j=0; j < cfg.njc; j++) {
                for (int i=0; i < cfg.nic; i++) {
                    FVCell& c = cells[cfg.activeCellIndex(i,j,k)];
                    c.face[Face::iminus] = cfg.iFaceIndex(i,j,k);
                    c.face[Face::iplus] = cfg.iFaceIndex(i+1,j,k);
                    c.face[Face::jminus] = cfg.jFaceIndex(i,j,k);
                    c.face[Face::jplus] = cfg.jFaceIndex(i,j+1,k);
                    c.face[Face::kminus] = cfg.kFaceIndex(i,j,k);
                    c.face[Face::kplus] = cfg.kFaceIndex(i,j,k+1);
                    c.vtx[0] = cfg.vtxIndex(i,j,k);
                    c.vtx[1] = cfg.vtxIndex(i+1,j,k);
                    c.vtx[2] = cfg.vtxIndex(i+1,j+1,k);
                    c.vtx[3] = cfg.vtxIndex(i,j+1,k);
                    c.vtx[4] = cfg.vtxIndex(i,j,k+1);
                    c.vtx[5] = cfg.vtxIndex(i+1,j,k+1);
                    c.vtx[6] = cfg.vtxIndex(i+1,j+1,k+1);
                    c.vtx[7] = cfg.vtxIndex(i,j+1,k+1);
                }
            }
        }
        //
        // Make connections from faces to cells and vertices.
        //
        // Face i  0     1     2     3     4
        //         +-----+-----+-----+-----+
        // Cell i  |  0  |  1  |  2  |  3  |
        //         +-----+-----+-----+-----+
        //
        // iFaces
        for (int k=0; k < cfg.nkc; k++) {
            for (int j=0; j < cfg.njc; j++) {
                for (int i=0; i < cfg.nic+1; i++) {
                    FVFace& f = faces[cfg.iFaceIndex(i,j,k)];
                    f.vtx[0] = cfg.vtxIndex(i,j,k);
                    f.vtx[1] = cfg.vtxIndex(i,j+1,k);
                    f.vtx[2] = cfg.vtxIndex(i,j+1,k+1);
                    f.vtx[3] = cfg.vtxIndex(i,j,k+1);
                    // Set neighbouring cells for convective fluxes.
                    if (i == 0) {
                        f.left_cells[1] = cfg.ghostCellIndex(Face::iminus,j,k,1);
                        f.left_cells[0] = cfg.ghostCellIndex(Face::iminus,j,k,0);
                        f.right_cells[0] = cfg.activeCellIndex(i,j,k);
                        f.right_cells[1] = cfg.activeCellIndex(i+1,j,k);
                        f.bcId = Face::iminus;
                        f.bcCode = cfg.bcCodes[Face::iminus];
                        if (f.bcCode == BCCode::inflow) f.inflowId = cfg.bc_fs[Face::iminus];
                        if (f.bcCode == BCCode::wall_no_slip_fixed_T) f.TWall = cfg.bc_TWall[Face::iminus];
                        if (f.bcCode == BCCode::inflow_function) f.bcFun = cfg.bc_fun[Face::iminus];
                    } else if (i == 1 && i == cfg.nic) {
                        cerr << "cfg.nic=" << cfg.nic << endl;
                        throw runtime_error("Too few cells in the i-index direction.");
                    } else if (i == 1 && i == cfg.nic-1) {
                        f.left_cells[1] = cfg.ghostCellIndex(Face::iminus,j,k,0);
                        f.left_cells[0] = cfg.activeCellIndex(i-1,j,k);
                        f.right_cells[0] = cfg.activeCellIndex(i,j,k);
                        f.right_cells[1] = cfg.ghostCellIndex(Face::iplus,j,k,0);
                    } else if (i == 1) {
                        f.left_cells[1] = cfg.ghostCellIndex(Face::iminus,j,k,0);
                        f.left_cells[0] = cfg.activeCellIndex(i-1,j,k);
                        f.right_cells[0] = cfg.activeCellIndex(i,j,k);
                        f.right_cells[1] = cfg.activeCellIndex(i+1,j,k);
                    } else if (i == cfg.nic-1) {
                        f.left_cells[1] = cfg.activeCellIndex(i-2,j,k);
                        f.left_cells[0] = cfg.activeCellIndex(i-1,j,k);
                        f.right_cells[0] = cfg.activeCellIndex(i,j,k);
                        f.right_cells[1] = cfg.ghostCellIndex(Face::iplus,j,k,0);
                    } else if (i == cfg.nic) {
                        f.left_cells[1] = cfg.activeCellIndex(i-2,j,k);
                        f.left_cells[0] = cfg.activeCellIndex(i-1,j,k);
                        f.right_cells[0] = cfg.ghostCellIndex(Face::iplus,j,k,0);
                        f.right_cells[1] = cfg.ghostCellIndex(Face::iplus,j,k,1);
                        f.bcId = Face::iplus;
                        f.bcCode = cfg.bcCodes[Face::iplus];
                        if (f.bcCode == BCCode::inflow) f.inflowId = cfg.bc_fs[Face::iplus];
                        if (f.bcCode == BCCode::wall_no_slip_fixed_T) f.TWall = cfg.bc_TWall[Face::iplus];
                        if (f.bcCode == BCCode::inflow_function) f.bcFun = cfg.bc_fun[Face::iplus];
                    } else {
                        // All interior cells.
                        f.left_cells[1] = cfg.activeCellIndex(i-2,j,k);
                        f.left_cells[0] = cfg.activeCellIndex(i-1,j,k);
                        f.right_cells[0] = cfg.activeCellIndex(i,j,k);
                        f.right_cells[1] = cfg.activeCellIndex(i+1,j,k);
                    }
                    // Set cloud of FlowStates for gradient calculations of viscous fluxes.
                    if (i == 0) {
                        f.bcCode = cfg.bcCodes[Face::iminus];
                        if (f.bcCode == BCCode::wall_with_slip ||
                            f.bcCode == BCCode::wall_no_slip_adiabatic ||
                            f.bcCode == BCCode::wall_no_slip_fixed_T) {
                            // Do not use ghost cell.
                            f.cells_in_cloud = {cfg.activeCellIndex(i,j,k), -1};
                            f.cloud_nc = 1;
                        } else {
                            f.cells_in_cloud = {cfg.ghostCellIndex(Face::iminus,j,k,0), cfg.activeCellIndex(i,j,k)};
                            f.cloud_nc = 2;
                        }
                        f.faces_in_cloud = {cfg.jFaceIndex(i,j,k), cfg.jFaceIndex(i,j+1,k),
                                            cfg.kFaceIndex(i,j,k), cfg.kFaceIndex(i,j,k+1),
                                            -1, -1, -1, -1};
                        f.cloud_nf = 4;
                    } else if (i == cfg.nic) {
                        f.bcCode = cfg.bcCodes[Face::iplus];
                        if (f.bcCode == BCCode::wall_with_slip ||
                            f.bcCode == BCCode::wall_no_slip_adiabatic ||
                            f.bcCode == BCCode::wall_no_slip_fixed_T) {
                            // Do not use ghost cell.
                            f.cells_in_cloud = {cfg.activeCellIndex(i-1,j,k), -1};
                            f.cloud_nc = 1;
                        } else {
                            f.cells_in_cloud = {cfg.activeCellIndex(i-1,j,k), cfg.ghostCellIndex(Face::iplus,j,k,0)};
                            f.cloud_nc = 2;
                        }
                        f.faces_in_cloud = {cfg.jFaceIndex(i-1,j,k), cfg.jFaceIndex(i-1,j+1,k),
                                            cfg.kFaceIndex(i-1,j,k), cfg.kFaceIndex(i-1,j,k+1),
                                            -1, -1, -1, -1};
                        f.cloud_nf = 4;
                    } else {
                        f.bcCode = -1; // Interior face.
                        f.cells_in_cloud = {cfg.activeCellIndex(i-1,j,k), cfg.activeCellIndex(i,j,k)};
                        f.cloud_nc = 2;
                        f.faces_in_cloud = {cfg.jFaceIndex(i-1,j,k), cfg.jFaceIndex(i-1,j+1,k),
                                            cfg.kFaceIndex(i-1,j,k), cfg.kFaceIndex(i-1,j,k+1),
                                            cfg.jFaceIndex(i,j,k), cfg.jFaceIndex(i,j+1,k),
                                            cfg.kFaceIndex(i,j,k), cfg.kFaceIndex(i,j,k+1)};
                        f.cloud_nf = 8;
                    }
                }
            }
        }
        // jFaces
        for (int k=0; k < cfg.nkc; k++) {
            for (int i=0; i < cfg.nic; i++) {
                for (int j=0; j < cfg.njc+1; j++) {
                    FVFace& f = faces[cfg.jFaceIndex(i,j,k)];
                    f.vtx[0] = cfg.vtxIndex(i+1,j,k);
                    f.vtx[1] = cfg.vtxIndex(i,j,k);
                    f.vtx[2] = cfg.vtxIndex(i,j,k+1);
                    f.vtx[3] = cfg.vtxIndex(i+1,j,k+1);
                    // Set neighbouring cells for convective fluxes.
                    if (j == 0) {
                        f.left_cells[1] = cfg.ghostCellIndex(Face::jminus,i,k,1);
                        f.left_cells[0] = cfg.ghostCellIndex(Face::jminus,i,k,0);
                        f.right_cells[0] = cfg.activeCellIndex(i,j,k);
                        f.right_cells[1] = cfg.activeCellIndex(i,j+1,k);
                        f.bcId = Face::jminus;
                        f.bcCode = cfg.bcCodes[Face::jminus];
                        if (f.bcCode == BCCode::inflow) f.inflowId = cfg.bc_fs[Face::jminus];
                        if (f.bcCode == BCCode::wall_no_slip_fixed_T) f.TWall = cfg.bc_TWall[Face::jminus];
                        if (f.bcCode == BCCode::inflow_function) f.bcFun = cfg.bc_fun[Face::jminus];
                    } else if (j == 1 && j == cfg.njc) {
                        cerr << "cfg.njc=" << cfg.njc << endl;
                        throw runtime_error("Too few cells in the j-index direction.");
                    } else if (j == 1 && j == cfg.njc-1) {
                        f.left_cells[1] = cfg.ghostCellIndex(Face::jminus,i,k,0);
                        f.left_cells[0] = cfg.activeCellIndex(i,j-1,k);
                        f.right_cells[0] = cfg.activeCellIndex(i,j,k);
                        f.right_cells[1] = cfg.ghostCellIndex(Face::jplus,i,k,0);
                    } else if (j == 1) {
                        f.left_cells[1] = cfg.ghostCellIndex(Face::jminus,i,k,0);
                        f.left_cells[0] = cfg.activeCellIndex(i,j-1,k);
                        f.right_cells[0] = cfg.activeCellIndex(i,j,k);
                        f.right_cells[1] = cfg.activeCellIndex(i,j+1,k);
                    } else if (j == cfg.njc-1) {
                        f.left_cells[1] = cfg.activeCellIndex(i,j-2,k);
                        f.left_cells[0] = cfg.activeCellIndex(i,j-1,k);
                        f.right_cells[0] = cfg.activeCellIndex(i,j,k);
                        f.right_cells[1] = cfg.ghostCellIndex(Face::jplus,i,k,0);
                    } else if (j == cfg.njc) {
                        f.left_cells[1] = cfg.activeCellIndex(i,j-2,k);
                        f.left_cells[0] = cfg.activeCellIndex(i,j-1,k);
                        f.right_cells[0] = cfg.ghostCellIndex(Face::jplus,i,k,0);
                        f.right_cells[1] = cfg.ghostCellIndex(Face::jplus,i,k,1);
                        f.bcId = Face::jplus;
                        f.bcCode = cfg.bcCodes[Face::jplus];
                        if (f.bcCode == BCCode::inflow) f.inflowId = cfg.bc_fs[Face::jplus];
                        if (f.bcCode == BCCode::wall_no_slip_fixed_T) f.TWall = cfg.bc_TWall[Face::jplus];
                        if (f.bcCode == BCCode::inflow_function) f.bcFun = cfg.bc_fun[Face::jplus];
                    } else {
                        // All interior cells.
                        f.left_cells[1] = cfg.activeCellIndex(i,j-2,k);
                        f.left_cells[0] = cfg.activeCellIndex(i,j-1,k);
                        f.right_cells[0] = cfg.activeCellIndex(i,j,k);
                        f.right_cells[1] = cfg.activeCellIndex(i,j+1,k);
                    }
                    // Set cloud of FlowStates for gradient calculations of viscous fluxes.
                    if (j == 0) {
                        f.bcCode = cfg.bcCodes[Face::jminus];
                        if (f.bcCode == BCCode::wall_with_slip ||
                            f.bcCode == BCCode::wall_no_slip_adiabatic ||
                            f.bcCode == BCCode::wall_no_slip_fixed_T) {
                            // Do not use ghost cell.
                            f.cells_in_cloud = {cfg.activeCellIndex(i,j,k), -1};
                            f.cloud_nc = 1;
                        } else {
                            f.cells_in_cloud = {cfg.ghostCellIndex(Face::jminus,i,k,0), cfg.activeCellIndex(i,j,k)};
                            f.cloud_nc = 2;
                        }
                        f.faces_in_cloud = {cfg.iFaceIndex(i,j,k), cfg.iFaceIndex(i+1,j,k),
                                            cfg.kFaceIndex(i,j,k), cfg.kFaceIndex(i,j,k+1),
                                            -1, -1, -1, -1};
                        f.cloud_nf = 4;
                    } else if (j == cfg.njc) {
                        f.bcCode = cfg.bcCodes[Face::jplus];
                        if (f.bcCode == BCCode::wall_with_slip ||
                            f.bcCode == BCCode::wall_no_slip_adiabatic ||
                            f.bcCode == BCCode::wall_no_slip_fixed_T) {
                            // Do not use ghost cell.
                            f.cells_in_cloud = {cfg.activeCellIndex(i,j-1,k), -1};
                            f.cloud_nc = 1;
                        } else {
                            f.cells_in_cloud = {cfg.activeCellIndex(i,j-1,k), cfg.ghostCellIndex(Face::jplus,i,k,0)};
                            f.cloud_nc = 2;
                        }
                        f.faces_in_cloud = {cfg.iFaceIndex(i,j-1,k), cfg.iFaceIndex(i+1,j-1,k),
                                            cfg.kFaceIndex(i,j-1,k), cfg.kFaceIndex(i,j-1,k+1),
                                            -1, -1, -1, -1};
                        f.cloud_nf = 4;
                    } else {
                        f.bcCode = -1; // Interior face.
                        f.cells_in_cloud = {cfg.activeCellIndex(i,j-1,k), cfg.activeCellIndex(i,j,k)};
                        f.cloud_nc = 2;
                        f.faces_in_cloud = {cfg.iFaceIndex(i,j-1,k), cfg.iFaceIndex(i+1,j-1,k),
                                            cfg.kFaceIndex(i,j-1,k), cfg.kFaceIndex(i,j-1,k+1),
                                            cfg.iFaceIndex(i,j,k), cfg.iFaceIndex(i+1,j,k),
                                            cfg.kFaceIndex(i,j,k), cfg.kFaceIndex(i,j,k+1)};
                        f.cloud_nf = 8;
                    }
                }
            }
        }
        // kFaces
        for (int j=0; j < cfg.njc; j++) {
            for (int i=0; i < cfg.nic; i++) {
                for (int k=0; k < cfg.nkc+1; k++) {
                    FVFace& f = faces[cfg.kFaceIndex(i,j,k)];
                    f.vtx[0] = cfg.vtxIndex(i,j,k);
                    f.vtx[1] = cfg.vtxIndex(i+1,j,k);
                    f.vtx[2] = cfg.vtxIndex(i+1,j+1,k);
                    f.vtx[3] = cfg.vtxIndex(i,j+1,k);
                    // Set neighbouring cells for convective fluxes.
                    if (k == 0) {
                        f.left_cells[1] = cfg.ghostCellIndex(Face::kminus,i,j,1);
                        f.left_cells[0] = cfg.ghostCellIndex(Face::kminus,i,j,0);
                        f.right_cells[0] = cfg.activeCellIndex(i,j,k);
                        f.right_cells[1] = cfg.activeCellIndex(i,j,k+1);
                        f.bcId = Face::kminus;
                        f.bcCode = cfg.bcCodes[Face::kminus];
                        if (f.bcCode == BCCode::inflow) f.inflowId = cfg.bc_fs[Face::kminus];
                        if (f.bcCode == BCCode::wall_no_slip_fixed_T) f.TWall = cfg.bc_TWall[Face::kminus];
                        if (f.bcCode == BCCode::inflow_function) f.bcFun = cfg.bc_fun[Face::kminus];
                    } else if (k == 1 && k == cfg.nkc) {
                        cerr << "cfg.nkc=" << cfg.nkc << endl;
                        throw runtime_error("Too few cells in the k-index direction.");
                    } else if (k == 1 && k == cfg.nkc-1) {
                        f.left_cells[1] = cfg.ghostCellIndex(Face::kminus,i,j,0);
                        f.left_cells[0] = cfg.activeCellIndex(i,j,k-1);
                        f.right_cells[0] = cfg.activeCellIndex(i,j,k);
                        f.right_cells[1] = cfg.ghostCellIndex(Face::kplus,i,j,0);
                    } else if (k == 1) {
                        f.left_cells[1] = cfg.ghostCellIndex(Face::kminus,i,j,0);
                        f.left_cells[0] = cfg.activeCellIndex(i,j,k-1);
                        f.right_cells[0] = cfg.activeCellIndex(i,j,k);
                        f.right_cells[1] = cfg.activeCellIndex(i,j,k+1);
                    } else if (k == cfg.nkc-1) {
                        f.left_cells[1] = cfg.activeCellIndex(i,j,k-2);
                        f.left_cells[0] = cfg.activeCellIndex(i,j,k-1);
                        f.right_cells[0] = cfg.activeCellIndex(i,j,k);
                        f.right_cells[1] = cfg.ghostCellIndex(Face::kplus,i,j,0);
                    } else if (k == cfg.nkc) {
                        f.left_cells[1] = cfg.activeCellIndex(i,j,k-2);
                        f.left_cells[0] = cfg.activeCellIndex(i,j,k-1);
                        f.right_cells[0] = cfg.ghostCellIndex(Face::kplus,i,j,0);
                        f.right_cells[1] = cfg.ghostCellIndex(Face::kplus,i,j,1);
                        f.bcId = Face::kplus;
                        f.bcCode = cfg.bcCodes[Face::kplus];
                        if (f.bcCode == BCCode::inflow) f.inflowId = cfg.bc_fs[Face::kplus];
                        if (f.bcCode == BCCode::wall_no_slip_fixed_T) f.TWall = cfg.bc_TWall[Face::kplus];
                        if (f.bcCode == BCCode::inflow_function) f.bcFun = cfg.bc_fun[Face::kplus];
                    } else {
                        // All interior cells.
                        f.left_cells[1] = cfg.activeCellIndex(i,j,k-2);
                        f.left_cells[0] = cfg.activeCellIndex(i,j,k-1);
                        f.right_cells[0] = cfg.activeCellIndex(i,j,k);
                        f.right_cells[1] = cfg.activeCellIndex(i,j,k+1);
                    }
                    // Set cloud of FlowStates for gradient calculations of viscous fluxes.
                    if (k == 0) {
                        f.bcCode = cfg.bcCodes[Face::kminus];
                        if (f.bcCode == BCCode::wall_with_slip ||
                            f.bcCode == BCCode::wall_no_slip_adiabatic ||
                            f.bcCode == BCCode::wall_no_slip_fixed_T) {
                            // Do not use ghost cell.
                            f.cells_in_cloud = {cfg.activeCellIndex(i,j,k), -1};
                            f.cloud_nc = 1;
                        } else {
                            f.cells_in_cloud = {cfg.ghostCellIndex(Face::kminus,i,j,0), cfg.activeCellIndex(i,j,k)};
                            f.cloud_nc = 2;
                        }
                        f.faces_in_cloud = {cfg.iFaceIndex(i,j,k), cfg.iFaceIndex(i+1,j,k),
                                            cfg.jFaceIndex(i,j,k), cfg.jFaceIndex(i,j+1,k),
                                            -1, -1, -1, -1};
                        f.cloud_nf = 4;
                    } else if (k == cfg.nkc) {
                        f.bcCode = cfg.bcCodes[Face::kplus];
                        if (f.bcCode == BCCode::wall_with_slip ||
                            f.bcCode == BCCode::wall_no_slip_adiabatic ||
                            f.bcCode == BCCode::wall_no_slip_fixed_T) {
                            // Do not use ghost cell.
                            f.cells_in_cloud = {cfg.activeCellIndex(i,j,k-1), -1};
                            f.cloud_nc = 1;
                        } else {
                            f.cells_in_cloud = {cfg.activeCellIndex(i,j,k-1), cfg.ghostCellIndex(Face::kplus,i,j,0)};
                            f.cloud_nc = 2;
                        }
                        f.faces_in_cloud = {cfg.iFaceIndex(i,j,k-1), cfg.iFaceIndex(i+1,j,k-1),
                                            cfg.jFaceIndex(i,j,k-1), cfg.jFaceIndex(i,j+1,k-1),
                                            -1, -1, -1, -1};
                        f.cloud_nf = 4;
                    } else {
                        f.bcCode = -1; // Interior face.
                        f.cells_in_cloud = {cfg.activeCellIndex(i,j,k-1), cfg.activeCellIndex(i,j,k)};
                        f.cloud_nc = 2;
                        f.faces_in_cloud = {cfg.iFaceIndex(i,j,k-1), cfg.iFaceIndex(i+1,j,k-1),
                                            cfg.jFaceIndex(i,j,k-1), cfg.jFaceIndex(i,j+1,k-1),
                                            cfg.iFaceIndex(i,j,k), cfg.iFaceIndex(i+1,j,k),
                                            cfg.jFaceIndex(i,j,k), cfg.jFaceIndex(i,j+1,k)};
                        f.cloud_nf = 8;
                    }
                }
            }
        }
        //
        // Now that the faces exist, we can overwrite the TWall values in the FixedWallT
        // boundary conditions were set with individual temperature values for each face.
        for (int bf=0; bf < 6; bf++) {
            if (cfg.bcCodes[bf] == BCCode::wall_no_slip_fixed_T && cfg.bc_TWall_form[bf] == TWallForm::fun) {
                char fileName[256];
                sprintf(fileName, "%s/TWall-%04d-%04d-%04d-%s.gz", Config::job.c_str(),
                        cfg.i, cfg.j, cfg.k, Face::names[bf].c_str());
                // Gzipped text file.
                auto fs = bxz::ifstream(fileName); // gzip file
                if (!fs) {
                    throw runtime_error("Did not open gzipped TWall file successfully: "+string(fileName));
                }
                switch (bf) {
                case Face::iminus: case Face::iplus: {
                    int i = (bf == Face::iminus) ? 0 : cfg.nic;
                    for (int k=0; k < cfg.nkc; k++) {
                        for (int j=0; j < cfg.njc; j++) {
                            fs >> faces[cfg.iFaceIndex(i,j,k)].TWall;
                        }
                    }
                    break;
                }
                case Face::jminus: case Face::jplus: {
                    int j = (bf == Face::jminus) ? 0 : cfg.njc;
                    for (int k=0; k < cfg.nkc; k++) {
                        for (int i=0; i < cfg.nic; i++) {
                            fs >> faces[cfg.jFaceIndex(i,j,k)].TWall;
                        }
                    }
                    break;
                }
                case Face::kminus: case Face::kplus: {
                    int k = (bf == Face::kminus) ? 0 : cfg.nkc;
                    for (int j=0; j < cfg.njc; j++) {
                        for (int i=0; i < cfg.nic; i++) {
                            fs >> faces[cfg.kFaceIndex(i,j,k)].TWall;
                        }
                    }
                    break;
                }
                default:
                    throw runtime_error("Oops.");
                }
                fs.close();
            }
        }
        return bytes_allocated;
    } // end configure()

    __host__
    void releaseMemory()
    {
        cells.resize(0);
        Q.resize(0);
        dQdt.resize(0);
        faces.resize(0);
        vertices.resize(0);
#ifdef CUDA
        if (cells_on_gpu) { cudaFree(&cells_on_gpu); cells_on_gpu = NULL; }
        if (Q_on_gpu) { cudaFree(&Q_on_gpu); Q_on_gpu = NULL; }
        if (dQdt_on_gpu) { cudaFree(&dQdt_on_gpu); dQdt_on_gpu = NULL; }
        if (faces_on_gpu) { cudaFree(&faces_on_gpu); faces_on_gpu = NULL; }
        if (vertices_on_gpu) { cudaFree(&vertices_on_gpu); vertices_on_gpu = NULL; }
#endif
        return;
    }

    __host__
    void computeGeometry(const BConfig& cfg)
    // Compute cell and face geometric data.
    // Do this after reading the grid and flow files because we need the vertex locations
    // and because cell positions and volumes are part of the flow data.
    // This function will overwrite them with (potentially) better values.
    {
        for (int ic=0; ic < cfg.nActiveCells; ic++) {
            FVCell& c = cells[ic];
            hex_cell_properties(vertices[c.vtx[0]], vertices[c.vtx[1]],
                                vertices[c.vtx[2]], vertices[c.vtx[3]],
                                vertices[c.vtx[4]], vertices[c.vtx[5]],
                                vertices[c.vtx[6]], vertices[c.vtx[7]],
                                false, c.pos, c.volume, c.iLength, c.jLength, c.kLength);
        }
        for (auto& f : faces) {
            quad_properties(vertices[f.vtx[0]], vertices[f.vtx[1]],
                            vertices[f.vtx[2]], vertices[f.vtx[3]],
                            f.pos, f.n, f.t1, f.t2, f.area);
        }
        //
        if (!cfg.active) return; // No ghost cells for an inactive block.
        //
        // Work around the boundaries and extrapolate cell positions and lengths
        // into the ghost cells.  We need this data for high-order reconstruction
        // for the inviscid fluxes and for computation of the flow-property gradients
        // for the viscous fluxes.
        //
        // Face::iminus
        for (int k=0; k < cfg.nkc; k++) {
            for (int j=0; j < cfg.njc; j++) {
                FVFace& f = faces[cfg.iFaceIndex(0,j,k)];
                FVCell& c0 = cells[f.right_cells[0]];
                FVCell& g0 = cells[f.left_cells[0]];
                g0.iLength = c0.iLength;
                g0.jLength = c0.jLength;
                g0.kLength = c0.kLength;
                Vector3 d = f.pos; d.sub(c0.pos);
                g0.pos = f.pos; g0.pos.add(d);
                //
                FVCell& g1 = cells[f.left_cells[1]];
                g1.iLength = c0.iLength;
                g1.jLength = c0.jLength;
                g1.kLength = c0.kLength;
                d.mul(3.0);
                g1.pos = f.pos; g1.pos.add(d);
            }
        }
        // Face::iplus
        for (int k=0; k < cfg.nkc; k++) {
            for (int j=0; j < cfg.njc; j++) {
                FVFace& f = faces[cfg.iFaceIndex(cfg.nic,j,k)];
                FVCell& c0 = cells[f.left_cells[0]];
                FVCell& g0 = cells[f.right_cells[0]];
                g0.iLength = c0.iLength;
                g0.jLength = c0.jLength;
                g0.kLength = c0.kLength;
                Vector3 d = f.pos; d.sub(c0.pos);
                g0.pos = f.pos; g0.pos.add(d);
                //
                FVCell& g1 = cells[f.right_cells[1]];
                g1.iLength = c0.iLength;
                g1.jLength = c0.jLength;
                g1.kLength = c0.kLength;
                d.mul(3.0);
                g1.pos = f.pos; g1.pos.add(d);
            }
        }
        // Face::jminus
        for (int k=0; k < cfg.nkc; k++) {
            for (int i=0; i < cfg.nic; i++) {
                FVFace& f = faces[cfg.jFaceIndex(i,0,k)];
                FVCell& c0 = cells[f.right_cells[0]];
                FVCell& g0 = cells[f.left_cells[0]];
                g0.iLength = c0.iLength;
                g0.jLength = c0.jLength;
                g0.kLength = c0.kLength;
                Vector3 d = f.pos; d.sub(c0.pos);
                g0.pos = f.pos; g0.pos.add(d);
                //
                FVCell& g1 = cells[f.left_cells[1]];
                g1.iLength = c0.iLength;
                g1.jLength = c0.jLength;
                g1.kLength = c0.kLength;
                d.mul(3.0);
                g1.pos = f.pos; g1.pos.add(d);
            }
        }
        // Face::jplus
        for (int k=0; k < cfg.nkc; k++) {
            for (int i=0; i < cfg.nic; i++) {
                FVFace& f = faces[cfg.jFaceIndex(i,cfg.njc,k)];
                FVCell& c0 = cells[f.left_cells[0]];
                FVCell& g0 = cells[f.right_cells[0]];
                g0.iLength = c0.iLength;
                g0.jLength = c0.jLength;
                g0.kLength = c0.kLength;
                Vector3 d = f.pos; d.sub(c0.pos);
                g0.pos = f.pos; g0.pos.add(d);
                //
                FVCell& g1 = cells[f.right_cells[1]];
                g1.iLength = c0.iLength;
                g1.jLength = c0.jLength;
                g1.kLength = c0.kLength;
                d.mul(3.0);
                g1.pos = f.pos; g1.pos.add(d);
            }
        }
        // Face::kminus
        for (int j=0; j < cfg.njc; j++) {
            for (int i=0; i < cfg.nic; i++) {
                FVFace& f = faces[cfg.kFaceIndex(i,j,0)];
                FVCell& c0 = cells[f.right_cells[0]];
                FVCell& g0 = cells[f.left_cells[0]];
                g0.iLength = c0.iLength;
                g0.jLength = c0.jLength;
                g0.kLength = c0.kLength;
                Vector3 d = f.pos; d.sub(c0.pos);
                g0.pos = f.pos; g0.pos.add(d);
                //
                FVCell& g1 = cells[f.left_cells[1]];
                g1.iLength = c0.iLength;
                g1.jLength = c0.jLength;
                g1.kLength = c0.kLength;
                d.mul(3.0);
                g1.pos = f.pos; g1.pos.add(d);
            }
        }
        // Face::kplus
        for (int j=0; j < cfg.njc; j++) {
            for (int i=0; i < cfg.nic; i++) {
                FVFace& f = faces[cfg.kFaceIndex(i,j,cfg.nkc)];
                FVCell& c0 = cells[f.left_cells[0]];
                FVCell& g0 = cells[f.right_cells[0]];
                g0.iLength = c0.iLength;
                g0.jLength = c0.jLength;
                g0.kLength = c0.kLength;
                Vector3 d = f.pos; d.sub(c0.pos);
                g0.pos = f.pos; g0.pos.add(d);
                //
                FVCell& g1 = cells[f.right_cells[1]];
                g1.iLength = c0.iLength;
                g1.jLength = c0.jLength;
                g1.kLength = c0.kLength;
                d.mul(3.0);
                g1.pos = f.pos; g1.pos.add(d);
            }
        }
        //
        return;
    } // end computeGeometry()

    __host__
    void readGrid(const BConfig& cfg, string fileName, bool binary_data, bool vtkHeader=false)
    // Reads the vertex locations from a compressed file, resizing storage as needed.
    // The numbers of cells are also checked.
    {
        if (vtkHeader && binary_data) {
            throw runtime_error("Do no ask to read vtk grid file with binary data.");
        }
        if (binary_data) {
            // Raw, binary data in the grid file.
            auto f = ifstream(fileName);
            if (!f) {
                throw runtime_error("Did not open binary grid file successfully: "+fileName);
            }
            number item;
            f.read(reinterpret_cast<char*>(&item), sizeof(number)); // dimensions
            f.read(reinterpret_cast<char*>(&item), sizeof(number)); // 0.0
            f.read(reinterpret_cast<char*>(&item), sizeof(number)); // 0.0
            f.read(reinterpret_cast<char*>(&item), sizeof(number)); int niv = int(item);
            f.read(reinterpret_cast<char*>(&item), sizeof(number)); int njv = int(item);
            f.read(reinterpret_cast<char*>(&item), sizeof(number)); int nkv = int(item);
            if ((cfg.nic != niv-1) || (cfg.njc != njv-1) || (cfg.nkc != nkv-1)) {
                throw runtime_error("Unexpected grid size: niv="+to_string(niv)+
                                    " njv="+to_string(njv)+ " nkv="+to_string(nkv));
            }
            if (vertices.size() != niv*njv*nkv) throw runtime_error("Incorrect size of vertices.");
            //
            // Standard order of vertices.
            for (int k=0; k < nkv; k++) {
                for (int j=0; j < njv; j++) {
                    for (int i=0; i < niv; i++) {
                        number x, y, z;
                        f.read(reinterpret_cast<char*>(&x), sizeof(number));
                        f.read(reinterpret_cast<char*>(&y), sizeof(number));
                        f.read(reinterpret_cast<char*>(&z), sizeof(number));
                        vertices[cfg.vtxIndex(i,j,k)].set(x, y, z);
                    } // for i
                } // for j
            } // for k
            f.close();
        } else {
            // Gzipped text file.
            auto f = bxz::ifstream(fileName); // gzip file
            if (!f) {
                throw runtime_error("Did not open gzipped grid file successfully: "+fileName);
            }
            constexpr int maxc = 256;
            char line[maxc];
            int niv, njv, nkv;
            if (vtkHeader) {
                f.getline(line, maxc); // expect "vtk"
                f.getline(line, maxc); // title line
                f.getline(line, maxc); // expect "ASCII"
                f.getline(line, maxc); // expect "STRUCTURED_GRID"
                f.getline(line, maxc); // DIMENSIONS line
                sscanf(line, "DIMENSIONS %d %d %d", &niv, &njv, &nkv);
            } else {
                f.getline(line, maxc); // expect "structured_grid 1.0"
                f.getline(line, maxc); // label:
                f.getline(line, maxc); // dimensions:
                f.getline(line, maxc);
                sscanf(line, "niv: %d", &niv);
                f.getline(line, maxc);
                sscanf(line, "njv: %d", &njv);
                f.getline(line, maxc);
                sscanf(line, "nkv: %d", &nkv);
            }
            if ((cfg.nic != niv-1) || (cfg.njc != njv-1) || (cfg.nkc != nkv-1)) {
                throw runtime_error("Unexpected grid size: niv="+to_string(niv)+
                                    " njv="+to_string(njv)+ " nkv="+to_string(nkv));
            }
            if (vertices.size() != niv*njv*nkv) throw runtime_error("Incorrect size of vertices.");
            //
            // Standard order of vertices.
            for (int k=0; k < nkv; k++) {
                for (int j=0; j < njv; j++) {
                    for (int i=0; i < niv; i++) {
                        number x, y, z;
                        f.getline(line, maxc);
#                       ifdef FLOAT_NUMBERS
                        sscanf(line "%f %f %f", &x, &y, &z);
#                       else
                        sscanf(line, "%lf %lf %lf", &x, &y, &z);
#                       endif
                        vertices[cfg.vtxIndex(i,j,k)].set(x, y, z);
                    } // for i
                } // for j
            } // for k
            f.close();
        } // end of gzipped text file
        return;
    } // end readGrid()

    __host__
    void readFlow(const BConfig& cfg, string fileName, bool binary_data)
    // Reads the flow data from the file, one IO-variable at a time.
    // The cell indexing is in the same as expected in a VTK file.
    {
        if (binary_data) {
            auto f = ifstream(fileName, ios::binary);
            if (!f) {
                throw runtime_error("Did not open binary flow file successfully: "+fileName);
            }
            for (int m=0; m < IOvar::n; m++) {
                string name = IOvar::names[m];
                for (int k=0; k < cfg.nkc; k++) {
                    for (int j=0; j < cfg.njc; j++) {
                        for (int i=0; i < cfg.nic; i++) {
                            FVCell& c = cells[cfg.activeCellIndex(i,j,k)];
                            number value;
                            f.read(reinterpret_cast<char*>(&value), sizeof(number));
                            c.iovar_set(m, value);
                        }
                    }
                }
            } // end for m...
            f.close();
        } else {
            // Gzipped text file.
            auto f = bxz::ifstream(fileName); // gzip file
            if (!f) {
                throw runtime_error("Did not open gzipped flow file successfully: "+fileName);
            }
            for (int m=0; m < IOvar::n; m++) {
                string name = IOvar::names[m];
                for (int k=0; k < cfg.nkc; k++) {
                    for (int j=0; j < cfg.njc; j++) {
                        for (int i=0; i < cfg.nic; i++) {
                            FVCell& c = cells[cfg.activeCellIndex(i,j,k)];
                            number value;
                            f >> value;
                            c.iovar_set(m, value);
                        }
                    }
                }
            } // end for m...
            f.close();
        }
        return;
    } // end readFlow()

    __host__
    void writeFlow(const BConfig& cfg, string fileName, bool binary_data)
    // Writes the flow data into a new file.
    // All IO-variables are written sequentially.
    // Any necessary directories are presumed to exist.
    {
        if (binary_data) {
            auto f = ofstream(fileName, ios::binary);
            if (!f) {
                throw runtime_error("Did not open binary flow file successfully for writing: "+fileName);
            }
            for (int m=0; m < IOvar::n; m++) {
                string name = IOvar::names[m];
                for (int k=0; k < cfg.nkc; k++) {
                    for (int j=0; j < cfg.njc; j++) {
                        for (int i=0; i < cfg.nic; i++) {
                            FVCell& c = cells[cfg.activeCellIndex(i,j,k)];
                            number item = c.iovar_get(m);
                            f.write(reinterpret_cast<char*>(&item), sizeof(number));
                        }
                    }
                }
            } // end for m...
            f.close();
        } else {
            // Gzipped text file.
            auto f = bxz::ofstream(fileName); // gzip file
            if (!f) {
                throw runtime_error("Did not open gzipped flow file successfully for writing: "+fileName);
            }
            for (int m=0; m < IOvar::n; m++) {
                string name = IOvar::names[m];
                for (int k=0; k < cfg.nkc; k++) {
                    for (int j=0; j < cfg.njc; j++) {
                        for (int i=0; i < cfg.nic; i++) {
                            FVCell& c = cells[cfg.activeCellIndex(i,j,k)];
                            f << c.iovar_get(m) << endl;
                        }
                    }
                }
            } // end for m...
            f.close();
        }
        return;
    } // end writeFlow()

    __host__
    number estimate_allowed_dt(const BConfig& cfg, number cfl)
    {
        number smallest_dt = numeric_limits<number>::max();
        for (int i=0; i < cfg.nActiveCells; i++) {
            FVCell& c = cells[i];
            Vector3 inorm = faces[c.face[Face::iminus]].n;
            Vector3 jnorm = faces[c.face[Face::jminus]].n;
            Vector3 knorm = faces[c.face[Face::kminus]].n;
            smallest_dt = fmin(smallest_dt, c.estimate_local_dt(inorm, jnorm, knorm, cfl));
        }
        return smallest_dt;
    } // end estimate_allowed_dt()

    __host__
    void encodeConserved(const BConfig& cfg, int level)
    {
        for (int i=0; i < cfg.nActiveCells; i++) {
            FlowState& fs = cells[i].fs;
            ConservedQuantities& U = Q[level*cfg.nActiveCells + i];
            fs.encode_conserved(U);
        }
    }

    __host__
    int decodeConserved(const BConfig& cfg, int level)
    {
        int bad_cell_count = 0;
        for (int i=0; i < cfg.nActiveCells; i++) {
            FVCell& c = cells[i];
            ConservedQuantities U = Q[level*cfg.nActiveCells + i];
            int bad_cell_flag = c.fs.decode_conserved(U);
            bad_cell_count += bad_cell_flag;
            if (bad_cell_flag) {
                cerr << "DEBUG-A Bad cell at pos=" << c.pos << endl;
            }
        }
        return bad_cell_count;
    }

    __host__
    void calculate_convective_fluxes(int flux_calc, int x_order)
    {
        for (auto& face : faces) {
            FlowState& fsL1 = cells[face.left_cells[1]].fs;
            FlowState& fsL0 = cells[face.left_cells[0]].fs;
            FlowState& fsR0 = cells[face.right_cells[0]].fs;
            FlowState& fsR1 = cells[face.right_cells[1]].fs;
            face.calculate_convective_flux(fsL1, fsL0, fsR0, fsR1, flux_calc, x_order);
        }
    } // end calculate_convective_fluxes()

    __host__
    void setup_LSQ_arrays()
    {
        int failures = 0;
        for (auto& face : faces) {
            int flag = setup_LSQ_arrays_at_face(face, cells.data(), faces.data());
            if (flag) {
                cerr << "Singular normal matrix at f.pos=" << face.pos << endl;
            }
        }
        if (failures > 0) {
            throw runtime_error("Singular matrices encountered while setting up LSQ weights.");
        }
    }

    __host__
    void add_viscous_fluxes()
    {
        for (auto& face : faces) {
            add_viscous_fluxes_at_face(face, cells.data(), faces.data());
        }
    }

    __host__
    void WRITE_DEBUG_DATA(string label, int cIndx)
    {
        cout << "DEBUG " << label << endl;
        FVCell c = cells[cIndx];
        cout << "DEBUG cells[" << cIndx << "]=" << c << endl;
        cout << "DEBUG faces=["; for (auto i : c.face) cout << faces[i] << ","; cout << "]" << endl;
        cout << "DEBUG FlowStates given to kminus face flux calculation:" << endl;
        FVFace& fkm = faces[c.face[Face::kminus]];
        FlowState& fsL1 = cells[fkm.left_cells[1]].fs; Vector3 posL1 = cells[fkm.left_cells[1]].pos;
        FlowState& fsL0 = cells[fkm.left_cells[0]].fs; Vector3 posL0 = cells[fkm.left_cells[0]].pos;
        FlowState& fsR0 = cells[fkm.right_cells[0]].fs; Vector3 posR0 = cells[fkm.right_cells[0]].pos;
        FlowState& fsR1 = cells[fkm.right_cells[1]].fs; Vector3 posR1 = cells[fkm.right_cells[1]].pos;
        cout << "DEBUG posL1=" << posL1 << " fsL1=" << fsL1 << endl;
        cout << "DEBUG posL0=" << posL0 << " fsL0=" << fsL0 << endl;
        cout << "DEBUG posR0=" << posR0 << " fsR0=" << fsR0 << endl;
        cout << "DEBUG posR1=" << posR1 << " fsR1=" << fsR1 << endl;
    }

    __host__
    void WRITE_DEBUG_DATA2(string label, int fIndx)
    {
        cout << "DEBUG " << label << endl;
        FVFace f = faces[fIndx];
        cout << "DEBUG faces[" << fIndx << "]=" << f << endl;
        cout << "DEBUG cells in cloud" << endl;
        for (int i=0; i < f.cloud_nc; i++) {
            cout << "DEBUG i=" << i << " pos=" << cells[f.cells_in_cloud[i]].pos
                 << " fs=" << cells[f.cells_in_cloud[i]].fs << endl;
        }
        cout << "DEBUG faces in cloud" << endl;
        for (int i=0; i < f.cloud_nf; i++) {
            cout << "DEBUG i=" << i << " pos=" << faces[f.faces_in_cloud[i]].pos
                 << " fs=" << faces[f.faces_in_cloud[i]].fs << endl;
        }
    }

    __host__
    int update_stage_1(const BConfig& cfg, number dt)
    // Stage 1 of the TVD-RK3 update scheme (predictor step).
    {
        int bad_cell_count = 0;
        for (int i=0; i < cfg.nActiveCells; i++) {
            FVCell& c = cells[i];
            ConservedQuantities& dUdt0 = dQdt[i];
            c.eval_dUdt(dUdt0, faces.data(), Config::source_terms);
            ConservedQuantities& U0 = Q[i];
            ConservedQuantities& U1 = Q[cfg.nActiveCells + i];
            for (int j=0; j < CQI::n; j++) {
                U1[j] = U0[j] + dt*dUdt0[j];
            }
            int bad_cell_flag = c.fs.decode_conserved(U1);
            bad_cell_count += bad_cell_flag;
            if (bad_cell_flag) {
                cerr << "Stage 1 update, Bad cell at pos=" << c.pos << endl;
            }
        }
        return bad_cell_count;
    } // end update_stage_1()

    __host__
    int update_stage_2(const BConfig& cfg, number dt)
    // Stage 2 of the TVD-RK3 update scheme.
    {
        int bad_cell_count = 0;
        for (int i=0; i < cfg.nActiveCells; i++) {
            FVCell& c = cells[i];
            ConservedQuantities& dUdt0 = dQdt[i];
            ConservedQuantities& dUdt1 = dQdt[cfg.nActiveCells + i];
            c.eval_dUdt(dUdt1, faces.data(), Config::source_terms);
            ConservedQuantities& U0 = Q[i];
            ConservedQuantities& U1 = Q[cfg.nActiveCells + i];
            for (int j=0; j < CQI::n; j++) {
                U1[j] = U0[j] + 0.25*dt*(dUdt0[j] + dUdt1[j]);
            }
            int bad_cell_flag = c.fs.decode_conserved(U1);
            bad_cell_count += bad_cell_flag;
            if (bad_cell_flag) {
                cerr << "Stage 2 update, Bad cell at pos=" << c.pos << endl;
            }
        }
        return bad_cell_count;
    } // end update_stage_2()

    __host__
    int update_stage_3(const BConfig& cfg, number dt)
    // Stage 3 of the TVD_RK3 update scheme.
    {
        int bad_cell_count = 0;
        for (int i=0; i < cfg.nActiveCells; i++) {
            FVCell& c = cells[i];
            ConservedQuantities& dUdt0 = dQdt[i];
            ConservedQuantities& dUdt1 = dQdt[cfg.nActiveCells + i];
            ConservedQuantities& dUdt2 = dQdt[2*cfg.nActiveCells + i];
            c.eval_dUdt(dUdt2, faces.data(), Config::source_terms);
            ConservedQuantities& U0 = Q[i];
            ConservedQuantities& U1 = Q[cfg.nActiveCells + i];
            for (int j=0; j < CQI::n; j++) {
                U1[j] = U0[j] + dt*(1.0/6.0*dUdt0[j] + 1.0/6.0*dUdt1[j] + 4.0/6.0*dUdt2[j]);
            }
            int bad_cell_flag = c.fs.decode_conserved(U1);
            bad_cell_count += bad_cell_flag;
            if (bad_cell_flag) {
                cerr << "Stage 3 update, Bad cell at pos=" << c.pos << endl;
            }
        }
        return bad_cell_count;
    } // end update_stage_3()

    __host__
    void copy_conserved_data(const BConfig& cfg, int from_level, int to_level)
    {
        for (auto i=0; i < cfg.nActiveCells; i++) {
            ConservedQuantities& U_from = Q[from_level*cfg.nActiveCells + i];
            ConservedQuantities& U_to = Q[to_level*cfg.nActiveCells + i];
            for (int j=0; j < CQI::n; j++) {
                U_to[j] = U_from[j];
            }
        }
    } // end copy_conserved_data()

    __host__
    void update_chemistry(const BConfig& cfg, number dt)
    {
        for (auto i=0; i < cfg.nActiveCells; i++) {
            FlowState& fs = cells[i].fs;
            GasState& gs = fs.gas;
            gs.update_chemistry(dt);
            ConservedQuantities& U = Q[i]; // presume level 0
            fs.encode_conserved(U);
        }
    } // end update_chemistry()

}; // end Block

#endif
