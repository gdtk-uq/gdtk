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
#include <zip.h>

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
int setup_LSQ_arrays_at_face(FVFace f, FVCell cells[], FVFace faces[])
// Prepare the inverse of the least-squares design matrix and use it to
// prepare the final weights of the points in the cloud at each face.
// The evaluation of each flow gradient becomes a matrix*vector product.
//
// Returns 0 if all successful, else 1 if a singular matrix is encountered.
//
// Adapted from the Eilmer4 code PJ, 2022-10-22.
// Placed in block.cu so that FVCell and FVFace struct definitions can be seen.
{
    // Get pointers to all of the cloud FlowStates and positions.
    FlowState* flows[cloud_nmax];
    Vector3* cloud_pos[cloud_nmax];
    for (int i=0; i < f.cloud_nc; i++) { cloud_pos[i] = &(cells[i].pos); }
    for (int i=0; i < f.cloud_nf; i++) { cloud_pos[f.cloud_nc+i] = &(faces[i].pos); }
    int cloud_n = f.cloud_nc + f.cloud_nf;
    //
    // Calculate the weights used in the least-squares gradient calculation.
    // These are the square of the weights on the original linear constraint eqns
    // and are calculated with the face centre as the reference point.
    number weights2[cloud_nmax];
    number x0 = f.pos.x; number y0 = f.pos.y; number z0 = f.pos.z;
    for (int i=0; i < cloud_n; i++) {
        number dx = cloud_pos[i]->x - x0;
        number dy = cloud_pos[i]->y - y0;
        number dz = cloud_pos[i]->z - z0;
        weights2[i] = 1.0/(dx*dx+dy*dy+dz*dz);
    }
    //
    number dx[cloud_nmax], dy[cloud_nmax], dz[cloud_nmax];
    number xx = 0.0; number xy = 0.0; number xz = 0.0;
    number yy = 0.0; number yz = 0.0; number zz = 0.0;
    for (int i=0; i < cloud_n; i++) {
        dx[i] = cloud_pos[i]->x - x0;
        dy[i] = cloud_pos[i]->y - y0;
        dz[i] = cloud_pos[i]->z - z0;
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
    if (0 != MInverse(xTx, xTxInv, very_small_value)) {
        return 1;
    }
    // Prepare final weights for later use in the reconstruction phase.
    for (int i=0; i < cloud_n; i++) {
        f.wx[i] = xTxInv[0][0]*dx[i] + xTxInv[0][1]*dy[i] + xTxInv[0][2]*dz[i];
        f.wx[i] *= weights2[i];
        f.wy[i] = xTxInv[1][0]*dx[i] + xTxInv[1][1]*dy[i] + xTxInv[1][2]*dz[i];
        f.wy[i] *= weights2[i];
        f.wz[i] = xTxInv[2][0]*dx[i] + xTxInv[2][1]*dy[i] + xTxInv[2][2]*dz[i];
        f.wz[i] *= weights2[i];
    }
    return 0; // All weights successfully computed.
} // end setup_LSQ_arrays()



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
                        if (f.bcCode == BCCode::wall_with_slip || f.bcCode == BCCode::wall_with_slip) {
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
                        if (f.bcCode == BCCode::wall_with_slip || f.bcCode == BCCode::wall_with_slip) {
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
                        if (f.bcCode == BCCode::wall_with_slip || f.bcCode == BCCode::wall_with_slip) {
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
                        if (f.bcCode == BCCode::wall_with_slip || f.bcCode == BCCode::wall_with_slip) {
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
                        if (f.bcCode == BCCode::wall_with_slip || f.bcCode == BCCode::wall_with_slip) {
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
                        if (f.bcCode == BCCode::wall_with_slip || f.bcCode == BCCode::wall_with_slip) {
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
    void readGrid(const BConfig& cfg, string fileName, bool vtkHeader=false)
    // Reads the vertex locations from a compressed file, resizing storage as needed.
    // The numbers of cells are also checked.
    {
        auto f = bxz::ifstream(fileName); // gzip file
        if (!f) {
            throw runtime_error("Did not open grid file successfully: "+fileName);
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
                    f.getline(line, maxc);
                    number x, y, z;
                    #ifdef FLOAT_NUMBERS
                    sscanf(line "%f %f %f", &x, &y, &z);
                    #else
                    sscanf(line, "%lf %lf %lf", &x, &y, &z);
                    #endif
                    vertices[cfg.vtxIndex(i,j,k)].set(x, y, z);
                } // for i
            } // for j
        } // for k
        f.close();
        return;
    } // end readGrid()

    __host__
    void readFlow(const BConfig& cfg, string fileName)
    // Reads the flow data archive from a ZIP file.
    // The correct data storage is presumed to exist.
    //
    // Code modelled on the simple example by Dodrigo Rivas Costa found at
    // https://stackoverflow.com/questions/10440113/simple-way-to-unzip-a-zip-file-using-zlib
    {
        int err = 0;
        zip *z = zip_open(fileName.c_str(), ZIP_RDONLY, &err);
        if (err) {
            cerr << "Failed to open zip archive for reading: " << fileName << endl;
        }
        if (z) {
            struct zip_stat st;
            for (int m=0; m < IOvar::n; m++) {
                string name = IOvar::names[m];
                // Search archive for a variable's data.
                zip_stat_init(&st);
                zip_stat(z, name.c_str(), 0, &st);
                // Allocate enough memory for the uncompressed content and read it.
                char* content = new char[st.size];
                zip_file* f = zip_fopen(z, name.c_str(), 0);
                if (f) {
                    zip_fread(f, content, st.size);
                    zip_fclose(f);
                    stringstream ss(content);
                    string item;
                    for (int k=0; k < cfg.nkc; k++) {
                        for (int j=0; j < cfg.njc; j++) {
                            for (int i=0; i < cfg.nic; i++) {
                                getline(ss, item, '\n');
                                FVCell& c = cells[cfg.activeCellIndex(i,j,k)];
                                c.iovar_set(m, stod(item));
                            }
                        }
                    }
                } else {
                    cerr << "Could not open file " << name << " in ZIP archive " << fileName << endl;
                }
                delete[] content;
            }
            zip_close(z);
        }
        return;
    } // end readFlow()

    __host__
    void writeFlow(const BConfig& cfg, string fileName)
    // Writes the flow data into a new ZIP archive file.
    // Any necessary directories are presumed to exist.
    {
        vector<string> data; // A place to retain the string data while the zip file is constructed.
        int err = 0;
        zip *z = zip_open(fileName.c_str(), ZIP_CREATE, &err);
        if (err) {
            cerr << "Failed to open zip archive for writing: " << fileName << endl;
        }
        if (z) {
            for (int m=0; m < IOvar::n; m++) {
                string name = IOvar::names[m];
                ostringstream ss;
                for (int k=0; k < cfg.nkc; k++) {
                    for (int j=0; j < cfg.njc; j++) {
                        for (int i=0; i < cfg.nic; i++) {
                            FVCell& c = cells[cfg.activeCellIndex(i,j,k)];
                            ss << c.iovar_get(m) << endl;
                        }
                    }
                }
                data.push_back(ss.str());
                int last = data.size()-1;
                // Add the data to the ZIP archive as a file.
                zip_source_t* zs = zip_source_buffer(z, data[last].c_str(), data[last].size(), 0);
                if (zs) {
                    int zindx = zip_file_add(z, name.c_str(), zs, ZIP_FL_OVERWRITE|ZIP_FL_ENC_UTF_8);
                    if (zindx < 0) {
                        cerr << "Could not add file " << name << " to ZIP archive " << fileName << endl;
                        zip_source_free(zs);
                    }
                } else {
                    cerr << "Error getting source to add file to zip: " << string(zip_strerror(z)) << endl;
                }
            }
            zip_close(z);
        }
        data.resize(0);
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
                cerr << "DEBUG-A Bad cell at pos=" << c.pos.toString() << endl;
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
                cerr << "Singular normal matrix at f.pos=" << face.pos.toString() << endl;
            }
        }
        if (failures > 0) {
            throw runtime_error("Singular matrices encountered while setting up LSQ weights.");
        }
    }

    __host__
    void apply_viscous_boundary_conditions()
    {
        for (auto& face : faces) {
            face.apply_viscous_boundary_condition();
        }
    }

    __host__
    void add_viscous_flux()
    {
        for (auto& face : faces) {
            face.add_viscous_flux();
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
            c.eval_dUdt(dUdt0, faces.data());
            ConservedQuantities& U0 = Q[i];
            ConservedQuantities& U1 = Q[cfg.nActiveCells + i];
            for (int j=0; j < CQI::n; j++) {
                U1[j] = U0[j] + dt*dUdt0[j];
            }
            int bad_cell_flag = c.fs.decode_conserved(U1);
            bad_cell_count += bad_cell_flag;
            if (bad_cell_flag) {
                cerr << "Stage 1 update, Bad cell at pos=" << c.pos.toString() << endl;
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
            c.eval_dUdt(dUdt1, faces.data());
            ConservedQuantities& U0 = Q[i];
            ConservedQuantities& U1 = Q[cfg.nActiveCells + i];
            for (int j=0; j < CQI::n; j++) {
                U1[j] = U0[j] + 0.25*dt*(dUdt0[j] + dUdt1[j]);
            }
            int bad_cell_flag = c.fs.decode_conserved(U1);
            bad_cell_count += bad_cell_flag;
            if (bad_cell_flag) {
                cerr << "Stage 2 update, Bad cell at pos=" << c.pos.toString() << endl;
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
            c.eval_dUdt(dUdt2, faces.data());
            ConservedQuantities& U0 = Q[i];
            ConservedQuantities& U1 = Q[cfg.nActiveCells + i];
            for (int j=0; j < CQI::n; j++) {
                U1[j] = U0[j] + dt*(1.0/6.0*dUdt0[j] + 1.0/6.0*dUdt1[j] + 4.0/6.0*dUdt2[j]);
            }
            int bad_cell_flag = c.fs.decode_conserved(U1);
            bad_cell_count += bad_cell_flag;
            if (bad_cell_flag) {
                cerr << "Stage 3 update, Bad cell at pos=" << c.pos.toString() << endl;
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

}; // end Block



// GPU global functions cannot be member functions of FluidBlock
// so we need to pass the FluidBlock reference into them and that
// Block struct also needs to be in the global memory of the GPU.

__global__
void estimate_allowed_dt_on_gpu(Block& blk, const BConfig& cfg, number cfl, long long int* smallest_dt_picos)
{
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i < cfg.nActiveCells) {
        FVCell& c = blk.cells_on_gpu[i];
        Vector3 inorm = blk.faces_on_gpu[c.face[Face::iminus]].n;
        Vector3 jnorm = blk.faces_on_gpu[c.face[Face::jminus]].n;
        Vector3 knorm = blk.faces_on_gpu[c.face[Face::kminus]].n;
        long long int dt_picos = trunc(c.estimate_local_dt(inorm, jnorm, knorm, cfl)*1.0e12);
        atomicMin(smallest_dt_picos, dt_picos);
    }
}

__global__
void encodeConserved_on_gpu(Block& blk, const BConfig& cfg, int level)
{
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i < cfg.nActiveCells) {
        FlowState& fs = blk.cells_on_gpu[i].fs;
        ConservedQuantities& U = blk.Q_on_gpu[level*cfg.nActiveCells + i];
        fs.encode_conserved(U);
    }
}

__global__
void copy_conserved_data_on_gpu(Block& blk, const BConfig& cfg, int from_level, int to_level)
{
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i < cfg.nActiveCells) {
        ConservedQuantities& U_from = blk.Q_on_gpu[from_level*cfg.nActiveCells + i];
        ConservedQuantities& U_to = blk.Q_on_gpu[to_level*cfg.nActiveCells + i];
        for (int j=0; j < CQI::n; j++) {
            U_to[j] = U_from[j];
        }
    }
}

__global__
void calculate_convective_fluxes_on_gpu(Block& blk, const BConfig& cfg, int flux_calc, int x_order)
{
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i < cfg.nFaces) {
        FVFace& face = blk.faces_on_gpu[i];
        FlowState& fsL1 = blk.cells_on_gpu[face.left_cells[1]].fs;
        FlowState& fsL0 = blk.cells_on_gpu[face.left_cells[0]].fs;
        FlowState& fsR0 = blk.cells_on_gpu[face.right_cells[0]].fs;
        FlowState& fsR1 = blk.cells_on_gpu[face.right_cells[1]].fs;
        face.calculate_convective_flux(fsL1, fsL0, fsR0, fsR1, flux_calc, x_order);
    }
}

__global__
void setup_LSQ_arrays_on_gpu(Block& blk, const BConfig& cfg, int* failures)
{
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i < cfg.nFaces) {
        FVFace& face = blk.faces_on_gpu[i];
        int flag = setup_LSQ_arrays_at_face(face, blk.cells_on_gpu, blk.faces_on_gpu);
        atomicAdd(failures, flag);
    }
}

__global__
void apply_viscous_boundary_conditions_on_gpu(Block& blk, const BConfig& cfg)
{
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i < cfg.nFaces) {
        FVFace& face = blk.faces_on_gpu[i];
        face.apply_viscous_boundary_condition();
    }
}

__global__
void add_viscous_flux_on_gpu(Block& blk, const BConfig& cfg)
{
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i < cfg.nFaces) {
        FVFace& face = blk.faces_on_gpu[i];
        face.add_viscous_flux();
    }
}

__global__
void update_stage_1_on_gpu(Block& blk, const BConfig& cfg, number dt, int* bad_cell_count)
{
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i < cfg.nActiveCells) {
        FVCell& c = blk.cells_on_gpu[i];
        ConservedQuantities& dUdt0 = blk.dQdt_on_gpu[i];
        c.eval_dUdt(dUdt0, blk.faces_on_gpu);
        ConservedQuantities& U0 = blk.Q_on_gpu[i];
        ConservedQuantities& U1 = blk.Q_on_gpu[cfg.nActiveCells + i];
        for (int j=0; j < CQI::n; j++) {
            U1[j] = U0[j] + dt*dUdt0[j];
        }
        int bad_cell_flag = c.fs.decode_conserved(U1);
        atomicAdd(bad_cell_count, bad_cell_flag);
        if (bad_cell_flag) {
            printf("Stage 1 update, Bad cell at pos x=%g y=%g z=%g\n", c.pos.x, c.pos.y, c.pos.z);
        }
    }
}

__global__
void update_stage_2_on_gpu(Block& blk, const BConfig& cfg, number dt, int* bad_cell_count)
{
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i < cfg.nActiveCells) {
        FVCell& c = blk.cells_on_gpu[i];
        ConservedQuantities& dUdt0 = blk.dQdt_on_gpu[i];
        ConservedQuantities& dUdt1 = blk.dQdt_on_gpu[cfg.nActiveCells + i];
        c.eval_dUdt(dUdt1, blk.faces_on_gpu);
        ConservedQuantities& U0 = blk.Q_on_gpu[i];
        ConservedQuantities& U1 = blk.Q_on_gpu[cfg.nActiveCells + i];
        for (int j=0; j < CQI::n; j++) {
            U1[j] = U0[j] + 0.25*dt*(dUdt0[j] + dUdt1[j]);
        }
        int bad_cell_flag = c.fs.decode_conserved(U1);
        atomicAdd(bad_cell_count, bad_cell_flag);
        if (bad_cell_flag) {
            printf("Stage 2 update, Bad cell at pos x=%g y=%g z=%g\n", c.pos.x, c.pos.y, c.pos.z);
        }
    }
}

__global__
void update_stage_3_on_gpu(Block& blk, const BConfig& cfg, number dt, int* bad_cell_count)
{
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i < cfg.nActiveCells) {
        FVCell& c = blk.cells_on_gpu[i];
        ConservedQuantities& dUdt0 = blk.dQdt_on_gpu[i];
        ConservedQuantities& dUdt1 = blk.dQdt_on_gpu[cfg.nActiveCells + i];
        ConservedQuantities& dUdt2 = blk.dQdt_on_gpu[2*cfg.nActiveCells + i];
        c.eval_dUdt(dUdt2, blk.faces_on_gpu);
        ConservedQuantities& U0 = blk.Q_on_gpu[i];
        ConservedQuantities& U1 = blk.Q_on_gpu[cfg.nActiveCells + i];
        for (int j=0; j < CQI::n; j++) {
            U1[j] = U0[j] + dt*(1.0/6.0*dUdt0[j] + 1.0/6.0*dUdt1[j] + 4.0/6.0*dUdt2[j]);
        }
        int bad_cell_flag = c.fs.decode_conserved(U1);
        atomicAdd(bad_cell_count, bad_cell_flag);
        if (bad_cell_flag) {
            printf("Stage 3 update, Bad cell at pos x=%g y=%g z=%g\n", c.pos.x, c.pos.y, c.pos.z);
        }
    }
}

#endif
