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
#include "vector3.cu"
#include "config.cu"
#include "gas.cu"
#include "vertex.cu"
#include "flow.cu"
#include "face.cu"
#include "cell.cu"
#include "flux.cu"

using namespace std;

struct Block {
    // Active cells are the "real" cells in the simulation.
    // We compute the evolution of the gas-dynamic flow properties within them.
    int nic; // Number of cells i-direction.
    int njc; // Number of cells j-direction.
    int nkc; // Number of cells k-direction.
    int nActiveCells; // Number of active cells (with conserved quantities) in the block.
    vector<FVCell> cells;
    //
    // Active cells have conserved quantities data, along with the time derivatives.
    vector<ConservedQuantities> Q;
    vector<ConservedQuantities> dQdt;
    //
    // Ghost cells are associated with each block boundary face and
    // will be stored at the end of the active cells collection.
    // The flux calculation functions will dip into this collection for
    // active cells and ghost cells without knowing the difference.
    // Also, for each boundary face, we store some config information
    // to make the boundary-condition code a bit more compact.
    array<int,6> n0c; // number of cells in first index direction for each face.
    array<int,6> n1c; // Number of cells in second index direction for each face.
    array<int,6> nGhostCells; // Number of ghost cells on each face.
    array<int,6> firstGhostCells; // Index of the first ghost cell for each face.
    //
    // Collections of faces which bound the active cells.
    // We compute fluxes of conserved flow properties across these faces.
    vector<FVFace> iFaces;
    vector<FVFace> jFaces;
    vector<FVFace> kFaces;
    //
    // The vertices are used to define the locations and geometric properties
    // of faces and cells.
    vector<Vector3> vertices;

    __host__
    string toString() {
        string repr = "Block(nic=" + to_string(nic) +
            ", njc=" + to_string(njc) + ", nkc=" + to_string(nkc) + ")";
        return repr;
    }

    // Methods to index the elements making up the block.

    __host__ __device__
    int activeCellIndex(int i, int j, int k)
    {
        return k*nic*njc + j*nic + i;
    }

    __host__ __device__
    int ghostCellIndex(int faceIndx, int i0, int i1, int depth)
    {
        int cellIndxOnFace = i1*n0c[faceIndx] + i0;
        int nCellsOnFace = n0c[faceIndx]*n1c[faceIndx];
        return firstGhostCells[faceIndx] + nCellsOnFace*depth + cellIndxOnFace;
    }

    __host__ __device__
    int iFaceIndex(int i, int j, int k)
    {
        return i*njc*nkc + k*njc + j;
    }

    __host__ __device__
    int jFaceIndex(int i, int j, int k)
    {
        return j*nic*nkc + k*nic + i;
    }

    __host__ __device__
    int kFaceIndex(int i, int j, int k)
    {
        return k*nic*njc + j*nic + i;
    }

    __host__ __device__
    int vtxIndex(int i, int j, int k)
    {
        return k*(nic+1)*(njc+1) + j*(nic+1) + i;
    }

    __host__
    void configure(int i, int j, int k)
    // Set up the block to hold the grid and flow data.
    // Do this before reading a grid or flow file.
    {
        nic = i;
        njc = j;
        nkc = k;
        int nActiveCells = nic*njc*nkc;
        //
        // For the moment assume that all boundary conditions require ghost cells.
        n0c[Face::iminus] = njc; n1c[Face::iminus] = nkc;
        n0c[Face::iplus] = njc; n1c[Face::iplus] = nkc;
        n0c[Face::jminus] = nic; n1c[Face::jminus] = nkc;
        n0c[Face::jplus] = nic; n1c[Face::jplus] = nkc;
        n0c[Face::kminus] = nic; n1c[Face::kminus] = njc;
        n0c[Face::kplus] = nic; n1c[Face::kplus] = njc;
        for (int ib=0; ib < 6; ib++) {
            nGhostCells[ib] = 2*n0c[ib]*n1c[ib];
            if (ib > 1) {
                firstGhostCells[ib] = firstGhostCells[ib-1] + nGhostCells[ib-1];
            } else {
                firstGhostCells[ib] = nActiveCells;
            }
        }
        //
        // Now that we know the numbers of cells, resize the data store to fit them all.
        cells.resize(firstGhostCells[5]+nGhostCells[5]);
        Q.resize(nActiveCells*TLevels);
        dQdt.resize(nActiveCells*TLevels);
        #ifdef CUDA
        // We need to allocate corresponding memory space on the GPU.
        // [TODO]
        #endif
        //
        // Each set of finite-volume faces is in the index-plane of the corresponding vertices.
        iFaces.resize((nic+1)*njc*nkc);
        jFaces.resize(nic*(njc+1)*nkc);
        kFaces.resize(nic*njc*(nkc+1));
        //
        // And the vertices.
        vertices.resize((nic+1)*(njc+1)*(nkc+1));
        //
        // Make connections from cells to faces and vertices.
        for (int k=0; k < nkc; k++) {
            for (int j=0; j < njc; j++) {
                for (int i=0; i < nic; i++) {
                    FVCell& c = cells[activeCellIndex(i,j,k)];
                    c.face[Face::iminus] = iFaceIndex(i,j,k);
                    c.face[Face::iplus] = iFaceIndex(i+1,j,k);
                    c.face[Face::jminus] = jFaceIndex(i,j,k);
                    c.face[Face::jplus] = jFaceIndex(i,j+1,k);
                    c.face[Face::kminus] = kFaceIndex(i,j,k);
                    c.face[Face::kplus] = kFaceIndex(i,j,k+1);
                    c.vtx[0] = vtxIndex(i,j,k);
                    c.vtx[1] = vtxIndex(i+1,j,k);
                    c.vtx[2] = vtxIndex(i+1,j+1,k);
                    c.vtx[3] = vtxIndex(i,j+1,k);
                    c.vtx[4] = vtxIndex(i,j,k+1);
                    c.vtx[5] = vtxIndex(i+1,j,k+1);
                    c.vtx[6] = vtxIndex(i+1,j+1,k+1);
                    c.vtx[7] = vtxIndex(i,j+1,k+1);
                }
            }
        }
        //
        // Make connections from faces to cells and vertices.
        // iFaces
        for (int k=0; k < nkc; k++) {
            for (int j=0; j < njc; j++) {
                for (int i=0; i < nic+1; i++) {
                    FVFace& f = iFaces[iFaceIndex(i,j,k)];
                    f.vtx[0] = vtxIndex(i,j,k);
                    f.vtx[1] = vtxIndex(i,j+1,k);
                    f.vtx[2] = vtxIndex(i,j+1,k+1);
                    f.vtx[3] = vtxIndex(i,j,k+1);
                    if (i == 0) {
                        f.left_cells[0] = ghostCellIndex(Face::iminus,j,k,1);
                        f.left_cells[1] = ghostCellIndex(Face::iminus,j,k,0);
                        f.right_cells[0] = activeCellIndex(i+1,j,k);
                        f.right_cells[1] = activeCellIndex(i+2,j,k);
                    } else if (i == 1) {
                        f.left_cells[0] = ghostCellIndex(Face::iminus,j,k,0);
                        f.left_cells[1] = activeCellIndex(i-1,j,k);
                        f.right_cells[0] = activeCellIndex(i+1,j,k);
                        f.right_cells[1] = activeCellIndex(i+2,j,k);
                    } else if (i == nic-1) {
                        f.left_cells[0] = activeCellIndex(i,j,k);
                        f.left_cells[1] = activeCellIndex(i-1,j,k);
                        f.right_cells[0] = activeCellIndex(i+1,j,k);
                        f.right_cells[1] = ghostCellIndex(Face::iplus,j,k,0);
                    } else if (i == nic) {
                        f.left_cells[0] = activeCellIndex(i,j,k);
                        f.left_cells[1] = activeCellIndex(i-1,j,k);
                        f.right_cells[0] = ghostCellIndex(Face::iplus,j,k,0);
                        f.right_cells[1] = ghostCellIndex(Face::iplus,j,k,1);
                    } else {
                        // Interior cell.
                        f.left_cells[0] = activeCellIndex(i,j,k);
                        f.left_cells[1] = activeCellIndex(i-1,j,k);
                        f.right_cells[0] = activeCellIndex(i+1,j,k);
                        f.right_cells[1] = activeCellIndex(i+2,j,k);
                    }
                }
            }
        }
        // jFaces
        for (int k=0; k < nkc; k++) {
            for (int i=0; i < nic; i++) {
                for (int j=0; j < njc+1; j++) {
                    FVFace& f = jFaces[jFaceIndex(i,j,k)];
                    f.vtx[0] = vtxIndex(i,j,k);
                    f.vtx[1] = vtxIndex(i+1,j,k);
                    f.vtx[2] = vtxIndex(i+1,j,k+1);
                    f.vtx[3] = vtxIndex(i,j,k+1);
                    if (j == 0) {
                        f.left_cells[0] = ghostCellIndex(Face::jminus,i,k,1);
                        f.left_cells[1] = ghostCellIndex(Face::jminus,i,k,0);
                        f.right_cells[0] = activeCellIndex(i,j+1,k);
                        f.right_cells[1] = activeCellIndex(i,j+2,k);
                    } else if (j == 1) {
                        f.left_cells[0] = ghostCellIndex(Face::jminus,i,k,0);
                        f.left_cells[1] = activeCellIndex(i,j-1,k);
                        f.right_cells[0] = activeCellIndex(i,j+1,k);
                        f.right_cells[1] = activeCellIndex(i,j+2,k);
                    } else if (j == njc-1) {
                        f.left_cells[0] = activeCellIndex(i,j,k);
                        f.left_cells[1] = activeCellIndex(i,j-1,k);
                        f.right_cells[0] = activeCellIndex(i,j+1,k);
                        f.right_cells[1] = ghostCellIndex(Face::jplus,i,k,0);
                    } else if (j == njc) {
                        f.left_cells[0] = activeCellIndex(i,j,k);
                        f.left_cells[1] = activeCellIndex(i,j-1,k);
                        f.right_cells[0] = ghostCellIndex(Face::jplus,i,k,0);
                        f.right_cells[1] = ghostCellIndex(Face::jplus,i,k,1);
                    } else {
                        // Interior cell.
                        f.left_cells[0] = activeCellIndex(i,j,k);
                        f.left_cells[1] = activeCellIndex(i,j-1,k);
                        f.right_cells[0] = activeCellIndex(i,j+1,k);
                        f.right_cells[1] = activeCellIndex(i,j+2,k);
                    }
                }
            }
        }
        // kFaces
        for (int j=0; j < njc; j++) {
            for (int i=0; i < nic; i++) {
                for (int k=0; k < nkc+1; k++) {
                    FVFace& f = kFaces[kFaceIndex(i,j,k)];
                    f.vtx[0] = vtxIndex(i,j,k);
                    f.vtx[1] = vtxIndex(i+1,j,k);
                    f.vtx[2] = vtxIndex(i+1,j+1,k);
                    f.vtx[3] = vtxIndex(i,j+1,k);
                    if (k == 0) {
                        f.left_cells[0] = ghostCellIndex(Face::kminus,i,j,1);
                        f.left_cells[1] = ghostCellIndex(Face::kminus,i,j,0);
                        f.right_cells[0] = activeCellIndex(i,j,k+1);
                        f.right_cells[1] = activeCellIndex(i,j,k+2);
                    } else if (k == 1) {
                        f.left_cells[0] = ghostCellIndex(Face::kminus,i,j,0);
                        f.left_cells[1] = activeCellIndex(i,j,k-1);
                        f.right_cells[0] = activeCellIndex(i,j,k+1);
                        f.right_cells[1] = activeCellIndex(i,j,k+2);
                    } else if (k == nkc-1) {
                        f.left_cells[0] = activeCellIndex(i,j,k);
                        f.left_cells[1] = activeCellIndex(i,j,k-1);
                        f.right_cells[0] = activeCellIndex(i,j,k+1);
                        f.right_cells[1] = ghostCellIndex(Face::kplus,i,j,0);
                    } else if (k == nkc) {
                        f.left_cells[0] = activeCellIndex(i,j,k);
                        f.left_cells[1] = activeCellIndex(i,j,k-1);
                        f.right_cells[0] = ghostCellIndex(Face::kplus,i,j,0);
                        f.right_cells[1] = ghostCellIndex(Face::kplus,i,j,1);
                    } else {
                        // Interior cell.
                        f.left_cells[0] = activeCellIndex(i,j,k);
                        f.left_cells[1] = activeCellIndex(i,j,k-1);
                        f.right_cells[0] = activeCellIndex(i,j,k+1);
                        f.right_cells[1] = activeCellIndex(i,j,k+2);
                    }
                }
            }
        }
        return;
    } // end configure()

    __host__
    void computeGeometry()
    // Compute cell and face geometric data.
    // Do this after reading the grid and flow files because we need the vertex locations
    // and because cell positions and volumes are part of the flow data.
    // This function will overwrite them with (potentially) better values.
    {
        for (int k=0; k < nkc; k++) {
            for (int j=0; j < njc; j++) {
                for (int i=0; i < nic; i++) {
                    FVCell& c = cells[activeCellIndex(i,j,k)];
                    hex_cell_properties(vertices[c.vtx[0]], vertices[c.vtx[1]],
                                        vertices[c.vtx[2]], vertices[c.vtx[3]],
                                        vertices[c.vtx[4]], vertices[c.vtx[5]],
                                        vertices[c.vtx[6]], vertices[c.vtx[7]],
                                        false, c.pos, c.volume, c.iLength, c.jLength, c.kLength);
                }
            }
        }
        // iFaces
        for (int k=0; k < nkc; k++) {
            for (int j=0; j < njc; j++) {
                for (int i=0; i < nic+1; i++) {
                    FVFace& f = iFaces[iFaceIndex(i,j,k)];
                    quad_properties(vertices[f.vtx[0]], vertices[f.vtx[1]],
                                    vertices[f.vtx[2]], vertices[f.vtx[3]],
                                    f.pos, f.n, f.t1, f.t2, f.area);
                }
            }
        }
        // jFaces
        for (int k=0; k < nkc; k++) {
            for (int i=0; i < nic; i++) {
                for (int j=0; j < njc+1; j++) {
                    FVFace& f = jFaces[jFaceIndex(i,j,k)];
                    quad_properties(vertices[f.vtx[0]], vertices[f.vtx[1]],
                                    vertices[f.vtx[2]], vertices[f.vtx[3]],
                                    f.pos, f.n, f.t1, f.t2, f.area);
                }
            }
        }
        // kFaces
        for (int j=0; j < njc; j++) {
            for (int i=0; i < nic; i++) {
                for (int k=0; k < nkc+1; k++) {
                    FVFace& f = kFaces[kFaceIndex(i,j,k)];
                    quad_properties(vertices[f.vtx[0]], vertices[f.vtx[1]],
                                    vertices[f.vtx[2]], vertices[f.vtx[3]],
                                    f.pos, f.n, f.t1, f.t2, f.area);
                }
            }
        }
        //
        // Work around the boundaries and extrapolate cell positions and lengths
        // into the ghost cells.  We need this data for high-order reconstruction
        // for the inviscid fluxes and for computation of the flow-property gradients
        // for the viscous fluxes.
        //
        // Face::iminus
        for (int k=0; k < nkc; k++) {
            for (int j=0; j < njc; j++) {
                FVFace& f = iFaces[iFaceIndex(0,j,k)];
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
        for (int k=0; k < nkc; k++) {
            for (int j=0; j < njc; j++) {
                FVFace& f = iFaces[iFaceIndex(nic,j,k)];
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
        for (int k=0; k < nkc; k++) {
            for (int i=0; i < nic; i++) {
                FVFace& f = jFaces[jFaceIndex(i,0,k)];
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
        for (int k=0; k < nkc; k++) {
            for (int i=0; i < nic; i++) {
                FVFace& f = jFaces[jFaceIndex(i,njc,k)];
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
        for (int j=0; j < njc; j++) {
            for (int i=0; i < nic; i++) {
                FVFace& f = kFaces[kFaceIndex(i,j,0)];
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
        for (int j=0; j < njc; j++) {
            for (int i=0; i < nic; i++) {
                FVFace& f = kFaces[kFaceIndex(i,j,nkc)];
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
    void readGrid(string fileName, bool vtkHeader=false)
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
        if ((nic != niv-1) || (njc != njv-1) || (nkc != nkv-1)) {
            throw runtime_error("Unexpected grid size: niv="+to_string(niv)+
                                " njv="+to_string(njv)+ " nkv="+to_string(nkv));
        }
        vertices.resize(niv*njv*nkv);
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
                    vertices[vtxIndex(i,j,k)].set(x, y, z);
                } // for i
            } // for j
        } // for k
        f.close();
        return;
    } // end readGrid()

    __host__
    void readFlow(string fileName)
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
                    for (int k=0; k < nkc; k++) {
                        for (int j=0; j < njc; j++) {
                            for (int i=0; i < nic; i++) {
                                getline(ss, item, '\n');
                                FVCell& c = cells[activeCellIndex(i,j,k)];
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
    void writeFlow(string fileName)
    // Writes the flow data into a new ZIP archive file.
    // Any necessary directories are presumed to exist.
    {
        int err = 0;
        zip *z = zip_open(fileName.c_str(), ZIP_CREATE, &err);
        if (err) {
            cerr << "Failed to open zip archive for writing: " << fileName << endl;
        }
        if (z) {
            for (int m=0; m < IOvar::n; m++) {
                string name = IOvar::names[m];
                ostringstream ss;
                for (int k=0; k < nkc; k++) {
                    for (int j=0; j < njc; j++) {
                        for (int i=0; i < nic; i++) {
                            FVCell& c = cells[activeCellIndex(i,j,k)];
                            ss << c.iovar_get(m) << endl;
                        }
                    }
                }
                string data = ss.str();
                // Add the data to the ZIP archive as a file.
                zip_source_t* zs = zip_source_buffer(z, data.c_str(), data.size(), 0);
                if (zs) {
                    int zindx = zip_file_add(z, name.c_str(), zs, ZIP_FL_OVERWRITE);
                    if (zindx < 0) {
                        cerr << "Could not add file " << name << " to ZIP archive " << fileName << endl;
                        zip_source_free(zs);
                    }
                } else {
                    cerr << "Error adding file to zip: " << string(zip_strerror(z)) << endl;
                    zip_source_free(zs);
                }
            }
            zip_close(z);
        }
        return;
    } // end writeFlow()

    __host__ __device__
    number estimate_local_dt(FVCell& c, Vector3 inorm, Vector3 jnorm, Vector3 knorm, number cfl)
    {
        // We assume that the cells are (roughly) hexagonal and work with
        // velocities normal to the faces.
        FlowState& fs = c.fs;
        number isignal = c.iLength/(fabs(fs.vel.dot(inorm))+fs.gas.a);
        number jsignal = c.jLength/(fabs(fs.vel.dot(jnorm))+fs.gas.a);
        number ksignal = c.kLength/(fabs(fs.vel.dot(knorm))+fs.gas.a);
        return cfl * fmin(fmin(isignal,jsignal),ksignal);
    } // end estimate_local_dt()

    __host__
    number estimate_allowed_dt(number cfl)
    {
        number smallest_dt = numeric_limits<number>::max();
        for (auto i=0; i < nActiveCells; i++) {
            FVCell& c = cells[i];
            Vector3 inorm = iFaces[c.face[Face::iminus]].n;
            Vector3 jnorm = jFaces[c.face[Face::jminus]].n;
            Vector3 knorm = kFaces[c.face[Face::kminus]].n;
            smallest_dt = fmin(smallest_dt, estimate_local_dt(c, inorm, jnorm, knorm, cfl));
        }
        return smallest_dt;
    } // end estimate_allowed_dt()

    __host__
    void encodeConserved(int level)
    {
        for (auto i=0; i < nActiveCells; i++) {
            FVCell& c = cells[i];
            ConservedQuantities& U = Q[level*nActiveCells + i];
            c.encode_conserved(U);
        }
        return;
    }

    __host__
    int decodeConserved(int level)
    {
        int bad_cell_count = 0;
        for (auto i=0; i < nActiveCells; i++) {
            FVCell& c = cells[i];
            ConservedQuantities U = Q[level*nActiveCells + i];
            int flag = c.decode_conserved(U);
            if (flag) { bad_cell_count += 1; }
        }
        return bad_cell_count;
    }

    __host__
    void calculate_fluxes(int x_order)
    {
        for (auto& face : iFaces) {
            FlowState fsL = cells[face.left_cells[0]].fs;
            FlowState fsR = cells[face.right_cells[0]].fs;
            ausmdv(face, fsL, fsR);
        }
        for (auto& face : jFaces) {
            FlowState fsL = cells[face.left_cells[0]].fs;
            FlowState fsR = cells[face.right_cells[0]].fs;
            ausmdv(face, fsL, fsR);
        }
        for (auto& face : kFaces) {
            FlowState& fsL = cells[face.left_cells[0]].fs;
            FlowState& fsR = cells[face.right_cells[0]].fs;
            ausmdv(face, fsL, fsR);
        }
        return;
    }

    __host__
    void eval_dUdt(FVCell& c, ConservedQuantities& dUdt)
    // These are the spatial (RHS) terms in the semi-discrete governing equations.
    {
        number vol_inv = 1.0/c.volume;
        auto& fim = iFaces[c.face[Face::iminus]];
        auto& fip = iFaces[c.face[Face::iplus]];
        auto& fjm = jFaces[c.face[Face::jminus]];
        auto& fjp = jFaces[c.face[Face::jplus]];
        auto& fkm = kFaces[c.face[Face::kminus]];
        auto& fkp = kFaces[c.face[Face::kplus]];
        //
        for (int i=0; i < CQI::n; i++) {
            // Integrate the fluxes across the interfaces that bound the cell.
            number surface_integral = 0.0;
            surface_integral = fim.area*fim.F[i] - fip.area*fip.F[i]
                + fjm.area*fjm.F[i] - fjp.area*fjp.F[i]
                + fkm.area*fkm.F[i] - fkp.area*fkp.F[i];
            // Then evaluate the derivatives of conserved quantity.
            // Note that conserved quantities are stored per-unit-volume.
            dUdt[i] = vol_inv*surface_integral;
        }
        return;
    } // end eval_dUdt()

    __host__
    int update_stage_1(number dt)
    // Predictor step.
    {
        cout << "Start update_state_1 for block " << endl;
        int bad_cell_count = 0;
        for (auto i=0; i < nActiveCells; i++) {
            FVCell& c = cells[i];
            ConservedQuantities dUdt = dQdt[i];
            eval_dUdt(c, dUdt);
            ConservedQuantities U0 = Q[i];
            ConservedQuantities U1 = Q[nActiveCells + i];
            for (int j=0; j < CQI::n; j++) {
                U1[j] = U0[j] + dt*dUdt[j];
            }
            int flag = c.decode_conserved(U1);
            if (flag) { bad_cell_count += 1; }
        }
        cout << "End update_state_1 for block: bad_cell_count=" << bad_cell_count << endl;
        return bad_cell_count;
    } // end update_stage_1()

    __host__
    void copy_conserved_data(int from_level, int to_level)
    {
        for (auto i=0; i < nActiveCells; i++) {
            ConservedQuantities U_from = Q[from_level*nActiveCells + i];
            ConservedQuantities U_to = Q[to_level*nActiveCells + i];
            for (int j=0; j < CQI::n; j++) {
                U_to[j] = U_from[j];
            }
        }
    } // end copy_conserved_data()

}; // end Block

#endif
