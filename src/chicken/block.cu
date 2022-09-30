// block.cu
// Include file for chicken.
// PJ 2022-09-11

#ifndef BLOCK_INCLUDED
#define BLOCK_INCLUDED

#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include "include/bxzstr/bxzstr.hpp"
#include <zip.h>

#include "number.cu"
#include "vector3.cu"
#include "gas.cu"
#include "vertex.cu"
#include "flow.cu"
#include "face.cu"
#include "cell.cu"

using namespace std;

namespace BCCode {
    // Boundary condition codes, to decide what to do for the ghost cells.
    // Periodic boundary conditions should just work if we wrap the index in each direction.
    // There's not enough information here to have arbitrary block connections.
    constexpr int wall_with_slip = 0;
    constexpr int wall_no_slip = 1;
    constexpr int exchange = 2;
    constexpr int inflow = 3;
    constexpr int outflow = 4;

    array<string,5> names{"wall_with_slip", "wall_no_slip", "exchange", "inflow", "outflow"};
};

int BC_code_from_name(string name)
{
    if (name == "wall_with_slip") return BCCode::wall_with_slip;
    if (name == "wall_no_slip") return BCCode::wall_no_slip;
    if (name == "exchange") return BCCode::exchange;
    if (name == "inflow") return BCCode::inflow;
    if (name == "outflow") return BCCode::outflow;
    return BCCode::wall_with_slip;
}

struct Block {
    int nic; // Number of cells i-direction.
    int njc; // Number of cells j-direction.
    int nkc; // Number of cells k-direction.
    int nActiveCells; // Number of active cells (with conserved quantities) in the block.
    // Ghost cells will be stored at the end of the active cells collection.
    array<int,6> nGhostCells; // Number of ghost cells on each face.
    array<int,6> firstGhostCells; // Index of the first ghost cell for each face.
    //
    vector<FVCell> cells;
    vector<FVFace> iFaces;
    vector<FVFace> jFaces;
    vector<FVFace> kFaces;
    vector<Vector3> vertices;
    //
    // Active cells have conserved quantities data, along with the time derivatives.
    vector<ConservedQuantities> Q;
    vector<ConservedQuantities> dQdt;
    //
    array<int,6> bcCodes;

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
        int cellIndxOnFace = 0;
        switch (faceIndx) {
        case Face::iminus:
        case Face::iplus:
            // jk face
            cellIndxOnFace = i1*njc + i0;
            break;
        case Face::jminus:
        case Face::jplus:
            // ik face
            cellIndxOnFace = i1*njc + i0;
            break;
        case Face::kminus:
        case Face::kplus:
            // ik face
            cellIndxOnFace = i1*njc + i0;
            break;
        }
        return firstGhostCells[faceIndx] + cellIndxOnFace;
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
    void configure(int i, int j, int k, int codes[])
    // Set up the block to hold the grid and flow data.
    // Do this before reading a grid or flow file.
    {
        nic = i;
        njc = j;
        nkc = k;
        for (int b=0; b < 6; b++) { bcCodes[b] = codes[b]; }
        //
        int nActiveCells = nic*njc*nkc;
        // For the moment assume that all boundary conditions require ghost cells.
        nGhostCells[Face::iminus] = 2*njc*nkc;
        nGhostCells[Face::iplus] = 2*njc*nkc;
        nGhostCells[Face::jminus] = 2*nic*nkc;
        nGhostCells[Face::jplus] = 2*nic*nkc;
        nGhostCells[Face::kminus] = 2*nic*njc;
        nGhostCells[Face::kplus] = 2*nic*njc;
        firstGhostCells[0] = nActiveCells;
        for (int f=1; f < 6; f++) firstGhostCells[f] = firstGhostCells[f-1] + nGhostCells[f-1];
        //
        // Now that we know the numbers of cells, resize the vector to fit them all.
        cells.resize(firstGhostCells[5]+nGhostCells[5]);
        Q.resize(nActiveCells*TLevels);
        dQdt.resize(nActiveCells*TLevels);
        // Each set of finite-volume faces is in the index-plane of the corresponding vertices.
        iFaces.resize((nic+1)*njc*nkc);
        jFaces.resize(nic*(njc+1)*nkc);
        kFaces.resize(nic*njc*(nkc+1));
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
            for (int j=0; j < njc+1; j++) {
                for (int i=0; i < nic; i++) {
                    // [TODO]
                }
            }
        }
        // kFaces
        for (int k=0; k < nkc+1; k++) {
            for (int j=0; j < njc; j++) {
                for (int i=0; i < nic; i++) {
                    // [TODO]
                }
            }
        }
        return;
    }

    __host__
    void readGrid(string fileName, bool vtkHeader=false)
    // Reads the vertex locations from a compressed file, resizing storage as needed.
    // The numbers of cells are also checked.
    {
        auto f = bxz::ifstream(fileName); // gzip file
        if (!f) {
            throw new runtime_error("Did not open grid file successfully: "+fileName);
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
            throw new runtime_error("Unexpected grid size: niv="+to_string(niv)+
                                    " njv="+to_string(njv)+
                                    " nkv="+to_string(nkv));
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

}; // end Block

#endif
