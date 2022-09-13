// block.cu
// Include file for chicken.
// PJ 2022-09-11

#ifndef BLOCK_INCLUDED
#define BLOCK_INCLUDED

#include <string>
#include <fstream>
#include <stdexcept>

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
    constexpr int wall = 0;
    constexpr int exchange = 1;
    constexpr int inflow = 2;
    constexpr int outflow = 3;
    // [TODO] need to consider periodic boundary conditions.
    // There's not enough information here to have arbitrary block connections.
};

struct Block {
    int nic; // Number of cells i-direction.
    int njc; // Number of cells j-direction.
    int nkc; // Number of cells k-direction.
    int firstGhostCells[6]; // Index of the first ghost cell for each face.
    // Ghost cells will be stored at the end of the active cells collection.
    //
    vector<Cell> cells;
    vector<Face> iFaces;
    vector<Face> jFaces;
    vector<Face> kFaces;
    vector<Vector3> vertices;
    //
    // Active cells have conserved quantities data, along with the time derivatives.
    vector<ConservedQuantities> Q;
    vector<ConservedQuantities> dQdt;
    //
    int bcCodes[6];

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
        case FaceNames::iminus:
        case FaceNames::iplus:
            // jk face
            cellIndxOnFace = i1*njc + i0;
            break;
        case FaceNames::jminus:
        case FaceNames::jplus:
            // ik face
            cellIndxOnFace = i1*njc + i0;
            break;
        case FaceNames::kminus:
        case FaceNames::kplus:
            // ik face
            cellIndxOnFace = i1*njc + i0;
            break;
        }
        return firstGhostCells[faceIndx] + cellIndxOnFace;
    }

    __host__ __device__
    int iFaceIndx(int i, int j, int k)
    {
        return i*njc*nkc + k*njc + j;
    }

    __host__ __device__
    int jFaceIndx(int i, int j, int k)
    {
        return j*nic*nkc + k*nic + i;
    }

    __host__ __device__
    int kFaceIndx(int i, int j, int k)
    {
        return k*nic*njc + j*nic + i;
    }

    __host__ __device__
    int vtxIndx(int i, int j, int k)
    {
        return k*(nic+1)*(njc+1) + j*(nic+1) + i;
    }

    __host__
    void configure(int i, int j, int k, int codes[])
    // Do this before reading a grid or flow file.
    {
        nic = i;
        njc = j;
        nkc = k;
        for (int b=0; b < 6; b++) { bcCodes[b] = codes[b]; }
        return;
    }

    __host__
    void allocateStorage()
    {
        // [TODO] maybe we should do this after reading a grid.
        return;
    }

    __host__
    void readGrid(string fileName, bool vtkHeader=true)
    // Reads the vertex locations from a file, resizing storage as needed.
    // The numbers of cells are also set.
    {
        auto f = fstream(fileName, fstream::in);
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
            f.getline(line, maxc);
            sscanf(line, "%d %d %d", &niv, &njv, &nkv);
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
                    scanf(line, "%lf %lf %lf", &x, &y, &z);
                    #endif
                    vertices[vtxIndx(i,j,k)].set(x, y, z);
                } // for i
            } // for j
        } // for k
        f.close();
        return;
    } // end readGrid()

}; // end Block

#endif
