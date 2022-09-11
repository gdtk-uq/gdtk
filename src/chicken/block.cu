// block.cu
// Include file for chicken.
// PJ 2022-09-11

#ifndef BLOCK_INCLUDED
#define BLOCK_INCLUDED

#include <string>
#include "number.cu"
#include "vector3.cu"
#include "gas.cu"
#include "flow.cu"
#include "face.cu"
#include "cell.cu"

using namespace std;

namespace BCCode {
    // Boundary condition codes, to decide what to do for the ghost cells.
    constexpr int wall = 0;
    constexpr int exchange = 1;
};

struct Block {
    int nic;
    int njc;
    int nkc;
    //
    vector<Cell> cells;
    vector<Face> ifaces;
    vector<Face> jfaces;
    vector<Face> kfaces;
    vector<Vector3> vertices;
    //
    int bc_codes[6];
    //
    string toString() {
        string repr = "Block(nic=" + to_string(nic) +
            ", njc=" + to_string(njc) + ", nkc=" + to_string(nkc) + ")";
        return repr;
    }

    // [TODO] PJ 2022-09-11 methods for allocation of the storage arrays, etc.

}; // end Block

#endif
