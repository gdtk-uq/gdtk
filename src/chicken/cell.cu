// cell.cu
// Include file for chicken.
// PJ 2022-09-11

#ifndef CELL_INCLUDED
#define CELL_INCLUDED

#include <string>
#include "number.cu"
#include "vector3.cu"
#include "gas.cu"
#include "flow.cu"
#include "face.cu"

using namespace std;

namespace FaceNames {
    // Symbolic names for the faces of the cell and of the block.
    constexpr int iminus = 0;
    constexpr int iplus = 1;
    constexpr int jminus = 2;
    constexpr int jplus = 3;
    constexpr int kminus = 4;
    constexpr int kplus = 5;
};

struct Cell {
    Vector3 pos; // position of centroid
    number volume;
    number iLength, jLength, kLength; // These lengths are used in the interpolation fns.
    FlowState fs;
    // We will keep connections to the pieces compising the cell
    // as indices into global arrays.
    int vtx[8]{0, 0, 0, 0};
    int face[6]{0, 0, 0, 0, 0, 0};

    string toString() {
        string repr = "Cell(pos=" + pos.toString() +
            ", volume=" + to_string(volume) + ", Q=[";
        // for (int i=0; i < TLevels; i++) {
        //     repr += "[";
        //     for (int j=0; j < CQI::n; j++) { repr += to_string(Q[i][j]) + " "; }
        //     repr += "] ";
        // }
        repr += "] )";
        return repr;
    }

    // [TODO] PJ 2022-09-11 construct cell details from the connected faces and vertices.
    // Cell update and IO methods, etc.

}; // end Cell

#endif
