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

constexpr size_t TLevels = 3;

struct Cell {
    Vector3 pos; // position of centroid
    number volume;
    FlowState fs;
    number Q[TLevels][CQI::n]; // Conserved quantities.
    number dQdt[TLevels][CQI::n]; // Time derivatives of conserved quantities.
    // We will keep connections to the pieces composing the cell
    // as indices into global arrays.
    int vtx[8]{0, 0, 0, 0};
    int face[6]{0, 0, 0, 0, 0, 0};

    string toString() {
        string repr = "Cell(pos=" + pos.toString() +
            ", volume=" + to_string(volume) + ", Q=[";
        for (int i=0; i < TLevels; i++) {
            repr += "[";
            for (int j=0; j < CQI::n; j++) { repr += to_string(Q[i][j]) + " "; }
            repr += "] ";
        }
        repr += "] )";
        return repr;
    }

    // [TODO] PJ 2022-09-11 construct cell details from the connected faces and vertices.
    // Cell update and IO methods, etc.

}; // end Cell

#endif
