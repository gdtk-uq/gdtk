// face.cu
// Include file for chicken.
// PJ 2022-09-11

#ifndef FACE_INCLUDED
#define FACE_INCLUDED

#include <string>
#include "number.cu"
#include "vector3.cu"
#include "vertex.cu"
#include "flow.cu"

using namespace std;

struct FVFace {
    Vector3 pos; // midpoint position in space
    number area;
    Vector3 n;  // unit normal
    Vector3 t1; // unit tangent 1
    Vector3 t2; // unit tangent 2
    number F[CQI::n]; // flux vector for conserved quantities
    // We will keep connections to the pieces composing the face
    // as indices into global arrays.
    int vtx[4]{0, 0, 0, 0};
    int left_cells[2]{0, 0};
    int right_cells[2]{0, 0};

    string toString() {
        return "FVFace(n=" + n.toString() + ", t1=" + t1.toString() +
            ", t2=" + t2.toString() + ")";
    }
}; // end FVFace

#endif
