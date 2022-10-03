// face.cu
// Include file for chicken.
// PJ 2022-09-11

#ifndef FACE_INCLUDED
#define FACE_INCLUDED

#include <string>
#include <sstream>
#include <vector>
#include <array>

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
    array<number,CQI::n> F; // flux vector for conserved quantities
    // We will keep connections to the pieces composing the face
    // as indices into global arrays.
    array<int,4> vtx{0, 0, 0, 0};
    array<int,2> left_cells{0, 0};
    array<int,2> right_cells{0, 0};

    string toString() {
        ostringstream repr;
        repr << "FVFace(pos=" << pos.toString() << ", n=" << n.toString()
             << ", t1=" << t1.toString() << ", t2=" << t2.toString() << ", area=" << area;
        repr << ", vtx=["; for (auto i : vtx) repr << i << ","; repr << "]";
        repr << ", left_cells=["; for (auto i : left_cells) repr << i << ","; repr << "]";
        repr << ", right_cells=["; for (auto i : right_cells) repr << i << ","; repr << "]";
        repr << ")";
        return repr.str();
    }

}; // end FVFace

#endif
