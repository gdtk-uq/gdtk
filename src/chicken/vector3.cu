// vector3.cu
// Include file for chicken, adapted from the Dlang module.
// PJ 2022-09-11

#ifndef VECTOR3_INCLUDED
#define VECTOR3_INCLUDED

#include <string>
#include "number.cu"

using namespace std;

struct Vector3 {
    number x;
    number y;
    number z;

    string toString() {
        return "[x=" + to_string(x) + ", y=" + to_string(y) + ", z=" + to_string(z) + "]";
    }

    __host__ __device__ number dot(const Vector3& other) const
    {
        return x*other.x + y*other.y + z*other.z;
    }

    // Transform functions used to reorient vector values in the CFD codes.

    __host__ __device__
    void transform_to_local_frame(const Vector3& n, const Vector3& t1, const Vector3& t2)
    // Rotate v from the global xyz coordinate system into the local frame
    // defined by the orthogonal unit vectors n,t1,t2.
    //
    // We assume, without checking, that these vectors do nicely define
    // such a local system.
    {
        number v_x = dot(n); // normal component
        number v_y = dot(t1); // tangential component 1
        number v_z = dot(t2); // tangential component 2
        x = v_x; y = v_y; z = v_z;
    }

    __host__ __device__
    void transform_to_global_frame(const Vector3& n, const Vector3& t1, const Vector3& t2)
    // Rotate v back into the global (xyz) coordinate system.
    {
        number v_x = x*n.x + y*t1.x + z*t2.x; // global-x
        number v_y = x*n.y + y*t1.y + z*t2.y; // global-y
        number v_z = x*n.z + y*t1.z + z*t2.z; // global-z
        x = v_x; y = v_y; z = v_z;
    }

    // 2D flavour for change of coordinate system functions.

    __host__ __device__
    void transform_to_local_frame(const Vector3& n, const Vector3& t1)
    {
        number v_x = x*n.x + y*n.y;   // normal component
        number v_y = x*t1.x + y*t1.y; // tangential component 1
        x = v_x; y = v_y; z = 0.0;
    }

    __host__ __device__
    void transform_to_global_frame(const Vector3& n, const Vector3& t1)
    // Rotate v back into the global (xy) coordinate system.
    {
        number v_x = x*n.x + y*t1.x; // global-x
        number v_y = x*n.y + y*t1.y; // global-y
        x = v_x; y = v_y; z = 0.0;
    }

}; // end Vector3

#endif
