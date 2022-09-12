// vertex.cu
// Include file for chicken.
// PJ 2022-09-12

#ifndef VERTEX_INCLUDED
#define VERTEX_INCLUDED

#include <string>
#include "number.cu"
#include "vector3.cu"

using namespace std;

struct Vertex {
    Vector3 pos; // position in space
    number dvelxdx, dvelxdy, dvelxdz;
    number dvelydx, dvelydy, dvelydz;
    number dvelzdx, dvelzdy, dvelzdz;
    number dTdx, dTdy, dTdz;

    string toString() {
        return "Vertex(pos=" + pos.toString() + ")";
    }
}; // end Vertex

#endif
