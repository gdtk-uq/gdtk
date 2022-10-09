// flow.cu
// Include file for chicken.
// PJ 2022-09-11

#ifndef FLOW_INCLUDED
#define FLOW_INCLUDED

#include <string>
#include <array>
#include "number.cu"
#include "vector3.cu"
#include "gas.cu"

using namespace std;

struct FlowState {
    GasState gas;
    Vector3 vel;

    string toString() {
        return "FlowState(" + gas.toString() + ", vel=" + vel.toString() + ")";
    }
}; // end FlowState

namespace CQI {
    // For sizing and indexing into the Flux vectors of conserved quantities.
    //
    // Whenever we introduce multiple species we can extent this namespace,
    // knowing that we are building the code for a particular number of species.
    // The plan is to have the gas model, number of species, etc,
    // known at compile time. [TODO] PJ 2022-09-11
    constexpr size_t n = 5;
    constexpr size_t mass = 0;
    constexpr size_t xMom = 1;
    constexpr size_t yMom = 2;
    constexpr size_t zMom = 3;
    constexpr size_t totEnergy = 4;
};

// We are going to store the vector of conserved quantities and its derivatives
// in a C++ arrays.
typedef array<number,CQI::n> ConservedQuantities;

#endif
