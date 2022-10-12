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


struct FlowState {
    GasState gas;
    Vector3 vel;

    string toString() {
        return "FlowState(" + gas.toString() + ", vel=" + vel.toString() + ")";
    }

    __host__ __device__
    void encode_conserved(ConservedQuantities& U)
    {
        number rho = gas.rho;
        // Mass per unit volume.
        U[CQI::mass] = rho;
        // Momentum per unit volume.
        U[CQI::xMom] = rho*vel.x;
        U[CQI::yMom] = rho*vel.y;
        U[CQI::zMom] = rho*vel.z;
        // Total Energy / unit volume
        number ke = 0.5*(vel.x*vel.x + vel.y*vel.y+vel.z*vel.z);
        U[CQI::totEnergy] = rho*(gas.e + ke);
        return;
    } // end encode_conserved()

    __host__ __device__
    int decode_conserved(ConservedQuantities& U)
    // Returns 0 for success, 1 for the presence of nans or negative density.
    {
        bool any_nans = false; for (auto v : U) { if (isnan(v)) { any_nans = true; } }
        if (any_nans) return 1;
        number rho = U[CQI::mass];
        if (rho <= 0.0) return 1;
        number dinv = 1.0/rho;
        vel.set(U[CQI::xMom]*dinv, U[CQI::yMom]*dinv, U[CQI::zMom]*dinv);
        // Split the total energy per unit volume.
        // Start with the total energy, then take out the other components.
        // Internal energy is what remains.
        number e = U[CQI::totEnergy] * dinv;
        // Remove kinetic energy for bulk flow.
        number ke = 0.5*(vel.x*vel.x + vel.y*vel.y + vel.z*vel.z);
        e -= ke;
        // Put data into cell's gas state.
        gas.rho = rho;
        gas.e = e;
        gas.update_from_rhoe();
        return 0;
    } // end decode_conserved()

}; // end FlowState

#endif
