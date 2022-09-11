// gas.cu
// Include file for chicken.
// PJ 2022-09-11

#ifndef GAS_INCLUDED
#define GAS_INCLUDED

#include <string>
#include "number.cu"

using namespace std;

namespace GasModel {
    constexpr number g = 1.4f;
    constexpr number R = 287.1f;
    constexpr number Cv = R / (g-1.0f);
};


struct GasState {
    number rho;
    number e;
    number p;
    number T;
    number a;

    string toString() {
        return "GasState(p=" + to_string(p) + ", rho=" + to_string(rho) +
            ", T=" + to_string(T) + ", e=" + to_string(e) + ", a=" + to_string(a) + ")";
    }

    __host__ __device__ void update_from_pT()
    {
        using namespace GasModel;
        e = T*Cv;
        rho = p/(T*R);
        a = sqrt(g*R*T);
    }

    __host__ __device__ void update_from_rhoe()
    {
        using namespace GasModel;
        T = e/Cv;
        p = rho*R*T;
        a = sqrt(g*R*T);
    }
}; // end GasState

#endif
