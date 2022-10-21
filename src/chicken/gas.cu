// gas.cu
// Include file for chicken.
// PJ 2022-09-11

#ifndef GAS_INCLUDED
#define GAS_INCLUDED

#include <cmath>
#include <string>
#include "number.cu"

using namespace std;

namespace GasModel {
    // Ideal gas, air.
    constexpr number g = 1.4;
    constexpr number R = 287.1;
    constexpr number Cv = R / (g-1.0);
    constexpr number Cp = Cv + R;
    // Sutherland viscosity model, air.
    constexpr number mu_ref = 1.716e-5;
    constexpr number T_ref = 273.0;
    constexpr number S = 111.0;
    //
    constexpr number Prandtl = 0.72;
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


    __host__ __device__
    void update_from_pT()
    {
        using namespace GasModel;
        e = T*Cv;
        rho = p/(T*R);
        a = sqrt(g*R*T);
    }

    __host__ __device__
    void update_from_rhoe()
    {
        using namespace GasModel;
        T = e/Cv;
        p = rho*R*T;
        a = sqrt(g*R*T);
    }

    __host__ __device__
    void set_as_average(const GasState& a, const GasState& b)
    {
        rho = 0.5*(a.rho+b.rho);
        e = 0.5*(a.e+b.e);
        update_from_rhoe();
    }

    __host__ __device__
    void trans_coeffs(number& mu, number& k)
    {
        using namespace GasModel;
        mu = mu_ref*sqrt(T/T_ref)*(T/T_ref)*(T_ref + S)/(T + S);
        k = Cp*mu/Prandtl;
    }

}; // end GasState

#endif
