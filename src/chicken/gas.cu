// gas.cu
// Include file for chicken.
// We implement two flavours of gas model:
//   an ideal air model with no reactions
//   an ideal gas with simple chemical reaction and parameters from Powers & Aslam (2006)
//   For notes on the implementation of this reacting-gas model,
//   see PJ's workbook notes Jan 2017.
//
// PJ 2022-09-11
//    2022-11-16 Added reacting-gas model

#ifndef GAS_INCLUDED
#define GAS_INCLUDED

#include <cmath>
#include <string>
#include "number.cu"

using namespace std;

namespace GasModel {
    // Our selection of gas models is limited to one of two.
    constexpr int ideal_air = 0; // Ideal non-reacting gas, air.
    constexpr int a_b_reacting_gas = 1; // Reacting gas, species A to species B
    //
    // We must make a selection at compile time, here.
    #ifdef IDEAL_AIR
    constexpr int model = ideal_air;
    #endif
    #ifdef AB_REACTING_GAS
    constexpr int model = a_b_reacting_gas;
    #endif
    //
    constexpr number g = (model == ideal_air) ? 1.4 : 6.0/5.0; // ratio of specific heats
    constexpr number R = (model == ideal_air) ? 287.1 : 287.0; // gas constant J/kg/K
    constexpr number q = (model == ideal_air) ? 0.0 : 300000.0; // heat of reaction, J/kg
    constexpr number alpha = (model == ideal_air) ? 0.0 : 1000.0; // reaction rate, 1/s
    constexpr number Ti = (model == ideal_air) ? 0.0 : 362.58; // ignition temperature, K
    //
    constexpr number Cv = R / (g-1.0);
    constexpr number Cp = Cv + R;
    //
    // Sutherland viscosity model, air.
    constexpr number mu_ref = 1.716e-5;
    constexpr number T_ref = 273.0;
    constexpr number S = 111.0;
    //
    constexpr number Prandtl = 0.72;
};


struct GasState {
    number rho;  // density, kg/m^3
    number e;    // internal energy, J/kg
    number p;    // pressure, Pa
    number T;    // temperature, K
    number YB;   // mass fraction of species B, rhoB/rho
    number a;    // sound speed, m/s

    string toString() {
        return "GasState(p=" + to_string(p) + ", rho=" + to_string(rho) +
            ", T=" + to_string(T) + ", e=" + to_string(e) +
            ", YB=" + to_string(YB) + ", a=" + to_string(a) + ")";
    }


    __host__ __device__
    void update_from_pT()
    {
        using namespace GasModel;
        e = T*Cv - YB*q;
        rho = p/(T*R);
        a = sqrt(g*R*T);
    }

    __host__ __device__
    void update_from_rhoe()
    {
        using namespace GasModel;
        T = (e + YB*q)/Cv;
        p = rho*R*T;
        a = sqrt(g*R*T);
    }

    __host__ __device__
    void set_as_average(const GasState& a, const GasState& b)
    {
        rho = 0.5*(a.rho+b.rho);
        e = 0.5*(a.e+b.e);
        YB = 0.5*(a.YB+b.YB);
        update_from_rhoe();
    }

    __host__ __device__
    void trans_coeffs(number& mu, number& k)
    {
        using namespace GasModel;
        mu = mu_ref*sqrt(T/T_ref)*(T/T_ref)*(T_ref + S)/(T + S);
        k = Cp*mu/Prandtl;
    }

    __host__ __device__
    void update_chemistry(number dt)
    // Update the gas data, reaction is active in an isolated blob of gas.
    {
        using namespace GasModel;
        number YA = 1.0 - YB;
        if (T > Ti) { YB = 1.0 - YA*exp(-alpha*dt); }
        // In an isolated reactor, density and energy remain constant.
        update_from_rhoe();
    }

}; // end GasState

#endif
