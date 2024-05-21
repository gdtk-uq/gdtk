/**
 * A module related to representing the state/properties in a solid cell.
 *
 * Author: RJG and KAD
 * Version: 2024-04-27
 */

module lmr.solid.solidstate;

import nm.number;

struct SolidState {
    number e; // energy, J/kg
    number T; // temperature, K
    number rho; // density, kg/m^3
    number k; // thermal conductivity, W/(m.K)
    number Cp; // specific heat at constant p, J/(kg.K)
}

@nogc
void harmonicAverage(ref SolidState lft, ref SolidState rght, ref SolidState iface)
{
    iface.rho = 2.0/(1.0/lft.rho + 1.0/rght.rho);
    iface.k = 2.0/(1.0/lft.k + 1.0/rght.k);
    iface.Cp = 2.0/(1.0/lft.Cp + 1.0/rght.Cp);
}


