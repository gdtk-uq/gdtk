/**
 * pvt_eos.d
 * Interface for all p-v-T equations of state.
 *
 * The PVT_EOS class is used to capture the
 * p-v-T behaviour of the gas, where v is the specific volume,
 * the inverse of density. This is sometimes referred to
 * as the thermal equation of state.
 *
 * The thermal equation of state is distinct
 * from the caloric equation of state. The latter
 * relates the internal energy of the gas to the
 * density and temperature. The two equations of
 * state together can be used to specify the 
 * complete thermodynamic state of the gas given
 * two state variables.
 *
 * Author: Rowan G. and Peter J.
 * Version: 2014-08-14 -- initial cut at building an interface
 */

module gas.thermo.pvt_eos;

import gas.gas_model;

/++
  PVT_EOS defines the services provided by a p-v-T equation
  of state.
+/
interface PVT_EOS {
    @nogc void update_pressure(ref GasState Q) const;
    @nogc void update_density(ref GasState Q) const;
    @nogc void update_temperature(ref GasState Q) const;
}
