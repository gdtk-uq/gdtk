/**
 * evt_eos.d
 * Interface for all e-v-T equations of state.
 *
 * The e-v-T equation of state relates the
 * internal energy of the gas to the density
 * and temperature. This is sometimes called the
 * caloric equation of state in the literature.
 *
 * See the counterpart interface file pvt_eos.d
 * for the thermal equation of state (p-v-T behaviour).
 *
 * Author: Rowan G. and Peter J.
 * Version: 2014-08-18 -- initial cut at building the interface
 */

module gas.thermo.evt_eos;

import gas.gas_model;

/++
  EVT_EOS defines the services provied by a caloric equation
  of state model.

  All caloric equations of state define a relationship between
  the internal energy, temperature and density of a gas.

+/
interface EVT_EOS {
    void update_energy(ref GasState Q);
    void update_temperature(ref GasState Q);
}

