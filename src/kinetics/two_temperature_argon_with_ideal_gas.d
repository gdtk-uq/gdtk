/**
 * two_temperature_argon_with_ideal_gas.d
 *
 * Composite model of two-temperature reacting argon and an Ideal gas.
 * The intention is to accommodate flows in which there are distinct regions
 * of ideal gas and reacting argon gas.
 * The reactions within the argon species proceed once they have been isolated
 * from the ideal gas component of the gas state.
 * This is different to the situation where the ideal gas participates as
 * a third-body in the argon reactions.
 *
 * Authors: Peter J.
 * Version: 2020-Feb-19
 */

module kinetics.two_temperature_argon_with_ideal_gas;

import std.stdio : writeln;
import std.format;
import std.math;
import std.conv;
import ntypes.complex;
import nm.number;
import nm.bbla;

import gas;
import gas.two_temperature_argon_plus_ideal;
import util.lua;
import util.lua_service;
import kinetics.thermochemical_reactor;
import kinetics.two_temperature_argon_kinetics;

final class UpdateArgonFracWithIdeal : ThermochemicalReactor {

    this(string fname, GasModel gmodel)
    {
        super(gmodel); // hang on to a reference to the gas model
        auto mygm = cast(TwoTemperatureArgonPlusIdealGas) gmodel;
        if (mygm !is null) {
            massf_tiny = mygm.massf_tiny;
            ideal_gas = mygm.ideal_gas;
            Q_ideal = new GasState(ideal_gas);
            argon_gas = mygm.argon_gas;
            Q_argon = new GasState(argon_gas);
            argon_reactor = new UpdateArgonFrac(fname, argon_gas);
        } else {
            throw new ThermochemicalReactorUpdateException("Need a 2T argon gas but did not get one.");
        }
    }

    @nogc
    override void opCall(ref GasState Q, double tInterval, ref double dtSuggest,
                         ref number[maxParams] params)
    {
        bool with_ideal = Q.massf[0] > massf_tiny;
        number argon_massf = Q.massf[1]+Q.massf[2]+Q.massf[3];
        bool with_argon = argon_massf > massf_tiny;
        number internal_energy = _gmodel.internal_energy(Q); // remains fixed
        if (with_argon) {
            // Abstract the reacting argon species and do the update on them in isolation.
            // Fixed volume, density and energy process
            Q_argon.rho = argon_massf*Q.rho;
            Q_argon.T = Q.T;
            Q_argon.T_modes[0] = Q.T_modes[0];
            Q_argon.massf[0] = Q.massf[1]/argon_massf;
            Q_argon.massf[1] = Q.massf[2]/argon_massf;
            Q_argon.massf[2] = Q.massf[3]/argon_massf;
            argon_gas.update_thermo_from_rhoT(*Q_argon);
            argon_reactor(*Q_argon, tInterval, dtSuggest, params);
            // Bring the updated argon species back into the composite model
            Q.massf[1] = Q_argon.massf[0]*argon_massf;
            Q.massf[2] = Q_argon.massf[1]*argon_massf;
            Q.massf[3] = Q_argon.massf[2]*argon_massf;
            Q.T_modes[0] = Q_argon.T_modes[0];
            Q.u_modes[0] = Q_argon.u_modes[0];
            // The changes in species mass fractions and electron energy
            // will have corresponding changes in thermal energy, pressure, etc.
            Q.u = internal_energy - argon_massf*Q.u_modes[0];
            _gmodel.update_thermo_from_rhou(Q);
        } else {
            // Just ideal gas is present, which does not react.
            // Do nothing.
        }
        return;
    }

    @nogc override void eval_source_terms(GasModel gmodel, ref GasState Q, ref number[] source) {
        string errMsg = "eval_source_terms not implemented for two_temperature_argon_with_ideal_gas.";
        throw new ThermochemicalReactorUpdateException(errMsg);
    }

private:
    double massf_tiny;
    GasModel ideal_gas;
    GasModel argon_gas;
    GasState* Q_ideal;
    GasState* Q_argon;
    ThermochemicalReactor argon_reactor;
} // end class UpdateArgonFracWithIdeal


version(two_temperature_argon_with_ideal_gas_test) {
    import std.stdio;
    import util.msg_service;
    import std.math : isClose;
    import gas.two_temperature_argon_plus_ideal;
    void main() {
        // [TODO] Put a meaningful test in this place.
        assert(isClose(1.0, 1.0, 1.0e-3), failedUnitTest());
    } // end main()
} // end two_temperature_argon_with_ideal_gas_test



