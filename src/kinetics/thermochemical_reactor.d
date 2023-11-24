/**
 * thermochemical_reactor.d
 * Base class for all of the kinetics package updaters.
 *
 * Author: Peter J and Rowan G.
 */

module kinetics.thermochemical_reactor;

import ntypes.complex;
import nm.number;
import gas;

class ThermochemicalReactorUpdateException : Exception {
    @nogc
    this(string message, string file=__FILE__, size_t line=__LINE__,
         Throwable next=null)
    {
        super(message, file, line, next);
    }
}

// The extra parameters are intended for some special gas and kinetics models,
// such as JJ Hoste's mixing-limited fuel-air model, that need context information
// from the flow solver.
// Most gas models and reaction schemes will just ignore the params array.
immutable size_t maxParams = 10;

class ThermochemicalReactor {
public:
    this(GasModel gmodel)
    {
        // We need a reference to the original gas model object
        // to update the GasState data at a later time.
        _gmodel = gmodel;
    }
    //
    // All the work happens when calling the concrete object
    // which updates the GasState over the (small) time, tInterval.
    //
    // The array params is there to allow extra information to be passed in.
    // For example, the mixing-limited combustion model by JJ Hoste needs
    // some information about the local flow state beyond the usual gas state.
    @nogc
    abstract void opCall(ref GasState Q, double tInterval, ref double dtSuggest,
                         ref number[maxParams] params);
    //
    // For the (fully-coupled) implicit solver we need public access to source terms.
    // The elements of this source-term array represent the rates of change of
    // the conserved quantities for the flow-solver.
    // For a multi-species gas, there will be first n_species elements for the chemical species
    // followed by n_modes of energies for a multi-temperature gas.
    @nogc
    abstract void eval_source_terms(GasModel gmodel, ref GasState Q, ref number[] source);
    //
    // We will need to access this referenced model from the Lua functions
    // so it needs to be public.
    GasModel _gmodel;
} // end class ThermochemicalReactor
