/**
 * equilibrium_update.d
 *
 * Thermally-perfect gas mix in equilibrium.
 *
 * This kinetics file accompanies the gas model in gas/therm_perf_gas_equil.d
 *
 * Authors: Nick Gibbons/PJ.
 * Version: 2020-xx-xx
 */

module kinetics.equilibrium_update;

import std.math;
import std.conv : to;
import std.stdio;
import std.algorithm;
import std.range : enumerate;
import std.string;
import nm.complex;
import nm.number;

import gas;
import gas.therm_perf_gas_equil;
import util.lua;
import util.lua_service;
import kinetics.thermochemical_reactor;
import ceq;

final class EquilibriumUpdate : ThermochemicalReactor {
    
    this(string fname, GasModel gmodel)
    {
        super(gmodel); // hang on to a reference to the gas model
        // We may need to pick a number of pieces out of the file.
        auto L = init_lua_State();
        doLuaFile(L, fname);

        // Copy some things we need from the gasmodel
        nsp = to!int(_gmodel.n_species);
        foreach(Mi; _gmodel.mol_masses) M ~= Mi;
        foreach(i; 0 .. nsp) species_names ~= _gmodel.species_name(i);

        // Build arrays to mimic pyeq memory management
        compile_lewis_array(L, lewis);
        compile_element_map(L, element_map);
        compile_element_set(element_map, element_set);
        compile_element_matrix(species_names, element_map, element_set, a);
        nel = to!int(element_set.length);
        X0.length = nsp; X1.length = nsp;

        aptr = a.ptr;
        X0ptr = X0.ptr;
        X1ptr = X1.ptr;
        lewisptr = lewis.ptr;
        Mptr = M.ptr;
    }
    
    override void opCall(GasState Q, double tInterval,
                         ref double dtChemSuggest, ref double dtThermSuggest, 
                         ref number[maxParams] params)
    {
        // Since the internal energy and density in the (isolated) reactor is fixed,
        // we need to evaluate the new temperature, pressure, etc.
        int error;

        version(complex_numbers){
            throw new Error("Do not use equilibrium_update with complex numbers.");
        } else {
            _gmodel.massf2molef(Q, X0); 
            error = ceq.rhou(Q.rho,Q.u,X0ptr,nsp,nel,lewisptr,Mptr,aptr,X1ptr,&Q.T,0);
            if (error!=0) throw new GasModelException("ceq.rhou convergence failure");
            _gmodel.molef2massf(X1, Q);

            // Note that T is already set from ceq.rhou
            // This slightly messes with the u, but it isn't by much
            _gmodel.update_thermo_from_rhoT(Q);
            _gmodel.update_sound_speed(Q);
        }
    }

private:
    int nel,nsp;
    double[string][] element_map;
    string[] element_set,species_names;
    double[] a,X0,X1,lewis,M;
    double* aptr, X0ptr, X1ptr, lewisptr, Mptr, Tptr;

    void compile_lewis_array(lua_State* L, ref double[] _lewis){
        /*
        ceq needs lewis thermo data coefficients packed into a single array
        */

        _lewis.length = 0;
        double[] lewis_s;
        double[] coeffs;

        foreach (species; species_names) {
            lewis_s.length = 0;

            lua_getglobal(L, "db");
            lua_getfield(L, -1, species.toStringz);
            lua_getfield(L, -1, "thermoCoeffs");
            auto nseg = getInt(L, -1, "nsegments");
            if (nseg<3) {
                throw new Error(format("Species %s wrong number of thermo segments: (%d)",species, nseg));
            }

            // TODO: ceq assumes fixed thermo tables with three segments in them.
            // if a species with more segments is present it only uses the first 3
            foreach (i; 0 .. 3) {
                auto key = format("segment%d", i);
                getArrayOfDoubles(L, -1, key, coeffs);
                lewis_s ~= coeffs;
            }
            lua_pop(L, 1);
            lua_pop(L, 1);
            lua_pop(L, 1);
            _lewis ~= lewis_s;
        }
        return;
    }

    void compile_element_map(lua_State* L, ref double[string][] _element_map){
        /*
        Get array of associative arrays that tells you how many of each element in each species
        */
        double[string] species_elements;

        foreach(species; species_names){
            species_elements.clear();
            lua_getglobal(L, "db");
            lua_getfield(L, -1, species.toStringz);
            lua_getfield(L, -1, "atomicConstituents");

            // For some reason this completely idiomatic operation take TEN lines ...
            // It also segfaults if atomicConstituents is not present :)
            lua_pushnil(L); // dummy first key
            while (lua_next(L, -2) != 0) { // -1 is the dummy key, -2 is the atomicConstituents table
                // lua_next removes the dummy key and puts the first key value pair on the stack
                // the key at -2, value at -1
                string key = to!string(lua_tostring(L, -2));
                int value = to!int(lua_tointeger(L, -1));
                species_elements[key] = to!double(value);
                lua_pop(L, 1); // discard value but keep key so that lua_next can remove it (?!)
            }
            lua_pop(L, 1); // remove atomicConstituents (lua_next removed the key when it broke loop)
            lua_pop(L, 1); // remove species table
            lua_pop(L, 1); // remove db table
            _element_map ~= species_elements.dup();
        }
        return;
    }

    void compile_element_set(double[string][] element_map, ref string[] _element_set){
        /*
        Get an alphabetical list of all the elements in the complete set of species
        Inputs: 
             element_map : array of associative arrays mapping species to elements
        Outputs:
            elements : list of elements in the entire system
        */
    
        foreach(e; element_map){
            foreach(key; e.keys()){
                if (!_element_set.canFind(key)) _element_set ~= key;
            }
        }
        _element_set.sort();
        return;
    }

    void compile_element_matrix(string[] species_names, double[string][] element_map,
                                string[] element_set,
                                ref double[] a){
        /*
        Get the number of each element in a species from its lewis table data
        Inputs: 
             species_names : array of strings of each species
             element_map   : array of associative arrays mapping species to elements
             element_set   : array of strings for each element (sorted alphabetically)
        Outputs:
            a : element_matrix mapping species to elemental composition 
        */
        size_t _nsp,_nel,j;
    
        _nsp = species_names.length;
        _nel = element_set.length; 
        a.length = _nel*_nsp;
        a[] = 0.0;
    
        foreach(i, atoms; enumerate(element_map)){
            foreach (key, value; atoms){
                j = countUntil(element_set, key);
                a[j*_nsp + i] = value;
            }
        }
        return;
    }

} // end class EquilibriumUpdate

version(equilibrium_update_test) {
    import std.stdio;
    import util.msg_service;
    import gas.therm_perf_gas;

    int main() {
        auto gm = new ThermallyPerfectGasEquilibrium("../gas/sample-data/therm-perf-equil-5-species-air.lua");
        auto reactor = new EquilibriumUpdate("../gas/sample-data/therm-perf-equil-5-species-air.lua", gm);

        auto gs2 = new GasState(5, 0);
        double rho_target = 0.0139638507337;
        double u_target = 2131154.032665843;
        double T_target = 2500.0;
        gs2.rho = rho_target;
        gs2.u = u_target;
        gs2.massf = [0.74311527, 0.25688473, 0.0, 0.0, 0.0];
        double tInterval = 0.0;
        double dtChemSuggest = 0.0;
        double dtThermSuggest = 0.0;
        double[maxParams] params;

        reactor(gs2, tInterval, dtChemSuggest, dtThermSuggest, params);
        assert(approxEqual(0.7321963 , gs2.massf[0], 1.0e-6)); 
        assert(approxEqual(0.23281198, gs2.massf[1], 1.0e-6));
        assert(approxEqual(0.0, gs2.massf[2], 1.0e-6));
        assert(approxEqual(0.01160037, gs2.massf[3], 1.0e-6));
        assert(approxEqual(0.02339135, gs2.massf[4], 1.0e-6));
        assert(approxEqual(T_target, gs2.T, 1.0e-6));
        assert(approxEqual(rho_target, gs2.rho, 1.0e-6));
        assert(approxEqual(u_target, gs2.u, 1.0e-1));
        return 0;
    }
}
