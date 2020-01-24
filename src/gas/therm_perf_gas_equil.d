// therm_perf_gas_equil.d
// Thermally-perfect gas with equilibrium chemistry.

module gas.therm_perf_gas_equil;

import std.math;
import std.stdio;
import std.string;
import std.conv : to;
import std.algorithm;
import std.range : enumerate;
import util.lua;
import util.lua_service;
import nm.complex;
import nm.number;

import gas.gas_model;
import gas.gas_state;
import gas.therm_perf_gas;
import ceq;

class ThermallyPerfectGasEquilibrium: ThermallyPerfectGas {
public:
    this(lua_State* L)
    // Construct the model from parameters that are contained in a Lua interpreter,
    // but delegate all of the hard work to Rowan's ThermallyPerfectGas.
    {
        super(L);

        // Build arrays to mimic pyeq memory management
        compile_lewis_array(L, lewis);
        compile_element_map(L, element_map);
        compile_element_set(element_map, element_set);
        compile_element_matrix(_species_names, element_map, element_set, a);
        nel = to!int(element_set.length);
        nsp = to!int(_n_species);
        X0.length = nsp; X1.length = nsp;

        aptr = a.ptr;
        X0ptr = X0.ptr;
        X1ptr = X1.ptr;
        lewisptr = lewis.ptr;
        Mptr = _mol_masses.ptr;

    } // end constructor using Lua interpreter

    this(in string fname)
    {
        auto L = init_lua_State();
        doLuaFile(L, fname);
        this(L);
        lua_close(L); // We no longer need the Lua interpreter.
    } // end constructor from a Lua file

version(complex_numbers){
    // This... Has turned into a difficult situation...
    override void update_thermo_from_pT(GasState Q) 
    {
        throw new Error("Do not use with complex numbers.");
    }
} else {
    override void update_thermo_from_pT(GasState Q) 
    {
        int error;
        massf2molef(Q, X0); 
        error = ceq.pt(Q.p,Q.T,X0ptr,nsp,nel,lewisptr,Mptr,aptr,X1ptr,0);
        molef2massf(X1, Q);
        _pgMixEOS.update_density(Q);
        _tpgMixEOS.update_energy(Q);
    }
}

private:
    int nel,nsp;
    double[string][] element_map;
    string[] element_set;
    double[] a,X0,X1,lewis;
    double* aptr, X0ptr, X1ptr, lewisptr, Mptr, Tptr;

    void compile_lewis_array(lua_State* L, ref double[] _lewis){
        /*
        ceq needs lewis thermo data coefficients packed into a single array
        */

        _lewis.length = 0;
        double[] lewis_s;
        double[] coeffs;

        foreach ( isp; 0.._n_species ) {
            lewis_s.length = 0;

            lua_getglobal(L, "db");
            lua_getfield(L, -1, _species_names[isp].toStringz);
            lua_getfield(L, -1, "thermoCoeffs");
            auto nseg = getInt(L, -1, "nsegments");
            assert(nseg==3);

            foreach (i; 0 .. nseg) {
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

        foreach(species; _species_names){
            species_elements.clear();
            lua_getglobal(L, "db");
            lua_getfield(L, -1, species.toStringz);
            lua_getfield(L, -1, "atomicConstituents");

            // For some reason this completely idiomatic operation take TEN lines ...
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

    void compile_element_matrix(string[] speciesList, double[string][] element_map,
                                string[] element_set,
                                ref double[] a){
        /*
        Get the number of each element in a species from its lewis table data
        Inputs: 
             speciesList : array of strings of each species
             element_map : array of associative arrays mapping species to elements
             element_set : array of strings for each element (sorted alphabetically)
        Outputs:
            a : element_matrix mapping species to elemental composition 
        */
        size_t _nsp,_nel,j;
    
        _nsp = speciesList.length;
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

} // end class ThermallyPerfectGasEquilibrium

version(therm_perf_gas_equil_test) {
    int main() {
        import util.msg_service;

        FloatingPointControl fpctrl;
        // Enable hardware exceptions for division by zero, overflow to infinity,
        // invalid operations, and uninitialized floating-point variables.
        // Copied from https://dlang.org/library/std/math/floating_point_control.html
        fpctrl.enableExceptions(FloatingPointControl.severeExceptions);
        //
        auto gm = new ThermallyPerfectGasEquilibrium("sample-data/therm-perf-equil-5-species-air.lua");
        auto gd = new GasState(5, 0);

        gd.p = 0.1*101.35e3;
        gd.T = 2500.0;
        gd.massf = [0.74311527, 0.25688473, 0.0, 0.0, 0.0];
        gm.update_thermo_from_pT(gd);
        assert(approxEqual(0.7321963 , gd.massf[0], 1.0e-6)); 
        assert(approxEqual(0.23281198, gd.massf[1], 1.0e-6));
        assert(approxEqual(0.0, gd.massf[2], 1.0e-6));
        assert(approxEqual(0.01160037, gd.massf[3], 1.0e-6));
        assert(approxEqual(0.02339135, gd.massf[4], 1.0e-6));

        return 0;
    } // end main()
}
