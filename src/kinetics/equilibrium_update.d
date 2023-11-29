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
import ntypes.complex;
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
        eqcalc = new EquilibriumCalculator(fname);
    }

    override void opCall(ref GasState Q, double tInterval, ref double dtSuggest,
                         ref number[maxParams] params)
    {
        // Since the internal energy and density in the (isolated) reactor is fixed,
        // we need to evaluate the new temperature, pressure, etc.
        int error;

        eqcalc.set_massf_and_T_from_rhou(Q);
        // Note that T is already set from ceq.rhou
        // This slightly messes with the u, but it isn't by much
        _gmodel.update_thermo_from_rhoT(Q);
        _gmodel.update_sound_speed(Q);
    }

    @nogc override void eval_source_terms(GasModel gmodel, ref GasState Q, ref number[] source) {
        string errMsg = "eval_source_terms not implemented for equilibrium_update.";
        throw new ThermochemicalReactorUpdateException(errMsg);
    }

private:
    EquilibriumCalculator eqcalc;

} // end class EquilibriumUpdate

class EquilibriumCalculator {
    /*
       D version of pyeq's equilibrium calculator.
       Rather than have the code here duplicated in separate reactor and gas model
       classes, all of the ceq machinery is encapulated in this class which can be
       instantiated whereever it is needed.

       @author: Nick Gibbons
    */
    this(string fname)
    {
        /*
           Constructor requires a thermally_perfect_gas or a thermally_perfect_gas_equilibrium file
           @author: Nick Gibbons
        */
        auto L = init_lua_State();
        doLuaFile(L, fname);

        getArrayOfStrings(L, "species", species_names);
        nsp = to!int(species_names.length);

        // Build arrays to mimic pyeq memory management
        compile_molar_masses(L, species_names, M);
        compile_lewis_array(L, lewis);
        compile_element_map(L, element_map);
        compile_element_set(element_map, element_set);
        compile_element_matrix(species_names, element_map, element_set, a);
        nel = to!int(element_set.length);
        X0.length = nsp; X1.length = nsp;

        // c routines interact with the arrays through these pointers
        aptr = a.ptr;
        X0ptr = X0.ptr;
        X1ptr = X1.ptr;
        lewisptr = lewis.ptr;
        Mptr = M.ptr;
    }

version(complex_numbers){

    @nogc void set_massf_from_pT(ref GasState Q)
    {
        throw new Error("Do not use with complex numbers.");
    }

    @nogc void set_massf_and_T_from_rhou(ref GasState Q)
    {
        throw new Error("Do not use with complex numbers.");
    }

    @nogc void set_massf_and_T_from_ps(ref GasState Q, double s)
    {
        throw new Error("Do not use with complex numbers.");
    }

    @nogc void set_massf_from_rhoT(ref GasState Q)
    {
        throw new Error("Do not use with complex numbers.");
    }

    @nogc double get_s(ref GasState Q)
    {
        throw new Error("Do not use with complex numbers.");
    }

} else {
    @nogc void set_massf_from_pT(ref GasState Q)
    {
        massf2molef(Q.massf, M, X0);
        int error = ceq.pt(Q.p,Q.T,X0ptr,nsp,nel,lewisptr,Mptr,aptr,X1ptr,0);
        if (error!=0) throw new GasModelException("ceq.pt convergence failure");
        molef2massf(X1, M, Q.massf);
    }

    @nogc void set_massf_and_T_from_rhou(ref GasState Q)
    {
        massf2molef(Q.massf, M, X0);
        int error = ceq.rhou(Q.rho,Q.u,X0ptr,nsp,nel,lewisptr,Mptr,aptr,X1ptr,&Q.T,0);
        if (error!=0) throw new GasModelException("ceq.rhou convergence failure");
        molef2massf(X1, M, Q.massf);
    }

    @nogc void set_massf_and_T_from_ps(ref GasState Q, double s)
    {
        massf2molef(Q.massf, M, X0);
        int error = ceq.ps(Q.p,s,X0ptr,nsp,nel,lewisptr,Mptr,aptr,X1ptr,&Q.T,0);
        if (error!=0) throw new GasModelException("ceq.ps convergence failure");
        molef2massf(X1, M, Q.massf);
    }

    @nogc void set_massf_from_rhoT(ref GasState Q)
    {
        massf2molef(Q.massf, M, X0);
        int error = ceq.rhot(Q.rho,Q.T,X0ptr,nsp,nel,lewisptr,Mptr,aptr,X1ptr,0);
        if (error!=0) throw new GasModelException("ceq.rhot convergence failure");
        molef2massf(X1, M, Q.massf);
    }

    @nogc double get_s(ref GasState Q)
    {
        massf2molef(Q.massf, M, X0);
        return ceq.get_s(Q.T, Q.p, X0ptr, nsp, lewisptr, Mptr);
    }
}
    @nogc int n_species() const { return nsp; }
    @nogc int n_elements() const { return nel; }

private:
    int nel,nsp;
    double[string][] element_map;
    string[] element_set,species_names;
    double[] a,X0,X1,lewis,M;
    double* aptr, X0ptr, X1ptr, lewisptr, Mptr, Tptr;

    void compile_molar_masses(lua_State* L, string[] _species_names, ref double[] _M){
        /*
        Extra an array of molar masses from a gas file
        */
        lua_getglobal(L, "db");
        foreach (species; _species_names) {
            lua_getfield(L, -1, species.toStringz);
            _M ~= getDouble(L, -1, "M");
            lua_pop(L, 1);
        }
        lua_pop(L, 1);
    }

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
            double Mi = getDouble(L, -1, "M");

            lua_getfield(L, -1, "thermoCoeffs");
            auto nseg = getInt(L, -1, "nsegments");
            if (nseg<2) {
                throw new Error(format("Species %s wrong number of thermo segments: (%d)",species, nseg));
            }

            foreach (i; 0 .. 2) {
                auto key = format("segment%d", i);
                getArrayOfDoubles(L, -1, key, coeffs);
                lewis_s ~= coeffs;
            }
            if (nseg==2) {
                lewis_s ~= fix_missing_thermo_segment(species, lewis_s, Mi);
            } else {
                auto key = format("segment%d", 2);
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

    double[] fix_missing_thermo_segment(string species, double[] lewis_s, double Mi){
    /*
        ceq's thermodynamic routines assume a fixed layout of curve-fit data with 3 sections
        for each species. Since some species in the Lewis database only have 2 segments,
        this routine can be used to make a make a fake segment that assumes constant Cp above
        6000 K.
    */
        if (lewis_s.length!=18)
            throw new Error(format("Cannot fix missing thermo segment of species (%s", species));

        double[27] lspa;
        foreach(i, li; lewis_s) lspa[i] = li;

        // We need to call the thermo machinery with the half-built lspa table
        int nsp = 1;
        double T = 5999.99999;
        double Ru = 8.3144621;
        double[1] X = [1.0];
        double[1] M = [Mi];

        double cp = ceq.get_cp(T, X.ptr, nsp, lspa.ptr, M.ptr);
        double h  = ceq.get_h(T, X.ptr, nsp, lspa.ptr, M.ptr);
        double s0 = ceq.get_s0(T, X.ptr, nsp, lspa.ptr, M.ptr);

        double a2 = cp*Mi/Ru;
        double b1 = h*Mi/Ru - a2*T;
        double b2 = s0*Mi/Ru - a2*log(T);
        double[] lsp2 = [0.0, 0.0, a2, 0.0, 0.0, 0.0, 0.0, b1, b2];
        return lsp2;
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

            // ceq deals with ionised species with a pretend Electron element named E
            int charge = getInt(L, -1, "charge");
            if (charge!=0) species_elements["E"] = -1*charge;

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

} // end class EquilibriumCalculator


version(equilibrium_update_test) {
    import std.stdio;
    import util.msg_service;
    import gas.therm_perf_gas;

    int main() {
        auto gm = new ThermallyPerfectGasEquilibrium("../gas/sample-data/therm-perf-equil-5-species-air.lua");
        auto reactor = new EquilibriumUpdate("../gas/sample-data/therm-perf-equil-5-species-air.lua", gm);

        auto gs2 = GasState(5, 0);
        double rho_target = 0.0139638507337;
        double u_target = 2131154.032665843;
        double T_target = 2500.0;
        gs2.rho = rho_target;
        gs2.u = u_target;
        gs2.massf = [0.74311527, 0.25688473, 0.0, 0.0, 0.0];
        gs2.T = 2000.0; // ceq doesn't guess temperature anymore
        double tInterval = 0.0;
        double dtSuggest = -1.0;
        double[maxParams] params;

        reactor(gs2, tInterval, dtSuggest, params);
        // writeln("T: ", gs2.T);
        assert(isClose(0.7321963 , gs2.massf[0], 1.0e-6));
        assert(isClose(0.23281198, gs2.massf[1], 1.0e-6));
        assert(isClose(0.0, gs2.massf[2], 1.0e-6));
        assert(isClose(0.01160037, gs2.massf[3], 1.0e-6));
        assert(isClose(0.02339135, gs2.massf[4], 1.0e-6));
        assert(isClose(T_target, gs2.T, 1.0e-6));
        assert(isClose(rho_target, gs2.rho, 1.0e-6));
        assert(isClose(u_target, gs2.u, 1.0e-1));
        return 0;
    }
}
