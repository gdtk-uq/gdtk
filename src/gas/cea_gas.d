/**
 * cea_gas.d
 * Chemical-equilibrium gas model for use in the CFD codes.
 * This gas model code delegates the computations to the NASA CEA2 code.
 *
 * It interfaces to the CEA code by writing a small input file,
 * running the CEA code as a child process and then reading the results
 * from the CEA output files.
 * You will need an executable file and its support database files.
 *
 * Author: Peter J. and Rowan G.
 * Version: 2017-03-04: port the essential parts of our Python module.
 */

module gas.cea_gas;

import gas.gas_model;
import gas.physical_constants;
import std.math;
import std.stdio;
import std.string;
import std.file;
import std.path;
import std.conv;
import std.algorithm;
import std.process;
import util.lua;
import util.lua_service;
import core.stdc.stdlib : exit;

class CEAGas: GasModel {
public:

    this(string mixtureName, string[] speciesList, double[string] reactants,
	 string inputUnits, string outputUnits, double trace, bool withIons)
    {
	_mixtureName = mixtureName;
	_n_modes = 0;
	_species_names = speciesList.dup();
	_n_species = to!uint(_species_names.length);
	_reactants = reactants.dup();
	_inputUnits = inputUnits;
	_outputUnits = outputUnits;
	_trace = trace;
	_withIons = withIons;
	//
	create_species_reverse_lookup();
	//
	_cea_exe_path = environment.get("CEA_EXE_PATH", "~/e3bin/cea2");
	_cea_cases_path = environment.get("CEA_CASES_PATH", "~/e3bin/cea-cases");
	string ceaExeFile = expandTilde(_cea_exe_path);
	writeln("ceaExeFile=", ceaExeFile);
	if (!exists(ceaExeFile)) {
	    throw new Exception("Cannot find cea2 exe file.");
	}
    } // end constructor
    
    this(lua_State *L) {
	// Construct from information in a Lua table.
	lua_getglobal(L, "CEAGas"); // Bring that table to TOS
	string name = getString(L, -1, "mixtureName");
	// We will specify the full set of species in a table.
	// Although CEA2 does not need to have this specified,
	// we will rely upon the table being correctly filled.
	string[] speciesList;
	getArrayOfStrings(L, -1, "speciesList", speciesList);
	// For the amounts of the reactants, we are expecting
	// to find only the non-zero components in the Lua file
	// but we will use whatever we find.
	lua_getfield(L, -1, "reactants");
	double[string] reactants;
	foreach (sname; speciesList) {
	    lua_getfield(L, -1, sname.toStringz);
	    if (lua_isnumber(L, -1)) {
		reactants[sname] = lua_tonumber(L, -1);
	    } else {
		reactants[sname] = 0.0;
	    }
	    lua_pop(L, 1);
	}
	lua_pop(L, 1); // dispose of reactants table
	string inputUnits = getString(L, -1, "inputUnits");
	string outputUnits = getString(L, -1, "outputUnits");
	double trace = getDouble(L, -1, "trace");
	bool withIons = getBool(L, -1, "withIons");
	// We have finished gathering the information from the Lua file.
	// If we have any ionic species already listed, we need to
	// be sure that the withIons flag is set true.
	foreach (sname; speciesList) {
	    if (canFind(sname, "+")) { withIons = true; }
	    if (canFind(sname, "-")) { withIons = true; }
	}	
	this(name, speciesList, reactants, inputUnits, outputUnits, trace, withIons);
    } // end constructor from a Lua file

    override string toString() const
    {
	char[] repr;
	repr ~= "CEAGas(";
	repr ~= "mixtureName=\"" ~ _mixtureName ~"\"";
	repr ~= ", speciesList=[";
	foreach (i, sname; _species_names) {
	    repr ~= format("\"%s\"", sname);
	    repr ~= (i+1 < _species_names.length) ? ", " : "]";
	}
	repr ~= ", reactants=[";
	string[] reactantNames = _reactants.keys;
	foreach (i, sname; reactantNames) {
	    repr ~= format("\"%s\"=%g", sname, _reactants[sname]);
	    repr ~= (i+1 < reactantNames.length) ? ", " : "]";
	}
	repr ~= ", inputUnits=\"" ~ _inputUnits ~ "\"";
	repr ~= ", outputUnits=\"" ~ _outputUnits ~ "\"";
	repr ~= ", trace=" ~ to!string(_trace);
	repr ~= ", withIons=" ~ to!string(_withIons);
	repr ~= ", cea_exe_path=\"" ~ _cea_exe_path ~ "\"";
	repr ~= ", cea_cases_path=\"" ~ _cea_cases_path ~ "\"";
	repr ~= ")";
	return to!string(repr);
    }

    override void update_thermo_from_pT(GasState Q) const 
    {
	// [TODO] replace
	Q.rho = Q.p/(Q.Ttr*_Rgas);
	Q.u = _Cv*Q.Ttr;
    }
    override void update_thermo_from_rhoe(GasState Q) const
    {
	// [TODO] replace
	Q.Ttr = Q.u/_Cv;
	Q.p = Q.rho*_Rgas*Q.Ttr;
    }
    override void update_thermo_from_rhoT(GasState Q) const
    {
	// [TODO] replace
	Q.p = Q.rho*_Rgas*Q.Ttr;
	Q.u = _Cv*Q.Ttr;
    }
    override void update_thermo_from_rhop(GasState Q) const
    {
	// [TODO] replace
	Q.Ttr = Q.p/(Q.rho*_Rgas);
	Q.u = _Cv*Q.Ttr;
    }
    
    override void update_thermo_from_ps(GasState Q, double s) const
    {
	// [TODO] replace
	Q.Ttr = _T1 * exp((1.0/_Cp)*((s - _s1) + _Rgas * log(Q.p/_p1)));
	update_thermo_from_pT(Q);
    }
    override void update_thermo_from_hs(GasState Q, double h, double s) const
    {
	// [TODO] replace
	Q.Ttr = h / _Cp;
	Q.p = _p1 * exp((1.0/_Rgas)*(_s1 - s + _Cp*log(Q.Ttr/_T1)));
	update_thermo_from_pT(Q);
    }
    override void update_sound_speed(GasState Q) const
    {
	// [TODO] replace
	Q.a = sqrt(_gamma*_Rgas*Q.Ttr);
    }
    override void update_trans_coeffs(GasState Q)
    {
	// [TODO] replace
	// need to do something, or not
    }
    override double dudT_const_v(in GasState Q) const
    {
	// [TODO] replace
	return _Cv;
    }
    override double dhdT_const_p(in GasState Q) const
    {
	// [TODO] replace
	return _Cp;
    }
    override double dpdrho_const_T(in GasState Q) const
    {
	// [TODO] replace
	double R = gas_constant(Q);
	return R*Q.Ttr;
    }
    override double gas_constant(in GasState Q) const
    {
	// [TODO] replace
	return _Rgas;
    }
    override double internal_energy(in GasState Q) const
    {
	// [TODO] replace
	return Q.u;
    }
    override double enthalpy(in GasState Q) const
    {
	// [TODO] replace
	return Q.u + Q.p/Q.rho;
    }
    override double entropy(in GasState Q) const
    {
	// [TODO] replace
	return _s1 + _Cp * log(Q.Ttr/_T1) - _Rgas * log(Q.p/_p1);
    }

private:
    string _mixtureName;
    string _inputUnits; // "moles" or "massf"
    string _outputUnits;
    double _trace; // fraction below which CEA ignores species
    double[string] _reactants;
    string[] _onlyList;
    bool _withIons;
    string _cea_exe_path; // expect to find cea2 executable file here
    string _cea_cases_path; // expect to find thermo.lib and trans.lib
    
    // [TODO] replace
    // Thermodynamic constants
    double _Rgas; // J/kg/K
    double _gamma;   // ratio of specific heats
    double _Cv; // J/kg/K
    double _Cp; // J/kg/K
    // Reference values for entropy
    double _s1;  // J/kg/K
    double _T1;  // K
    double _p1;  // Pa

    // [TODO] private helper functions to write and scan cea2 data files
    // and to run the cea2 program.
} // end class CEAGgas

version(cea_gas_test) {
    import std.stdio;
    import util.msg_service;

    int main() {
	lua_State* L = init_lua_State("sample-data/cea-air5species-gas-model.lua");
	auto gm = new CEAGas(L);
	lua_close(L);
	writeln("gm=", gm); // for debug
	auto gd = new GasState(1, 0);
	gd.p = 1.0e5;
	gd.Ttr = 300.0;
	// FIXME gd.massf[0] = 1.0;
	assert(approxEqual(gm.R(gd), 287.086, 1.0e-4), failedUnitTest());
	assert(gm.n_modes == 0, failedUnitTest());
	assert(gm.n_species == 1, failedUnitTest());
	assert(approxEqual(gd.p, 1.0e5, 1.0e-6), failedUnitTest());
	assert(approxEqual(gd.Ttr, 300.0, 1.0e-6), failedUnitTest());
	assert(approxEqual(gd.massf[0], 1.0, 1.0e-6), failedUnitTest());

	gm.update_thermo_from_pT(gd);
	gm.update_sound_speed(gd);
	assert(approxEqual(gd.rho, 1.16109, 1.0e-4), failedUnitTest());
	assert(approxEqual(gd.u, 215314.0, 1.0e-4), failedUnitTest());
	assert(approxEqual(gd.a, 347.241, 1.0e-4), failedUnitTest());
	/+
	gm.update_trans_coeffs(gd);
	assert(approxEqual(gd.mu, 1.84691e-05, 1.0e-6), failedUnitTest());
	assert(approxEqual(gd.k, 0.0262449, 1.0e-6), failedUnitTest());
	+/

	return 0;
    }
}
