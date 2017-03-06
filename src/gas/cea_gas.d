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
 *
 * Note: Do not use this gas model in a parallel calculation.
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
import std.array;
import std.format;
import util.lua;
import util.lua_service;
import core.stdc.stdlib : exit;


struct CEASavedData {
    // This holds some parameters of the gas state that are scanned
    // from the CEA output file and may be needed at a later time.
    // We will put them into the GasState object only when the
    // GasModel is a CEAGas.
public:
    double p; // Pa
    double rho; // kg/m**3
    double u; // J/kg
    double T; // K
    double a; // m/s
    double Rgas; // J/kg/K
    double gamma; // ratio of specific heats
    double Cv; // J/kg/K
    double Cp; // J/kg/K
    double s;  // J/kg/K
    double[] massf;
    
    this(const CEASavedData other) {
	this.p = other.p;
	this.rho = other.rho;
	this.u = other.u;
	this.T = other.T;
	this.a = other.a;
	this.Rgas = other.Rgas;
	this.gamma = other.gamma;
	this.Cv = other.Cv;
	this.Cp = other.Cp;
	this.s = other.s;
	this.massf = other.massf.dup();
    }
} // end CEASavedData


class CEAGas: GasModel {
public:

    this(string mixtureName, string[] speciesList, double[string] reactants,
	 string inputUnits, double trace, bool withIons)
    {
	_mixtureName = mixtureName;
	_n_modes = 0;
	_species_names = speciesList.dup();
	_n_species = to!uint(_species_names.length);
	_reactants = reactants.dup();
	_inputUnits = inputUnits;
	_trace = trace;
	_withIons = withIons;
	//
	create_species_reverse_lookup();
	//
	// Make sure that all of the pieces are in place to run CEA calculations.
	_cea_exe_path = expandTilde(environment.get("CEA_EXE_PATH", "~/e3bin/cea2"));
	_cea_cases_path = expandTilde(environment.get("CEA_CASES_PATH", "~/e3bin/cea-cases"));
	// writeln("_cea_exe_path=", _cea_exe_path); // DEBUG
	if (!exists(_cea_exe_path)) {
	    throw new Exception("Cannot find cea2 exe file.");
	}
	if (!exists("thermo.lib") || getSize("thermo.lib") == 0) {
	    string ceaThermoFile = buildNormalizedPath(_cea_cases_path, "thermo.inp");
	    writeln("ceaThermoFile=", ceaThermoFile); // DEBUG
	    if (!exists(ceaThermoFile)) {
		throw new Exception("Cannot find cea2 thermo.inp file.");
	    }
	    auto copy_thermo = executeShell("cp " ~ ceaThermoFile ~ " .");
	    runCEAProgram("thermo", false);
	}
	if (!exists("trans.lib") || getSize("trans.lib") == 0) {
	    string ceaTransFile = buildNormalizedPath(_cea_cases_path, "trans.inp");
	    writeln("ceaTransFile=", ceaTransFile); // DEBUG
	    if (!exists(ceaTransFile)) {
		throw new Exception("Cannot find cea2 trans.inp file.");
	    }
	    auto copy_trans = executeShell("cp " ~ ceaTransFile ~ " .");
	    runCEAProgram("trans", false);
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
	double trace = getDouble(L, -1, "trace");
	bool withIons = getBool(L, -1, "withIons");
	// We have finished gathering the information from the Lua file.
	// If we have any ionic species already listed, we need to
	// be sure that the withIons flag is set true.
	foreach (sname; speciesList) {
	    if (canFind(sname, "+")) { withIons = true; }
	    if (canFind(sname, "-")) { withIons = true; }
	}	
	this(name, speciesList, reactants, inputUnits, trace, withIons);
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
	repr ~= ", trace=" ~ to!string(_trace);
	repr ~= ", withIons=" ~ to!string(_withIons);
	repr ~= ", cea_exe_path=\"" ~ _cea_exe_path ~ "\"";
	repr ~= ", cea_cases_path=\"" ~ _cea_cases_path ~ "\"";
	repr ~= ")";
	return to!string(repr);
    }

    // All of the update functions call up CEA behind the scene
    // to do the real calculations.

    override void update_thermo_from_pT(GasState Q) const 
    {
	callCEA(Q, 0.0, 0.0, "pT", false);
    }
    override void update_thermo_from_rhoe(GasState Q) const
    {
	callCEA(Q, 0.0, 0.0, "rhoe", false);
    }
    override void update_thermo_from_rhoT(GasState Q) const
    {
	callCEA(Q, 0.0, 0.0, "rhoT", false);
    }
    override void update_thermo_from_rhop(GasState Q) const
    {
	callCEA(Q, 0.0, 0.0, "rhop", false);
    }
    
    override void update_thermo_from_ps(GasState Q, double s) const
    {
	callCEA(Q, 0.0, s, "ps", false);
    }
    override void update_thermo_from_hs(GasState Q, double h, double s) const
    {
	callCEA(Q, h, 0.0, "hs", false);
    }
    override void update_sound_speed(GasState Q) const
    {
	callCEA(Q, 0.0, 0.0, "sound_speed", false);
    }
    override void update_trans_coeffs(GasState Q)
    {
	callCEA(Q, 0.0, 0.0, "pT", true);
    }

    // The following functions return the most-recently-computed values.
    // It is assumed that you have previously called an update function
    // so that the values sitting in the CEAGas object are relevant.

    override double dudT_const_v(in GasState Q) const
    {
	return Q.ceaSavedData.Cv;
    }
    override double dhdT_const_p(in GasState Q) const
    {
	return Q.ceaSavedData.Cp;
    }
    override double dpdrho_const_T(in GasState Q) const
    {
	return Q.ceaSavedData.Rgas * Q.Ttr;
    }
    override double gas_constant(in GasState Q) const
    {
	return Q.ceaSavedData.Rgas;
    }
    override double internal_energy(in GasState Q) const
    {
	return Q.u;
    }
    override double enthalpy(in GasState Q) const
    {
	return Q.u + Q.p/Q.rho;
    }
    override double entropy(in GasState Q) const
    {
	return Q.ceaSavedData.s;
    }

private:
    string _mixtureName;
    string _inputUnits; // "moles" or "massf"
    string _outputUnits;
    double _trace; // fraction below which CEA ignores species
    double[string] _reactants;
    bool _withIons;
    string _cea_exe_path; // expect to find cea2 executable file here
    string _cea_cases_path; // expect to find thermo.lib and trans.lib

    void runCEAProgram(string jobName, bool checkTableHeader=true) const
    {
	string inpFile = jobName ~ ".inp";
	string outFile = jobName ~ ".out";
	string pltFile = jobName ~ ".plt";
	if (!exists(inpFile)) {
	    throw new Exception(format("CEAGas cannot find %s", inpFile));
	}
	if (exists(outFile)) { remove(outFile); }
	if (exists(pltFile)) { remove(pltFile); }
	auto pp = pipeProcess(_cea_exe_path, Redirect.stdin);
	pp.stdin.writeln(jobName);
	pp.stdin.flush();
	pp.stdin.close();
	auto returnCode = wait(pp.pid);
	if (returnCode) {
	    throw new Exception(format("CEAGas: cea program return code nonzero %d",
				       returnCode));
	}
	if (checkTableHeader) {
	    // scan the outFile looking for summary table header.
	    auto outFileText = readText(outFile);
	    if (!canFind(outFileText, "THERMODYNAMIC PROPERTIES")) {
		throw new Exception("CEAGas: the cea output file seems incomplete.");
	    }
	}
    } // end runCEAProgram()

    string cleanFloat(string[] tokens) const
    // Clean up the CEA2 short-hand notation for exponential format.
    // CEA2 seems to write exponential-format numbers in a number of ways:
    // 1.023-2
    // 1.023+2
    // 1.023 2
    {
	string valueStr;
	switch (tokens.length) {
	case 0:
	    valueStr = "0.0";
	    break;
	case 1:
	    valueStr = tokens[0];
	    if (canFind(valueStr, "****")) { valueStr = "nil"; }
	    if (canFind(valueStr, "-")) { valueStr = valueStr.replace("-", "e-"); }
	    if (canFind(valueStr, "+")) { valueStr = valueStr.replace("+", "e+"); }
	    break;
	case 2:
	    valueStr = tokens[0] ~ "e+" ~ tokens[1];
	    break;
	default:
	    throw new Exception("CEAGas: cleanFloat received too many tokens.");
	}
	return valueStr;
    } // end cleanFloat()

    void callCEA(GasState Q, double h, double s,
		 string problemType="pT", bool transProps=true) const
    {
	// Write input file for CEA that is specific to the problemType.
	string inpFileName = "tmp.inp";
	auto writer = appender!string();
	writer.put(format("# %s generated by CEAGas\n", inpFileName));
	switch (problemType) {
	case "pT":
	    writer.put("problem case=CEAGas tp");
            if (_withIons) { writer.put(" ions"); }
	    writer.put("\n");
            assert(Q.p > 0.0 && Q.Ttr > 0.0, "CEAGas: Invalid pT");
            writer.put(format("   p(bar)      %e\n", Q.p / 1.0e5));
            writer.put(format("   t(k)        %e\n", Q.Ttr));
	    break;
	case "rhoe":
	    // [TODO]
	    break;
	case "rhoT":
	    // [TODO]
	    break;
	case "rhop":
	    // [TODO]
	    break;
	case "ps":
	    // [TODO]
	    break;
 	case "hs":
	    // [TODO]
	    break;
	case "sound_speed":
	    // [TODO]
	    break;
	default: 
	    throw new Exception("Unknown problemType for CEA.");
	}
	// Select the gas components for CEA.
	// Note that the speciesList contains all of the reactants,
	// and that the constructor would have made a reactants entry
	// for each species.
	writer.put("react\n");
	foreach (name; _species_names) {
            double frac = _reactants[name];
            if (frac > 0.0) {
		writer.put(format("   name= %s  %s=%g", name,
				  ((_inputUnits == "moles") ? "moles" : "wtf"),
				  frac));
                if (canFind(["ph", "rhoe"], problemType)) { writer.put(" t=300"); }
                writer.put("\n");
	    }
	}
	writer.put("only");
	foreach (name; _species_names) { writer.put(format(" %s", name)); }
	writer.put("\n");
	writer.put("output massf");
	writer.put(format(" trace=%e", _trace));
	if (transProps) { writer.put(" trans"); }
	writer.put("\n");
	writer.put("end\n");
	std.file.write("tmp.inp", writer.data);
	//
	// Run the CEA program and check that the summary header is present
	// in the output file.  This is a crude test for success.
	runCEAProgram("tmp", true);
	//
	// Scan the output file for all of the bits that we need.
	// [TODO]
	//
	switch (problemType) {
	case "pT":
	    // [TODO]
	    Q.rho = Q.p/(Q.Ttr*Q.ceaSavedData.Rgas);
	    Q.u = Q.ceaSavedData.Cv*Q.Ttr;
	    break;
	case "rhoe":
	    // [TODO]
	    Q.Ttr = Q.u/Q.ceaSavedData.Cv;
	    Q.p = Q.rho*Q.ceaSavedData.Rgas*Q.Ttr;
	    break;
	case "rhoT":
	    // [TODO]
	    Q.p = Q.rho*Q.ceaSavedData.Rgas*Q.Ttr;
	    Q.u = Q.ceaSavedData.Cv*Q.Ttr;
	    break;
	case "rhop":
	    // [TODO]
	    Q.Ttr = Q.p/(Q.rho*Q.ceaSavedData.Rgas);
	    Q.u = Q.ceaSavedData.Cv*Q.Ttr;
	    break;
	case "ps":
	    // [TODO]
	    // Q.Ttr = _T1 * exp((1.0/_Cp)*((s - _s1) + _Rgas * log(Q.p/_p1)));
	    update_thermo_from_pT(Q);
	    break;
 	case "hs":
	    // [TODO]
	    Q.Ttr = h / Q.ceaSavedData.Cp;
	    // Q.p = _p1 * exp((1.0/_Rgas)*(_s1 - s + _Cp*log(Q.Ttr/_T1)));
	    update_thermo_from_pT(Q);
	    break;
	case "sound_speed":
	    // [TODO]
	    Q.a = sqrt(Q.ceaSavedData.gamma*Q.ceaSavedData.Rgas*Q.Ttr);
	    break;
	default: 
	    throw new Exception("Unknown problemType for CEA.");
	}
    } // end callCEA()
} // end class CEAGgas

version(cea_gas_test) {
    import std.stdio;
    import util.msg_service;

    int main() {
	lua_State* L = init_lua_State("sample-data/cea-air5species-gas-model.lua");
	auto gm = new CEAGas(L);
	lua_close(L);
	writeln("gm=", gm); // for debug
	auto gd = new GasState(1, 0, true);
	gd.p = 1.0e5;
	gd.Ttr = 300.0;
	gm.update_thermo_from_pT(gd);
	assert(approxEqual(gm.R(gd), 287.1, 0.1), failedUnitTest());
	assert(gm.n_modes == 0, failedUnitTest());
	assert(gm.n_species == 1, failedUnitTest());
	assert(approxEqual(gd.p, 1.0e5, 1.0e-6), failedUnitTest());
	assert(approxEqual(gd.Ttr, 300.0, 1.0e-6), failedUnitTest());
	assert(approxEqual(gd.massf[0], 1.0, 1.0e-6), failedUnitTest());

	gm.update_sound_speed(gd);
	assert(approxEqual(gd.rho, 1.16109, 1.0e-4), failedUnitTest());
	assert(approxEqual(gd.u, 215314.0, 1.0e-4), failedUnitTest());
	assert(approxEqual(gd.a, 347.241, 1.0e-4), failedUnitTest());

	gm.update_trans_coeffs(gd);
	assert(approxEqual(gd.mu, 1.84691e-05, 1.0e-6), failedUnitTest());
	assert(approxEqual(gd.k, 0.0262449, 1.0e-6), failedUnitTest());

	return 0;
    }
}
