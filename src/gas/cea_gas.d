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
import ntypes.complex;
import nm.number;

import gas.gas_model;
import gas.gas_state;
import gas.physical_constants;


struct CEASavedData {
    // This holds some parameters of the gas state that are scanned
    // from the CEA output file and may be needed at a later time.
    // We will put them into the GasState object only when the
    // GasModel is a CEAGas.
public:
    double p; // Pa
    double rho; // kg/m**3
    double u; // J/kg
    double h; // J/kg
    double T; // K
    double a; // m/s
    double Mmass; // average molecular mass, kg/mole
    double Rgas; // J/kg/K
    double gamma; // ratio of specific heats
    double Cp; // J/kg/K
    double s;  // J/kg/K
    double[string] massf;
    double k; // thermal conductivity, W/m/K
    double mu; // viscosity, Pa.s
    string previous_problemType = ""; // "pT", "rhoT", etc

    this(const CEASavedData other) {
        this.copy_values_from(other);
    }

    void copy_values_from(const CEASavedData other)
    {
        this.p = other.p;
        this.rho = other.rho;
        this.u = other.u;
        this.h = other.h;
        this.T = other.T;
        this.a = other.a;
        this.Mmass = other.Mmass;
        this.Rgas = other.Rgas;
        this.gamma = other.gamma;
        this.Cp = other.Cp;
        this.s = other.s;
        foreach (key; other.massf.byKey()) { this.massf[key] = other.massf[key]; }
        this.mu = other.mu;
        this.k = other.k;
        this.previous_problemType = other.previous_problemType;
    }
} // end CEASavedData


class CEAGas: GasModel {
public:

    this(string mixtureName, string[] speciesList, double[string] reactants,
         string inputUnits, double trace, bool withIons)
    {
        type_str = "CEAGas";
        // In the context of Rowan's perfect gas mix, the CEA gas is a strange beast.
        // We will hide it's internal species behind a single pseudo-species that
        // we will call by the mixtureName.  It will not have a fixed molecular mass.
        _n_modes = 0;
        _n_species = 1;
        _species_names.length =  _n_species;
        _species_names[0] = mixtureName;
        _mol_masses.length = 1;
        _mol_masses[0] = 0.0; // dummy value; we shouldn't use it
        _cea_species_names = speciesList.dup();
        _reactants = reactants.dup();
        _inputUnits = inputUnits;
        _trace = trace;
        _withIons = withIons;
        //
        create_species_reverse_lookup();
        //
        // Make sure that all of the pieces are in place to run CEA calculations.
        string dgdinst = expandTilde(environment.get("DGD", "~/gdtkinst"));
        _cea_exe_path = buildPath(dgdinst, "bin", "cea2");
        _cea_cases_path = buildPath(dgdinst, "share", "cea-cases");
        // writeln("_cea_exe_path=", _cea_exe_path); // DEBUG
        if (!exists(_cea_exe_path)) {
            throw new Exception("Cannot find cea2 exe file.");
        }
        if (!exists("thermo.lib") || getSize("thermo.lib") == 0) {
            string ceaThermoFile = buildPath(_cea_cases_path, "thermo.inp");
            if (!exists(ceaThermoFile)) {
                throw new Exception("Cannot find cea2 thermo.inp file.");
            }
            auto copy_thermo = executeShell("cp " ~ ceaThermoFile ~ " .");
            runCEAProgram("thermo", false);
        }
        if (!exists("trans.lib") || getSize("trans.lib") == 0) {
            string ceaTransFile = buildPath(_cea_cases_path, "trans.inp");
            if (!exists(ceaTransFile)) {
                throw new Exception("Cannot find cea2 trans.inp file.");
            }
            auto copy_trans = executeShell("cp " ~ ceaTransFile ~ " .");
            runCEAProgram("trans", false);
        }
        //
        version(complex_numbers) {
            throw new Error("Do not use with complex numbers.");
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
        repr ~= "mixtureName=\"" ~ _species_names[0] ~"\"";
        repr ~= ", speciesList=[";
        foreach (i, sname; _cea_species_names) {
            repr ~= format("\"%s\"", sname);
            repr ~= (i+1 < _cea_species_names.length) ? ", " : "]";
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

    override void update_thermo_from_pT(ref GasState Q) const
    {
        debug {
            callCEA(Q, 0.0, 0.0, "pT", false);
        } else {
            throw new Exception("not implemented for @nogc. PJ 2018-09-23");
        }
    }
    override void update_thermo_from_rhou(ref GasState Q) const
    {
        debug {
            callCEA(Q, 0.0, 0.0, "rhoe", false);
        } else {
            throw new Exception("not implemented for @nogc. PJ 2018-09-23");
        }
    }
    override void update_thermo_from_rhoT(ref GasState Q) const
    {
        debug {
            callCEA(Q, 0.0, 0.0, "rhoT", false);
        } else {
            throw new Exception("not implemented for @nogc. PJ 2018-09-23");
        }
    }
    override void update_thermo_from_rhop(ref GasState Q) const
    {
        throw new Exception("CEAGas update_thermo_from_rhop not implemented.");
    }

    override void update_thermo_from_ps(ref GasState Q, number s) const
    {
        debug {
            callCEA(Q, 0.0, s.re, "ps", false);
        } else {
            throw new Exception("not implemented for @nogc. PJ 2018-09-23");
        }
    }
    override void update_thermo_from_hs(ref GasState Q, number h, number s) const
    {
        throw new Exception("CEAGas update_thermo_from_hs not implemented.");
    }
    override void update_sound_speed(ref GasState Q) const
    {
        // We presume that the thermo update has just (previously) been done.
        // It's not really a separate operation since all thermo updates will
        // also get a new estimate of sound speed, but we do it again.
        if (Q.ceaSavedData.previous_problemType.length == 0) {
            throw new Exception("Seems that we have not called CEA before.");
        }
        debug {
            string problemType = Q.ceaSavedData.previous_problemType;
            if (problemType == "ps") { problemType = "pT"; }
            callCEA(Q, 0.0, 0.0, problemType, false);
        } else {
            throw new Exception("not implemented for @nogc. PJ 2018-09-23");
        }
    }
    override void update_trans_coeffs(ref GasState Q)
    {
        // We presume that the thermo update has just (previously) been done.
        // Repeat it but do the transport coefficients this time.
        if (Q.ceaSavedData.previous_problemType.length == 0) {
            throw new Exception("Seems that we have not called CEA before.");
        }
        debug {
            string problemType = Q.ceaSavedData.previous_problemType;
            if (problemType == "ps") { problemType = "pT"; }
            callCEA(Q, 0.0, 0.0, problemType, true);
        } else {
            throw new Exception("not implemented for @nogc. PJ 2018-09-23");
        }
    }

    // The following functions return the most-recently-computed values.
    // It is assumed that you have previously called an update function
    // so that the values sitting in the CEAGas object are relevant.

    override number dudT_const_v(in GasState Q) const
    {
        return to!number(Q.ceaSavedData.Cp - Q.ceaSavedData.Rgas);
    }
    override number dhdT_const_p(in GasState Q) const
    {
        return to!number(Q.ceaSavedData.Cp);
    }
    override number dpdrho_const_T(in GasState Q) const
    {
        return Q.ceaSavedData.Rgas * Q.T;
    }
    override number gas_constant(in GasState Q) const
    {
        return to!number(Q.ceaSavedData.Rgas);
    }
    override number internal_energy(in GasState Q) const
    {
        return Q.u;
    }
    override number enthalpy(in GasState Q) const
    {
        return to!number(Q.ceaSavedData.h);
    }
    override number entropy(in GasState Q) const
    {
        return to!number(Q.ceaSavedData.s);
    }

private:
    string[] _cea_species_names;
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
        if (exists(outFile)) { core.stdc.stdio.remove(outFile.toStringz); }
        if (exists(pltFile)) { core.stdc.stdio.remove(pltFile.toStringz); }
	int returnCode = 0;
	try {
            // debug { import std.stdio; writeln("DEBUG Z2: About to create pipeProcess."); }
	    // Putting a pause here improves the chances of the call working in a Ruby environment
	    // but we still get an eventual failure with the message:
	    // Could not read from pipe to get child status (Interrupted system call).
	    // The problem does not occur when being called in a Python environment.
	    // import core.thread.osthread; import core.time; Thread.sleep(dur!("msecs")(100));
            auto pp = pipeProcess(_cea_exe_path, Redirect.all);
            // debug { import std.stdio; writeln("DEBUG Z3: About to send jobName."); }
            pp.stdin.writeln(jobName);
            pp.stdin.flush();
            pp.stdin.close();
            // debug { import std.stdio; writeln("DEBUG Z4: About to wait for process."); }
            returnCode = wait(pp.pid);
            // debug { import std.stdio; writeln("DEBUG Z5: After waiting for process."); }
	} catch (Exception e) {
            throw new Exception(format("Exception while trying to run CEA2 as a pipeProcess: %s", e.msg));
	}
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

    double ceaFloat(char[][] tokens) const
    // Clean up the CEA2 short-hand notation for exponential format.
    // CEA2 seems to write exponential-format numbers in a number of ways:
    // -1.023
    //  1.023-2
    //  1.023+2
    //  1.023 2
    {
        char[] valueStr;
        char signChr;
        switch (tokens.length) {
        case 0:
            valueStr = to!(char[])("0.0");
            break;
        case 1:
            // Preserve the sign.
            if (tokens[0][0] == '+' || tokens[0][0] == '-') {
                signChr = tokens[0][0];
                valueStr = tokens[0][1..$];
            } else {
                signChr = 0;
                valueStr = tokens[0][];
            }
            if (canFind(valueStr, "****")) {
                // Occasionally, we get dodgy strings like ******e-3 because
                // CEA2 seems to write them for values like 0.009998.
                // Assume the we shoul replace the stars with 1.0.
                char[] exponent = find(valueStr, "e");
                valueStr = to!(char[])("1.0") ~ exponent;
            }
            // Fix the exponent notation, if necessary.
            if (canFind(valueStr, "-") && !canFind(valueStr, "e")) {
                valueStr = valueStr.replace("-", "e-");
            }
            if (canFind(valueStr, "+") && !canFind(valueStr, "e")) {
                valueStr = valueStr.replace("+", "e+");
            }
            // Restore the sign.
            if (signChr) { valueStr = signChr ~ valueStr; }
            break;
        case 2:
            valueStr = tokens[0] ~ "e+" ~ tokens[1];
            break;
        default:
            throw new Exception("CEAGas: ceaFloat received too many tokens.");
        }
        return to!double(valueStr);
    } // end ceaFloat()

    void callCEA(ref GasState Q, double h, double s,
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
            if (!(Q.p > 0.0 && Q.T > 0.0)) { throw new Exception("CEAGas: Invalid pT"); }
            writer.put(format("   p(bar)      %e\n", Q.p / 1.0e5));
            writer.put(format("   t(k)        %e\n", Q.T));
            break;
        case "rhoe":
            writer.put("problem case=CEAGas vu");
            if (_withIons) { writer.put(" ions"); }
            writer.put("\n");
            if (!(Q.rho > 0.0)) { throw new Exception("CEAGas: Invalid rho"); }
            writer.put(format("   rho,kg/m**3 %e\n", Q.rho));
            // R_universal is in J/mol/K and u from flow solver is in J/kg
            writer.put(format("   u/r         %e\n", Q.u/R_universal/1000.0));
            break;
        case "rhoT":
            writer.put("problem case=CEAGas tv");
            if (_withIons) { writer.put(" ions"); }
            writer.put("\n");
            if (!(Q.rho > 0.0 && Q.T > 0.0)) { throw new Exception("CEAGas: Invalid rhoT"); }
            writer.put(format("   rho,kg/m**3 %e\n", Q.rho));
            writer.put(format("   t(k)        %e\n", Q.T));
            break;
        case "ps":
            writer.put("problem case=CEAGas ps");
            if (_withIons) { writer.put(" ions"); }
            writer.put("\n");
            if (!(Q.p > 0.0)) { throw new Exception("CEAGas: Invalid p"); }
            writer.put(format("   p(bar)      %e\n", Q.p / 1.0e5));
            // R_universal is in J/mol/K and s from flow solver is in J/kg/K
            writer.put(format("   s/r         %e\n", s/R_universal/1000.0));
            break;
        default:
            throw new Exception("Unknown problemType for CEA.");
        }
        // Select the gas components for CEA.
        // Note that the speciesList contains all of the reactants,
        // and that the constructor would have made a reactants entry
        // for each species.
        writer.put("react\n");
        foreach (name; _cea_species_names) {
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
        foreach (name; _cea_species_names) { writer.put(format(" %s", name)); }
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
        // Scan the output file for all of the data that we need.
        auto lines = File("tmp.out", "r").byLine();
        bool thermo_props_found = false;
        bool conductivity_found = false;
        foreach (line; lines) {
            if (line.length == 0) continue;
            if (line.canFind("PRODUCTS WHICH WERE CONSIDERED BUT WHOSE")) break;
            if (line.canFind("THERMODYNAMIC EQUILIBRIUM PROPERTIES AT ASSIGNED") ||
                line.canFind("THERMODYNAMIC EQUILIBRIUM COMBUSTION PROPERTIES AT ASSIGNED")) {
                thermo_props_found = true;
            } else if (thermo_props_found) {
                char[][] tokens = line.split();
                // Scan for the thermodynamic properties.
                if (line.canFind("H, KJ/KG")) { Q.ceaSavedData.h = ceaFloat(tokens[2..$])*1.0e3; }
                if (line.canFind("U, KJ/KG")) { Q.ceaSavedData.u = ceaFloat(tokens[2..$])*1.0e3; }
                if (line.canFind("S, KJ/(KG)(K)")) { Q.ceaSavedData.s = ceaFloat(tokens[2..$])*1.0e3; }
                if (line.canFind("Cp, KJ/(KG)(K)")) { Q.ceaSavedData.Cp = ceaFloat(tokens[2..$])*1.0e3; }
                if (line.canFind("GAMMAs")) { Q.ceaSavedData.gamma = ceaFloat(tokens[1..$]); }
                if (line.canFind("M, (1/n)")) {
                    Q.ceaSavedData.Mmass = ceaFloat(tokens[2..$])/1.0e3; // kg/mole
                    Q.ceaSavedData.Rgas = R_universal / Q.ceaSavedData.Mmass; // kJ/kg/K
                }
                if (line.canFind("SON VEL,M/SEC")) { Q.ceaSavedData.a = ceaFloat(tokens[2..$]); }
                if (line.canFind("P, BAR")) { Q.ceaSavedData.p = ceaFloat(tokens[2..$])*1.0e5; }
                if (line.canFind("T, K")) { Q.ceaSavedData.T = ceaFloat(tokens[2..$]); }
                if (line.canFind("RHO, KG/CU M")) { Q.ceaSavedData.rho = ceaFloat(tokens[3..$]); }
                //
                if (transProps) {
                    // Scan for transport properties.
                    if (line.canFind("VISC,MILLIPOISE")) {
                        Q.ceaSavedData.mu = ceaFloat(tokens[1..$])*1.0e-4;
                    } else if (!conductivity_found && line.canFind("CONDUCTIVITY") && tokens.length == 2) {
                        Q.ceaSavedData.k = ceaFloat(tokens[1..$])*1.0e-1;
                        conductivity_found = true;
                    }
                } else {
                    Q.ceaSavedData.mu = 0.0;
                    Q.ceaSavedData.k = 0.0;
                }
            } // end thermo_props_found
        } // end foreach line
        //
        // Scan the output file again, looking for the mass fractions of species.
        // Trace species will not have been listed by CEA but the caller
        // might like to assume they are available. As such, we'll loop
        // over the complete _cea_species_names and set the mass fractions to zero.
        foreach (sname; _cea_species_names) { Q.ceaSavedData.massf[sname] = 0.0; }
        lines = File("tmp.out", "r").byLine();
        bool species_fractions_found = false;
        foreach (line; lines) {
            if (line.length == 0) continue;
            if (line.strip().startsWith("MASS FRACTIONS")) {
                species_fractions_found = true;
                continue;
            }
            if (line.canFind("* THERMODYNAMIC PROPERTIES FITTED")) { break; }
            if (species_fractions_found) {
                char[][] tokens = line.split();
                string speciesName = to!string(tokens[0]).replace("*", "");
                Q.ceaSavedData.massf[speciesName] = ceaFloat(tokens[1..$]);
            }
        } // end foreach line
        Q.ceaSavedData.previous_problemType = problemType;
        //
        // Put the relevant pieces of the scanned data into the GasState object.
        switch (problemType) {
        case "pT":
            Q.rho = Q.ceaSavedData.rho;
            Q.u = Q.ceaSavedData.u;
            Q.a = Q.ceaSavedData.a;
            break;
        case "rhoe":
            Q.T = Q.ceaSavedData.T;
            Q.p = Q.ceaSavedData.p;
            Q.a = Q.ceaSavedData.a;
            break;
        case "rhoT":
            Q.p = Q.ceaSavedData.p;
            Q.u = Q.ceaSavedData.u;
            Q.a = Q.ceaSavedData.a;
            break;
        case "ps":
            Q.T = Q.ceaSavedData.T;
            Q.rho = Q.ceaSavedData.rho;
            Q.u = Q.ceaSavedData.u;
            Q.a = Q.ceaSavedData.a;
            break;
        default:
            throw new Exception("Unknown problemType for CEA.");
        }
        Q.mu = Q.ceaSavedData.mu;
        Q.k = Q.ceaSavedData.k;
    } // end callCEA()
} // end class CEAGgas

version(cea_gas_test) {
    import std.stdio;
    import util.msg_service;

    int main() {
        lua_State* L = init_lua_State();
        doLuaFile(L, "sample-data/cea-air5species-gas-model.lua");
        auto gm = new CEAGas(L);
        lua_close(L);

        auto gd = GasState(1, 0, true);
        gd.p = 1.0e5;
        gd.T = 300.0;
        gd.massf[0] = 1.0;
        gm.update_thermo_from_pT(gd);
        assert(isClose(gm.R(gd), 288.198, 0.01), failedUnitTest());
        assert(gm.n_modes == 0, failedUnitTest());
        assert(gm.n_species == 1, failedUnitTest());
        assert(isClose(gd.p, 1.0e5, 1.0e-6), failedUnitTest());
        assert(isClose(gd.T, 300.0, 1.0e-6), failedUnitTest());
        assert(isClose(gd.massf[0], 1.0, 1.0e-6), failedUnitTest());
        assert(isClose(gd.rho, 1.1566, 1.0e-4), failedUnitTest());
        assert(isClose(gd.u, -84587.0, 0.1), failedUnitTest());
        assert(isClose(gd.ceaSavedData.massf["N2"], 0.76708, 1.0e-5), failedUnitTest());
        assert(isClose(gd.ceaSavedData.massf["O2"], 0.23292, 1.0e-5), failedUnitTest());

        gm.update_sound_speed(gd);
        assert(isClose(gd.a, 347.7, 0.1), failedUnitTest());

        gm.update_trans_coeffs(gd);
        assert(isClose(gd.mu, 1.87e-05, 0.01), failedUnitTest());
        assert(isClose(gd.k, 0.02647, 1.0e-5), failedUnitTest());

        gm.update_thermo_from_ps(gd, gd.ceaSavedData.s);
        assert(isClose(gd.p, 1.0e5, 1.0), failedUnitTest());
        assert(isClose(gd.T, 300.0, 0.1), failedUnitTest());

        gm.update_thermo_from_rhou(gd);
        assert(isClose(gd.p, 1.0e5, 1.0), failedUnitTest());
        assert(isClose(gd.T, 300.0, 0.1), failedUnitTest());

        return 0;
    }
}
