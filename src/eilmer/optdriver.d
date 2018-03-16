/** optdriver.d
 * 
 * Eilmer4 optimisation driver code.
 *
 * This is the core driver for performing CFD-based optimisation with Eilmer4. 
 *
 * It is written to interface with the open-source DAKOTA optimisation software,
 * found here, https://dakota.sandia.gov/
 *
 * Author: Kyle D.
 * Date: 2018-03-16
 *
**/

import core.stdc.stdlib : exit;
import core.memory;
import std.stdio;
import std.file;
import std.format;
import std.conv;
import std.parallelism;
import std.algorithm;
import std.getopt;
import std.string;
import std.array;
import std.math;
import std.datetime;
import std.process;

void main(string[] args) {
    writeln("Eilmer optimisation driver:");
    writeln("Revision: PUT_REVISION_STRING_HERE");

    // some internal settings
    bool computeObjFns = false;
    bool computeObjGrads = false;
    string runFlowSolverCmd = "run_flow_solver.sh";
    string runGradSolverCmd = "run_shape_sensitivity_calc.sh";
    int nfns; // number of functions (objective/constraints)
    int nrsps; // number of expected responses
    int ndvars; // number of design variables
    double[] fneval; // function evaluation values
    double[] dvars; // design variable values
    double[] gvars; // gradients of design variables
    string dakotaInFileName = "params.in";
    string dakotaResultFileName = "results.out";
    
    // read DAKOTA input file
    assert(exists(dakotaInFileName), "DAKOTA input file not present");
    auto f = File(dakotaInFileName, "r");
    // first line holds ndvars
    {
        auto lineContent = f.readln().strip();
        auto tokens = lineContent.split();
        ndvars = to!int(tokens[0]);
    }
    // the following ndvar lines holds the design variable values
    foreach ( i; 0..ndvars) {
        auto lineContent = f.readln().strip();
        auto tokens = lineContent.split();
        dvars ~= to!double(tokens[0]);
    }
    // the following line holds nfns
    {
        auto lineContent = f.readln().strip();
        auto tokens = lineContent.split();
        nfns = to!int(tokens[0]);
    }
    // the following nfns lines holds the number of expected response data for each fns
    // NB. for now I assume that each variable will asks for the same number of response data (so just read first variable)
    // also, I assume this value will either be:
    //          1          (return the ojective function evaluation), or
    //          nvars      (return objective function gradients), or
    //          1 + nvars  (return objective function evaluation, and gradients)
    {
        auto lineContent = f.readln().strip();
        auto tokens = lineContent.split();
        nrsps = to!int(tokens[0]);
        if ( nrsps == 1) { computeObjFns = true; computeObjGrads = false; }
        else if ( nrsps == ndvars) { computeObjFns = false; computeObjGrads = true; }
        else if ( nrsps == 1+ndvars) { computeObjFns = true; computeObjGrads = true; }
        else assert(0, "No requested response data in DAKOTA input file");
    }
    // finished reading DAKOTA input file
    
    // execute flow solver if objective function evaluations are needed
    if (computeObjFns) {
        assert(exists(runFlowSolverCmd), "e4sss execution bash file not present");
        string cmd = "bash " ~ runFlowSolverCmd;
        auto output = executeShell(cmd);
    }
    // execute shape sensitivity calculator if gradients are needed
    if (computeObjGrads) {
        assert(exists(runGradSolverCmd), "e4ssc execution bash file not present");
        string cmd = "bash " ~ runGradSolverCmd;
        auto output = executeShell(cmd);
    }

    // check DAKOTA results file exists
    assert(exists(dakotaResultFileName), "DAKOTA results.out file not found");
}
