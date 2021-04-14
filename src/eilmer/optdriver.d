/** optdriver.d
 *
 * This is the core driver for performing CFD-based optimisation with Eilmer4.
 * It is currently written to interface with the open-source DAKOTA optimisation software (https://dakota.sandia.gov/).
 *
 * Author: Kyle Damm
 * First code: 16-03-2018
 * Refactor:   13-04-2021
 *
 */

import std.algorithm;
import std.conv;
import std.file;
import std.getopt;
import std.process;
import std.stdio;
import std.string;

void main(string[] args) {
    writeln("Eilmer 4.0 optimisation driver.");
    writeln("Revision: PUT_REVISION_STRING_HERE");
    writeln("Compiler-name: PUT_COMPILER_NAME_HERE");

    string jobName;
    getopt(args, "job", &jobName);

    // some internal settings
    bool eval_objective = false;
    bool eval_gradients = false;
    read_dakota_input_file(eval_objective, eval_gradients);

    // prep the simulation
    string prep_sim_cmd = "e4shared --prep --job=" ~ jobName;
    auto prep_sim_output = executeShell(prep_sim_cmd);
    if (prep_sim_output.status == 1) { writeln(prep_sim_output.output); }

    // perturb grid
    string perturb_grid_cmd = "e4ssc --update-grid --job=" ~ jobName;
    auto perturb_grid_output = executeShell(perturb_grid_cmd);
    if (perturb_grid_output.status == 1) { writeln(perturb_grid_output.output); }

    // execute flow solver
    string run_flow_solver_cmd = "e4-nk-shared --job=" ~ jobName;
    auto run_flow_solver_output = executeShell(run_flow_solver_cmd);
    if (run_flow_solver_output.status == 1) { writeln(run_flow_solver_output.output); }

    if (eval_objective) {
        // evaluate objective function
        string eval_objective_cmd = "e4ssc --return-objective-function=true --job=" ~ jobName;
        auto eval_objective_output = executeShell(eval_objective_cmd);
        if (eval_objective_output.status == 1) { writeln(eval_objective_output.output); }
    }

    if (eval_gradients) {
        // execute shape sensitivity calculator
        string eval_grads_cmd = "e4ssc --job=" ~ jobName;
        auto eval_grads_output = executeShell(eval_grads_cmd);
        if (eval_grads_output.status == 1) { writeln(eval_grads_output.output); }
    }

    // check DAKOTA results file exists
    string dakotaResultFileName = "results.out";
    assert(exists(dakotaResultFileName), "DAKOTA results.out file not found");
}

void read_dakota_input_file(ref bool eval_objective, ref bool eval_gradients) {

    int nfns; // number of functions (objective/constraints)
    int nrsps; // number of expected responses
    int ndvars; // number of design variables
    double[] dvars; // design variable values

    string dakotaInFileName = "params.in";
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
        if ( nrsps == 1) { eval_objective = true; eval_gradients = false; }
        else if ( nrsps == 2) { eval_objective = false; eval_gradients = true; }
        else if ( nrsps == 3) { eval_objective = true; eval_gradients = true; }
        else assert(0, "No requested response data in DAKOTA input file");
    }

    return;
} // end read_dakota_input_file
