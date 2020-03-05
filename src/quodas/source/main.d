/* 
 QUODAS: Quasi-One Dimensional Adjoint Solver
 
 author:   Kyle A. Damm
 date:     July 2017 (refactored February 2020)
 location: Seoul, South Korea
*/

import std.stdio;
import std.string;
import std.file;
import std.conv;
import std.exception;
import config;
import flowsolve;
import optimizer;
import block;

void main() {
    Params params;
    read_config_file(params);
        
    if (params.solver == "flow_solver") {
        Block block;
        if (params.simulation_type == "nozzle") {
            block = cast(Nozzle) new Nozzle(params);
        } else if (params.simulation_type == "shocktube") {
            block = cast(Shocktube) new Shocktube(params);
        } else {
            string msg = text("simulation_type is not set to a valid option in config file. Select either 'shocktube' or 'nozzle'.");
            throw new Error(msg);
        }
        
        Solver solver = new Solver(params);
        solver.integrate_in_time(block);
        block.write_flow_solution_to_file();    
    } else if (params.solver == "flow_optimizer") {
        Nozzle block = new Nozzle(params);
        assert(block.label == "nozzle", "flow_optimizer only compatible with nozzle.");
        Solver solver = new Solver(params);
        Optimizer optimizer = new Optimizer(params, block, solver);
        block.read_flow_solution_from_file(optimizer.target);
        optimizer.optimize();
        solver.integrate_in_time(block);
        block.write_flow_solution_to_file();
    } else {
        string msg = text("solver is not set to a valid option in config file. Select either 'flow_solver' or 'flow_optimizer'.");
        throw new Error(msg);
    }
} // end main
