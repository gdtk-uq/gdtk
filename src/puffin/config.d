// config.d -- Part of the Puffin steady-flow calculator.
//
// A place to keep the configuration details.
// PA Jacobs
// 2022-01-21
//
module config;

import std.conv;
import std.format;
import nm.schedule;

enum FluxCalcCode {hanel=0, riemann=1}; // To choose a flux calculator.


final class Config {
public:
    shared static int verbosity_level = 1;
    // Messages have a hierarchy:
    // 0 : only error messages will be omitted
    // 1 : emit messages that are useful for a long-running job (default)
    // 2 : plus verbose init messages
    // 3 : plus verbose stepping messages
    //
    shared static string job_name = "job"; // Change this to suit at run time.
    shared static string title = "";
    static string gas_model_file;
    static string reaction_file_1;
    static string reaction_file_2;
    shared static bool reacting = false;
    shared static double T_frozen;
    shared static bool axisymmetric;
    shared static double max_x;
    shared static int max_step;
    shared static double dx;
    shared static int max_step_relax;
    shared static double cfl;
    shared static int print_count;
    shared static double plot_dx;
    shared static int x_order;
    shared static int t_order;
    shared static FluxCalcCode flux_calc = FluxCalcCode.hanel;
    shared static int n_streams;
}
