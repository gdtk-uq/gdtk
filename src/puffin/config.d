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

// To choose a flux calculator.
enum FluxCalcCode {ausmdv=0, hanel=1, riemann=2, ausmdv_plus_hanel=3, riemann_plus_hanel=4};


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
    shared static FluxCalcCode flux_calc;
    shared static double compression_tol;
    shared static double shear_tol;
    shared static int n_streams;
}
