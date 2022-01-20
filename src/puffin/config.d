// config.d
// A place to keep the configuration details.

module config;

import std.conv;
import std.format;
import nm.schedule;


final class PuffinConfig {
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
    static string[] gas_model_files;
    static string[] reaction_files_1;
    static string[] reaction_files_2;
    shared static bool reacting = false;
    shared static double T_frozen;
    shared static double max_x;
    shared static int max_step;
    shared static double dx;
    shared static int cfl_count;
    shared static int print_count;
    shared static int plot_count;
    shared static int x_order;
    shared static int nstreams;
}
