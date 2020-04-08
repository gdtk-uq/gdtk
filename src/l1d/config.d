// config.d
// A place to keep the configuration details.

module config;

final class L1dConfig {
public:
    shared static int verbosity_level = 1;
    // Messages have a hierarchy:
    // 0 : only error messages will be omitted
    // 1 : emit messages that are useful for a long-running job (default)
    // 2 : plus verbose init messages
    // 3 : plus verbose boundary condition messages
    //
    shared static string job_name = "job"; // Change this to suit at run time.
    shared static string title = "";
    static string[] gas_model_files;
    static string[] reaction_files_1;
    static string[] reaction_files_2;
    shared static bool reacting = false;
    shared static double max_time;
    shared static int max_step;
    shared static double dt_init;
    shared static double cfl_value;
    shared static int x_order;
    shared static int t_order;
    shared static int n_dt_plot;
    shared static double[] t_change;
    shared static double[] dt_plot;
    shared static double[] dt_hist;
    shared static int hloc_n;
    shared static double[] hloc_x;
    shared static int nslugs;
    shared static int npistons;
    shared static int ndiaphragms;
    shared static int necs;
}
