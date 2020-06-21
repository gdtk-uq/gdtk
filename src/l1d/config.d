// config.d
// A place to keep the configuration details.

module config;

import std.conv;
import std.format;


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
    shared static double T_frozen;
    shared static double max_time;
    shared static int max_step;
    shared static double dt_init;
    shared static double cfl_value;
    shared static int cfl_count;
    shared static int print_count;
    shared static int x_order;
    shared static int t_order;
    static Schedule dt_plot;
    static Schedule dt_hist;
    shared static int hloc_n;
    shared static double[] hloc_x;
    shared static int nslugs;
    shared static int npistons;
    shared static int ndiaphragms;
    shared static int necs;
}

class Schedule {
public:
    this(double[] t_change, double[] values)
    {
        auto n = t_change.length;
        assert(n > 0, "Need at least one value in the array.");
        assert(n == values.length, "Inconsistent array lengths.");
        this.t_change.length = n; this.t_change[] = t_change[];
        this.values.length = n; this.values[] = values[];
    }

    override string toString()
    {
        return format("Schedule(t_change=%s, values=%s)", t_change, values);
    }

    double get_value(double t)
    {
        int i = to!int(t_change.length) - 1;
        while ((i > 0) && (t < t_change[i])) { i--; }
        return values[i];
    }

private:
    shared double[] t_change;
    shared double[] values;
} // end class Schedule
