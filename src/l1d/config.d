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
    static Schedule cfl_schedule;
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
    shared static int nvalves;
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

    @nogc
    double get_value(double t)
    {
        // Select one of our tabulated schedule of values.
        int i = cast(int)t_change.length - 1;
        while ((i > 0) && (t < t_change[i])) { i--; }
        return values[i];
    }

    @nogc
    double interpolate_value(double t)
    {
        // Attempt an interpolation of the tabulated schedule of values.
        if (t <= t_change[0]) { return values[0]; }
        if (t >= t_change[$-1]) { return values[$-1]; }
        // If we get to this point, we must have at least 2 values in our schedule
        // and we can interpolate between a pair of them.
        int i = cast(int)t_change.length - 1;
        while ((i > 0) && (t < t_change[i])) { i--; }
        double frac = (t-t_change[i])/(t_change[i+1]-t_change[i]);
        return (1.0-frac)*values[i] + frac*values[i+1];
    }

private:
    shared double[] t_change;
    shared double[] values;
} // end class Schedule
