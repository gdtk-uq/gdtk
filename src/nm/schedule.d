// schedule.d
// Select or interpolate a value from a schedule of values.
//
// PJ, 2020-10-10 extracted from l1d/config.d
//

module nm.schedule;

import std.format;
import std.math;


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


version(schedule_test) {
    import util.msg_service;
    import std.conv;
    int main() {
        auto mysched = new Schedule([0.0, 1.0, 1.5, 2.0], [1.0, 1.0, 0.5, 0.0]);
        assert(approxEqual(mysched.get_value(-1.0), 1.0), failedUnitTest());
        assert(approxEqual(mysched.interpolate_value(-1.2), 1.0), failedUnitTest());
        assert(approxEqual(mysched.get_value(1.2), 1.0), failedUnitTest());
        assert(approxEqual(mysched.interpolate_value(1.2), 0.8), failedUnitTest());
        assert(approxEqual(mysched.get_value(2.2), 0.0), failedUnitTest());
        assert(approxEqual(mysched.interpolate_value(2.2), 0.0), failedUnitTest());
        return 0;
    }
}
