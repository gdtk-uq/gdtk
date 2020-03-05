// config.d
module config;
import std.json;
import std.conv;
import std.stdio;
import std.file;
import core.stdc.stdlib : exit;

// define some global values
immutable size_t nprimitive = 3; // number of primitive flow variables
immutable size_t ndesign = 3;    // number of design variables
immutable size_t nghost = 2;     // number of ghost cells attached to boundaries

// data structure to store configuration file parameters
struct Params {
    string simulation_type;
    string solver;
    int ncells;
    double length;
    string flux_calc;
    double gamma;
    int interpolation_order;
    double nozzle_back_pressure_factor;
    double cfl;
    int max_flow_steps;
    double max_flow_time;
    bool adjoint_method;
    double eps;
    int max_opt_iters;
    double opt_tol;
    double b;
    double c;
    double d;
    double yo;
    double scale;
} // end Params

// Utility functions to condense the following code.
//
string getJSONstring(JSONValue jsonData, string key, string defaultValue)
{
    string value;
    try {
        value = to!string(jsonData[key].str);
    } catch (Exception e) {
        value = defaultValue;
    }
    return value;
} // end getJSONstring()

int getJSONint(JSONValue jsonData, string key, int defaultValue)
{
    int value;
    try {
        value = to!int(jsonData[key].integer);
    } catch (Exception e) {
        value = defaultValue;
    }
    return value;
} // end getJSONint()

double getJSONdouble(JSONValue jsonData, string key, double defaultValue)
{
    double value;
    try {
        auto json_val = jsonData[key];
        // We wish to accept value like 0.0 or 0
        if (json_val.type() == JSONType.float_) {
            value = json_val.floating;
        } else {
            value = to!double(json_val.str);
        }
    } catch (Exception e) {
        value = defaultValue;
    }
    return value;
} // end getJSONdouble()

bool getJSONbool(JSONValue jsonData, string key, bool defaultValue)
{
    bool value;
    try {
        value = jsonData[key].type is JSONType.true_;
    } catch (Exception e) {
        value = defaultValue;
    }
    return value;
} // end getJSONbool()
//
string update_string(string key, string field)
{
    return "params."~field~" = getJSONstring(jsonData, \""~key~"\", params."~field~");";
}
string update_bool(string key, string field)
{
    return "params."~field~" = getJSONbool(jsonData, \""~key~"\", params."~field~");";
}
string update_int(string key, string field)
{
    return "params."~field~" = getJSONint(jsonData, \""~key~"\", params."~field~");";
}
string update_double(string key, string field)
{
    return "params."~field~" = getJSONdouble(jsonData, \""~key~"\", params."~field~");";
}
//

void read_config_file(ref Params params) {

    writeln("Read config file.");
    string fileName = "solver.config";
    string content;
    try {
        content = readText(fileName);
    } catch (Exception e) {
        writeln("Failed to read config file: ", fileName);
        writeln("Message is: ", e.msg);
        exit(1);
    }
    JSONValue jsonData;
    try {
        jsonData = parseJSON!string(content);
    } catch (Exception e) {
        writeln("Failed to parse JSON from config file: ", fileName);
        writeln("Message is: ", e.msg);
        exit(1);
    }
    
    mixin(update_string("simulation_type", "simulation_type"));
    mixin(update_string("solver", "solver"));
    mixin(update_int("ncells", "ncells"));
    mixin(update_double("length", "length"));
    mixin(update_string("flux_calc", "flux_calc"));
    mixin(update_double("gamma", "gamma"));
    mixin(update_int("interpolation_order", "interpolation_order"));
    mixin(update_double("nozzle_back_pressure_factor", "nozzle_back_pressure_factor"));
    mixin(update_double("cfl", "cfl"));
    mixin(update_int("max_flow_steps", "max_flow_steps"));
    mixin(update_double("max_flow_time", "max_flow_time"));
    mixin(update_bool("adjoint_method", "adjoint_method"));   
    mixin(update_double("eps", "eps"));
    mixin(update_int("max_opt_iters", "max_opt_iters"));
    mixin(update_double("opt_tol", "opt_tol"));
    mixin(update_double("b", "b"));
    mixin(update_double("c", "c"));
    mixin(update_double("d", "d"));
    mixin(update_double("yo", "yo"));
    mixin(update_double("scale", "scale"));
}

void remove_old_files(string filename) {
    if (exists(filename)) {
        remove(filename);
        writef("WARNING: previous %s file removed.\n", filename);
    }
} // end remove_old_files()
