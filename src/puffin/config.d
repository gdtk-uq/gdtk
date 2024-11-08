// config.d -- Part of the Puffin and Lorikeet flow calculators.
//
// A place to keep the configuration details.
// PA Jacobs
// 2022-01-21
//
module config;

import std.conv;
import std.stdio;
import std.format;
import std.json;

import util.json_helper;
import nm.schedule;


// To choose a flux calculator.
enum FluxCalcCode {ausmdv=0, hanel=1, riemann=2, ausmdv_plus_hanel=3, riemann_plus_hanel=4};


final class Config {
public:
    // Parameters common to both space- and time-marching codes.
    shared static int verbosity_level = 1;
    // Messages have a hierarchy:
    // 0 : only error messages will be omitted
    // 1 : emit messages that are useful for a long-running job (default)
    // 2 : plus verbose init messages
    // 3 : plus verbose stepping messages
    //
    shared static maxCPUs = 1;
    shared static string job_name = ""; // Change this to suit at run time.
    shared static string title = "";
    static string gas_model_file;
    static string[] iovar_names;
    static string reaction_file_1;
    static string reaction_file_2;
    shared static bool reacting = false;
    shared static double T_frozen;
    shared static bool axisymmetric;
    shared static bool add_user_supplied_source_terms;
    shared static double cfl;
    shared static int max_step;
    shared static int print_count;
    shared static int x_order;
    shared static int t_order;
    shared static FluxCalcCode flux_calc;
    shared static double compression_tol;
    shared static double shear_tol;
    //
    // Parameters specific to space-marching code.
    shared static double max_x;
    shared static double dx;
    shared static int max_step_relax;
    shared static double plot_dx;
    shared static int n_streams;
    //
    // Parameters specific to time-marching code.
    shared static double dt_init;
    shared static int cfl_count = 10;
    shared static double max_t;
    shared static double plot_dt;
    shared static int n_fluid_blocks;
    shared static int nib;
    shared static int njb;
    shared static int[][] blk_ids;
    shared static int[] nics;
    shared static int[] njcs;
}


void parse_config_data_for_transient_solver(JSONValue configData)
{
    Config.title = getJSONstring(configData, "title", "");
    Config.gas_model_file = getJSONstring(configData, "gas_model_file", "");
    Config.iovar_names = getJSONstringarray(configData, "iovar_names", [""]);
    Config.reaction_file_1 = getJSONstring(configData, "reaction_file_1", "");
    Config.reaction_file_2 = getJSONstring(configData, "reaction_file_2", "");
    Config.reacting = getJSONbool(configData, "reacting", false);
    Config.T_frozen = getJSONdouble(configData, "T_frozen", 300.0);
    Config.axisymmetric = getJSONbool(configData, "axisymmetric", false);
    Config.add_user_supplied_source_terms = getJSONbool(configData, "add_user_supplied_source_terms", false);
    Config.max_t = getJSONdouble(configData, "max_time", 0.0);
    Config.max_step = getJSONint(configData, "max_step", 0);
    Config.dt_init = getJSONdouble(configData, "dt_init", 1.0e-6);
    Config.cfl = getJSONdouble(configData, "cfl", 0.5);
    Config.cfl_count = getJSONint(configData, "cfl_count", 0);
    Config.print_count = getJSONint(configData, "print_count", 50);
    Config.plot_dt = getJSONdouble(configData, "plot_dt", 1.0e-2);
    Config.x_order = getJSONint(configData, "x_order", 2);
    Config.t_order = getJSONint(configData, "t_order", 2);
    Config.flux_calc = to!FluxCalcCode(getJSONint(configData, "flux_calc", 0));
    Config.compression_tol = getJSONdouble(configData, "compression_tol", -0.3);
    Config.shear_tol = getJSONdouble(configData, "shear_tol", 0.2);
    Config.n_fluid_blocks = getJSONint(configData, "n_fluid_blocks", 0);
    Config.nib = getJSONint(configData, "nib", 0);
    Config.njb = getJSONint(configData, "njb", 0);
    JSONValue jsonIds = configData["blk_ids"];
    Config.blk_ids.length = Config.nib;
    foreach (i; 0 .. Config.nib) {
        Config.blk_ids[i].length = Config.njb;
        JSONValue jsonRow = jsonIds[i];
        foreach (j; 0 .. Config.njb) { Config.blk_ids[i][j] = to!int(jsonRow[j].integer); }
    }
    int[] nics = getJSONintarray(configData, "nics", [0]);
    foreach (n; nics) { Config.nics ~= n; }
    int[] njcs = getJSONintarray(configData, "njcs", [0]);
    foreach (n; njcs) { Config.njcs ~= n; }
    if (Config.verbosity_level > 1) {
        writeln("Config:");
        writefln("  title= \"%s\"", Config.title);
        writeln("  gas_model_files= ", Config.gas_model_file);
        writeln("  iovar_names= ", Config.iovar_names);
        writeln("  reaction_files_1= ", Config.reaction_file_1);
        writeln("  reaction_files_2= ", Config.reaction_file_2);
        writeln("  reacting= ", Config.reacting);
        writeln("  T_frozen= ", Config.T_frozen);
        writeln("  axisymmetric= ", Config.axisymmetric);
        writeln("  max_time= ", Config.max_t);
        writeln("  max_step= ", Config.max_step);
        writeln("  dt_init= ", Config.dt_init);
        writeln("  cfl= ", Config.cfl);
        writeln("  cfl_count= ", Config.cfl_count);
        writeln("  print_count= ", Config.print_count);
        writeln("  plot_dt= ", Config.plot_dt);
        writeln("  x_order= ", Config.x_order);
        writeln("  t_order= ", Config.t_order);
        writeln("  flux_calc= ", Config.flux_calc);
        writeln("  compression_tol= ", Config.compression_tol);
        writeln("  shear_tol= ", Config.shear_tol);
        writeln("  n_fluid_blocks= ", Config.n_fluid_blocks);
        writeln("  nib= ", Config.nib);
        writeln("  njb= ", Config.njb);
        writeln("  blk_ids= ", Config.blk_ids);
        writeln("  nics= ", Config.nics);
        writeln("  njcs= ", Config.njcs);
    }
} // end parse_config_data_for_transient_solver()


void parse_config_data_for_marching_solver(JSONValue configData)
{
    Config.title = getJSONstring(configData, "title", "");
    Config.gas_model_file = getJSONstring(configData, "gas_model_file", "");
    Config.reaction_file_1 = getJSONstring(configData, "reaction_file_1", "");
    Config.reaction_file_2 = getJSONstring(configData, "reaction_file_2", "");
    Config.reacting = getJSONbool(configData, "reacting", false);
    Config.T_frozen = getJSONdouble(configData, "T_frozen", 300.0);
    Config.axisymmetric = getJSONbool(configData, "axisymmetric", false);
    Config.add_user_supplied_source_terms = getJSONbool(configData, "add_user_supplied_source_terms", false);
    Config.max_x = getJSONdouble(configData, "max_x", 0.0);
    Config.max_step = getJSONint(configData, "max_step", 0);
    Config.dx = getJSONdouble(configData, "dx", 0.0);
    Config.max_step_relax = getJSONint(configData, "max_step_relax", 100);
    Config.cfl = getJSONdouble(configData, "cfl", 0.5);
    Config.print_count = getJSONint(configData, "print_count", 50);
    Config.plot_dx = getJSONdouble(configData, "plot_dx", 1.0e-2);
    Config.x_order = getJSONint(configData, "x_order", 2);
    Config.t_order = getJSONint(configData, "t_order", 2);
    Config.flux_calc = to!FluxCalcCode(getJSONint(configData, "flux_calc", 0));
    Config.compression_tol = getJSONdouble(configData, "compression_tol", -0.3);
    Config.shear_tol = getJSONdouble(configData, "shear_tol", 0.2);
    Config.n_streams = getJSONint(configData, "n_streams", 1);
    if (Config.verbosity_level > 1) {
        writeln("Config:");
        writefln("  title= \"%s\"", Config.title);
        writeln("  gas_model_files= ", Config.gas_model_file);
        writeln("  reaction_files_1= ", Config.reaction_file_1);
        writeln("  reaction_files_2= ", Config.reaction_file_2);
        writeln("  reacting= ", Config.reacting);
        writeln("  T_frozen= ", Config.T_frozen);
        writeln("  axisymmetric= ", Config.axisymmetric);
        writeln("  max_x= ", Config.max_x);
        writeln("  max_step= ", Config.max_step);
        writeln("  dx= ", Config.dx);
        writeln("  max_step_relax= ", Config.max_step_relax);
        writeln("  cfl= ", Config.cfl);
        writeln("  print_count= ", Config.print_count);
        writeln("  plot_dx= ", Config.plot_dx);
        writeln("  x_order= ", Config.x_order);
        writeln("  t_order= ", Config.t_order);
        writeln("  flux_calc= ", Config.flux_calc);
        writeln("  compression_tol= ", Config.compression_tol);
        writeln("  shear_tol= ", Config.shear_tol);
        writeln("  n_streams= ", Config.n_streams);
    }
} // end parse_config_data_for_marching_solver()
