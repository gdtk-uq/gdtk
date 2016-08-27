/**
 * gasmodelutil.d
 * Utility functions that make use of the gasmodel class and its derived classes.
 *
 * Author: Peter J. and Rowan G.
 * Version: 2014-06-22, first cut, exploring the options.
 */

module gas.gas_model_util;

import gas.gas_model;
import gas.ideal_gas;
import gas.therm_perf_gas;
import gas.very_viscous_air;
import gas.co2gas;
import gas.co2gas_sw;
import gas.sf6virial;
import gas.uniform_lut;
import gas.adaptive_lut_CEA;
import std.file;
import std.stdio;
import std.string;
import util.lua;
import util.lua_service;

