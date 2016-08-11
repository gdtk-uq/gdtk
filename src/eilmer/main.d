/** main.d
 * Eilmer4 compressible-flow simulation code, top-level function.
 *
 * Author: Peter J. and Rowan G. 
 * First code: 2015-02-05
 */

import core.memory;
import core.stdc.stdlib : exit;
import std.stdio;
import std.string;
import std.file;
import std.path;
import std.getopt;
import std.conv;
import std.parallelism;
import std.algorithm;

import geom;
import gas;
import gas.luagas_model;
import kinetics.luareaction_mechanism;
import kinetics.luachemistry_update;
import fvcore: FlowSolverException;
import globalconfig;
import simcore;
import util.lua;
import luaglobalconfig;
import luaflowstate;
import luaflowsolution;
import luageom;
import luagpath;
import luasurface;
import luavolume;
import luaunifunction;
import luasgrid;
import luausgrid;
import luasketch;
import luasolidprops;
import postprocess;
import luaflowsolution;

void main(string[] args)
{
    writeln("Eilmer4 compressible-flow simulation code.");
    writeln("Revision: PUT_REVISION_STRING_HERE");

    string msg = "Usage:                               Comment:\n";
    msg       ~= "e4shared [--job=<string>]            file names built from this string\n";
    msg       ~= "         [--verbosity=<int>]         defaults to 0\n";
    msg       ~= "\n";
    msg       ~= "         [--prep]                    prepare config, grid and flow files\n";
    msg       ~= "\n";
    msg       ~= "         [--run]                     run the simulation over time\n";
    msg       ~= "         [--tindx-start=<int>|last|9999]  defaults to 0\n";
    msg       ~= "         [--max-cpus=<int>]          defaults to ";
    msg       ~= to!string(totalCPUs) ~" on this machine\n";
    msg       ~= "         [--max-wall-clock=<int>]    in seconds\n";
    msg       ~= "\n";
    msg       ~= "         [--post]                    post-process simulation data\n";
    msg       ~= "         [--list-info]               report some details of this simulation\n";
    msg       ~= "         [--tindx-plot=<int>|all|last|9999]  default to last\n";
    msg       ~= "         [--add-vars=\"mach,pitot,total-h,total-p\"]\n";
    msg       ~= "         [--ref-soln=<filename>]     Lua file for reference solution\n";
    msg       ~= "         [--vtk-xml]                 produce XML VTK-format plot files\n";
    msg       ~= "         [--binary-format]           use binary within the VTK-XML\n";
    msg       ~= "         [--tecplot]                 write an ASCII file for Tecplot\n";
    msg       ~= "         [--plot-dir=<string>]       defaults to plot\n";
    msg       ~= "         [--output-file=<string>]    defaults to stdout\n";
    msg       ~= "         [--slice-list=\"blk-range,i-range,j-range,k-range;...\"]\n";
    msg       ~= "         [--extract-streamline=\"x,y,z;...\"]        streamline locus points\n";
    msg       ~= "         [--extract-line=\"x0,y0,z0,x1,y1,z1,n;...\"]    sample along a line\n";
    msg       ~= "         [--probe=\"x,y,z;...\"]       locations to sample flow data\n";
    msg       ~= "         [--norms=\"varName,varName,...\"] report L1,L2,Linf norms\n";
    msg       ~= "         [--region=\"x0,y0,z0,x1,y1,z1\"]  limit norms calculation to a box\n";
    msg       ~= "\n";
    msg       ~= "         [--custom-post]             run custom post-processing script\n";
    msg       ~= "         [--script-file=<string>]    defaults to post.lua\n";
    msg       ~= "\n";
    msg       ~= "         [--help]                    writes this message\n";
    if ( args.length < 2 ) {
	writeln("Too few arguments.");
	write(msg);
	exit(1);
    }
    string jobName = "";
    int verbosityLevel = 0;
    bool prepFlag = false;
    bool runFlag = false;
    string tindxStartStr = "0";
    int tindxStart = 0;
    int maxCPUs = totalCPUs;
    int maxWallClock = 5*24*3600; // 5 days default
    bool postFlag = false;
    bool listInfoFlag = false;
    string tindxPlot = "last";
    string addVarsStr = "";
    string luaRefSoln = "";
    bool vtkxmlFlag = false;
    bool binaryFormat = false;
    bool tecplotFlag = false;
    string plotDir = "plot";
    string outputFileName = "";
    string sliceListStr = "";
    string extractStreamStr = "";
    string extractLineStr = "";
    string probeStr = "";
    string normsStr = "";
    string regionStr = "";
    bool customPostFlag = false;
    string scriptFile = "post.lua";
    bool helpWanted = false;
    try {
	getopt(args,
	       "job", &jobName,
	       "verbosity", &verbosityLevel,
	       "prep", &prepFlag,
	       "run", &runFlag,
	       "tindx-start", &tindxStartStr,
	       "max-cpus", &maxCPUs,
	       "max-wall-clock", &maxWallClock,
	       "post", &postFlag,
	       "list-info", &listInfoFlag,
	       "tindx-plot", &tindxPlot,
	       "add-vars", &addVarsStr,
               "ref-soln", &luaRefSoln,
	       "vtk-xml", &vtkxmlFlag,
	       "binary-format", &binaryFormat,
	       "tecplot", &tecplotFlag,
	       "plot-dir", &plotDir,
	       "output-file", &outputFileName,
	       "slice-list", &sliceListStr,
	       "extract-streamline", &extractStreamStr,
	       "extract-line", &extractLineStr,
	       "probe", &probeStr,
	       "norms", &normsStr,
	       "region", &regionStr,
	       "custom-post", &customPostFlag,
	       "script-file", &scriptFile,
	       "help", &helpWanted
	       );
    } catch (Exception e) {
	writeln("Problem parsing command-line options.");
	writeln("Arguments not processed: ");
	args = args[1 .. $]; // Dispose of program name in first argument.
	foreach (myarg; args) writeln("    arg: ", myarg);
	write(msg);
	exit(1);
    }
    if (helpWanted) {
	write(msg);
	exit(0);
    }

    if (prepFlag) {
	writeln("Begin preparation stage for a simulation.");
	if (jobName.length == 0) {
	    writeln("Need to specify a job name.");
	    write(msg);
	    exit(1);
	}
	writeln("Start lua connection.");
	auto L = luaL_newstate();
	luaL_openlibs(L);
	registerVector3(L);
	registerGlobalConfig(L);
	registerFlowSolution(L);
	registerFlowState(L);
	registerPaths(L);
	registerSurfaces(L);
	registerVolumes(L);
	registerUnivariateFunctions(L);
	registerStructuredGrid(L);
	registerUnstructuredGrid(L);
	registerSketch(L);
	registerSolidProps(L);
	registerGasModel(L, LUA_GLOBALSINDEX);
	if ( luaL_dofile(L, toStringz(dirName(thisExePath())~"/prep.lua")) != 0 ) {
	    writeln("There was a problem in the Eilmer Lua code: prep.lua");
	    string errMsg = to!string(lua_tostring(L, -1));
	    throw new FlowSolverException(errMsg);
	}
	if ( luaL_dofile(L, toStringz(jobName~".lua")) != 0 ) {
	    writeln("There was a problem in the user-supplied input lua script: ", jobName~".lua");
	    string errMsg = to!string(lua_tostring(L, -1));
	    throw new FlowSolverException(errMsg);
	}
	if ( luaL_dostring(L, toStringz("build_job_files(\""~jobName~"\")")) != 0 ) {
	    writeln("There was a problem in the Eilmer build function build_job_files() in prep.lua");
	    string errMsg = to!string(lua_tostring(L, -1));
	    throw new FlowSolverException(errMsg);
	}
	writeln("Done preparation.");
    } // end if prepFlag

    if (runFlag) {
	if (jobName.length == 0) {
	    writeln("Need to specify a job name.");
	    write(msg);
	    exit(1);
	}
	GlobalConfig.base_file_name = jobName;
	GlobalConfig.verbosity_level = verbosityLevel;
	maxCPUs = min(max(maxCPUs, 1), totalCPUs); // don't ask for more than available
	switch (tindxStartStr) {
	case "9999":
	case "last":
	    auto times_dict = readTimesFile(jobName);
            auto tindx_list = times_dict.keys;
	    sort(tindx_list);
	    tindxStart = tindx_list[$-1];
	    break;
	default:
	    // We assume that the command-line argument was an integer.
	    tindxStart = to!int(tindxStartStr);
	} // end switch
	if (verbosityLevel > 0) {
	    writeln("Begin simulation with command-line arguments.");
	    writeln("  jobName: ", jobName);
	    writeln("  tindxStart: ", tindxStart);
	    writeln("  maxWallClock: ", maxWallClock);
	    writeln("  verbosityLevel: ", verbosityLevel);
	    writeln("  maxCPUs: ", maxCPUs);
	}
	
	init_simulation(tindxStart, maxCPUs, maxWallClock);
	writeln("starting simulation time= ", simcore.sim_time);
	if (GlobalConfig.block_marching) {
	    march_over_blocks();
	} else {
	    integrate_in_time(GlobalConfig.max_time);
	}
	finalize_simulation();
	writeln("Done simulation.");
    } // end if runFlag

    if (postFlag) {
	if (jobName.length == 0) {
	    writeln("Need to specify a job name.");
	    write(msg);
	    exit(1);
	}
	GlobalConfig.base_file_name = jobName;
	GlobalConfig.verbosity_level = verbosityLevel;
	if (verbosityLevel > 0) {
	    writeln("Begin post-processing with command-line arguments.");
	    writeln("  jobName: ", jobName);
	    writeln("  listInfoFlag: ", listInfoFlag);
	    writeln("  tindxPlot: ", tindxPlot);
	    writeln("  addVarsStr: ", addVarsStr);
	    writeln("  luaRefSoln: ", luaRefSoln);
	    writeln("  vtkxmlFlag: ", vtkxmlFlag);
	    writeln("  binaryFormat: ", binaryFormat);
	    writeln("  tecplotFlag: ", tecplotFlag);
	    writeln("  plotDir: ", plotDir);
	    writeln("  outputFileName: ", outputFileName);
	    writeln("  sliceListStr: ", sliceListStr);
	    writeln("  extractStreamStr: ", extractStreamStr);
	    writeln("  extractLineStr: ", extractLineStr);
	    writeln("  probeStr: ", probeStr);
	    writeln("  normsStr: ", normsStr);
	    writeln("  regionStr: ", regionStr);
	    writeln("  verbosityLevel: ", verbosityLevel);
	}
	post_process(plotDir, listInfoFlag, tindxPlot,
		     addVarsStr, luaRefSoln,
		     vtkxmlFlag, binaryFormat, tecplotFlag,
		     outputFileName, sliceListStr, extractStreamStr,
		     extractLineStr, probeStr,
		     normsStr, regionStr);
	writeln("Done postprocessing.");
    } // end if postFlag

    if (customPostFlag) {
	writeln("Begin custom post-processing using user-supplied script.");
	// For this case, there is very little job context loaded and
	// after loading all of the libraries, we pretty much hand over
	// to a Lua file to do everything.
	writeln("Start lua connection.");
	auto L = luaL_newstate();
	luaL_openlibs(L);
	registerVector3(L);
	registerGlobalConfig(L);
	registerFlowSolution(L);
	registerFlowState(L);
	registerPaths(L);
	registerSurfaces(L);
	registerVolumes(L);
	registerUnivariateFunctions(L);
	registerStructuredGrid(L);
	registerUnstructuredGrid(L);
	registerSketch(L);
	registerSolidProps(L);
	registerGasModel(L, LUA_GLOBALSINDEX);
	registerReactionMechanism(L, LUA_GLOBALSINDEX);
	registerChemistryUpdate(L, LUA_GLOBALSINDEX);
	if ( luaL_dofile(L, toStringz(scriptFile)) != 0 ) {
	    writeln("There was a problem in the user-supplied Lua script: ", scriptFile);
	    string errMsg = to!string(lua_tostring(L, -1));
	    throw new FlowSolverException(errMsg);
	}
	writeln("Done custom postprocessing.");
    } // end if customPostFlag
} // end main()


