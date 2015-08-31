import std.stdio;
import std.string;
import std.conv;
import std.math;
import std.algorithm;
import luad.all;
import util.lua_service;

import fileutil;
import readconfig;
import globalconfig;
import globaldata;

void main(string[] args)
{
    string jobName = args[1];
    string luafname = args[2];
    string luafname2 = args[3];
    int tindx = to!int(args[4]);
    writeln("jobName= ", jobName, " luafname= ", luafname, " luafname2= ", luafname2);
    
    GlobalConfig.base_file_name = jobName;
    read_config_file();

    double LInf_norm = 0.0;
    double sum = 0.0;
    double total_volume = 0;

    // First work on solid block
    if ( solidBlocks.length > 0 ) {
	auto blk = solidBlocks[0];
	blk.assembleArrays();
	blk.bindFacesAndVerticesToCells();
	blk.readGrid(make_file_name!"solid-grid"(jobName, blk.id, 0));
	blk.readSolution(make_file_name!"solid"(jobName, blk.id, tindx));

	auto lua = initLuaState(luafname);
	auto lfunc = lua.get!LuaFunction("Ts");


	foreach (cell; blk.activeCells) {
	    double T_ex = T_analytical(lfunc, cell.pos.x, cell.pos.y);
	    double abs_diff = fabs(cell.T[0] - T_ex);
	    double volume = cell.volume;
	    LInf_norm = max(LInf_norm, abs_diff);
	    sum += volume*abs_diff^^2;
	    total_volume += volume;
	}
    }

    // Second, work on gas block
    if ( gasBlocks.length > 0 ) {
	auto gblk = gasBlocks[0];
	gblk.assemble_arrays();
	gblk.bind_interfaces_and_vertices_to_cells();
	gblk.bind_vertices_and_cells_to_interfaces();
	gblk.read_grid(make_file_name!"grid"(jobName, gblk.id, 0), 0);
	auto sim_time = gblk.read_solution(make_file_name!"flow"(jobName, gblk.id, tindx));
	auto lua2 = initLuaState(luafname2);
	auto lfunc2 = lua2.get!LuaFunction("T");
	foreach (cell; gblk.active_cells) {
	    double T_ex = T_analytical(lfunc2, cell.pos[0].x, cell.pos[0].y);
	    double abs_diff = fabs(cell.fs.gas.T[0] - T_ex);
	    double volume = cell.volume[0];
	    LInf_norm = max(LInf_norm, abs_diff);
	    sum += volume*abs_diff^^2;
	    total_volume += volume;
	}
    }

    double L2_norm = sqrt(sum/total_volume);

    writeln(format("L2-norm = %20.12e", L2_norm));
    writeln(format("L-inf-norm = %20.12e", LInf_norm));
}

double T_analytical(LuaFunction lfunc, double x, double y)
{
    LuaObject[] ret = lfunc(x, y);
    return ret[0].to!double();
}
