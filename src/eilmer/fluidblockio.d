module fluidblockio;

import std.conv;


import fluidblockio_old;
import fluidblockio_new;

import fluidblock;
import sfluidblock;
import ufluidblock;
import fvcell;
import globalconfig;
import geom;
import geom.luawrap;
import flowstate;
import luaflowstate;
import util.lua;
import util.lua_service;
import nm.number;
import flowsolution;

// These functions provide switching between the new and legacy data formats


bool is_legacy_format(string flow_format)
{
    switch(flow_format) {
        case "gziptext" : return true;
        case "rawbinary" : return true;
        case "eilmer4text" : return false;
        case "eilmer4binary" : return false;
        default : throw new Error("How did we get here?");
    }
}

string flow_format_ext(string flow_format)
{
    switch(flow_format) {
        case "gziptext" : return "gz";
        case "rawbinary" : return "bin";
        case "eilmer4text" : return "zip";
        case "eilmer4binary" : return "zip";
        default : throw new Error("How did we get here?");
    }
}

double read_solution(FluidBlock blk, string filename, bool overwrite_geometry_data)
{

    if (is_legacy_format(GlobalConfig.flow_format)) {
        auto sblk = cast(SFluidBlock) blk;
        if (sblk) { return read_legacy_solution(sblk, filename, overwrite_geometry_data); }
        auto ublk = cast(UFluidBlock) blk;
        if (ublk) { return read_legacy_solution(ublk, filename, overwrite_geometry_data); }
        throw new Error("Oops, unknown type of FluidBlock");
    } else {
        return read_zip_solution(blk, filename);
    }
}

extern(C) int write_initial_sg_flow_file_from_lua(lua_State* L)
{
    if (is_legacy_format(GlobalConfig.flow_format)) {
        return write_initial_sg_flow_legacy_file_from_lua(L);
    } else {
        return write_initial_sg_flow_zip_file_from_lua(L);
    }
}

extern(C) int write_initial_usg_flow_file_from_lua(lua_State* L)
{
    if (is_legacy_format(GlobalConfig.flow_format)) {
        return write_initial_usg_flow_legacy_file_from_lua(L);
    } else {
        return write_initial_usg_flow_zip_file_from_lua(L);
    }
}

void read_block_data_from_file(FluidBlockLite blk, string filename, Grid_t gridType, string flow_format)
{
    if (is_legacy_format(flow_format)) {
        return read_block_data_from_legacy_file(blk, filename, gridType, flow_format);
    } else {
        return read_block_data_from_zip_file(blk, filename, gridType, flow_format);
    }
}

