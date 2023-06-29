/**
 * Classes to handle block-level I/O.
 *
 * Author: Rowan J. Gollan
 * Date: 2023-06-21
 */

module blockio;

immutable string BLK_IO_VERSION = "1.0";

import std.stdio;
import std.conv : to;
import std.format : format;

import gzip : GzipOut;

import util.lua;
import nm.number;
import geom.luawrap;

import lmrexceptions : LmrException;
import globalconfig : GlobalConfig;
import lmrconfig;
import fvcell : FVCell;
import fvcellio : FVCellIO, buildFlowVariables;
import fluidblock : FluidBlock;
import sfluidblock : SFluidBlock;
import ufluidblock : UFluidBlock;

BlockIO blkIO;

class BlockIO {
public:

    this()
    {
	mVariables = buildFlowVariables();
	mCIO = new FVCellIO(mVariables);
    }
    /**
     * Writes metadata for field data files as YAML
     *
     * Authors: Rowan J. Gollan
     */
    final
    void writeMetadataToFile(string fname)
    {
	auto f = File(fname, "w");
	f.writefln("---");
	f.writefln("  version: %s", BLK_IO_VERSION);
	f.writeln("  revision-id: PUT_REVISION_STRING_HERE");
	f.writefln("  variables:");
	foreach (var; mVariables) {
	    f.writefln("   - %s", var);
	}
	f.close();
    }
    
    abstract
    void writeVariablesToFile(string fname, FVCell[] cells);

private:
    string[] mVariables;
    FVCellIO mCIO;
}

class BinaryBlockIO : BlockIO {
public:

    this() { super(); }

    override
    void writeVariablesToFile(string fname, FVCell[] cells)
    {
	double[1] dbl1;
	File outfile = File(fname, "wb");
	foreach (var; mVariables) {
	    foreach (cell; cells) {
		dbl1[0] = mCIO[cell,var];
		outfile.rawWrite(dbl1);
	    }
	}
	outfile.close();
    }
}

class GzipBlockIO : BlockIO {
public:

    this() {
	super();
	mDblVarFmt = lmrCfg.dblVarFmt;
    }

    override
    void writeVariablesToFile(string fname, FVCell[] cells)
    {
	auto outfile = new GzipOut(fname);
	foreach (var; mVariables) {
	    foreach (cell; cells) {
		outfile.compress(format(mDblVarFmt ~ "\n", mCIO[cell,var]));
	    }
	}
	outfile.finish();
    }

private:
    string mDblVarFmt;
    
}

extern(C) int luafn_writeFlowMetadata(lua_State *L)
{
    alias cfg = GlobalConfig;

    if (cfg.flow_format == "rawbinary") {
	blkIO = new BinaryBlockIO();
    }
    else if (cfg.flow_format == "gziptext") {
	blkIO = new GzipBlockIO();
    }
    else {
	throw new LmrException(format("Flow format type '%s' unknown", cfg.flow_format));
    }

    blkIO.writeMetadataToFile(lmrCfg.flowMetadataFile);

    return 0;
}

extern(C) int luafn_writeInitialFlowFile(lua_State *L)
{
    auto blkId = to!int(luaL_checkinteger(L, 1));
    auto fname = flowFilename(lmrCfg.initialFieldDir, blkId);
    auto grid = checkStructuredGrid(L, 2);
    FluidBlock blk;
    if (grid) { // We do have a structured grid
	blk = new SFluidBlock(L);
    }
    else { // Assume an unstructured grid
	blk = new UFluidBlock(L);
    }
    blkIO.writeVariablesToFile(fname, blk.cells);
    return 0;
}

 
