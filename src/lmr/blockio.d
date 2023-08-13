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
import std.format : format, formattedRead;

import gzip : GzipOut, GzipByLine;
import dyaml;

import util.lua;
import nm.number;
import geom.luawrap;

import lmrexceptions : LmrException;
import globalconfig : GlobalConfig;
import lmrconfig;
import fvcell : FVCell;
import fvcellio;
import fluidblock : FluidBlock;
import sfluidblock : SFluidBlock;
import ufluidblock : UFluidBlock;
import flowsolution : FluidBlockLite;


BlockIO blkIO;

class BlockIO {
public:
    this() {}
    this(FVCellIO cio)
    {
	mCIO = cio.dup;
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
	f.writefln("  version: \"%s\"", BLK_IO_VERSION);
	f.writeln("  revision-id: ", lmrCfg.revisionId);
	f.writefln("  field-type: \"%s\"", fieldVarsTypeName(mCIO.FVT));
	f.writefln("  variables:");
	foreach (var; mCIO.variables) {
	    f.writefln("   - %s", var);
	}
	f.close();
    }

    final readMetadataFromFile(string fname)
    {
	Node metadata = dyaml.Loader.fromFile(fname).load();
	string blkIOversion = metadata["version"].as!string;
	if (blkIOversion != BLK_IO_VERSION) {
	    string errMsg = "ERROR: The flow metadata version is unknown.\n";
	    errMsg ~= format("Expected: %s   Received: %s", BLK_IO_VERSION, blkIOversion);
	    throw new LmrException(errMsg);
	}

	string fvtString = metadata["field-type"].as!string;
	FieldVarsType fvt = fieldVarsTypeFromName(fvtString);

	string[] variables;
	variables.length = 0;
	foreach (node; metadata["variables"].sequence) {
	    variables ~= node.as!string;
	}
	mCIO = createFVCellIO(fvt, variables);
    }

abstract:
    void writeVariablesToFile(string fname, FVCell[] cells);
    void readVariablesFromFile(string fname, FVCell[] cells);

private:
    FVCellIO mCIO;
}

class BinaryBlockIO : BlockIO {
public:

    this() {}

    this(FVCellIO cio)
    {
	super(cio);
    }

    override
    void writeVariablesToFile(string fname, FVCell[] cells)
    {
	double[1] dbl1;
	auto outfile = File(fname, "wb");
	foreach (var; mCIO.variables) {
	    foreach (cell; cells) {
		dbl1[0] = mCIO[cell,var];
		outfile.rawWrite(dbl1);
	    }
	}
	outfile.close();
    }

    override
    void readVariablesFromFile(string fname, FVCell[] cells)
    {
	double[1] dbl1;
	auto infile = File(fname, "rb");
	foreach (var; mCIO.variables) {
	    foreach (cell; cells) {
		infile.rawRead(dbl1);
		mCIO[cell,var] = dbl1[0];
	    }
	}
	infile.close();
    }

}

class GzipBlockIO : BlockIO {
public:
    this() {}

    this(FVCellIO cio)
    {
	super(cio);
    }

    override
    void writeVariablesToFile(string fname, FVCell[] cells)
    {
	auto outfile = new GzipOut(fname);
	foreach (var; mCIO.variables) {
	    foreach (cell; cells) {
		outfile.compress(format(lmrCfg.dblVarFmt ~ "\n", mCIO[cell,var]));
	    }
	}
	outfile.finish();
    }

    override
    void readVariablesFromFile(string fname, FVCell[] cells)
    {
	auto byLine = new GzipByLine(fname);
	string line;
	double dbl;
	foreach (var; mCIO.variables) {
	    foreach (cell; cells) {
		line = byLine.front; byLine.popFront;
		formattedRead(line, "%e", &dbl);
		mCIO[cell,var] = dbl;
	    }
	}
    }
}


//-------------------------------------------------------------
// Functions intended for use from Lua at prep stage:
//
//     1. luafn_writeFlowMetadata
//     2. luafn_writeInitialFlowFile
//
//-------------------------------------------------------------

extern(C) int luafn_writeFlowMetadata(lua_State *L)
{
    alias cfg = GlobalConfig;
    FVCellIO cio = new FVCellFlowIO(buildFlowVariables());

    if (cfg.flow_format == "rawbinary") {
	blkIO = new BinaryBlockIO(cio);
    }
    else if (cfg.flow_format == "gziptext") {
	blkIO = new GzipBlockIO(cio);
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


//-------------------------------------------------------------
// Functions intended for use at post-processing stage:
//
//     1. readFlowVariablesFromFlowMetadata
//     2. readFlowVariablesFromFile
//
//-------------------------------------------------------------
string[] readFlowVariablesFromFlowMetadata()
{
    Node metadata = dyaml.Loader.fromFile(lmrCfg.flowMetadataFile).load();
    string[] variables;
    foreach (node; metadata["variables"].sequence) {
	variables ~= node.as!string;
    }
    return variables;
}

void readFlowVariablesFromFile(FluidBlockLite blk, string fname, string[] variables, int ncells)
{
    size_t nvars = variables.length;
    blk._data.length = ncells;
    foreach (ref d; blk._data) d.length = nvars;

    final switch (blk.flow_format) {
    case "rawbinary":
	double[1] dbl1;
	auto infile = File(fname, "rb");
	foreach (ivar; 0 .. nvars) {
	    foreach (icell; 0 .. ncells) {
		infile.rawRead(dbl1);
		blk._data[icell][ivar] = dbl1[0];
	    }
	}
	infile.close();
	break;
    case "gziptext":
	auto byLine = new GzipByLine(fname);
	string line;
	double dbl;
	foreach (ivar; 0 .. nvars) {
	    foreach (icell; 0 .. ncells) {
		line = byLine.front; byLine.popFront;
		formattedRead(line, "%e", &dbl);
		blk._data[icell][ivar] = dbl;
	    }
	}
	break;
    }
}
