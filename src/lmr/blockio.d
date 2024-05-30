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
import std.algorithm: canFind;

import gzip : GzipOut, GzipByLine;
import dyaml;

import util.lua;
import nm.number;
import geom.luawrap;

import lmrexceptions : LmrException;
import globalconfig : GlobalConfig;
import lmrconfig;
import lmr.fvcell : FVCell;
import fvcellio;
import fluidblock : FluidBlock;
import sfluidblock : SFluidBlock;
import ufluidblock : UFluidBlock;
import flowsolution : FluidBlockLite;
import solidblock : SolidBlock;
import ssolidblock : SSolidBlock;
import solidsolution : SolidBlockLite;
import globaldata : fluidBlkIO, solidBlkIO;

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
    final void writeMetadataToFile(string fname)
    {
    	auto f = File(fname, "w");
        f.writefln("---");
        f.writefln("  version: \"%s\"", BLK_IO_VERSION);
        f.writeln("  revision-id: ", lmrCfg.revisionId);
        f.writefln("  field-type: \"%s\"", fieldVarsTypeName(mCIO.FVT));
        f.writefln("  variables:");
        foreach (var; mCIO.variables) {
            f.writefln("   - '%s'", var);
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
    void readVariablesFromFile(string fname, FVCell[] cells, string[] skipList = []);

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
    void readVariablesFromFile(string fname, FVCell[] cells, string[] skipList=[])
    {
        double[1] dbl1;
        auto infile = File(fname, "rb");
        foreach (var; mCIO.variables) {
            bool skip = (canFind(skipList, var)) ? true : false;
            foreach (cell; cells) {
                infile.rawRead(dbl1);
                if (skip) continue;
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
    void readVariablesFromFile(string fname, FVCell[] cells, string[] skipList=[])
    {
        auto byLine = new GzipByLine(fname);
        string line;
        double dbl;
        foreach (var; mCIO.variables) {
            bool skip = (canFind(skipList, var)) ? true : false;
            foreach (cell; cells) {
                if (skip) { // pop the line and move on, we don't need to read this data into a cell
                    byLine.popFront;
                    continue;
                }
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
//     1. luafn_writeFluidMetadata
//     2. luafn_writeInitialFluidFile
//     3. luafn_writeSolidMetadata
//     4. luafn_writeInitialSolidFile
//
//-------------------------------------------------------------

/**
 * Lua-exposed function to write fluid metadata.
 *
 * Author: RJG
 * Date: 2023
 */
extern(C) int luafn_writeFluidMetadata(lua_State *L)
{
    alias cfg = GlobalConfig;
    FVCellIO cio = new FluidFVCellIO();

    if (cfg.field_format == "rawbinary") {
    	fluidBlkIO = new BinaryBlockIO(cio);
    }
    else if (cfg.field_format == "gziptext") {
	    fluidBlkIO = new GzipBlockIO(cio);
    }
    else {
	    throw new LmrException(format("Fluid format type '%s' unknown", cfg.field_format));
    }

    fluidBlkIO.writeMetadataToFile(lmrCfg.fluidMetadataFile);

    return 0;
}

/**
 * Lua-exposed function to write an initial fluid field file.
 *
 * Author: RJG
 * Date: 2023
 */
extern(C) int luafn_writeInitialFluidFile(lua_State *L)
{
    auto blkId = to!int(luaL_checkinteger(L, 1));
    auto fname = fluidFilename(lmrCfg.initialFieldDir, blkId);
    auto grid = checkStructuredGrid(L, 2);
    FluidBlock blk;
    if (grid) { // We do have a structured grid
    	blk = new SFluidBlock(L);
    }
    else { // Assume an unstructured grid
        blk = new UFluidBlock(L);
    }
    // Make temporary holder because we can't pass
    // an array of type FluidFVCells, so promote
    FVCell[] cells;
    cells.length = blk.cells.length;
    foreach (i, ref c; cells) c = blk.cells[i];
    fluidBlkIO.writeVariablesToFile(fname, cells);
    return 0;
}

/**
 * Lua-exposed function to write solid metadata.
 *
 * Author: RJG
 * Date: 2024-03-02
 */
extern(C) int luafn_writeSolidMetadata(lua_State *L)
{
    FVCellIO cio = new SolidFVCellIO();

    if (GlobalConfig.field_format == "rawbinary") {
    	solidBlkIO = new BinaryBlockIO(cio);
    }
    else if (GlobalConfig.field_format == "gziptext") {
	    solidBlkIO = new GzipBlockIO(cio);
    }
    else {
	    throw new LmrException(format("Solid format type '%s' unknown", GlobalConfig.field_format));
    }

    solidBlkIO.writeMetadataToFile(lmrCfg.solidMetadataFile);

    return 0;
}

/**
 * Lua-exposed function to write an initial solid field file.
 *
 * Author: RJG
 * Date: 2024-03-02
 */
extern(C) int luafn_writeInitialSolidFile(lua_State *L)
{
    auto blkId = to!int(luaL_checkinteger(L, 1));
    auto fname = solidFilename(lmrCfg.initialFieldDir, blkId);
    auto grid = checkStructuredGrid(L, 2);
    SolidBlock solidBlk;
    if (grid) { // We do have a structured grid
    	solidBlk = new SSolidBlock(L);
    }
    else { // Assume an unstructured grid
        throw new LmrException("Unstuctured solid blocks are not implemented.");
    }
    // Make temporary holder because we can't pass
    // an array of type SolidFVCells, so promote
    FVCell[] cells;
    cells.length = solidBlk.cells.length;
    foreach (i, ref c; cells) c = solidBlk.cells[i];
    solidBlkIO.writeVariablesToFile(fname, cells);
    return 0;
}

//-------------------------------------------------------------
// Functions intended for use at post-processing stage:
//
//     1. readFluidVariablesFromFlowMetadata
//     2. readFluidVariablesFromFile
//     3. readVariablesFromMetadata
//     4. readValuesFromFile
//
//  1. and 2. are special cases because we may want to do further
//  manipulations with fluid data, so we load a FluidBlockLite.
//  3. and 4. are general functions. We simply want to pick up
//  from Eilmer-native format and write out in some other format.
//-------------------------------------------------------------
string[] readFluidVariablesFromFluidMetadata(string fluidMetadataFile="")
{
    fluidMetadataFile = (fluidMetadataFile == "") ? lmrCfg.fluidMetadataFile : fluidMetadataFile;
    Node metadata = dyaml.Loader.fromFile(fluidMetadataFile).load();
    string[] variables;
    foreach (node; metadata["variables"].sequence) {
        variables ~= node.as!string;
    }
    return variables;
}

/*
 NOTE: RJG, 2024-03-04
 It would be possible to unify the next two functions,
 but it would require creation of a BlockLite superclass.
 Perhaps that overhead isn't worth it for the small
 amount of duplication here.
*/

void readFluidVariablesFromFile(FluidBlockLite blk, string fname, string[] variables, int ncells)
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

void readSolidVariablesFromFile(SolidBlockLite blk, string fname, string[] variables, int ncells)
{
    size_t nvars = variables.length;
    blk.data.length = ncells;
    foreach (ref d; blk.data) d.length = nvars;

    final switch (blk.fieldFmt) {
    case "rawbinary":
        double[1] dbl1;
        auto infile = File(fname, "rb");
        foreach (ivar; 0 .. nvars) {
            foreach (icell; 0 .. ncells) {
                infile.rawRead(dbl1);
                blk.data[icell][ivar] = dbl1[0];
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
                blk.data[icell][ivar] = dbl;
            }
        }
	    break;
    }
}

string[] readVariablesFromMetadata(string metadataFile)
{
    Node metadata = dyaml.Loader.fromFile(metadataFile).load();
    string[] variables;
    foreach (node; metadata["variables"].sequence) {
        variables ~= node.as!string;
    }
    return variables;
}

void readValuesFromFile(ref double[][] data, string fname, string[] variables, size_t ncells, string fileFmt)
{
    size_t nvars = variables.length;
    data.length = ncells;
    foreach (ref d; data) d.length = nvars;

    final switch (fileFmt) {
    case "rawbinary":
	double[1] dbl1;
	auto infile = File(fname, "rb");
	foreach (ivar; 0 .. nvars) {
        foreach (icell; 0 .. ncells) {
            infile.rawRead(dbl1);
            data[icell][ivar] = dbl1[0];
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
            data[icell][ivar] = dbl;
	    }
	}
	break;
    }
}
