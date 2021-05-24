module fluidblockio_new;

import std.json;
import std.string;
import std.conv;

import std.range : iota;
import std.array;
import std.traits : isNumeric;
import std.file: write, read;
import core.stdc.string : strlen;
import std.stdio;
import std.typecons;
import std.datetime.systime;
import std.format;

import zip;
import fluidblock;
import sfluidblock;
import ufluidblock;
import fvcell;
import globalconfig;
import fluidblockio_old : k_modesName, massfName, u_modesName, T_modesName;
import sfluidblock : cell_index_to_logical_coordinates;
import geom;
import geom.luawrap;
import flowstate;
import luaflowstate;
import util.lua;
import util.lua_service;
import nm.number;
import flowsolution;
import gas;
import globaldata;
import simcore;
import celldata;

/*

The new I/O for fluid variables is defined in two layers.

Layer 1: Each cell holds a list of auxiliary data objects which provide lists of
accessor objects. These accessors are a 'pipe' to data that we want to get or
set for I/O purposes. The auxiliary data objects may be a pass-through to data
that is held elsewhere, such as the parent FVCell, or will actually hold some
data themselves. If they hold data there is an update() function that allows for
that data to be modified at any iteration in the main loop (if desired). Note
that this layer could be skipped altogether as we can simply define all of the
accessors directly in the next layer!

Layer 2: Each FluidBlock holds a list of FluidBlockIO instances, each of which
hold lists of accessors obtained from any number of auxiliary data objects. Note
that these accessors are singletons (only one exists) and must be passed a cell
instance to actually access any data. The FluidBlockIO then provides all of the
functionality for reading/writing from file and getting/setting data through the
accessor pipes. The file structure is defined as one zip file per block with a
single text header file in JSON format and an arbitrary number of binary or text
blobs of cell data. The header thus provides all of the metadata which allows
the blobs to be correctly interpreted. Note that each FluidBlockIO instance will
write to a separate folder within the overall file structure.

The flow of data could look like this

FVCell <-> AuxCellData <-> VariableAccess <-> FluidBlockIO <-> Zip

or this

FVCell <-> VariableAccess <-> AuxCellData <-> VariableAccess <-> FluidBlockIO
<-> Zip

or this

FVCell <-> VariableAccess <-> FluidBlockIO <-> Zip

I will now describe how flow averages were added to I/O:
- Define a new AuxCellData sub-class (FlowAverage) that accumulated the flow
  averages
- Make sure that it is added to the parent AuxCellData indexing functions!!
- Define VariableAccess sub-classes to access all of the named member variables
  within FlowAverage
- Implement a static function within FlowAverage that serves up all pertinent
  accessors
- Define a new FluidBlockIO sub-class (TimeAverageIO) that holds all of the
  accessors provided from FlowAverage
- Make sure that it is added to the get_fluid_block_io() function!!
- Done!

To post-process any of the generated flow files use the new option --plotTag

Additional notes for developers:
- To add extra items to a FluidBlockIO instance we simply append more acccessors
- AuxCellData was implemented so that we had a list where auxiliary data within
  each cell could be allocated on demand rather than always being present. The
  downside is that it is a dynamic allocation.
- Slicing is something that can be implemented in the FluidBlockIO instances as
  each one carries around a list of cell indexs (potential optimization, if
  index list is empty assume all present)
- The string format within each blob is a little more verbose due to the
  requirement of knowing the size of the binary represntation. This is required
  as we are loading directly into the zip file buffer rather than into a dynamic
  array and then into a zip buffer. This is to reduce memory overhead. It may be
  preferable to write directly to disk, unfortunately this would generate a
  proliferation of files due to the multiple blob structure, hence the zip
  format.

flow format definition:
- The output zip file contains a header file, a JSON data object which contains
  metadata about the block (ncells, sim_time, etc) with a list of the recorded
  quantities with some metadata about the type of each quantity, and one data
  file for each recorded variable. Each data file is basically a data blob, a
  series of numbers corresponding to the value of the quantity in each cell.
- Information in the header file looks like the following where options are
  given in [...]. See the functions save_to_file() and read_from_file() for
  implementation details.

    {
        # version number of the file format
        "version": "1.0",

        # the format of the data held in the blobs
        "data_type": ["string","binary"],

        # delimiter that is used when writing to string
        "delimiter": " ",

        # how many cells total
        "n_cells": 44,

        # number of cells in i-index direction
        # if it is unstructured then nic = n_cells and njc=nkc=0
        "nic": 11,

        # number of cells in j-index direction
        "njc": 4,

        # number of cells in k-index direction
        "nkc": 1,

        # simulation time
        "sim_time": 0,

        # structured or unstructured data
        "structured": [true,false],

        # unique string identifier for this particular grouping of data
        "tag": "field",

        # list of all data blobs in this zip file
        "variables": {

            # the name of the data
            "B.x": {

                # the name of the data blob,
                "data": "B.x.dat",

                # a text description of the data
                "description": "[cell.fs.B._p[0].re]",

                # dimensionality of this particular data blob, the first dimension
                # corresponds to the linearised cell index. This allows for 1D arrays
                # of data to be stored for each cell.
                "dimension": [
                    44,
                    1
                ]
            },

            ...

        }
    }

*/


ubyte[] bytes(T)(T num) if (isNumeric!T)
// convert any number to bytes
{
    auto buf = new ubyte[T.sizeof];
    (*cast(T*)(buf.ptr)) = num;
    return buf;
}

T value(T)(ubyte[] buf) if (isNumeric!T)
// convert bytes to number
{
    return (*cast(T*)(buf.ptr));
}

ubyte[] to_binary(T)(T[] buffer, const bool binary, const string delimiter=" ")
{
    const int data_size = T.sizeof;
    if (binary) {
        ubyte[] bbuffer;
        bbuffer.length = buffer.length*data_size;
        foreach (idx, val; buffer) {
            bbuffer[idx*data_size .. (idx+1)*data_size] = val.bytes[];
        }
        return bbuffer;
    } else {
        string sbuffer;
        string str_format;
        if (is(T == int) || is(T == size_t)) {
            str_format = FluidBlockIO.int_fmt;
        } else {
            str_format = FluidBlockIO.float_fmt;
        }
        foreach (val; buffer) {
            sbuffer ~= format(str_format,val,delimiter);
        }
        return sbuffer.dup.representation;
    }
}

// sub-function to package data for inclusion in a zip archive
ArchiveMember archive_data(T)(const string name, T[] buffer, const bool binary, const string delimiter=" ")
{
    const int data_size = T.sizeof;
    // add the data to the archive
    ArchiveMember dat = new ArchiveMember();
    dat.name = name;
    dat.time(Clock.currTime);
    dat.expandedData(to_binary(buffer, binary, delimiter));
    dat.compressionMethod = CompressionMethod.deflate;
    return dat;
}

T[] from_binary(T)(ubyte[] dat_bytes, const size_t size, const bool binary, const string delimiter=" ")
{
    T[] buffer;
    buffer.length = size;
    if (binary) {
        const size_t data_size = T.sizeof;
        foreach (idx, ref val; buffer) {
            val = dat_bytes[idx*data_size .. (idx+1)*data_size].value!T;
        }
    } else {
        char* dat_cstr = cast(char*) dat_bytes;
        string dat_str = cast(string) dat_cstr[0..strlen(dat_cstr)];
        auto dat_items = split(dat_str, delimiter);
        foreach (idx, ref val; buffer) {
            formattedRead(dat_items[idx], "%s", val);
        }
    }
    return buffer;
}

T[] retrieve_data(T)(const string name, ref ZipArchive zip, const int size,
    const bool binary, const string delimiter=" ")
{
    ArchiveMember dat = zip.directory[name];
    ubyte[] dat_bytes = zip.expand(dat);

    return from_binary!T(dat_bytes, size, binary, delimiter);
}

class FluidBlockIO {
    // sits in-between the block and the file on disk
    // acts as buffer and communicator to get and set data in the block
    public:

    string tag = "null";
    VariableAccess[] accessors;
    string[] names;
    size_t[] index;
    bool success;
    bool binary = false;
    string delimiter = " ";
    size_t data_element_size;
    bool initialized = false;

    static const string data_ext = ".dat";
    static const string float_fmt = "%026.18e%s";
    static const string int_fmt = "%d%s";

    FluidBlock block; // the parent of this IO instance

    this(){}

    void set_binary(bool is_binary=false)
    {
        binary = is_binary;
        // get the size in bytes of a single element of data
        if (binary) {
            data_element_size = double.sizeof;
        } else {
            data_element_size = format(float_fmt,0.0, delimiter).representation.length;
        }
    }

    void add_accessor(string name, VariableAccess accessor)
    {
        names ~= name;
        accessors ~= accessor;
    }

    void make_buffers(size_t n_cells)
    {
        if (n_cells <= 0) return;
        foreach (i, accessor ; accessors) {
            auto arc = new ArchiveMember();
            arc.name = names[i]~data_ext;
            arc.compressionMethod = CompressionMethod.deflate;
            arc.setExpandedDataSize(accessor.num_return*n_cells*data_element_size);
            archive_members ~= arc;
        }
        initialized = true;
    }

    void get_cell_data(FVCell cell, size_t cell_count)
    {
        // retrieve stuff from a cell
        double[] dat;
        foreach (access_count, accessor; accessors) {
            size_t size = accessor.num_return;
            dat.length = size;
            dat = accessor.get(cell);
            // convert to binary
            size *= data_element_size;
            size_t start_idx = cell_count*size;
            // size_t stop_idx = (cell_count+1)*size;
            auto arc = archive_members[access_count];
            auto buffer = to_binary(dat, binary, delimiter);
            arc.expandedData(buffer, start_idx);
        }
    }

    void set_cell_data(FVCell cell, size_t cell_count)
    {
        // set stuff in a cell
        double[] dat;
        ubyte[] buffer;
        foreach (access_count, accessor; accessors) {
            size_t size = accessor.num_return;
            dat.length = size;
            buffer.length = size*data_element_size;
            size_t start_idx = cell_count*size*data_element_size;
            size_t stop_idx = (cell_count+1)*size*data_element_size;
            buffer = archive_members[access_count].expandedData()[start_idx .. stop_idx];
            dat = from_binary!double(buffer, size, binary, delimiter);
            accessor.set(cell, dat);
        }
    }

    void get_block_data()
    {
        foreach (cell_count, cell_idx; index) {
            FVCell cell = block.cells[cell_idx];
            get_cell_data(cell, cell_count);
        }
    }

    void set_block_data()
    {
        foreach (cell_count, cell_idx; index) {
            FVCell cell = block.cells[cell_idx];
            set_cell_data(cell, cell_count);
        }
    }

    void save_to_file(const string fname, double time=0.0)
    // utility function for writing block data
    // assumes data is already loaded into this object by get_block_data
    // any additional header data needs to be passed in through info
    // i.e. label, dimensions, etc
    {
        // harvest the data
        get_block_data();
        string delimiter = " ";
        ZipArchive save = new ZipArchive();
        // header
        JSONValue header = block.get_header();
        header["version"] = "1.0";
        header["tag"] = tag;
        header["sim_time"] = time;
        // the cell indexes held by this block of data
        size_t n_cells = index.length;
        header["n_cells"] = to!int(n_cells);
        save.addMember(archive_data("cell_id", index, binary, delimiter));
        // now all of the data
        JSONValue variables = [names[0]: "place holder"]; // need to initialise with something
        // add an entry to the header file, get the data and add it to the archive
        foreach(i, ref member; archive_members) {
            string varname = names[i];
            VariableAccess accessor = accessors[i];
            const size_t size = accessor.num_return;
            // describe the data that is being saved
            JSONValue jj = ["data":varname~data_ext]; // where is the data
            jj.object["dimension"] = JSONValue([n_cells, size]); // how much data
            jj.object["description"] = JSONValue(accessor.description());
            // add the data to the archive
            member.time(Clock.currTime);
            save.addMember(member);
            variables.object[varname] = jj;
        }
        header["variables"] = variables;
        if (binary) {
            header["data_type"] = JSONValue("binary");
        } else {
            header["data_type"] = JSONValue("string");
            header["delimiter"] = JSONValue(delimiter);
        }
        // add the header to the archive
        ArchiveMember hdr = new ArchiveMember();
        hdr.name = "header.txt";
        hdr.time(Clock.currTime);
        hdr.expandedData(header.toPrettyString.dup.representation);
        hdr.compressionMethod = CompressionMethod.deflate;
        save.addMember(hdr);
        // write to file
        write(fname, save.build());
    }

    JSONValue read_from_file(const string fname)
    // utility function for reading block data from zip file
    {

        // get the header and convert to string
        ZipArchive zip = new ZipArchive(read(fname));
        // parse the header
        ArchiveMember hdr = zip.directory["header.txt"];
        ubyte[] hdr_bytes = zip.expand(hdr);
        char* hdr_cstr = cast(char*) hdr_bytes;
        string hdr_str = cast(string) hdr_cstr[0..strlen(hdr_cstr)];
        JSONValue header = parseJSON(hdr_str);
        string version_str = header["version"].get!string;
        // check the version
        if (version_str != "1.0") {
            string msg = text("File format version found: " ~ version_str);
            throw new FlowSolverException(msg);
        }
        string tag_read = header["tag"].to!string;
        tag_read = tag_read.strip("\"");
        // check the tag
        if (tag != tag_read) {
            success = false;
            return header;
        }
        // get the data type
        const string data_type = header["data_type"].get!string;
        if (data_type == "string") {
            delimiter = header["delimiter"].get!string;
        }
        // retrive the cell indexing
        int size = header["n_cells"].get!int;
        index = retrieve_data!size_t("cell_id", zip, size, data_type!="string", delimiter);
        //
        if (!initialized)
            make_buffers(size);
        //
        foreach (i, accessor ; accessors) {
            string varname = names[i];
            if (!(varname in header["variables"])) {
                throw new Exception("Could not find variable: " ~ varname);
            }
            JSONValue data_hdr = header["variables"][varname];
            const string dat_name = data_hdr["data"].get!string;
            JSONValue[] jdim = data_hdr["dimension"].get!(JSONValue[]);
            size = 1;
            int[2] dimension;
            foreach (size_t idx, dim; jdim) {
                dimension[idx] = dim.get!int;
                size *= dimension[idx];
            }
            archive_members[i] = zip.directory[dat_name];
            zip.expand(archive_members[i]);
        }
        success = true;
        return header;
    }

    bool do_save()
    {
        return true;
    }

    private:
    ArchiveMember[] archive_members;
} // end FluidBlockIO

FluidBlockIO[] get_fluid_block_io(FluidBlock blk=null)
{
    FluidBlockIO[] io_list;
    if (!GlobalConfig.new_flow_format) return io_list;
    // output cell variables by default
    io_list ~= new CellVariablesIO(blk);
    if (GlobalConfig.do_flow_average) {
        io_list ~= new TimeAverageIO(blk);
    }
    if (GlobalConfig.viscous && GlobalConfig.save_viscous_gradients) {
        io_list ~= new CellGradIO(blk);
    }
    if (GlobalConfig.do_temporal_DFT) {
        io_list ~= new DFTIO(blk);
    }
    if (GlobalConfig.save_limiter_values) {
        io_list ~= new CellLimiterIO(blk);
    }
    if (GlobalConfig.solve_electric_field) {
        io_list ~= new FieldIO(blk);
    }

    return io_list;
}

extern(C) int write_initial_sg_flow_zip_file_from_lua(lua_State* L)
{
    auto fname = to!string(luaL_checkstring(L, 1));
    SFluidBlock block = new SFluidBlock(L);
    foreach(io; block.block_io) {
        if (io.tag == CellData.tag)
            io.save_to_file(fname);
    }
    return 0;
} // end write_initial_sg_flow_file_from_lua()

extern(C) int write_initial_usg_flow_zip_file_from_lua(lua_State* L)
{
    auto fname = to!string(luaL_checkstring(L, 1));
    UFluidBlock block = new UFluidBlock(L);
    foreach(io; block.block_io) {
        if (io.tag == CellData.tag)
            io.save_to_file(fname);
    }
    return 0;
} // end write_initial_usg_flow_file_from_lua()

double read_zip_solution(FluidBlock blk, string filename)
{
    double sim_time = -1;
    foreach(i, io; blk.block_io) {
        JSONValue header = io.read_from_file(filename);
        if (!io.success) continue;
        io.set_block_data();
        if (i == 0) {
            sim_time = header["sim_time"].get!double;
        } else {
            assert(sim_time == header["sim_time"].get!double);
        }
    }
    return sim_time;
}

void read_block_data_from_zip_file(FluidBlockLite blk, string filename, Grid_t gridType, string flow_format)
{
    // get the header and convert to string
    ZipArchive zip = new ZipArchive(read(filename));
    // parse the header
    ArchiveMember hdr = zip.directory["header.txt"];
    ubyte[] hdr_bytes = zip.expand(hdr);
    char* hdr_cstr = cast(char*) hdr_bytes;
    string hdr_str = cast(string) hdr_cstr[0..strlen(hdr_cstr)];
    JSONValue header = parseJSON(hdr_str);
    string version_str = header["version"].get!string;
    // check the version
    if (version_str != "1.0") { // <-- version should be a variable
        string msg = text("File format version found: " ~ version_str);
        throw new FlowSolverException(msg);
    }
    // get the data type
    const string data_type = header["data_type"].get!string;
    bool is_binary = data_type!="string";
    string delimiter = " ";
    if (!is_binary) delimiter = header["delimiter"].get!string;
    // retrive the cell indexing
    const int n_cells = header["n_cells"].get!int;
    size_t[] cell_index = retrieve_data!size_t("cell_id", zip, n_cells, is_binary, delimiter);
    blk.sim_time = header["sim_time"].get!double;
    blk.nic = header["nic"].get!int;
    blk.njc = header["njc"].get!int;
    blk.nkc = header["nkc"].get!int;
    blk.ncells = n_cells;
    blk._data.length = n_cells;
    size_t count_names = 0;
    foreach (string varname, ref JSONValue data; header["variables"]) {
        // get the dimensionality of the data (does it have multiple values for each cell?)
        const JSONValue[] jdim = data["dimension"].get!(JSONValue[]);
        const int n_per_cell = jdim[1].get!int;
        const int n_total = jdim[0].get!int*n_per_cell;
        // names for the data in FluidBlockLite
        if (n_per_cell == 1) {
            blk.variableNames ~= varname;
        } else {
            foreach (idx; 0 .. n_per_cell) {
                blk.variableNames ~= format("%s_%d",varname,idx);
            }
        }
        // allocate the data array, extending as necessary
        // could do this upfront...
        foreach (j; 0 .. n_cells) {
            blk._data[j].length = blk.variableNames.length;
        }
        // now read the data
        const string dat_name = data["data"].get!string;
        double[] dat = retrieve_data!double(dat_name, zip, n_total, is_binary, delimiter);
        foreach (k; 0 .. n_per_cell) {
            const size_t name_idx = count_names + k;
            foreach (j; cell_index) {
                blk._data[j][name_idx] = dat[j*n_per_cell + k];
            }
            blk.variableIndex[blk.variableNames[name_idx]] = name_idx;
        }
        count_names += n_per_cell;
    }
    return;
}

class CellVariablesIO : FluidBlockIO
{
    public:

    this(FluidBlock blk=null)
    {
        tag = CellData.tag;
        size_t n_cells = 0;
        LocalConfig myConfig;
        if (blk !is null) {
            block = blk;
            myConfig = block.myConfig;
            n_cells = block.cells.length;
        } else {
            myConfig = new LocalConfig(-1);
            myConfig.gmodel = GlobalConfig.gmodel_master;
        }
        // get all of the cells data
        index = iota(0, n_cells).array();
        set_binary(myConfig.flow_format == "eilmer4binary");
        // get the accessors
        foreach (key, acc ; CellData.get_accessors(myConfig)) {
            add_accessor(key, acc);
        }
        make_buffers(n_cells);
    }
}

class TimeAverageIO : FluidBlockIO
{
    public:

    this(FluidBlock blk)
    {
        tag = FlowAverage.tag;
        size_t n_cells = 0;
        LocalConfig myConfig;
        if (blk !is null) {
            block = blk;
            myConfig = block.myConfig;
            n_cells = block.cells.length;
        } else {
            myConfig = new LocalConfig(-1);
            myConfig.gmodel = GlobalConfig.gmodel_master;
        }
        // get all of the cells data
        index = iota(0, n_cells).array();
        set_binary(myConfig.flow_format == "eilmer4binary");
        // get the accessors
        foreach (key, acc ; FlowAverage.get_accessors(myConfig)) {
            add_accessor(key, acc);
        }
        make_buffers(n_cells);
    }
}

class CellGradIO : FluidBlockIO
{
    public:

    this(FluidBlock blk=null)
    {
        tag = CellGradientData.tag;
        size_t n_cells = 0;
        LocalConfig myConfig;
        if (blk !is null) {
            block = blk;
            myConfig = block.myConfig;
            n_cells = block.cells.length;
        } else {
            myConfig = new LocalConfig(-1);
            myConfig.gmodel = GlobalConfig.gmodel_master;
        }
        // get all of the cells data
        index = iota(0, n_cells).array();
        set_binary(myConfig.flow_format == "eilmer4binary");
        // get the accessors
        foreach (key, acc ; CellGradientData.get_accessors(myConfig)) {
            add_accessor(key, acc);
        }
        make_buffers(n_cells);
    }
}

class DFTIO : FluidBlockIO
{
    public:

    this(FluidBlock blk)
    {
        tag = GeneralDFT.tag;
        size_t n_cells = 0;
        if (blk !is null) {
            block = blk;
            n_cells = block.cells.length;
        }
        // get all of the cells data
        index = iota(0, n_cells).array();
        set_binary(GlobalConfig.flow_format == "eilmer4binary");
        foreach (key, acc ; GeneralDFT.get_accessors()) {
            add_accessor(key, acc);
        }
        make_buffers(n_cells);
    }

    override bool do_save()
    {
        if (SimState.time >= SimState.target_time) { return true; }
        if (SimState.step >= GlobalConfig.max_step) { return true; }
        if (GlobalConfig.halt_now == 1) { return true; }
        auto wall_clock_elapsed = (Clock.currTime() - SimState.wall_clock_start).total!"seconds"();
        if (SimState.maxWallClockSeconds > 0 && (wall_clock_elapsed > SimState.maxWallClockSeconds)) {
            return true;
        }
        return false;
    }
}

class CellLimiterIO : FluidBlockIO
{
    public:

    this(FluidBlock blk=null)
    {
        tag = CellLimiterData.tag;
        size_t n_cells = 0;
        LocalConfig myConfig;
        if (blk !is null) {
            block = blk;
            myConfig = block.myConfig;
            n_cells = block.cells.length;
        } else {
            myConfig = new LocalConfig(-1);
            myConfig.gmodel = GlobalConfig.gmodel_master;
        }
        // get all of the cells data
        index = iota(0, n_cells).array();
        set_binary(myConfig.flow_format == "eilmer4binary");
        // get the accessors
        foreach (key, acc ; CellLimiterData.get_accessors(myConfig)) {
            add_accessor(key, acc);
        }
        make_buffers(n_cells);
    }
}

class FieldIO : FluidBlockIO
{
    public:

    this(FluidBlock blk)
    {
        tag = FieldData.tag;

        size_t n_cells = 0;
        LocalConfig myConfig;

        if (blk !is null) {
            block = blk;
            myConfig = block.myConfig;
            n_cells = block.cells.length;
        } else {
            myConfig = new LocalConfig(-1);
            myConfig.gmodel = GlobalConfig.gmodel_master;
        }

        // get all of the cells data
        index = iota(0, n_cells).array();

        set_binary(myConfig.flow_format == "eilmer4binary");

        // get the accessors
        foreach (key, acc ; FieldData.get_accessors(myConfig)) {
            add_accessor(key, acc);
        }

        make_buffers(n_cells);
    }
}

