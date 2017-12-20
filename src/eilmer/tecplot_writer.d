/**
 * tecplot_writer.d
 *
 * Author: Pierpaolo Toniato
 * Date: 2017-05-14
 *
 * History:
 *  2017-06-21 -- Eilmer-specific code moved into its own module by RJG
 *
 */

module tecplot_writer;

import std.stdio;
import std.conv;
import std.string;

import tecio;
import geom;
import flowsolution;
import fvcore: FlowSolverException;

string[] zones =
    [
    "ORDERED",
    "FELINESEG",
    "FETRIANGLE",
    "FEQUADRILATERAL",
    "FETETRAHEDRON",
    "FEBRICK",
    "FEPOLYGON",
    "FEPOLYHEDRON"
     ];


// Transform the vtk connectivity to be suitable to use all elements 
// as quads/hexa in tecplot zones (where every single zone have single type)
size_t[] transformCellConnectivity(size_t[] ptOrder, int cellType) 
{
    size_t[] conn;
    switch (cellType) {
    default: // valid: ends with 'throw'
        throw new Exception("unknown cell type in transformCellConnectivity()");
    case -1: //structured grid to quad/bricks,all good
        conn = ptOrder.dup;
        break;
    case 1: //triangle 0122
        conn = [ptOrder[0], ptOrder[1], ptOrder[2], ptOrder[2]];
        break;
    case 2: //quad
        conn = ptOrder.dup;
        break;
    case 3: //polygonal data
        throw new FlowSolverException("Tecplot output for polygon data not avaiable");
    case 4: //tetra 01223333
        conn = [ptOrder[0], ptOrder[1], ptOrder[2], ptOrder[2],
		ptOrder[3], ptOrder[3], ptOrder[3], ptOrder[3]];
        break;
    case 5: //wedge 01223455
        conn = [ptOrder[0], ptOrder[1], ptOrder[2], ptOrder[2],
		ptOrder[3], ptOrder[4], ptOrder[5], ptOrder[5]];
        break;
    case 6: //hexa
        conn = ptOrder.dup;
        break;
    case 7: //pyramid 01234444
        conn = [ptOrder[0], ptOrder[1], ptOrder[2], ptOrder[3],
		ptOrder[4], ptOrder[4], ptOrder[4], ptOrder[4]];
        break;
    }
    return conn;
}

/**
 * Iterate through all the grid of a cell of a block and do the following:
 * determines whether the zone is only made of tetra/tria
 * if yes, return the unmodified list of of connectivity
 * if it is a mixed zone, or anything else, return a quad/hexa equivalent 
 * list of connectivity, ie make tecplot believe everything is a quad hexa
 * Necessary because tecplot does not support mixed elements in the same zone.
 *
 * Caller should pass an empty connList container.
 *
 * Returns:
 *  Amends zoneType and connList.
 */

void prepareGridConnectivity(Grid grid, ref int zoneType, ref size_t[][] connList)
{
    bool specialCase = true; // if true, only tris or tets
    int prevCellType;
    //Special cases: triangle and tetrahedra have a fast path
    //Assume that we have a zone full of tria/tetra
    zoneType = (grid.dimensions == 2) ? 2 : 4;
    prevCellType = (grid.dimensions == 2) ? 1 : 4;

    foreach (i; 0 .. grid.ncells)
    {
        auto cellType = grid.get_cell_type(i);
        if (specialCase) {
	    connList ~= grid.get_vtx_id_list_for_cell(i);
	    if (cellType != prevCellType ) {
		// We won't be visiting this branch of code anymore.
		specialCase = false;
		// writeln("Non tria/tetra cell detected: switching to quad/hexa elements");
		zoneType = (grid.dimensions == 2) ? 3 : 5;
		// We need to convert all of our previous connList
		foreach (j; 0 .. connList.length) {
		    if ( grid.dimensions == 2 ) {
			connList[j] = transformCellConnectivity(connList[j], 1);
		    }
		    else {
			connList[j] = transformCellConnectivity(connList[j], 4);
		    }
		}
            }
            prevCellType = cellType;
        }
	else {
	    connList ~= transformCellConnectivity(grid.get_vtx_id_list_for_cell(i), cellType);
	}
    }
}

int writeTecplotBinaryHeader(string jobName, int tindx, string fileName, string[] varNames)
{   
    string varstr;
    foreach (var; varNames) {
	varstr ~= " ";
	varstr ~= var;
    }
    auto title = format("%s-t%04d", jobName, tindx);
    return dtecini142(title, varstr, fileName);
}

File writeTecplotAsciiHeader(string jobName, int tindx, string fileName, string[] varNames)
{
    auto fp = File(fileName, "w");
    fp.writef("TITLE= %s-t%04d \n", jobName, tindx);
    fp.write("VARIABLES = ");
    foreach (var; varNames)
    {
        fp.writef("\"%s\", ", var);
    }
    fp.write("\n");
    return fp;
}

int writeTecplotBinaryZoneHeader(BlockFlow flow, Grid grid, size_t idx,
				 string[] varNames, double timestamp, int zoneType)
{
    if (flow.ncells != grid.ncells) {
        string msg = format("Mismatch between grid and flow: grid.ncells= %d flow.ncells= %d\n",
			    grid.ncells, flow.ncells);
        throw new FlowSolverException(msg);
    }
    int[] ValueLocation;
    auto zonetitle = format("block-%d",idx);
    foreach (var; varNames)
    {
        if ( (var == "pos.x") || (var == "pos.y") || (var == "pos.z") )
            ValueLocation ~= 1;
        else
            ValueLocation ~= 0;
    }
    //k=0 and StrandId=0
    return dteczne142(zonetitle, zoneType, to!int(grid.nvertices), to!int(grid.ncells), 0,
		      timestamp, 0, ValueLocation);
}

void writeTecplotAsciiZoneHeader(BlockFlow flow, Grid grid, size_t idx, File fp,
				 string[] varNames, double timestamp, int zoneType)
{
    size_t NumberOfPoints = grid.nvertices;
    size_t NumberOfCells = flow.ncells;
    int n_centered_vars = 0;

    if (flow.ncells != grid.ncells) {
        string msg = format("Mismatch between grid and flow: grid.ncells= %d flow.ncells= %d\n",
			    grid.ncells, flow.ncells);
        throw new FlowSolverException(msg);
    }
    fp.writef("ZONE T=\"block-%d\",", idx);
    fp.writef("SOLUTIONTIME=%.18e, ", timestamp);
    fp.writef("Nodes=%d , Elements=%d ", grid.nvertices, flow.ncells);
    fp.write(" DATAPACKING=BLOCK, ");
    fp.writef("ZONETYPE=%s,", zones[zoneType]);

    //Specifying which are cell centered variables
    foreach (var; varNames)
    {
        if (!((var == "pos.x") || (var == "pos.y") || (var == "pos.z")))
        {
            n_centered_vars += 1;
        }
    }
    fp.write("VARLOCATION =([");
    foreach (i; 0 .. n_centered_vars)
    {
        fp.writef("%d,", i+4); // +4 to account for pos.x, pos.y and pos.z in locations 1, 2, and 3.
    }
    fp.write("] = CELLCENTERED) \n");
}

void writeTecplotBinaryZoneData(BlockFlow flow, Grid grid, string[] varNames, size_t[][] conn)
{
    //writing nodal information
    double[] x, y, z;
    foreach (i; 0 .. grid.nvertices) {
        x ~= grid[i].x;
        y ~= grid[i].y;
        z ~= grid[i].z;
    }
    if (dtecdat142(x) != 0)
	throw new FlowSolverException("Error writing data in writeTecplotBinaryZoneData");
    if (dtecdat142(y) != 0)
	throw new FlowSolverException("Error writing data in writeTecplotBinaryZoneData");
    if (dtecdat142(z) != 0)    
	throw new FlowSolverException("Error writing data in writeTecplotBinaryZoneData");

    foreach (var; varNames) {
        if ((var == "pos.x") || (var == "pos.y") || (var == "pos.z")) {
            continue;
        }
        else {
	    double[] v;
	    foreach (i; 0 .. grid.ncells) v ~= flow[var,i];
	    if (dtecdat142(v) != 0) 
		throw new FlowSolverException("Error writing data in writeTecplotBinaryZoneData");
        }
    }
    // Writing connectivity information
    foreach (connElem; conn) {
	size_t[] c;
	c.length = connElem.length;
	c[] = connElem[] + 1; // Tecplot indices are 1-based
        if (dtecnode142(c) != 0)
	    throw new FlowSolverException("Error writing data in writeTecplotBinaryZoneData");
    }
}

void writeTecplotAsciiZoneData(BlockFlow flow, Grid grid, File fp,
			       string[] varNames, size_t[][] conn)
{
    // Write the cell-centred flow data from a single block (index jb)
    // as an unstructured grid of finite-volume cells.
    int j = 0;
    int values_per_row = 10; // 10 values per lines. Free to change to any number
    //Tecplot default is SINGLE precision
    //If more precision is necessary, make sure you modify ZONE header accordingly
    string values_format = " %.8e";
    // Writing nodal data for x,y,z 
    foreach (i; 0 .. grid.nvertices) {
        fp.writef(values_format, uflowz(grid[i].x));
        j += 1;
        if (j == values_per_row) {
            fp.write("\n");
            j = 0;
        }
    }
    fp.write("\n");
    j = 0;
    foreach (i; 0 .. grid.nvertices)
    {
        fp.writef(values_format, uflowz(grid[i].y));
        j += 1;
        if (j == values_per_row) {
            fp.write("\n");
            j = 0;
        }
    }
    fp.write("\n");
    j = 0;
    foreach (i; 0 .. grid.nvertices)
    {
        fp.writef(values_format, uflowz(grid[i].z));
        j += 1;
        if (j == values_per_row) {
            fp.write("\n");
            j = 0;
        }
    }
    fp.write("\n");
    j = 0;
    // Writing element data for all vars
    foreach (var; varNames) {
        if ((var == "pos.x") || (var == "pos.y") || (var == "pos.z"))
            continue;
        foreach (i; 0 .. grid.ncells) {
            fp.writef(values_format, uflowz(flow[var, i]));
            j += 1;
            if (j == values_per_row) {
                fp.write("\n");
                j = 0;
            }
        }
        fp.write("\n");
    }
    fp.write("\n");
    // Writing connectivity information
    foreach (connElem; conn) {
        size_t[] c;
	c.length = connElem.length;
	c[] = connElem[] + 1; // Add 1 for 1-based indexing in Tecplot
        foreach (id; c) fp.writef(" %d", id);
        fp.write("\n");
    }
}

int closeTecplotBinaryFile()
{
    return dtecend142();
}
