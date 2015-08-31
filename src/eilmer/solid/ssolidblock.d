/**
 * ssolidblock.d
 *
 * Author: Rowan G. and Peter J.
 * Version: 2015-22-04
 */

module ssolidblock;

import std.stdio;
import std.conv;
import std.format;
import std.array;
import std.math;
import std.json;
import util.lua;
import json_helper;
import gzip;
import geom;
import globalconfig;
import block;
import solidblock;
import solidfvcell;
import solidfvinterface;
import solidfvvertex;
import solidbc;
import solidprops;

class SSolidBlock : SolidBlock {
public:
    size_t nicell;
    size_t njcell;
    size_t nkcell;
    size_t imin, imax;
    size_t jmin, jmax;
    size_t kmin, kmax;
    size_t[] hicell, hjcell, hkcell; // locations of sample cells for history record
    SolidBoundaryCondition[6] bc;

private:
    size_t _nidim;
    size_t _njdim;
    size_t _nkdim;

    SolidFVCell[] _ctr;
    SolidFVInterface[] _ifi;
    SolidFVInterface[] _ifj;
    SolidFVInterface[] _ifk;
    SolidFVVertex[] _vtx;
    SolidFVInterface[] _sifi;
    SolidFVInterface[] _sifj;
    SolidFVInterface[] _sifk;

public:
    // Constructors
    this(int id, size_t nicell, size_t njcell, size_t nkcell, string label)
    {
	super(id, label);
	this.nicell = nicell;
	this.njcell = njcell;
	this.nkcell = nkcell;
	fillInOtherSizeData();
    }

    this(int id, JSONValue jsonData)
    {
	this.id = id;
	nicell = getJSONint(jsonData, "nic", 0);
	njcell = getJSONint(jsonData, "njc", 0);
	nkcell = getJSONint(jsonData, "nkc", 0);
	label = getJSONstring(jsonData, "label", "");
	this(id, nicell, njcell, nkcell, label);
	active = getJSONbool(jsonData, "active", true);
	sp = makeSolidPropsFromJson(jsonData["properties"]);
    }

    override void initLuaGlobals()
    {
	lua_pushinteger(myL, nicell); lua_setglobal(myL, "nicell");
	lua_pushinteger(myL, njcell); lua_setglobal(myL, "njcell");
	lua_pushinteger(myL, nkcell); lua_setglobal(myL, "nkcell");
	lua_pushinteger(myL, Face.north); lua_setglobal(myL, "north");
	lua_pushinteger(myL, Face.east); lua_setglobal(myL, "east");
	lua_pushinteger(myL, Face.south); lua_setglobal(myL, "south");
	lua_pushinteger(myL, Face.west); lua_setglobal(myL, "west");
	lua_pushinteger(myL, Face.top); lua_setglobal(myL, "top");
	lua_pushinteger(myL, Face.bottom); lua_setglobal(myL, "bottom");
    }

    override void initBoundaryConditions(JSONValue jsonData)
    {
	foreach (boundary; 0 .. (GlobalConfig.dimensions == 3 ? 6 : 4)) {
	    string jsonKey = "face_" ~ face_name[boundary];
	    auto bcJsonData = jsonData[jsonKey];
	    bc[boundary] = makeSolidBCFromJson(bcJsonData, id, boundary,
					       nicell, njcell, nkcell);
	}
    }

    void fillInOtherSizeData()
    // Helper function for the constructor
    {
	_nidim = nicell + 2 * nghost;
	_njdim = njcell + 2 * nghost;
	// Indices, in each grid direction for the active cells.
	// These limits are inclusive. The mincell and max cell
	// are both within the active set of cells.
	imin = nghost; imax = imin + nicell - 1;
	jmin = nghost; jmax = jmin + njcell - 1;
	if ( GlobalConfig.dimensions == 2 ) {
	    // In 2D simulations, the k range is from 0 to 0 for the
	    // storage arrays of cells and relevant faces.
	    if ( nkcell != 1 ) {
		writeln("Warning: inconsistent dimensions nkcell set to 1 for 2D");
		nkcell = 1;
	    }
	    _nkdim = 1;
	    kmin = 0; kmax = 0;
	} else {
	    // In 3D simulations the k index is just like the i and j indices.
	    _nkdim = nkcell + 2 * nghost;
	    kmin = nghost; kmax = kmin + nkcell - 1;
	}
    } // end fill_in_other_size_data()

    // -- Service methods
    @nogc
    size_t toGlobalIndex(size_t i, size_t j, size_t k) const
    in {
	assert(i < _nidim && j < _njdim && k < _nkdim, "Index out of bounds.");
    }
    body {
	return k * (_njdim * _nidim) + j * _nidim + i; 
    }

    size_t[] toIJKIndices(size_t gid) const
    {
	size_t k = gid / (_njdim * _nidim);
	size_t j = (gid - k * (_njdim * _nidim)) / _nidim;
	size_t i = gid - k * (_njdim * _nidim) - j * _nidim;
	return [i, j, k];
    }

    @nogc ref SolidFVCell getCell(size_t i, size_t j, size_t k=0) { return _ctr[toGlobalIndex(i,j,k)]; }
    @nogc ref SolidFVInterface getIfi(size_t i, size_t j, size_t k=0) { return _ifi[toGlobalIndex(i,j,k)]; }
    @nogc ref SolidFVInterface getIfj(size_t i, size_t j, size_t k=0) { return _ifj[toGlobalIndex(i,j,k)]; }
    @nogc ref SolidFVInterface getIfk(size_t i, size_t j, size_t k=0) { return _ifk[toGlobalIndex(i,j,k)]; }
    @nogc ref SolidFVVertex getVtx(size_t i, size_t j, size_t k=0) { return _vtx[toGlobalIndex(i,j,k)]; }
    @nogc ref SolidFVInterface getSifi(size_t i, size_t j, size_t k=0) { return _sifi[toGlobalIndex(i,j,k)]; }
    @nogc ref SolidFVInterface getSifj(size_t i, size_t j, size_t k=0) { return _sifj[toGlobalIndex(i,j,k)]; }
    @nogc ref SolidFVInterface getSifk(size_t i, size_t j, size_t k=0) { return _sifk[toGlobalIndex(i,j,k)]; }



    override void assembleArrays()
    {
	if ( GlobalConfig.verbosity_level >= 2 ) 
	    writefln("SSolidBlock.assembleArrays(): Begin for solid_block %d", id);
	// Check for obvious errors.
	if ( _nidim <= 0 || _njdim <= 0 || _nkdim <= 0 ) {
	    throw new Error(text("SSolidBlock.assembleArrays(): invalid dimensions nidim=",
				 _nidim, " njdim=", _njdim, " nkdim=", _nkdim));
	}
	size_t ntot = _nidim * _njdim * _nkdim;
	try {
	    // Create the cell and interface objects for the entire block.
	    foreach (gid; 0 .. ntot) {
		_ctr ~= new SolidFVCell(); _ctr[gid].id = to!int(gid);
		auto ijk = toIJKIndices(gid);
		if ( ijk[0] >= imin && ijk[0] <= imax && 
		     ijk[1] >= jmin && ijk[1] <= jmax && 
		     ijk[2] >= kmin && ijk[2] <= kmax ) {
		    activeCells ~= _ctr[gid];
		}
		_ifi ~= new SolidFVInterface(); _ifi[gid].id = gid;
		_ifj ~= new SolidFVInterface(); _ifj[gid].id = gid;
		if ( GlobalConfig.dimensions == 3 ) {
		    _ifk ~= new SolidFVInterface(); _ifk[gid].id = gid;
		}
		_vtx ~= new SolidFVVertex(); _vtx[gid].id = gid;
		_sifi ~= new SolidFVInterface(); _sifi[gid].id = gid;
		_sifj ~= new SolidFVInterface(); _sifj[gid].id = gid;
		if ( GlobalConfig.dimensions == 3 ) {
		    _sifk ~= new SolidFVInterface(); _sifk[gid].id = gid;
		}
	    } // gid loop
	} catch (Error err) {
	    writeln("Crapped out while assembling solid_block arrays.");
	    writefln("nicell=%d njcell=%d nkcell=%d", nicell, njcell, nkcell);
	    writefln("nidim=%d njdim=%d nkdim=%d", _nidim, _njdim, _nkdim);
	    writeln("Probably ran out of memory.");
	    writeln("Be a little less ambitious and try a smaller grid next time.");
	    writefln("System message: %s", err.msg);
	    throw new Error("SolidBlock.assembleArrays() failed.");
	}
	if ( GlobalConfig.verbosity_level >= 2 ) {
	    writefln("Done assembling arrays for %d solid cells.", ntot);
	}
    }
    override void bindFacesAndVerticesToCells()
    {
	size_t kstart, kend;
	if ( GlobalConfig.dimensions == 3 ) {
	    kstart = kmin - 1;
	    kend = kmax + 1;
	} else {
	    kstart = 0;
	    kend = 0;
	}
	// With these ranges, we also do the first layer of ghost cells.
	for ( size_t k = kstart; k <= kend; ++k ) {
	    for ( size_t j = jmin-1; j <= jmax+1; ++j ) {
		for ( size_t i = imin-1; i <= imax+1; ++i ) {
		    SolidFVCell cell = getCell(i,j,k);
		    cell.iface ~= getIfj(i,j+1,k); // north
		    cell.iface ~= getIfi(i+1,j,k); // east
		    cell.iface ~= getIfj(i,j,k); // south
		    cell.iface ~= getIfi(i,j,k); // west
		    cell.vtx ~= getVtx(i,j,k);
		    cell.vtx ~= getVtx(i+1,j,k);
		    cell.vtx ~= getVtx(i+1,j+1,k);
		    cell.vtx ~= getVtx(i,j+1,k);
		    if ( GlobalConfig.dimensions == 3 ) {
			cell.iface ~= getIfk(i,j,k+1); // top
			cell.iface ~= getIfk(i,j,k); // bottom
			cell.vtx ~= getVtx(i,j,k+1);
			cell.vtx ~= getVtx(i+1,j,k+1);
			cell.vtx ~= getVtx(i+1,j+1,k+1);
			cell.vtx ~= getVtx(i,j+1,k+1);
		    } // end if
		} // for i
	    } // for j
	} // for k
    }

    override void readGrid(string filename)
    {
	size_t nivtx, njvtx, nkvtx;
	double x, y, z;
	if ( GlobalConfig.verbosity_level >= 1 && id == 0 ) {
	    writeln("SSolidBlock.readGrid(): Start block ", id);
	}
	auto byLine = new GzipByLine(filename);
	auto line = byLine.front; byLine.popFront();
	formattedRead(line, "%d %d %d", &nivtx, &njvtx, &nkvtx);
	if ( GlobalConfig.dimensions == 3 ) {
	    if ( nivtx-1 != nicell || njvtx-1 != njcell || nkvtx-1 != nkcell ) {
		throw new Error(text("For solid_block[", id, "] we have a mismatch in 3D grid size.",
                                     " Have read nivtx=", nivtx, " njvtx=", njvtx,
				     " nkvtx=", nkvtx));
	    }
	    for ( size_t k = kmin; k <= kmax+1; ++k ) {
		for ( size_t j = jmin; j <= jmax+1; ++j ) {
		    for ( size_t i = imin; i <= imax+1; ++i ) {
			line = byLine.front; byLine.popFront();
			// Note that the line starts with whitespace.
			formattedRead(line, " %g %g %g", &x, &y, &z);
			auto vtx = getVtx(i,j,k);
			vtx.pos.refx = x;
			vtx.pos.refy = y;
			vtx.pos.refz = z;
		    } // for i
		} // for j
	    } // for k
	} else { // 2D case
	    if ( nivtx-1 != nicell || njvtx-1 != njcell || nkvtx != 1 ) {
		throw new Error(text("For solid_block[", id, "] we have a mismatch in 2D grid size.",
				     " Have read nivtx=", nivtx, " njvtx=", njvtx,
				     " nkvtx=", nkvtx));
	    }
	    for ( size_t j = jmin; j <= jmax+1; ++j ) {
		for ( size_t i = imin; i <= imax+1; ++i ) {
		    line = byLine.front; byLine.popFront();
		    // Note that the line starts with whitespace.
		    formattedRead(line, " %g %g", &x, &y);
		    auto vtx = getVtx(i,j);
		    vtx.pos.refx = x;
		    vtx.pos.refy = y;
		    vtx.pos.refz = 0.0;
		} // for i
	    } // for j
	}
    } // end read_grid()

    override void writeGrid(string filename, double sim_time)
    {
	if ( GlobalConfig.verbosity_level >= 1 && id == 0 ) {
	    writeln("SSolidBlock.writeGrid(): Start block ", id);
	}
	size_t kmaxrange;
	auto outfile = new GzipOut(filename);
	auto writer = appender!string();
	if ( GlobalConfig.dimensions == 3 ) {
	    formattedWrite(writer, "%d %d %d  # ni nj nk\n", nicell+1, njcell+1, nkcell+1);
	    kmaxrange = kmax + 1;
	} else { // 2D case
	    formattedWrite(writer, "%d %d %d  # ni nj nk\n", nicell+1, njcell+1, nkcell);
	    kmaxrange = kmax;
	}
	outfile.compress(writer.data);
	for ( size_t k = kmin; k <= kmaxrange; ++k ) {
	    for ( size_t j = jmin; j <= jmax+1; ++j ) {
		for ( size_t i = imin; i <= imax+1; ++i ) {
		    auto vtx = getVtx(i,j,k);
		    writer = appender!string();
		    formattedWrite(writer, "%20.12e %20.12e %20.12e\n",
				   vtx.pos.x, vtx.pos.y, vtx.pos.z);
		    outfile.compress(writer.data);
		} // for i
	    } // for j
	} // for k
	outfile.finish();
    } // end write_grid()

    override void readSolution(string filename)
    {
	size_t ni, nj, nk;
	double sim_time;
	if ( GlobalConfig.verbosity_level >= 1 && id == 0 ) {
	    writeln("read_solution(): Start solid_block ", id);
	}
	auto byLine = new GzipByLine(filename);
	auto line = byLine.front; byLine.popFront();
	formattedRead(line, " %g", &sim_time);
	line = byLine.front; byLine.popFront();
	// ignore second line; it should be just the names of the variables
	// [TODO] We should test the incoming strings against the current variable names.
	line = byLine.front; byLine.popFront();
	formattedRead(line, "%d %d %d", &ni, &nj, &nk);
	if ( ni != nicell || nj != njcell || 
	     nk != ((GlobalConfig.dimensions == 3) ? nkcell : 1) ) {
	    throw new Error(text("For solid_block[", id, "] we have a mismatch in solution size.",
				 " Have read ni=", ni, " nj=", nj, " nk=", nk));
	}	
	for ( size_t k = kmin; k <= kmax; ++k ) {
	    for ( size_t j = jmin; j <= jmax; ++j ) {
		for ( size_t i = imin; i <= imax; ++i ) {
		    line = byLine.front; byLine.popFront();
		    getCell(i,j,k).scanValuesFromString(line);
		} // for i
	    } // for j
	} // for k
    } // end read_solution()

    override void writeSolution(string fileName, double simTime)
    {
	if ( GlobalConfig.verbosity_level >= 1 && id == 0 ) {
	    writeln("write_solution(): Start solid_block ", id);
	}
	auto outfile = new GzipOut(fileName);
	auto writer = appender!string();
	formattedWrite(writer, "%20.12e\n", simTime);
	outfile.compress(writer.data);
	writer = appender!string();
	foreach(varname; varListForSolidCell()) {
	    formattedWrite(writer, " \"%s\"", varname);
	}
	formattedWrite(writer, "\n");
	outfile.compress(writer.data);
	writer = appender!string();
	formattedWrite(writer, "%d %d %d\n", nicell, njcell, nkcell);
	outfile.compress(writer.data);
	for ( size_t k = kmin; k <= kmax; ++k ) {
	    for ( size_t j = jmin; j <= jmax; ++j ) {
		for ( size_t i = imin; i <= imax; ++i ) {
		    outfile.compress(" " ~ getCell(i,j,k).writeValuesToString() ~ "\n");
		} // for i
	    } // for j
	} // for k
	outfile.finish();
    }
    

    override void computePrimaryCellGeometricData()
    {
	if ( GlobalConfig.dimensions == 2 ) {
	    calcVolumes2D();
	    calcFaces2D();
	    return;
	}
	throw new Error("SSolidBlock.computePrimaryCellGeometryData() not implemented for 3D.");

    } // end compute_primary_cell_geometric_data()

    override void computeSecondaryCellGeometricData()
    {
	if ( GlobalConfig.dimensions == 2 ) {
	    secondaryAreas2D();
	    return;
	}
	throw new Error("SSolidBlock.computeSecondaryCellGeometryData() not implemented for 3D.");

    }

    void calcVolumes2D()
    {
	size_t i, j;
	double xA, yA, xB, yB, xC, yC, xD, yD;
	double xN, yN, xS, yS, xE, yE, xW, yW;
	double vol, max_vol, min_vol, xyarea;
	double dx, dy, dxN, dyN, dxE, dyE;
	double lengthN, lengthE, length_max, length_min, length_cross;
	double max_aspect, aspect_ratio;
	SolidFVCell cell, source_cell, target_cell;

	// Cell layout
	// C-----B     3-----2
	// |     |     |     |
	// |  c  |     |  c  |
	// |     |     |     |
	// D-----A     0-----1
    
	max_vol = 0.0;
	min_vol = 1.0e30;    /* arbitrarily large */
	max_aspect = 0.0;
	for ( i = imin; i <= imax; ++i ) {
	    for ( j = jmin; j <= jmax; ++j ) {
		cell = getCell(i,j);
		// These are the corners.
		xA = cell.vtx[1].pos.x;
		yA = cell.vtx[1].pos.y;
		xB = cell.vtx[2].pos.x;
		yB = cell.vtx[2].pos.y;
		xC = cell.vtx[3].pos.x;
		yC = cell.vtx[3].pos.y;
		xD = cell.vtx[0].pos.x;
		yD = cell.vtx[0].pos.y;
		// Cell area in the (x,y)-plane.
		xyarea = 0.5 * ((xB + xA) * (yB - yA) + (xC + xB) * (yC - yB) +
				(xD + xC) * (yD - yC) + (xA + xD) * (yA - yD));
		// Cell Centroid.
		cell.pos.refx = 1.0 / (xyarea * 6.0) * 
		    ((yB - yA) * (xA * xA + xA * xB + xB * xB) + 
		     (yC - yB) * (xB * xB + xB * xC + xC * xC) +
		     (yD - yC) * (xC * xC + xC * xD + xD * xD) + 
		     (yA - yD) * (xD * xD + xD * xA + xA * xA));
		cell.pos.refy = -1.0 / (xyarea * 6.0) * 
		    ((xB - xA) * (yA * yA + yA * yB + yB * yB) + 
		     (xC - xB) * (yB * yB + yB * yC + yC * yC) +
		     (xD - xC) * (yC * yC + yC * yD + yD * yD) + 
		     (xA - xD) * (yD * yD + yD * yA + yA * yA));
		cell.pos.refz = 0.0;
		// Cell Volume.
		if ( GlobalConfig.axisymmetric ) {
		    // Volume per radian = centroid y-ordinate * cell area
		    vol = xyarea * cell.pos.y;
		} else {
		    // Assume unit depth in the z-direction.
		    vol = xyarea;
		}
		if (vol < 0.0) {
		    throw new Error(text("Negative cell volume: solid_block ", id,
					 " vol[", i, " ,", j, "]= ", vol));
		}
		if (vol > max_vol) max_vol = vol;
		if (vol < min_vol) min_vol = vol;
		cell.volume = vol;
		cell.areaxy = xyarea;
	    } // j loop
	} // i loop
    } // end calc_volumes_2D()

    void calcFaces2D()
    {
	SolidFVInterface iface;
	size_t i, j;
	double xA, xB, yA, yB, xC, yC;
	double LAB, LBC;

	// East-facing interfaces.
	for (i = imin; i <= imax+1; ++i) {
	    for (j = jmin; j <= jmax; ++j) {
		iface = getIfi(i,j);
		// These are the corners.
		xA = getVtx(i,j).pos.x; 
		yA = getVtx(i,j).pos.y;
		xB = getVtx(i,j+1).pos.x; 
		yB = getVtx(i,j+1).pos.y;
		LAB = sqrt((xB - xA) * (xB - xA) + (yB - yA) * (yB - yA));
		if (LAB < 1.0e-9) {
		    writefln("Zero length ifi[%d,%d]: %e", i, j, LAB);
		}
		// Direction cosines for the unit normal.
		iface.n.refx = (yB - yA) / LAB;
		iface.n.refy = -(xB - xA) / LAB;
		iface.n.refz = 0.0;  // 2D plane
		iface.t2 = Vector3(0.0, 0.0, 1.0);
		iface.t1 = cross(iface.n, iface.t2);
		// Length in the XY-plane.
		iface.length = LAB;
		// Mid-point and area.
		iface.Ybar = 0.5 * (yA + yB);
		if ( GlobalConfig.axisymmetric ) {
		    // Interface area per radian.
		    iface.area = LAB * iface.Ybar;
		} else {
		    // Assume unit depth in the Z-direction.
		    iface.area = LAB;
		}
		iface.pos = (getVtx(i,j).pos + getVtx(i,j+1).pos)/2.0;
	    } // j loop
	} // i loop
    
	// North-facing interfaces.
	for (i = imin; i <= imax; ++i) {
	    for (j = jmin; j <= jmax+1; ++j) {
		iface = getIfj(i,j);
		// These are the corners.
		xB = getVtx(i+1,j).pos.x;
		yB = getVtx(i+1,j).pos.y;
		xC = getVtx(i,j).pos.x;
		yC = getVtx(i,j).pos.y;
		LBC = sqrt((xC - xB) * (xC - xB) + (yC - yB) * (yC - yB));
		if (LBC < 1.0e-9) {
		    writefln("Zero length ifj[%d,%d]: %e", i, j, LBC);
		}
		// Direction cosines for the unit normal.
		iface.n.refx = (yC - yB) / LBC;
		iface.n.refy = -(xC - xB) / LBC;
		iface.n.refz = 0.0;  // 2D plane
		iface.t2 = Vector3(0.0, 0.0, 1.0);
		iface.t1 = cross(iface.n, iface.t2);
		// Length in the XY-plane.
		iface.length = LBC;
		// Mid-point and area.
		iface.Ybar = 0.5 * (yC + yB);
		if ( GlobalConfig.axisymmetric ) {
		    // Interface area per radian.
		    iface.area = LBC * iface.Ybar;
		} else {
		    // Assume unit depth in the Z-direction.
		    iface.area = LBC;
		}
		iface.pos = (getVtx(i+1,j).pos + getVtx(i,j).pos)/2.0;
	    } // j loop
	} // i loop
    } // end calc_faces_2D()

    void secondaryAreas2D()
    // Compute the secondary cell cell areas in the (x,y)-plane.
    //
    // The secondary cells are centred on the vertices of the 
    // primary cells and have primary cell centres as their corners.
    // For this particular secondary cell, centred on a vertex v (i,j),
    // the near-by primary cell i,j is centred on B'.
    //
    //          +-----+
    //          |     |
    //       C'-+--B' |
    //       |  |  |  |
    //       |  v--+--+
    //       |     |
    //       D'----A'
    //
    {
	size_t i, j;
	double xA, yA, xB, yB, xC, yC, xD, yD;
	double xyarea, max_area, min_area;

	max_area = 0.0;
	min_area = 1.0e6;   // arbitrarily large
	// First, do all of the internal secondary cells.
	// i.e. The ones centred on primary vertices which 
	// are not on a boundary.
	for (i = imin+1; i <= imax; ++i) {
	    for (j = jmin+1; j <= jmax; ++j) {
		// These are the corners.
		xA = getCell(i,j-1).pos.x;
		yA = getCell(i,j-1).pos.y;
		xB = getCell(i,j).pos.x;
		yB = getCell(i,j).pos.y;
		xC = getCell(i-1,j).pos.x;
		yC = getCell(i-1,j).pos.y;
		xD = getCell(i-1,j-1).pos.x;
		yD = getCell(i-1,j-1).pos.y;
		// Cell area in the (x,y)-plane.
		xyarea = 0.5 * ((xB + xA) * (yB - yA) + (xC + xB) * (yC - yB) +
				(xD + xC) * (yD - yC) + (xA + xD) * (yA - yD));
		if (xyarea < 0.0) {
		    throw new Error(text("Negative secondary-cell area: Block ", id,
					 " vtx[", i, " ,", j, "]= ", xyarea));
		}
		if (xyarea > max_area) max_area = xyarea;
		if (xyarea < min_area) min_area = xyarea;
		getVtx(i,j).areaxy = xyarea;
	    } // j loop
	} // i loop

	// Note that the secondary cells along block boundaries are HALF cells.
	//
	// East boundary.
	i = imax+1;
	for (j = jmin+1; j <= jmax; ++j) {
	    xA = getIfi(i,j-1).pos.x;
	    yA = getIfi(i,j-1).pos.y;
	    xB = getIfi(i,j).pos.x;
	    yB = getIfi(i,j).pos.y;
	    xC = getCell(i-1,j).pos.x;
	    yC = getCell(i-1,j).pos.y;
	    xD = getCell(i-1,j-1).pos.x;
	    yD = getCell(i-1,j-1).pos.y;
	    // Cell area in the (x,y)-plane.
	    xyarea = 0.5 * ((xB + xA) * (yB - yA) + (xC + xB) * (yC - yB) +
			    (xD + xC) * (yD - yC) + (xA + xD) * (yA - yD));
	    if (xyarea < 0.0) {
		throw new Error(text("Negative secondary-cell area: Block ", id,
				     " vtx[", i, " ,", j, "]= ", xyarea));
	    }
	    if (xyarea > max_area) max_area = xyarea;
	    if (xyarea < min_area) min_area = xyarea;
	    getVtx(i,j).areaxy = xyarea;
	} // j loop 

	// Fudge corners -- not expecting to use this data.
	getVtx(i,jmin).areaxy = 0.5 * getVtx(i,jmin+1).areaxy;
	getVtx(i,jmax+1).areaxy = 0.5 * getVtx(i,jmax).areaxy;
    
	// West boundary.
	i = imin;
	for (j = jmin+1; j <= jmax; ++j) {
	    xA = getCell(i,j-1).pos.x;
	    yA = getCell(i,j-1).pos.y;
	    xB = getCell(i,j).pos.x;
	    yB = getCell(i,j).pos.y;
	    xC = getIfi(i,j).pos.x;
	    yC = getIfi(i,j).pos.y;
	    xD = getIfi(i,j-1).pos.x;
	    yD = getIfi(i,j-1).pos.y;
	    // Cell area in the (x,y)-plane.
	    xyarea = 0.5 * ((xB + xA) * (yB - yA) + (xC + xB) * (yC - yB) +
			    (xD + xC) * (yD - yC) + (xA + xD) * (yA - yD));
	    if (xyarea < 0.0) {
		throw new Error(text("Negative secondary-cell area: Block ", id,
				     " vtx[", i, " ,", j, "]= ", xyarea));
	    }
	    if (xyarea > max_area) max_area = xyarea;
	    if (xyarea < min_area) min_area = xyarea;
	    getVtx(i,j).areaxy = xyarea;
	} // j loop 

	// Fudge corners.
	getVtx(i,jmin).areaxy = 0.5 * getVtx(i,jmin+1).areaxy;
	getVtx(i,jmax+1).areaxy = 0.5 * getVtx(i,jmax).areaxy;

	// North boundary.
	j = jmax+1;
	for (i = imin+1; i <= imax; ++i) {
	    // These are the corners.
	    xA = getCell(i,j-1).pos.x;
	    yA = getCell(i,j-1).pos.y;
	    xB = getIfj(i,j).pos.x;
	    yB = getIfj(i,j).pos.y;
	    xC = getIfj(i-1,j).pos.x;
	    yC = getIfj(i-1,j).pos.y;
	    xD = getCell(i-1,j-1).pos.x;
	    yD = getCell(i-1,j-1).pos.y;
	    // Cell area in the (x,y)-plane.
	    xyarea = 0.5 * ((xB + xA) * (yB - yA) + (xC + xB) * (yC - yB) +
			    (xD + xC) * (yD - yC) + (xA + xD) * (yA - yD));
	    if (xyarea < 0.0) {
		throw new Error(text("Negative secondary-cell area: Block ", id,
				     " vtx[", i, " ,", j, "]= ", xyarea));
	    }
	    if (xyarea > max_area) max_area = xyarea;
	    if (xyarea < min_area) min_area = xyarea;
	    getVtx(i,j).areaxy = xyarea;
	} // i loop 

	// Fudge corners.
	getVtx(imin,j).areaxy = 0.5 * getVtx(imin+1,j).areaxy;
	getVtx(imax+1,j).areaxy = 0.5 * getVtx(imax,j).areaxy;

	// South boundary.
	j = jmin;
	for (i = imin+1; i <= imax; ++i) {
	    xA = getIfj(i,j).pos.x;
	    yA = getIfj(i,j).pos.y;
	    xB = getCell(i,j).pos.x;
	    yB = getCell(i,j).pos.y;
	    xC = getCell(i-1,j).pos.x;
	    yC = getCell(i-1,j).pos.y;
	    xD = getIfj(i-1,j).pos.x;
	    yD = getIfj(i-1,j).pos.y;
	    // Cell area in the (x,y)-plane.
	    xyarea = 0.5 * ((xB + xA) * (yB - yA) + (xC + xB) * (yC - yB) +
			    (xD + xC) * (yD - yC) + (xA + xD) * (yA - yD));
	    if (xyarea < 0.0) {
		throw new Error(text("Negative secondary-cell area: Block ", id,
				     " vtx[", i, " ,", j, "]= ", xyarea));
	    }
	    if (xyarea > max_area) max_area = xyarea;
	    if (xyarea < min_area) min_area = xyarea;
	    getVtx(i,j).areaxy = xyarea;
	} // i loop

	// Fudge corners.
	getVtx(imin,j).areaxy = 0.5 * getVtx(imin+1,j).areaxy;
	getVtx(imax+1,j).areaxy = 0.5 * getVtx(imax,j).areaxy;

    } // end secondary_areas_2D()


    override void applyPreSpatialDerivAction(double t, int tLevel)
    {
	bc[Face.north].applyPreSpatialDerivAction(t, tLevel);
	bc[Face.east].applyPreSpatialDerivAction(t, tLevel);
	bc[Face.south].applyPreSpatialDerivAction(t, tLevel);
	bc[Face.west].applyPreSpatialDerivAction(t, tLevel);
	if ( GlobalConfig.dimensions == 3 ) {
	    bc[Face.top].applyPreSpatialDerivAction(t, tLevel);
	    bc[Face.bottom].applyPreSpatialDerivAction(t, tLevel);
	}
    }

    override void applyPostFluxAction(double t, int tLevel)
    {
	bc[Face.north].applyPostFluxAction(t, tLevel);
	bc[Face.east].applyPostFluxAction(t, tLevel);
	bc[Face.south].applyPostFluxAction(t, tLevel);
	bc[Face.west].applyPostFluxAction(t, tLevel);
	if ( GlobalConfig.dimensions == 3 ) {
	    bc[Face.top].applyPostFluxAction(t, tLevel);
	    bc[Face.bottom].applyPostFluxAction(t, tLevel);
	}
    }

    override void computeSpatialDerivatives(int ftl)
    {
	// NOTE: This presently uses the Eilmer3 2D formulation for
	// computing spatial derivatives. We might move this to the
	// 3D formulation at some point in the future.
	// [2015-25-04]

	if ( GlobalConfig.dimensions == 3 ) {
	    throw new Error("computeSpatialDerivatives() not implemented for 3D yet.");
	}
	size_t i, j;
	SolidFVVertex vtx;
	SolidFVCell cell;
	SolidFVInterface a, b;
	double xA, xB, xC, xD;
	double yA, yB, yC, yD;
	double TA, TB, TC, TD;
	double areaInv;

	// Work on all internal secondary cells.
	for ( j = jmin+1; j <= jmax; ++j ) {
	    for ( i = imin+1; i <= imax; ++i ) {
		vtx = getVtx(i,j);
		areaInv = 1.0/vtx.areaxy;
		// Corners of secondary cells
		xA = getCell(i,j-1).pos.x;
		yA = getCell(i,j-1).pos.y;
		xB = getCell(i,j).pos.x;
		yB = getCell(i,j).pos.y;
		xC = getCell(i-1,j).pos.x;
		yC = getCell(i-1,j).pos.y;
		xD = getCell(i-1,j-1).pos.x;
		yD = getCell(i-1,j-1).pos.y;
		// Temperature at the corners of the secondary cells
		TA = getCell(i,j-1).T[ftl];
		TB = getCell(i,j).T[ftl];
		TC = getCell(i-1,j).T[ftl];
		TD = getCell(i-1,j-1).T[ftl];
		// Compute derivative using Gauss-Green theorem
		vtx.dTdx = 0.5 * areaInv * 
		    ((TB + TA) * (yB - yA) + (TC + TB) * (yC - yB) +
		     (TD + TC) * (yD - yC) + (TA + TD) * (yA - yD)); 
		vtx.dTdy = -0.5 * areaInv *
		    ((TB + TA) * (xB - xA) + (TC + TB) * (xC - xB) +
		     (TD + TC) * (xD - xC) + (TA + TD) * (xA - xD));
	    } // end for j
	} // end for i

	// Next, work on the edges.
	// EAST boundary.
	i = imax + 1;
	for ( j = jmin+1; j <= jmax; ++j ) {
	    vtx = getVtx(i,j);
	    areaInv = 1.0 / vtx.areaxy;
	    // Corners for the secondary cell
	    xA = getIfi(i, j-1).pos.x;
	    yA = getIfi(i, j-1).pos.y;
	    xB = getIfi(i,j).pos.x;
	    yB = getIfi(i,j).pos.y;
	    xC = getCell(i-1,j).pos.x;
	    yC = getCell(i-1,j).pos.y;
	    xD = getCell(i-1,j-1).pos.x;
	    yD = getCell(i-1,j-1).pos.y;
	    // Temperatures
	    TA = getIfi(i, j-1).T;
	    TB = getIfi(i, j).T;
	    TC = getCell(i-1,j).T[ftl];
	    TD = getCell(i-1,j-1).T[ftl];
	    // Compute derivative using Gauss-Green theorem
	    vtx.dTdx = 0.5 * areaInv * 
		((TB + TA) * (yB - yA) + (TC + TB) * (yC - yB) +
		 (TD + TC) * (yD - yC) + (TA + TD) * (yA - yD)); 
	    vtx.dTdy = -0.5 * areaInv *
		((TB + TA) * (xB - xA) + (TC + TB) * (xC - xB) +
		 (TD + TC) * (xD - xC) + (TA + TD) * (xA - xD));

	}
	// WEST boundary
	i = imin;
	for ( j = jmin+1; j <= jmax; ++j ) {
	    vtx = getVtx(i, j);
	    areaInv = 1.0 / vtx.areaxy;
	    // Corners of the secondary cell
	    xA = getCell(i, j-1).pos.x;
	    yA = getCell(i, j-1).pos.y;
	    xB = getCell(i, j).pos.x;
	    yB = getCell(i, j).pos.y;
	    xC = getIfi(i, j).pos.x;
	    yC = getIfi(i, j).pos.y;
	    xD = getIfi(i, j-1).pos.x;
	    yD = getIfi(i, j-1).pos.y;
	    // Temperatures
	    TA = getCell(i, j-1).T[ftl];
	    TB = getCell(i, j).T[ftl];
	    TC = getIfi(i, j).T;
	    TD = getIfi(i, j-1).T;
	    // Compute derivative using Gauss-Green theorem
	    vtx.dTdx = 0.5 * areaInv * 
		    ((TB + TA) * (yB - yA) + (TC + TB) * (yC - yB) +
		     (TD + TC) * (yD - yC) + (TA + TD) * (yA - yD)); 
	    vtx.dTdy = -0.5 * areaInv *
		((TB + TA) * (xB - xA) + (TC + TB) * (xC - xB) +
		 (TD + TC) * (xD - xC) + (TA + TD) * (xA - xD));

	}
	// NORTH boundary
	j = jmax + 1;
	for ( i = imin+1; i <= imax; ++i ) {
	    vtx = getVtx(i, j);
	    areaInv = 1.0 / vtx.areaxy;
	    // Corners of the secondary cell
	    xA = getCell(i, j-1).pos.x;
	    yA = getCell(i, j-1).pos.y;
	    xB = getIfj(i, j).pos.x;
	    yB = getIfj(i, j).pos.y;
	    xC = getIfj(i-1, j).pos.x;
	    yC = getIfj(i-1, j).pos.y;
	    xD = getCell(i-1, j-1).pos.x;
	    yD = getCell(i-1, j-1).pos.y;
	    // Temperatures
	    TA = getCell(i, j-1).T[ftl];
	    TB = getIfj(i, j).T;
	    TC = getIfj(i-1, j).T;
	    TD = getCell(i-1, j-1).T[ftl];
	    // Compute derivative using Gauss-Green theorem
	    vtx.dTdx = 0.5 * areaInv * 
		((TB + TA) * (yB - yA) + (TC + TB) * (yC - yB) +
		 (TD + TC) * (yD - yC) + (TA + TD) * (yA - yD)); 
	    vtx.dTdy = -0.5 * areaInv *
		((TB + TA) * (xB - xA) + (TC + TB) * (xC - xB) +
		 (TD + TC) * (xD - xC) + (TA + TD) * (xA - xD));
	}
	// SOUTH boundary
	j = jmin;
	for ( i = imin+1; i <= imax; ++i ) {
	    vtx = getVtx(i, j);
	    areaInv = 1.0 / vtx.areaxy;
	    // Corners of the secondary cell
	    xA = getIfj(i, j).pos.x;
	    yA = getIfj(i, j).pos.y;
	    xB = getCell(i, j).pos.x;
	    yB = getCell(i, j).pos.y;
	    xC = getCell(i-1, j).pos.x;
	    yC = getCell(i-1, j).pos.y;
	    xD = getIfj(i-1, j).pos.x;
	    yD = getIfj(i-1, j).pos.y;
	    // Temperatures
	    TA = getIfj(i, j).T;
	    TB = getCell(i, j).T[ftl];
	    TC = getCell(i-1, j).T[ftl];
	    TD = getIfj(i-1, j).T;
	    // Compute derivative using Gauss-Green theorem
	    vtx.dTdx = 0.5 * areaInv * 
		    ((TB + TA) * (yB - yA) + (TC + TB) * (yC - yB) +
		     (TD + TC) * (yD - yC) + (TA + TD) * (yA - yD));
	    vtx.dTdy = -0.5 * areaInv *
		((TB + TA) * (xB - xA) + (TC + TB) * (xC - xB) +
		 (TD + TC) * (xD - xC) + (TA + TD) * (xA - xD));
	}

	// Finally, derivatives at corners
	// NORTH-EAST corner
	i = imax;
	j = jmax;
	vtx = getVtx(i+1, j+1);
	cell = getCell(i, j);
	a = getIfj(i, j+1);
	b = getIfi(i+1, j);
	xA = a.pos.x; yA = a.pos.y;
	xB = b.pos.x; yB = b.pos.y;
	xC = cell.pos.x; yC = cell.pos.y;
	double denom = (xC - xA)*(yB - yA) - (xB - xA)*(yC - yA);
	TA = a.T;
	TB = b.T;
	TC = cell.T[ftl];
	vtx.dTdx = ((TC-TA)*(yB-yA) - (TB-TA)*(yC-yA))/denom;
	vtx.dTdy = ((TB-TA)*(xC-xA) - (TC-TA)*(xB-xA))/denom;
	// SOUTH-EAST corner
	i = imax;
	j = jmin;
	vtx = getVtx(i+1, j);
	cell = getCell(i, j);
	a = getIfj(i, j);
	b = getIfi(i+1, j);
	xA = a.pos.x; yA = a.pos.y;
	xB = b.pos.x; yB = b.pos.y;
	xC = cell.pos.x; yC = cell.pos.y;
	denom = (xC - xA)*(yB - yA) - (xB - xA)*(yC - yA);
	TA = a.T;
	TB = b.T;
	TC = cell.T[ftl];
	vtx.dTdx = ((TC-TA)*(yB-yA) - (TB-TA)*(yC-yA))/denom;
	vtx.dTdy = ((TB-TA)*(xC-xA) - (TC-TA)*(xB-xA))/denom;
	// SOUTH-WEST corner
	i = imin;
	j = jmin;
	vtx = getVtx(i, j);
	cell = getCell(i, j);
	a = getIfj(i, j);
	b = getIfi(i, j);
	xA = a.pos.x; yA = a.pos.y;
	xB = b.pos.x; yB = b.pos.y;
	xC = cell.pos.x; yC = cell.pos.y;
	denom = (xC - xA)*(yB - yA) - (xB - xA)*(yC - yA);
	TA = a.T;
	TB = b.T;
	TC = cell.T[ftl];
	vtx.dTdx = ((TC-TA)*(yB-yA) - (TB-TA)*(yC-yA))/denom;
	vtx.dTdy = ((TB-TA)*(xC-xA) - (TC-TA)*(xB-xA))/denom;
	// NORTH-WEST corner
	i = imin;
	j = jmax;
	vtx = getVtx(i, j+1);
	cell = getCell(i, j);
	a = getIfj(i, j+1);
	b = getIfi(i, j);
	xA = a.pos.x; yA = a.pos.y;
	xB = b.pos.x; yB = b.pos.y;
	xC = cell.pos.x; yC = cell.pos.y;
	denom = (xC - xA)*(yB - yA) - (xB - xA)*(yC - yA);
	TA = a.T;
	TB = b.T;
	TC = cell.T[ftl];
	vtx.dTdx = ((TC-TA)*(yB-yA) - (TB-TA)*(yC-yA))/denom;
	vtx.dTdy = ((TB-TA)*(xC-xA) - (TC-TA)*(xB-xA))/denom;
    }

    override void computeFluxes()
    {
	size_t i, j, k;
	SolidFVInterface IFace;
	SolidFVVertex vtx1, vtx2;
	double dTdx, dTdy;
	double qx, qy;
	// East-facing interfaces
	for ( j = jmin; j <= jmax; ++j ) {
	    for ( i = imin; i <= imax + 1; ++i ) {
		if ( i == imin && bc[Face.west].setsFluxDirectly )
		    continue;
		if ( i == imax && bc[Face.east].setsFluxDirectly )
		    continue;
		IFace = getIfi(i, j);
		vtx1 = getVtx(i, j+1);
		vtx2 = getVtx(i, j);
		dTdx = 0.5*(vtx1.dTdx + vtx2.dTdx);
		dTdy = 0.5*(vtx1.dTdy + vtx2.dTdy);
		qx = -sp.k * dTdx;
		qy = -sp.k * dTdy;
		IFace.flux = qx * IFace.n.x + qy * IFace.n.y;
	    }
	}
	// North-facing interfaces
	for ( j = jmin; j <= jmax + 1; ++j ) {
	    for ( i = imin; i <= imax; ++i ) {
		if ( j == jmin && bc[Face.south].setsFluxDirectly ) 
		    continue;
		if ( j == jmax && bc[Face.north].setsFluxDirectly )
		    continue;
		IFace = getIfj(i, j);
		vtx1 = getVtx(i, j);
		vtx2 = getVtx(i+1, j);
		dTdx = 0.5*(vtx1.dTdx + vtx2.dTdx);
		dTdy = 0.5*(vtx1.dTdy + vtx2.dTdy);
		qx = -sp.k * dTdx;
		qy = -sp.k * dTdy;
		IFace.flux = qx * IFace.n.x + qy * IFace.n.y;
	    }
	}
    }

    override void clearSources()
    {
	foreach ( cell; activeCells ) cell.Q = 0.0;
    }
}
