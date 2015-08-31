// sblock.d
// Class for structured blocks of cells, for use within Eilmer4.
// This is the "classic" block within the mbcns/Eilmer series 
// of flow simulation codes.

// Peter J. 2014-07-20 first cut.

module sblock;

import std.conv;
import std.file;
import std.json;
import std.stdio;
import std.format;
import std.string;
import std.array;
import std.math;
import std.json;

import util.lua;
import json_helper;
import lua_helper;
import gzip;
import geom;
import gas;
import kinetics;
import globalconfig;
import globaldata;
import flowstate;
import fluxcalc;
import viscousflux;
import fvcore;
import fvvertex;
import fvinterface;
import fvcell;
import onedinterp;
import block;
import bc;

class SBlock: Block {
public:
    size_t nicell;
    size_t njcell;
    size_t nkcell;
    size_t imin, imax;
    size_t jmin, jmax;
    size_t kmin, kmax;
    size_t[] hicell, hjcell, hkcell; // locations of sample cells for history record
    size_t[] micell, mjcell, mkcell; // locations of monitor cells
    BoundaryCondition[6] bc;

private:
    // Total number of cells in each direction for this block.
    // these will be used in the array allocation routines.
    size_t _nidim;
    size_t _njdim;
    size_t _nkdim;
    // Most of the data is stored in the following arrays.
    // ctr = cell center values
    // ifi = east-facing face properties and fluxes (unit normal in the i-index direction)
    // ifj = north-facing face properties and fluxes (normal in the j-index direction)
    // ifk = top-facing 
    // vtx = cell vertex values (used for the viscous terms, mostly)
    // sifi, sifj and sifk are secondary-cell faces (whose corner nodes are the
    //                     the primary-cell centres.
    FVCell[] _ctr;
    FVInterface[] _ifi;
    FVInterface[] _ifj;
    FVInterface[] _ifk;
    FVVertex[] _vtx;
    FVInterface[] _sifi;
    FVInterface[] _sifj;
    FVInterface[] _sifk;
    // Work-space that gets reused.
    // The following objects are used in the convective_flux method.
    FlowState Lft;
    FlowState Rght;
    OneDInterpolator one_d;

public:
    this(int id, size_t nicell, size_t njcell, size_t nkcell, string label)
    {
	super(id, label);
	this.nicell = nicell;
	this.njcell = njcell;
	this.nkcell = nkcell;
	// Fill in other data sizes.
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
	// Workspace for flux_calc method.
	Lft = new FlowState(dedicatedConfig[id].gmodel);
	Rght = new FlowState(dedicatedConfig[id].gmodel);
	one_d = new OneDInterpolator(dedicatedConfig[id]);
    } // end constructor

    this(in int id, JSONValue json_data)
    {
	nicell = getJSONint(json_data, "nic", 0);
	njcell = getJSONint(json_data, "njc", 0);
	nkcell = getJSONint(json_data, "nkc", 0);
	label = getJSONstring(json_data, "label", "");
	this(id, nicell, njcell, nkcell, label);
	active = getJSONbool(json_data, "active", true);
	omegaz = getJSONdouble(json_data, "omegaz", 0.0);
    } // end constructor from json

    override void init_lua_globals()
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
	lua_pushcfunction(myL, &luafn_sampleFlow); lua_setglobal(myL, "sampleFlow"); // will work only in serial mode
    } // end init_lua_globals()

    override void init_boundary_conditions(JSONValue json_data)
    // Initialize boundary conditions after the blocks are fully constructed,
    // because we want access to the full collection of valid block references.
    {
	foreach (boundary; 0 .. (myConfig.dimensions == 3 ? 6 : 4)) {
	    string json_key = "face_" ~ face_name[boundary];
	    auto bc_json_data = json_data[json_key];
	    bc[boundary] = make_BC_from_json(bc_json_data, id, boundary,
					     nicell, njcell, nkcell);
	}
    } // end init_boundary_conditions()

    override string toString() const
    {
	char[] repr;
	repr ~= "SBlock(";
	repr ~= "id=" ~ to!string(id);
	repr ~= " label=\"" ~ label ~ "\"";
	repr ~= ", active=" ~ to!string(active);
	repr ~= ", omegaz=" ~ to!string(omegaz);
	repr ~= ", nicell=" ~ to!string(nicell);
	repr ~= ", njcell=" ~ to!string(njcell);
	repr ~= ", nkcell=" ~ to!string(nkcell);
	repr ~= ", \n    bc=["~ face_name[0] ~ "=" ~ to!string(bc[0]);
	foreach (i; 1 .. (myConfig.dimensions == 3 ? 6 : 4)) {
	    repr ~= ",\n        " ~ face_name[i] ~ "=" ~ to!string(bc[i]);
	}
	repr ~= "\n       ]"; // end bc list
	repr ~= ")";
	return to!string(repr);
    } // end toString()

    @nogc
    size_t to_global_index(size_t i, size_t j, size_t k) const
    in {
	assert(i < _nidim && j < _njdim && k < _nkdim, "Index out of bounds.");
    }
    body {
	return k * (_njdim * _nidim) + j * _nidim + i; 
    } // end to_global_index()

    size_t[] to_ijk_indices(size_t gid) const
    {
	size_t k = gid / (_njdim * _nidim);
	size_t j = (gid - k * (_njdim * _nidim)) / _nidim;
	size_t i = gid - k * (_njdim * _nidim) - j * _nidim;
	return [i, j, k];
    } // end to_ijk_indices()

    @nogc ref FVCell get_cell(size_t i, size_t j, size_t k=0) { return _ctr[to_global_index(i,j,k)]; }
    @nogc ref FVInterface get_ifi(size_t i, size_t j, size_t k=0) { return _ifi[to_global_index(i,j,k)]; }
    @nogc ref FVInterface get_ifj(size_t i, size_t j, size_t k=0) { return _ifj[to_global_index(i,j,k)]; }
    @nogc ref FVInterface get_ifk(size_t i, size_t j, size_t k=0) { return _ifk[to_global_index(i,j,k)]; }
    @nogc ref FVVertex get_vtx(size_t i, size_t j, size_t k=0) { return _vtx[to_global_index(i,j,k)]; }

    override void assemble_arrays()
    // We shouldn't be calling this until the essential bits of the GlobalConfig
    // have been set up.
    {
	if ( myConfig.verbosity_level >= 2 ) 
	    writefln("assemble_arrays(): Begin for block %d", id);
	// Check for obvious errors.
	if ( _nidim <= 0 || _njdim <= 0 || _nkdim <= 0 ) {
	    throw new Error(text("resize_arrays(): invalid dimensions nidim=",
				 _nidim, " njdim=", _njdim, " nkdim=", _nkdim));
	}
	size_t ntot = _nidim * _njdim * _nkdim;
	try {
	    // Create the cell and interface objects for the entire block.
	    foreach (gid; 0 .. ntot) {
		_ctr ~= new FVCell(myConfig); _ctr[gid].id = to!int(gid);
		auto ijk = to_ijk_indices(gid);
		if ( ijk[0] >= imin && ijk[0] <= imax && 
		     ijk[1] >= jmin && ijk[1] <= jmax && 
		     ijk[2] >= kmin && ijk[2] <= kmax ) {
		    active_cells ~= _ctr[gid];
		}
		_ifi ~= new FVInterface(myConfig.gmodel); _ifi[gid].id = gid;
		if ( ijk[0] >= imin && ijk[0] <= imax+1 && 
		     ijk[1] >= jmin && ijk[1] <= jmax && 
		     ijk[2] >= kmin && ijk[2] <= kmax ) {
		    active_ifaces ~= _ifi[gid];
		}
		_ifj ~= new FVInterface(myConfig.gmodel); _ifj[gid].id = gid;
		if ( ijk[0] >= imin && ijk[0] <= imax && 
		     ijk[1] >= jmin && ijk[1] <= jmax+1 && 
		     ijk[2] >= kmin && ijk[2] <= kmax ) {
		    active_ifaces ~= _ifj[gid];
		}
		if ( myConfig.dimensions == 3 ) {
		    _ifk ~= new FVInterface(myConfig.gmodel); _ifk[gid].id = gid;
		    if ( ijk[0] >= imin && ijk[0] <= imax && 
			 ijk[1] >= jmin && ijk[1] <= jmax && 
			 ijk[2] >= kmin && ijk[2] <= kmax+1 ) {
			active_ifaces ~= _ifk[gid];
		    }
		}
		_vtx ~= new FVVertex(myConfig.gmodel); _vtx[gid].id = gid;
	    } // gid loop
	} catch (Error err) {
	    writeln("Crapped out while assembling block arrays.");
	    writefln("nicell=%d njcell=%d nkcell=%d", nicell, njcell, nkcell);
	    writefln("nidim=%d njdim=%d nkdim=%d", _nidim, _njdim, _nkdim);
	    writeln("Probably ran out of memory.");
	    writeln("Be a little less ambitious and try a smaller grid next time.");
	    writefln("System message: %s", err.msg);
	    throw new Error("Block.assemble_arrays() failed.");
	}
	if ( myConfig.verbosity_level >= 2 ) {
	    writefln("Done assembling arrays for %d cells.", ntot);
	}
    } // end of assemble_arrays()

    override void bind_interfaces_and_vertices_to_cells()
    // There is a fixed order of faces and vertices for each cell.
    // Refer to fvcore.d
    {
	size_t kstart, kend;
	if ( myConfig.dimensions == 3 ) {
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
		    FVCell cell = get_cell(i,j,k);
		    cell.iface ~= get_ifj(i,j+1,k); // north
		    cell.iface ~= get_ifi(i+1,j,k); // east
		    cell.iface ~= get_ifj(i,j,k); // south
		    cell.iface ~= get_ifi(i,j,k); // west
		    cell.vtx ~= get_vtx(i,j,k);
		    cell.vtx ~= get_vtx(i+1,j,k);
		    cell.vtx ~= get_vtx(i+1,j+1,k);
		    cell.vtx ~= get_vtx(i,j+1,k);
		    if ( myConfig.dimensions == 3 ) {
			cell.iface ~= get_ifk(i,j,k+1); // top
			cell.iface ~= get_ifk(i,j,k); // bottom
			cell.vtx ~= get_vtx(i,j,k+1);
			cell.vtx ~= get_vtx(i+1,j,k+1);
			cell.vtx ~= get_vtx(i+1,j+1,k+1);
			cell.vtx ~= get_vtx(i,j+1,k+1);
		    } // end if
		} // for i
	    } // for j
	} // for k
    } // end bind_faces_and_vertices_to_cells()

    override void bind_vertices_and_cells_to_interfaces()
    // Sometimes it is convenient for an interface to come complete 
    // with information about the vertices that define it and the cells
    // that adjoin it.
    {
	// ifi interfaces are East-facing interfaces.
	// In 2D, vtx0==p11, vtx1==p10.
	// In 3D, the cycle [vtx0,vtx1,vtx2,vtx3] progresses counter-clockwise around 
	// the periphery of the face when the normal unit vector is pointing toward you.
	// t1 vector aligned with j-index direction
	// t2 vector aligned with k-index direction
	// The i,j,k indices are effectively cell indices in the following loops.
	for ( size_t k = kmin; k <= kmax; ++k ) {
	    for ( size_t j = jmin; j <= jmax; ++j ) {
		for ( size_t i = imin-1; i <= imax; ++i ) {
		    auto IFace = get_ifi(i+1,j,k);
		    if (myConfig.dimensions == 3) {
			IFace.vtx ~= get_vtx(i+1,j,k);
			IFace.vtx ~= get_vtx(i+1,j+1,k);
			IFace.vtx ~= get_vtx(i+1,j+1,k+1);
			IFace.vtx ~= get_vtx(i+1,j,k+1);
			IFace.left_cell = get_cell(i,j,k);
			IFace.right_cell = get_cell(i+1,j,k);
		    } else {
			IFace.vtx ~= get_vtx(i+1,j+1);
			IFace.vtx ~= get_vtx(i+1,j);
			IFace.left_cell = get_cell(i,j);
			IFace.right_cell = get_cell(i+1,j);
		    }
		} // i loop
	    } // j loop
	} // for k
	// ifj interfaces are North-facing interfaces.
	// In 2D, vtx0==p01, vtx1==p11.
	// t1 vector aligned with k-index direction
	// t2 vector aligned with i-index direction
	for ( size_t k = kmin; k <= kmax; ++k ) {
	    for ( size_t i = imin; i <= imax; ++i ) {
		for ( size_t j = jmin-1; j <= jmax; ++j ) {
		    auto IFace = get_ifj(i,j+1,k);
		    if (myConfig.dimensions == 3) {
			IFace.vtx ~= get_vtx(i,j+1,k);
			IFace.vtx ~= get_vtx(i,j+1,k+1);
			IFace.vtx ~= get_vtx(i+1,j+1,k+1);
			IFace.vtx ~= get_vtx(i+1,j+1,k);
			IFace.left_cell = get_cell(i,j,k);
			IFace.right_cell = get_cell(i,j+1,k);
		    } else {
			IFace.vtx ~= get_vtx(i,j+1);
			IFace.vtx ~= get_vtx(i+1,j+1);
			IFace.left_cell = get_cell(i,j);
			IFace.right_cell = get_cell(i,j+1);
		    }
		} // j loop
	    } // i loop
	} // for k
	if (myConfig.dimensions == 2) return;
	// ifk interfaces are Top-facing interfaces.
	// t1 vector aligned with i-index direction
	// t2 vector aligned with j-index direction
	for ( size_t i = imin; i <= imax; ++i ) {
	    for ( size_t j = jmin; j <= jmax; ++j ) {
		for ( size_t k = kmin-1; k <= kmax; ++k ) {
		    auto IFace = get_ifk(i,j,k+1);
		    IFace.vtx ~= get_vtx(i,j,k+1);
		    IFace.vtx ~= get_vtx(i+1,j,k+1);
		    IFace.vtx ~= get_vtx(i+1,j+1,k+1);
		    IFace.vtx ~= get_vtx(i,j+1,k+1);
		    IFace.left_cell = get_cell(i,j,k);
		    IFace.right_cell = get_cell(i,j,k+1);
		} // for k 
	    } // j loop
	} // i loop
	return;
    } // end bind_vertices_and_cells_to_interfaces()

    override int count_invalid_cells(int gtl)
    // Returns the number of cells that contain invalid data,
    // optionally patching bad cell data as it goes.
    //
    // This data can be identified by the density of internal energy 
    // being on the minimum limit or the velocity being very large.
    {
	int number_of_invalid_cells = 0;
	foreach(FVCell cell; active_cells) {
	    if ( cell.check_flow_data() == false ) {
		++number_of_invalid_cells;
		auto ijk = to_ijk_indices(cell.id);
		size_t i = ijk[0]; size_t j = ijk[1]; size_t k = ijk[2];
		writefln("count_invalid_cells: block_id = %d, cell[%d,%d,%d]\n", id, i, j, k);
		writeln(cell);
		if ( myConfig.adjust_invalid_cell_data ) {
		    // We shall set the cell data to something that
		    // is valid (and self consistent).
		    FVCell other_cell;
		    FVCell[] neighbours;
		    printf( "Adjusting cell data to a local average.\n" );
		    other_cell = get_cell(i-1,j,k);
		    if ( other_cell.check_flow_data() ) neighbours ~= other_cell;
		    other_cell = get_cell(i+1,j,k);
		    if ( other_cell.check_flow_data() ) neighbours ~= other_cell;
		    other_cell = get_cell(i,j-1,k);
		    if ( other_cell.check_flow_data() ) neighbours ~= other_cell;
		    other_cell = get_cell(i,j+1,k);
		    if ( other_cell.check_flow_data() ) neighbours ~= other_cell;
		    if ( myConfig.dimensions == 3 ) {
			other_cell = get_cell(i,j,k-1);
			if ( other_cell.check_flow_data() ) neighbours ~= other_cell;
			other_cell = get_cell(i,j,k+1);
			if ( other_cell.check_flow_data() ) neighbours ~= other_cell;
		    }
		    if ( neighbours.length == 0 ) {
			throw new Error(text("Block::count_invalid_cells(): "
					     "There were no valid neighbours to replace cell data."));
		    }
		    cell.replace_flow_data_with_average(neighbours);
		    cell.encode_conserved(gtl, 0, omegaz);
		    cell.decode_conserved(gtl, 0, omegaz);
		    writefln("after flow-data replacement: block_id = %d, cell[%d,%d,%d]\n",
			     id, i, j, k);
		    writeln(cell);
		} // end adjust_invalid_cell_data 
	    } // end of if invalid data...
	} // foreach cell
	return number_of_invalid_cells;
    } // end count_invalid_cells()

    override void compute_primary_cell_geometric_data(int gtl)
    // Compute cell and interface geometric properties.
    {
	size_t i, j, k;
	Vector3 dummy;
	Vector3 ds;
	if ( myConfig.dimensions == 2 ) {
	    calc_volumes_2D(gtl);
	    calc_faces_2D(gtl);
	    calc_ghost_cell_geom_2D(gtl);
	    return;
	}
	// Cell properties of volume and position.
	// Estimates of cross-cell distances for use in high-order reconstruction.
	for ( i = imin; i <= imax; ++i ) {
	    for ( j = jmin; j <= jmax; ++j ) {
		for ( k = kmin; k <= kmax; ++k ) {
		    auto cell = get_cell(i,j,k);
		    auto p0 = get_vtx(i,j,k).pos[gtl];
		    auto p1 = get_vtx(i+1,j,k).pos[gtl];
		    auto p2 = get_vtx(i+1,j+1,k).pos[gtl];
		    auto p3 = get_vtx(i,j+1,k).pos[gtl];
		    auto p4 = get_vtx(i,j,k+1).pos[gtl];
		    auto p5 = get_vtx(i+1,j,k+1).pos[gtl];
		    auto p6 = get_vtx(i+1,j+1,k+1).pos[gtl];
		    auto p7 = get_vtx(i,j+1,k+1).pos[gtl];
		    hex_cell_properties(p0, p1, p2, p3, p4, p5, p6, p7, 
					cell.pos[gtl], cell.volume[gtl], cell.iLength,
					cell.jLength, cell.kLength);
		    cell.L_min = cell.iLength;
		    if ( cell.jLength < cell.L_min ) cell.L_min = cell.jLength;
		    if ( cell.kLength < cell.L_min ) cell.L_min = cell.kLength;
		}
	    }
	}
	// work on ifi face as a WEST face
	// t1 in the j-ordinate direction
	// t2 in the k-ordinate direction
	for ( i = imin; i <= imax + 1; ++i ) {
	    for ( j = jmin; j <= jmax; ++j ) {
		for ( k = kmin; k <= kmax; ++k ) {
		    auto iface = get_ifi(i,j,k);
		    auto p0 = get_vtx(i,j,k).pos[gtl];
		    auto p3 = get_vtx(i,j+1,k).pos[gtl];
		    auto p7 = get_vtx(i,j+1,k+1).pos[gtl];
		    auto p4 = get_vtx(i,j,k+1).pos[gtl];
		    quad_properties(p0, p3, p7, p4,
				    iface.pos, iface.n, iface.t1, iface.t2,
				    iface.area[gtl]);
		}
	    }
	}
	// work on ifj face as a SOUTH face
	// t1 in the k-ordinate direction
	// t2 in the i-ordinate direction
	for ( i = imin; i <= imax; ++i ) {
	    for ( j = jmin; j <= jmax + 1; ++j ) {
		for ( k = kmin; k <= kmax; ++k ) {
		    auto iface = get_ifj(i,j,k);
		    auto p0 = get_vtx(i,j,k).pos[gtl];
		    auto p4 = get_vtx(i,j,k+1).pos[gtl];
		    auto p5 = get_vtx(i+1,j,k+1).pos[gtl];
		    auto p1 = get_vtx(i+1,j,k).pos[gtl];
		    quad_properties(p0, p4, p5, p1,
				    iface.pos, iface.n, iface.t1, iface.t2,
				    iface.area[gtl]);
		}
	    }
	}
	// work on ifk face as a BOTTOM face
	// t1 in the i-ordinate direction
	// t2 in the j-ordinate direction
	for ( i = imin; i <= imax; ++i ) {
	    for ( j = jmin; j <= jmax; ++j ) {
		for ( k = kmin; k <= kmax + 1; ++k ) {
		    auto iface = get_ifk(i,j,k);
		    auto p0 = get_vtx(i,j,k).pos[gtl];
		    auto p1 = get_vtx(i+1,j,k).pos[gtl];
		    auto p2 = get_vtx(i+1,j+1,k).pos[gtl];
		    auto p3 = get_vtx(i,j+1,k).pos[gtl];
		    quad_properties(p0, p1, p2, p3,
				    iface.pos, iface.n, iface.t1, iface.t2,
				    iface.area[gtl]);
		}
	    }
	}
	// Propagate cross-cell lengths into the ghost cells.
	// 25-Feb-2014
	// Jason Qin and Paul Petrie-Repar have identified the lack of exact symmetry in
	// the reconstruction process at the wall as being a cause of the leaky wall
	// boundary conditions.  Note that the symmetry is not consistent with the 
	// linear extrapolation used for the positions and volumes in the next section.
	// [TODO] -- think about this carefully.
	auto option = CopyDataOption.cell_lengths_only;
	for ( j = jmin; j <= jmax; ++j ) {
	    for ( k = kmin; k <= kmax; ++k ) {
		i = imin;
		get_cell(i-1,j,k).copy_values_from(get_cell(i,j,k), option);
		get_cell(i-2,j,k).copy_values_from(get_cell(i+1,j,k), option);
		i = imax;
		get_cell(i+1,j,k).copy_values_from(get_cell(i,j,k), option);
		get_cell(i+2,j,k).copy_values_from(get_cell(i-1,j,k), option);
	    }
	}
	for ( i = imin; i <= imax; ++i ) {
	    for ( k = kmin; k <= kmax; ++k ) {
		j = jmin;
		get_cell(i,j-1,k).copy_values_from(get_cell(i,j,k), option);
		get_cell(i,j-2,k).copy_values_from(get_cell(i,j+1,k), option);
		j = jmax;
		get_cell(i,j+1,k).copy_values_from(get_cell(i,j,k), option);
		get_cell(i,j+2,k).copy_values_from(get_cell(i,j-1,k), option);
	    }
	}
	for ( i = imin; i <= imax; ++i ) {
	    for ( j = jmin; j <= jmax; ++j ) {
		k = kmin;
		get_cell(i,j,k-1).copy_values_from(get_cell(i,j,k), option);
		get_cell(i,j,k-2).copy_values_from(get_cell(i,j,k+1), option);
		k = kmax;
		get_cell(i,j,k+1).copy_values_from(get_cell(i,j,k), option);
		get_cell(i,j,k+2).copy_values_from(get_cell(i,j,k-1), option);
	    }
	}
	/* Extrapolate (with first-order) cell positions and volumes to ghost cells. */
	// TODO -- think about how to make these things consistent.
	for ( j = jmin; j <= jmax; ++j ) {
	    for ( k = kmin; k <= kmax; ++k ) {
		i = imin;
		auto cell_1 = get_cell(i,j,k);
		auto cell_2 = get_cell(i+1,j,k);
		auto ghost_cell = get_cell(i-1,j,k);
		ghost_cell.pos[gtl] = 2.0*cell_1.pos[gtl] - cell_2.pos[gtl];
		ghost_cell.volume[gtl] = 2.0*cell_1.volume[gtl] - cell_2.volume[gtl];
		cell_2 = cell_1;
		cell_1 = ghost_cell;
		ghost_cell = get_cell(i-2,j,k);
		ghost_cell.pos[gtl] = 2.0*cell_1.pos[gtl] - cell_2.pos[gtl];
		ghost_cell.volume[gtl] = 2.0*cell_2.volume[gtl] - cell_2.volume[gtl];
		i = imax;
		cell_1 = get_cell(i,j,k);
		cell_2 = get_cell(i-1,j,k);
		ghost_cell = get_cell(i+1,j,k);
		ghost_cell.pos[gtl] = 2.0*cell_1.pos[gtl] - cell_2.pos[gtl];
		ghost_cell.volume[gtl] = 2.0*cell_1.volume[gtl] - cell_2.volume[gtl];
		cell_2 = cell_1;
		cell_1 = ghost_cell;
		ghost_cell = get_cell(i+2,j,k);
		ghost_cell.pos[gtl] = 2.0*cell_1.pos[gtl] - cell_2.pos[gtl];
		ghost_cell.volume[gtl] = 2.0*cell_1.volume[gtl] - cell_2.volume[gtl];
	    }
	}
	for ( i = imin; i <= imax; ++i ) {
	    for ( k = kmin; k <= kmax; ++k ) {
		j = jmin;
		auto cell_1 = get_cell(i,j,k);
		auto cell_2 = get_cell(i,j+1,k);
		auto ghost_cell = get_cell(i,j-1,k);
		ghost_cell.pos[gtl] = 2.0*cell_1.pos[gtl] - cell_2.pos[gtl];
		ghost_cell.volume[gtl] = 2.0*cell_1.volume[gtl] - cell_2.volume[gtl];
		cell_2 = cell_1;
		cell_1 = ghost_cell;
		ghost_cell = get_cell(i,j-2,k);
		ghost_cell.pos[gtl] = 2.0*cell_1.pos[gtl] - cell_2.pos[gtl];
		ghost_cell.volume[gtl] = 2.0*cell_1.volume[gtl] - cell_2.volume[gtl];
		j = jmax;
		cell_1 = get_cell(i,j,k);
		cell_2 = get_cell(i,j-1,k);
		ghost_cell = get_cell(i,j+1,k);
		ghost_cell.pos[gtl] = 2.0*cell_1.pos[gtl] - cell_2.pos[gtl];
		ghost_cell.volume[gtl] = 2.0*cell_1.volume[gtl] - cell_2.volume[gtl];
		cell_2 = cell_1;
		cell_1 = ghost_cell;
		ghost_cell = get_cell(i,j+2,k);
		ghost_cell.pos[gtl] = 2.0*cell_1.pos[gtl] - cell_2.pos[gtl];
		ghost_cell.volume[gtl] = 2.0*cell_1.volume[gtl] - cell_2.volume[gtl];
	    }
	}
	for ( i = imin; i <= imax; ++i ) {
	    for ( j = jmin; j <= jmax; ++j ) {
		k = kmin;
		auto cell_1 = get_cell(i,j,k);
		auto cell_2 = get_cell(i,j,k+1);
		auto ghost_cell = get_cell(i,j,k-1);
		ghost_cell.pos[gtl] = 2.0*cell_1.pos[gtl] - cell_2.pos[gtl];
		ghost_cell.volume[gtl] = 2.0*cell_1.volume[gtl] - cell_2.volume[gtl];
		cell_2 = cell_1;
		cell_1 = ghost_cell;
		ghost_cell = get_cell(i,j,k-2);
		ghost_cell.pos[gtl] = 2.0*cell_1.pos[gtl] - cell_2.pos[gtl];
		ghost_cell.volume[gtl] = 2.0*cell_1.volume[gtl] - cell_2.volume[gtl];
		k = kmax;
		cell_1 = get_cell(i,j,k);
		cell_2 = get_cell(i,j,k-1);
		ghost_cell = get_cell(i,j,k+1);
		ghost_cell.pos[gtl] = 2.0*cell_1.pos[gtl] - cell_2.pos[gtl];
		ghost_cell.volume[gtl] = 2.0*cell_1.volume[gtl] - cell_2.volume[gtl];
		cell_2 = cell_1;
		cell_1 = ghost_cell;
		ghost_cell = get_cell(i,j,k+2);
		ghost_cell.pos[gtl] = 2.0*cell_1.pos[gtl] - cell_2.pos[gtl];
		ghost_cell.volume[gtl] = 2.0*cell_1.volume[gtl] - cell_2.volume[gtl];
	    }
	}
    } // end compute_primary_cell_geometric_data()

    override void compute_distance_to_nearest_wall_for_all_cells(int gtl)
    // Used for the turbulence modelling.
    {
	FVCell[6] cell_at_wall;
	double[6] dist, half_width;
	FVInterface face_at_wall;

	foreach(ref FVCell cell; active_cells) {
	    auto ijk = to_ijk_indices(cell.id);
	    size_t i = ijk[0]; size_t j = ijk[1]; size_t k = ijk[2];
	    // Step 1: get distances to all boundaries along all index directions.
	    // If the block is not too distorted, these directions should take us
	    // straight to the bounding walls.
	    // North
	    face_at_wall = get_ifj(i,jmax+1,k);
	    dist[Face.north] = abs(cell.pos[gtl] - face_at_wall.pos);
	    cell_at_wall[Face.north] = get_cell(i,jmax,k);
	    half_width[Face.north] = abs(cell_at_wall[Face.north].pos[gtl] - face_at_wall.pos);
	    // East
	    face_at_wall = get_ifi(imax+1,j,k);
	    dist[Face.east] = abs(cell.pos[gtl] - face_at_wall.pos);
	    cell_at_wall[Face.east] = get_cell(imax,j,k);
	    half_width[Face.east] = abs(cell_at_wall[Face.east].pos[gtl] - face_at_wall.pos);
	    // South
	    face_at_wall = get_ifj(i,jmin,k);
	    dist[Face.south] = abs(cell.pos[gtl] - face_at_wall.pos);
	    cell_at_wall[Face.south] = get_cell(i,jmin,k);
	    half_width[Face.south] = abs(cell_at_wall[Face.south].pos[gtl] - face_at_wall.pos);
	    // West
	    face_at_wall = get_ifi(imin,j,k);
	    dist[Face.west] = abs(cell.pos[gtl] - face_at_wall.pos);
	    cell_at_wall[Face.west] = get_cell(imin,j,k);
	    half_width[Face.west] = abs(cell_at_wall[Face.west].pos[gtl] - face_at_wall.pos);
	    if ( myConfig.dimensions == 3 ) {
		// Top
		face_at_wall = get_ifk(i,j,kmax+1);
		dist[Face.top] = abs(cell.pos[gtl] - face_at_wall.pos);
		cell_at_wall[Face.top] = get_cell(i,j,kmax);
		half_width[Face.top] = abs(cell_at_wall[Face.top].pos[gtl] - face_at_wall.pos);
		// Bottom
		face_at_wall = get_ifk(i,j,kmin);
		dist[Face.bottom] = abs(cell.pos[gtl] - face_at_wall.pos);
		cell_at_wall[Face.bottom] = get_cell(i,j,kmin);
		half_width[Face.bottom] = abs(cell_at_wall[Face.bottom].pos[gtl] - face_at_wall.pos);
	    }

	    // Step 2: Just in case there are no real walls for this block...
	    //
	    // We'll start by storing the largest distance and 
	    // corresponding wall-cell half-width so that we have 
	    // a relatively large distance in case there are no walls
	    // on the boundary of the block.
	    size_t num_faces = (myConfig.dimensions == 3) ? 6 : 4;
	    cell.distance_to_nearest_wall = dist[0];
	    cell.half_cell_width_at_wall = half_width[0];
	    for ( size_t iface = 1; iface < num_faces; ++iface ) {
		if (dist[iface] > cell.distance_to_nearest_wall ) {
		    cell.distance_to_nearest_wall = dist[iface];
		    cell.half_cell_width_at_wall = half_width[iface];
		    cell.cell_at_nearest_wall = cell_at_wall[iface];
		}
	    }
		
	    // Step 3: find the closest real wall.
	    for ( size_t iface = 0; iface < num_faces; ++iface ) {
		if ( bc[iface].is_wall && 
		     // bc[iface].type_code == BCCode.slip_wall &&
		     // Isn't it enough that is_wal is true?
		     dist[iface] < cell.distance_to_nearest_wall ) {
		    cell.distance_to_nearest_wall = dist[iface];
		    cell.half_cell_width_at_wall = half_width[iface];
		    cell.cell_at_nearest_wall = cell_at_wall[iface];
		}
	    }
	} // end foreach cell
    } // end compute_distance_to_nearest_wall_for_all_cells()

    void calc_volumes_2D(int gtl)
    // Compute the PRIMARY cell volumes, areas, and centers 
    //  from the vertex positions.
    //
    // For 2D planar, assume unit length in the Z-direction.
    // For axisymmetry, compute a volume per radian.
    //
    // Determine minimum length and aspect ratio, also.
    {
	size_t i, j;
	double xA, yA, xB, yB, xC, yC, xD, yD;
	double xN, yN, xS, yS, xE, yE, xW, yW;
	double vol, max_vol, min_vol, xyarea;
	double dx, dy, dxN, dyN, dxE, dyE;
	double lengthN, lengthE, length_max, length_min, length_cross;
	double max_aspect, aspect_ratio;
	FVCell cell, source_cell, target_cell;

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
		cell = get_cell(i,j);
		// These are the corners.
		xA = cell.vtx[1].pos[gtl].x;
		yA = cell.vtx[1].pos[gtl].y;
		xB = cell.vtx[2].pos[gtl].x;
		yB = cell.vtx[2].pos[gtl].y;
		xC = cell.vtx[3].pos[gtl].x;
		yC = cell.vtx[3].pos[gtl].y;
		xD = cell.vtx[0].pos[gtl].x;
		yD = cell.vtx[0].pos[gtl].y;
		// Cell area in the (x,y)-plane.
		xyarea = 0.5 * ((xB + xA) * (yB - yA) + (xC + xB) * (yC - yB) +
				(xD + xC) * (yD - yC) + (xA + xD) * (yA - yD));
		// Cell Centroid.
		cell.pos[gtl].refx = 1.0 / (xyarea * 6.0) * 
		    ((yB - yA) * (xA * xA + xA * xB + xB * xB) + 
		     (yC - yB) * (xB * xB + xB * xC + xC * xC) +
		     (yD - yC) * (xC * xC + xC * xD + xD * xD) + 
		     (yA - yD) * (xD * xD + xD * xA + xA * xA));
		cell.pos[gtl].refy = -1.0 / (xyarea * 6.0) * 
		    ((xB - xA) * (yA * yA + yA * yB + yB * yB) + 
		     (xC - xB) * (yB * yB + yB * yC + yC * yC) +
		     (xD - xC) * (yC * yC + yC * yD + yD * yD) + 
		     (xA - xD) * (yD * yD + yD * yA + yA * yA));
		cell.pos[gtl].refz = 0.0;
		// Cell Volume.
		if ( myConfig.axisymmetric ) {
		    // Volume per radian = centroid y-ordinate * cell area
		    vol = xyarea * cell.pos[gtl].y;
		} else {
		    // Assume unit depth in the z-direction.
		    vol = xyarea;
		}
		if (vol < 0.0) {
		    throw new Error(text("Negative cell volume: Block ", id,
					 " vol[", i, " ,", j, "]= ", vol));
		}
		if (vol > max_vol) max_vol = vol;
		if (vol < min_vol) min_vol = vol;
		cell.volume[gtl] = vol;
		cell.areaxy[gtl] = xyarea;
		// Check cell length scale using North and East boundaries.
		// Also, save the minimum length for later use in the CFL
		// checking routine.  Record max aspect ratio over all cells.
		dxN = xC - xB;
		dyN = yC - yB;
		dxE = xA - xB;
		dyE = yA - yB;
		lengthN = sqrt(dxN * dxN + dyN * dyN);
		lengthE = sqrt(dxE * dxE + dyE * dyE);

		length_max = lengthN;
		if (lengthE > length_max) length_max = lengthE;
		length_cross = xyarea / length_max; 
		// estimate of minimum width of cell
		length_min = lengthN;
		if (lengthE < length_min) length_min = lengthE;
		if (length_cross < length_min) length_min = length_cross;
		cell.L_min = length_min;
		aspect_ratio = length_max / length_min;
		if (aspect_ratio > max_aspect) max_aspect = aspect_ratio;

		// Record the cell widths in the i- and j-index directions.
		// The widths are measured between corresponding midpoints of
		// the bounding interfaces.
		// This data is later used by the high-order reconstruction.
		xN = 0.5 * (xC + xB);
		yN = 0.5 * (yC + yB);
		xS = 0.5 * (xD + xA);
		yS = 0.5 * (yD + yA);
		xE = 0.5 * (xA + xB);
		yE = 0.5 * (yA + yB);
		xW = 0.5 * (xD + xC);
		yW = 0.5 * (yD + yC);
		dx = xN - xS;
		dy = yN - yS;
		cell.jLength = sqrt(dx * dx + dy * dy);
		dx = xE - xW;
		dy = yE - yW;
		cell.iLength = sqrt(dx * dx + dy * dy);
		cell.kLength = 0.0;
	    } // j loop
	} // i loop

	// We now need to mirror the cell iLength and jLength
	// around the boundaries.
	// Those boundaries that are adjacent to another block
	// will be updated later with the other-block's cell lengths.
	for ( i = imin; i <= imax; ++i ) {
	    // North boundary
	    j = jmax;
	    source_cell = get_cell(i,j);
	    target_cell = get_cell(i,j+1);
	    target_cell.iLength = source_cell.iLength;
	    target_cell.jLength = source_cell.jLength;
	    target_cell.kLength = 0.0;
	    source_cell = get_cell(i,j-1);
	    target_cell = get_cell(i,j+2);
	    target_cell.iLength = source_cell.iLength;
	    target_cell.jLength = source_cell.jLength;
	    target_cell.kLength = 0.0;
	    // South boundary
	    j = jmin;
	    source_cell = get_cell(i,j);
	    target_cell = get_cell(i,j-1);
	    target_cell.iLength = source_cell.iLength;
	    target_cell.jLength = source_cell.jLength;
	    target_cell.kLength = 0.0;
	    source_cell = get_cell(i,j+1);
	    target_cell = get_cell(i,j-2);
	    target_cell.iLength = source_cell.iLength;
	    target_cell.jLength = source_cell.jLength;
	    target_cell.kLength = 0.0;
	} // end for i

	for ( j = jmin; j <= jmax; ++j ) {
	    // East boundary
	    i = imax;
	    source_cell = get_cell(i,j);
	    target_cell = get_cell(i+1,j);
	    target_cell.iLength = source_cell.iLength;
	    target_cell.jLength = source_cell.jLength;
	    target_cell.kLength = 0.0;
	    source_cell = get_cell(i-1,j);
	    target_cell = get_cell(i+2,j);
	    target_cell.iLength = source_cell.iLength;
	    target_cell.jLength = source_cell.jLength;
	    target_cell.kLength = 0.0;
	    // West boundary
	    i = imin;
	    source_cell = get_cell(i,j);
	    target_cell = get_cell(i-1,j);
	    target_cell.iLength = source_cell.iLength;
	    target_cell.jLength = source_cell.jLength;
	    target_cell.kLength = 0.0;
	    source_cell = get_cell(i+1,j);
	    target_cell = get_cell(i-2,j);
	    target_cell.iLength = source_cell.iLength;
	    target_cell.jLength = source_cell.jLength;
	    target_cell.kLength = 0.0;
	} // end for j
    } // end calc_volumes_2D()

    void calc_faces_2D(int gtl)
    {
	FVInterface iface;
	size_t i, j;
	double xA, xB, yA, yB, xC, yC;
	double LAB, LBC;

	// East-facing interfaces.
	for (i = imin; i <= imax+1; ++i) {
	    for (j = jmin; j <= jmax; ++j) {
		iface = get_ifi(i,j);
		// These are the corners.
		xA = get_vtx(i,j).pos[gtl].x; 
		yA = get_vtx(i,j).pos[gtl].y;
		xB = get_vtx(i,j+1).pos[gtl].x; 
		yB = get_vtx(i,j+1).pos[gtl].y;
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
		if ( myConfig.axisymmetric ) {
		    // Interface area per radian.
		    iface.area[gtl] = LAB * iface.Ybar;
		} else {
		    // Assume unit depth in the Z-direction.
		    iface.area[gtl] = LAB;
		}
		iface.pos = (get_vtx(i,j).pos[gtl] + get_vtx(i,j+1).pos[gtl])/2.0;
	    
	    } // j loop
	} // i loop
    
	// North-facing interfaces.
	for (i = imin; i <= imax; ++i) {
	    for (j = jmin; j <= jmax+1; ++j) {
		iface = get_ifj(i,j);
		// These are the corners.
		xB = get_vtx(i+1,j).pos[gtl].x;
		yB = get_vtx(i+1,j).pos[gtl].y;
		xC = get_vtx(i,j).pos[gtl].x;
		yC = get_vtx(i,j).pos[gtl].y;
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
		if ( myConfig.axisymmetric ) {
		    // Interface area per radian.
		    iface.area[gtl] = LBC * iface.Ybar;
		} else {
		    // Assume unit depth in the Z-direction.
		    iface.area[gtl] = LBC;
		}
		iface.pos = (get_vtx(i+1,j).pos[gtl] + get_vtx(i,j).pos[gtl])/2.0;
	    } // j loop
	} // i loop
    } // end calc_faces_2D()

    void calc_ghost_cell_geom_2D(int gtl)
    // Compute the ghost cell positions and volumes.
    //
    // 'Compute' is a bit too strong to describe what we do here.
    //  Rather this is a first-order extrapolation
    // from interior cells to estimate the position
    // and volume of the ghost cells.
    {
	size_t i, j;
	FVCell cell_1, cell_2, ghost_cell;
	// East boundary
	i = imax;
	for ( j = jmin; j <= jmax; ++j ) {
	    cell_1 = get_cell(i,j);
	    cell_2 = get_cell(i-1,j);
	    ghost_cell = get_cell(i+1,j);
	    ghost_cell.pos[gtl] = 2.0*cell_1.pos[gtl] - cell_2.pos[gtl];
	    ghost_cell.volume[gtl] = 2.0*cell_1.volume[gtl] - cell_2.volume[gtl];
	    cell_2 = cell_1;
	    cell_1 = ghost_cell;
	    ghost_cell = get_cell(i+2,j);
	    ghost_cell.pos[gtl] = 2.0*cell_1.pos[gtl] - cell_2.pos[gtl];
	    ghost_cell.volume[gtl] = 2.0*cell_1.volume[gtl] - cell_2.volume[gtl];
	}
	// West boundary
	i = imin;
	for ( j = jmin; j <= jmax; ++j ) {
	    cell_1 = get_cell(i,j);
	    cell_2 = get_cell(i+1,j);
	    ghost_cell = get_cell(i-1,j);
	    ghost_cell.pos[gtl] = 2.0*cell_1.pos[gtl] - cell_2.pos[gtl];
	    ghost_cell.volume[gtl] = 2.0*cell_1.volume[gtl] - cell_2.volume[gtl];
	    cell_2 = cell_1;
	    cell_1 = ghost_cell;
	    ghost_cell = get_cell(i-2,j);
	    ghost_cell.pos[gtl] = 2.0*cell_1.pos[gtl] - cell_2.pos[gtl];
	    ghost_cell.volume[gtl] = 2.0*cell_1.volume[gtl] - cell_2.volume[gtl];
	}
	// North boundary
	j = jmax;
	for ( i = imin; i <= imax; ++i ) {
	    cell_1 = get_cell(i,j);
	    cell_2 = get_cell(i,j-1);
	    ghost_cell = get_cell(i,j+1);
	    ghost_cell.pos[gtl] = 2.0*cell_1.pos[gtl] - cell_2.pos[gtl];
	    ghost_cell.volume[gtl] = 2.0*cell_1.volume[gtl] - cell_2.volume[gtl];
	    cell_2 = cell_1;
	    cell_1 = ghost_cell;
	    ghost_cell = get_cell(i,j+2);
	    ghost_cell.pos[gtl] = 2.0*cell_1.pos[gtl] - cell_2.pos[gtl];
	    ghost_cell.volume[gtl] = 2.0*cell_1.volume[gtl] - cell_2.volume[gtl];
	}
	// South boundary
	j = jmin;
	for ( i = imin; i <= imax; ++i ) {
	    cell_1 = get_cell(i,j);
	    cell_2 = get_cell(i,j+1);
	    ghost_cell = get_cell(i,j-1);
	    ghost_cell.pos[gtl] = 2.0*cell_1.pos[gtl] - cell_2.pos[gtl];
	    ghost_cell.volume[gtl] = 2.0*cell_1.volume[gtl] - cell_2.volume[gtl];
	    cell_2 = cell_1;
	    cell_1 = ghost_cell;
	    ghost_cell = get_cell(i,j-2);
	    ghost_cell.pos[gtl] = 2.0*cell_1.pos[gtl] - cell_2.pos[gtl];
	    ghost_cell.volume[gtl] = 2.0*cell_1.volume[gtl] - cell_2.volume[gtl];
	}
    } // end calc_ghost_cell_geom_2D()

    override void assign_flow_locations_for_derivative_calc(size_t gtl)
    {
	size_t i, j, k;
	if (myConfig.dimensions == 2) {
	    // First, do all of the internal secondary cells.
	    // i.e. Those not on a boundary.
	    for ( i = imin+1; i <= imax; ++i ) {
		for ( j = jmin+1; j <= jmax; ++j ) {
		    // Secondary-cell centre is a primary-cell vertex.
		    FVVertex vtx = get_vtx(i,j);
		    // These are the corners of the secondary cell.
		    FVCell A = get_cell(i,j-1);
		    FVCell B = get_cell(i,j);
		    FVCell C = get_cell(i-1,j);
		    FVCell D = get_cell(i-1,j-1);
		    // Retain locations and references to flow states for later.
		    vtx.cloud_pos = [A.pos[gtl], B.pos[gtl], C.pos[gtl], D.pos[gtl]];
		    vtx.cloud_fs = [A.fs, B.fs, C.fs, D.fs];
		} // j loop
	    } // i loop
	    // Half-cells along the edges of the block.
	    // East boundary
	    i = imax+1;
	    for (j = jmin+1; j <= jmax; ++j) {
		FVVertex vtx = get_vtx(i,j);
		FVInterface A = get_ifi(i,j-1);
		FVInterface B = get_ifi(i,j);
		FVCell C = get_cell(i-1,j);
		FVCell D = get_cell(i-1,j-1);
		vtx.cloud_pos = [A.pos, B.pos, C.pos[gtl], D.pos[gtl]];
		vtx.cloud_fs = [A.fs, B.fs, C.fs, D.fs];
	    } // j loop
	    // West boundary
	    i = imin;
	    for (j = jmin+1; j <= jmax; ++j) {
		FVVertex vtx = get_vtx(i,j);
		// These are the corners of the secondary cell.
		FVCell A = get_cell(i,j-1);
		FVCell B = get_cell(i,j);
		FVInterface C = get_ifi(i,j);
		FVInterface D = get_ifi(i,j-1);
		vtx.cloud_pos = [A.pos[gtl], B.pos[gtl], C.pos, D.pos];
		vtx.cloud_fs = [A.fs, B.fs, C.fs, D.fs];
	    } // j loop
	    // North boundary
	    j = jmax+1;
	    for (i = imin+1; i <= imax; ++i) {
		FVVertex vtx = get_vtx(i,j);
		FVCell A = get_cell(i,j-1);
		FVInterface B = get_ifj(i,j);
		FVInterface C = get_ifj(i-1,j);
		FVCell D = get_cell(i-1,j-1);
		vtx.cloud_pos = [A.pos[gtl], B.pos, C.pos, D.pos[gtl]];
		vtx.cloud_fs = [A.fs, B.fs, C.fs, D.fs];
	    } // i loop
	    // South boundary
	    j = jmin;
	    for (i = imin+1; i <= imax; ++i) {
		FVVertex vtx = get_vtx(i,j);
		FVInterface A = get_ifj(i,j);
		FVCell B = get_cell(i,j);
		FVCell C = get_cell(i-1,j);
		FVInterface D = get_ifj(i-1,j);
		vtx.cloud_pos = [A.pos, B.pos[gtl], C.pos[gtl], D.pos];
		vtx.cloud_fs = [A.fs, B.fs, C.fs, D.fs];
	    } // i loop
	    // For the corners, we are going to use the same divergence-theorem-based
	    // gradient calculator and let one edge collapse to a point, thus giving
	    // it a triangle to compute over.  This should be fine. 
	    // North-east corner
	    {
		i = imax+1; j = jmax+1;
		FVVertex vtx = get_vtx(i,j);
		FVInterface A = get_ifi(i,j-1);
		FVInterface B = get_ifj(i-1,j);
		FVInterface C = get_ifj(i-1,j);
		FVCell D = get_cell(i-1,j-1);
		vtx.cloud_pos = [A.pos, B.pos, C.pos, D.pos[gtl]];
		vtx.cloud_fs = [A.fs, B.fs, C.fs, D.fs];
	    }
	    // South-east corner
	    {
		i = imax+1; j = jmin;
		FVVertex vtx = get_vtx(i,j);
		FVInterface A = get_ifi(i,j);
		FVInterface B = get_ifi(i,j);
		FVCell C = get_cell(i-1,j);
		FVInterface D = get_ifj(i-1,j);
		vtx.cloud_pos = [A.pos, B.pos, C.pos[gtl], D.pos];
		vtx.cloud_fs = [A.fs, B.fs, C.fs, D.fs];
	    }
	    // South-west corner
	    {
		i = imin; j = jmin;
		FVVertex vtx = get_vtx(i,j);
		FVInterface A = get_ifj(i,j);
		FVCell B = get_cell(i,j);
		FVInterface C = get_ifi(i,j);
		FVInterface D = A;
		vtx.cloud_pos = [A.pos, B.pos[gtl], C.pos, D.pos];
		vtx.cloud_fs = [A.fs, B.fs, C.fs, D.fs];
	    }
	    // North-west corner
	    {
		i = imin; j = jmax+1;
		FVVertex vtx = get_vtx(i,j);
		FVCell A = get_cell(i,j-1);
		FVInterface B = get_ifj(i,j);
		FVInterface C = B;
		FVInterface D = get_ifi(i,j-1);
		vtx.cloud_pos = [A.pos[gtl], B.pos, C.pos, D.pos];
		vtx.cloud_fs = [A.fs, B.fs, C.fs, D.fs];
	    }
	} else { // Flow quantity derivatives for 3D.
	    // Internal secondary cell geometry information
	    for ( i = imin; i <= imax-1; ++i ) {
		for ( j = jmin; j <= jmax-1; ++j ) {
		    for ( k = kmin; k <= kmax-1; ++k ) {
			FVVertex vtx = get_vtx(i+1,j+1,k+1);
			FVCell c0 = get_cell(i,j,k);
			FVCell c1 = get_cell(i+1,j,k);
			FVCell c2 = get_cell(i+1,j+1,k);
			FVCell c3 = get_cell(i,j+1,k);
			FVCell c4 = get_cell(i,j,k+1);
			FVCell c5 = get_cell(i+1,j,k+1);
			FVCell c6 = get_cell(i+1,j+1,k+1);
			FVCell c7 = get_cell(i,j+1,k+1);
			vtx.cloud_pos = [c0.pos[gtl], c1.pos[gtl], c2.pos[gtl], c3.pos[gtl],
					 c4.pos[gtl], c5.pos[gtl], c6.pos[gtl], c7.pos[gtl]];
			vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs, c6.fs, c7.fs];
		    }
		}
	    }
	    // East boundary secondary cell geometry information
	    i = imax;
	    for ( j = jmin; j <= jmax-1; ++j ) {
		for ( k = kmin; k <= kmax-1; ++k ) {
		    FVVertex vtx = get_vtx(i+1,j+1,k+1);
		    FVCell c0 = get_cell(i,j,k);
		    FVInterface c1 = get_ifi(i+1,j,k);
		    FVInterface c2 = get_ifi(i+1,j+1,k);
		    FVCell c3 = get_cell(i,j+1,k);
		    FVCell c4 = get_cell(i,j,k+1);
		    FVInterface c5 = get_ifi(i+1,j,k+1);
		    FVInterface c6 = get_ifi(i+1,j+1,k+1);
		    FVCell c7 = get_cell(i,j+1,k+1);
		    vtx.cloud_pos = [c0.pos[gtl], c1.pos, c2.pos, c3.pos[gtl],
				     c4.pos[gtl], c5.pos, c6.pos, c7.pos[gtl]];
		    vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs, c6.fs, c7.fs];
		}
	    }
	    // West boundary secondary cell geometry information
	    i = imin - 1;
	    for ( j = jmin; j <= jmax-1; ++j ) {
		for ( k = kmin; k <= kmax-1; ++k ) {
		    FVVertex vtx = get_vtx(i+1,j+1,k+1);
		    FVInterface c0 = get_ifi(i+1,j,k);
		    FVCell c1 = get_cell(i+1,j,k);
		    FVCell c2 = get_cell(i+1,j+1,k);
		    FVInterface c3 = get_ifi(i+1,j+1,k);
		    FVInterface c4 = get_ifi(i+1,j,k+1);
		    FVCell c5 = get_cell(i+1,j,k+1);
		    FVCell c6 = get_cell(i+1,j+1,k+1);
		    FVInterface c7 = get_ifi(i+1,j+1,k+1);
		    vtx.cloud_pos = [c0.pos, c1.pos[gtl], c2.pos[gtl], c3.pos,
				     c4.pos, c5.pos[gtl], c6.pos[gtl], c7.pos];
		    vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs, c6.fs, c7.fs];
		}
	    }
	    // North boundary secondary cell geometry information
	    j = jmax;
	    for ( i = imin; i <= imax-1; ++i ) {
		for ( k = kmin; k <= kmax-1; ++k ) {
		    FVVertex vtx = get_vtx(i+1,j+1,k+1);
		    FVCell c0 = get_cell(i,j,k);
		    FVCell c1 = get_cell(i+1,j,k);
		    FVInterface c2 = get_ifj(i+1,j+1,k);
		    FVInterface c3 = get_ifj(i,j+1,k);
		    FVCell c4 = get_cell(i,j,k+1);
		    FVCell c5 = get_cell(i+1,j,k+1);
		    FVInterface c6 = get_ifj(i+1,j+1,k+1);
		    FVInterface c7 = get_ifj(i,j+1,k+1);
		    vtx.cloud_pos = [c0.pos[gtl], c1.pos[gtl], c2.pos, c3.pos,
				     c4.pos[gtl], c5.pos[gtl], c6.pos, c7.pos];
		    vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs, c6.fs, c7.fs];
		}
	    }
	    // South boundary secondary cell geometry information
	    j = jmin - 1;
	    for ( i = imin; i <= imax-1; ++i ) {
		for ( k = kmin; k <= kmax-1; ++k ) {
		    FVVertex vtx = get_vtx(i+1,j+1,k+1);
		    FVInterface c0 = get_ifj(i,j+1,k);
		    FVInterface c1 = get_ifj(i+1,j+1,k);
		    FVCell c2 = get_cell(i+1,j+1,k);
		    FVCell c3 = get_cell(i,j+1,k);
		    FVInterface c4 = get_ifj(i,j+1,k+1);
		    FVInterface c5 = get_ifj(i+1,j+1,k+1);
		    FVCell c6 = get_cell(i+1,j+1,k+1);
		    FVCell c7 = get_cell(i,j+1,k+1);
		    vtx.cloud_pos = [c0.pos, c1.pos, c2.pos[gtl], c3.pos[gtl],
				     c4.pos, c5.pos, c6.pos[gtl], c7.pos[gtl]];
		    vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs, c6.fs, c7.fs];
		}
	    }
	    // Top boundary secondary cell geometry information
	    k = kmax;
	    for ( i = imin; i <= imax-1; ++i ) {
		for ( j = jmin; j <= jmax-1; ++j ) {
		    FVVertex vtx = get_vtx(i+1,j+1,k+1);
		    FVCell c0 = get_cell(i,j,k);
		    FVCell c1 = get_cell(i+1,j,k);
		    FVCell c2 = get_cell(i+1,j+1,k);
		    FVCell c3 = get_cell(i,j+1,k);
		    FVInterface c4 = get_ifk(i,j,k+1);
		    FVInterface c5 = get_ifk(i+1,j,k+1);
		    FVInterface c6 = get_ifk(i+1,j+1,k+1);
		    FVInterface c7 = get_ifk(i,j+1,k+1);
		    vtx.cloud_pos = [c0.pos[gtl], c1.pos[gtl], c2.pos[gtl], c3.pos[gtl],
				     c4.pos, c5.pos, c6.pos, c7.pos];
		    vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs, c6.fs, c7.fs];
		}
	    }
	    // Bottom boundary secondary cell geometry information
	    k = kmin - 1;
	    for ( i = imin; i <= imax-1; ++i ) {
		for ( j = jmin; j <= jmax-1; ++j ) {
		    FVVertex vtx = get_vtx(i+1,j+1,k+1);
		    FVInterface c0 = get_ifk(i,j,k+1);
		    FVInterface c1 = get_ifk(i+1,j,k+1);
		    FVInterface c2 = get_ifk(i+1,j+1,k+1);
		    FVInterface c3 = get_ifk(i,j+1,k+1);
		    FVCell c4 = get_cell(i,j,k+1);
		    FVCell c5 = get_cell(i+1,j,k+1);
		    FVCell c6 = get_cell(i+1,j+1,k+1);
		    FVCell c7 = get_cell(i,j+1,k+1);
		    vtx.cloud_pos = [c0.pos, c1.pos, c2.pos, c3.pos,
				     c4.pos[gtl], c5.pos[gtl], c6.pos[gtl], c7.pos[gtl]];
		    vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs, c6.fs, c7.fs];
		}
	    }
	    // Now, do the 4 edges around the bottom face.
	    // Bottom-South edge [0]-->[1]
	    j = jmin; k = kmin;    	
	    for ( i = imin+1; i <= imax; ++i ) {
		FVVertex vtx = get_vtx(i,j,k);
		FVCell c0 = get_cell(i-1,j,k);
		FVCell c1 = get_cell(i,j,k);
		FVInterface c2 = get_ifj(i-1,j,k);
		FVInterface c3 = get_ifk(i-1,j,k);
		FVInterface c4 = get_ifj(i,j,k);
		FVInterface c5 = get_ifk(i,j,k);
		vtx.cloud_pos = [c0.pos[gtl], c1.pos[gtl], c2.pos, c3.pos, c4.pos, c5.pos];
		vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs];
	    }
	    // Bottom-North edge [3]-->[2]
	    j = jmax; k = kmin;
	    for ( i = imin+1; i <= imax; ++i ) {
		FVVertex vtx = get_vtx(i,j+1,k);
		FVCell c0 = get_cell(i-1,j,k);
		FVCell c1 = get_cell(i,j,k);
		FVInterface c2 = get_ifj(i-1,j+1,k);
		FVInterface c3 = get_ifk(i-1,j,k);
		FVInterface c4 = get_ifj(i,j+1,k);
		FVInterface c5 = get_ifk(i,j,k);
		vtx.cloud_pos = [c0.pos[gtl], c1.pos[gtl], c2.pos, c3.pos, c4.pos, c5.pos];
		vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs];
	    }
	    // Bottom-West edge [0]-->[3]
	    i = imin; k = kmin;
	    for ( j = jmin+1; j <= jmax; ++j ) {
		FVVertex vtx = get_vtx(i,j,k);
		FVCell c0 = get_cell(i,j-1,k);
		FVCell c1 = get_cell(i,j,k);
		FVInterface c2 = get_ifi(i,j-1,k);
		FVInterface c3 = get_ifk(i,j-1,k);
		FVInterface c4 = get_ifi(i,j,k);
		FVInterface c5 = get_ifk(i,j,k);
		vtx.cloud_pos = [c0.pos[gtl], c1.pos[gtl], c2.pos, c3.pos, c4.pos, c5.pos];
		vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs];
	    }
	    // Bottom-East edge [1]-->[2]
	    i = imax; k = kmin;
	    for ( j = jmin+1; j <= jmax; ++j ) {
		FVVertex vtx = get_vtx(i+1,j,k);
		FVCell c0 = get_cell(i,j-1,k);
		FVCell c1 = get_cell(i,j,k);
		FVInterface c2 = get_ifi(i+1,j-1,k);
		FVInterface c3 = get_ifk(i,j-1,k);
		FVInterface c4 = get_ifi(i+1,j,k);
		FVInterface c5 = get_ifk(i,j,k);
		vtx.cloud_pos = [c0.pos[gtl], c1.pos[gtl], c2.pos, c3.pos, c4.pos, c5.pos];
		vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs];
	    }
	    // 4 edges around the top face.
	    // Top-South edge [4]-->[5]
	    j = jmin; k = kmax;
	    for ( i = imin+1; i <= imax; ++i ) {
		FVVertex vtx = get_vtx(i,j,k+1);
		FVCell c0 = get_cell(i-1,j,k);
		FVCell c1 = get_cell(i,j,k);
		FVInterface c2 = get_ifj(i-1,j,k);
		FVInterface c3 = get_ifk(i-1,j,k+1);
		FVInterface c4 = get_ifj(i,j,k);
		FVInterface c5 = get_ifk(i,j,k+1);
		vtx.cloud_pos = [c0.pos[gtl], c1.pos[gtl], c2.pos, c3.pos, c4.pos, c5.pos];
		vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs];
	    }
	    // Top-North edge [7]-->[6]
	    j = jmax; k = kmax;
	    for ( i = imin+1; i <= imax; ++i ) {
		FVVertex vtx = get_vtx(i,j+1,k+1);
		FVCell c0 = get_cell(i-1,j,k);
		FVCell c1 = get_cell(i,j,k);
		FVInterface c2 = get_ifj(i-1,j+1,k);
		FVInterface c3 = get_ifk(i-1,j,k+1);
		FVInterface c4 = get_ifj(i,j+1,k);
		FVInterface c5 = get_ifk(i,j,k+1);
		vtx.cloud_pos = [c0.pos[gtl], c1.pos[gtl], c2.pos, c3.pos, c4.pos, c5.pos];
		vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs];
	    }
	    // Top-West edge [4]-->[7]
	    i = imin; k = kmax;
	    for ( j = jmin+1; j <= jmax; ++j ) {
		FVVertex vtx = get_vtx(i,j,k+1);
		FVCell c0 = get_cell(i,j-1,k);
		FVCell c1 = get_cell(i,j,k);
		FVInterface c2 = get_ifi(i,j-1,k);
		FVInterface c3 = get_ifk(i,j-1,k+1);
		FVInterface c4 = get_ifi(i,j,k);
		FVInterface c5 = get_ifk(i,j,k+1);
		vtx.cloud_pos = [c0.pos[gtl], c1.pos[gtl], c2.pos, c3.pos, c4.pos, c5.pos];
		vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs];
	    }
	    // Top-East edge [5]-->[6]
	    i = imax; k = kmax;
	    for ( j = jmin+1; j <= jmax; ++j ) {
		FVVertex vtx = get_vtx(i+1,j,k+1);
		FVCell c0 = get_cell(i,j-1,k);
		FVCell c1 = get_cell(i,j,k);
		FVInterface c2 = get_ifi(i+1,j-1,k);
		FVInterface c3 = get_ifk(i,j-1,k+1);
		FVInterface c4 = get_ifi(i+1,j,k);
		FVInterface c5 = get_ifk(i,j,k+1);
		vtx.cloud_pos = [c0.pos[gtl], c1.pos[gtl], c2.pos, c3.pos, c4.pos, c5.pos];
		vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs];
	    }
	    // 4 edges running from bottom to top.
	    // South-West edge [0]-->[4]
	    i = imin; j = jmin;
	    for ( k = kmin+1; k <= kmax; ++k ) {
		FVVertex vtx = get_vtx(i,j,k);
		FVCell c0 = get_cell(i,j,k-1);
		FVCell c1 = get_cell(i,j,k);
		FVInterface c2 = get_ifi(i,j,k-1);
		FVInterface c3 = get_ifj(i,j,k-1);
		FVInterface c4 = get_ifi(i,j,k);
		FVInterface c5 = get_ifj(i,j,k);
		vtx.cloud_pos = [c0.pos[gtl], c1.pos[gtl], c2.pos, c3.pos, c4.pos, c5.pos];
		vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs];
	    }
	    // South-East edge [1]-->[5]
	    i = imax; j = jmin;
	    for ( k = kmin+1; k <= kmax; ++k ) {
		FVVertex vtx = get_vtx(i+1,j,k);
		FVCell c0 = get_cell(i,j,k-1);
		FVCell c1 = get_cell(i,j,k);
		FVInterface c2 = get_ifi(i+1,j,k-1);
		FVInterface c3 = get_ifj(i,j,k-1);
		FVInterface c4 = get_ifi(i+1,j,k);
		FVInterface c5 = get_ifj(i,j,k);
		vtx.cloud_pos = [c0.pos[gtl], c1.pos[gtl], c2.pos, c3.pos, c4.pos, c5.pos];
		vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs];
	    }
	    // North-East edge [2]-->[6]
	    i = imax; j = jmax;
	    for ( k = kmin+1; k <= kmax; ++k ) {
		FVVertex vtx = get_vtx(i+1,j+1,k);
		FVCell c0 = get_cell(i,j,k-1);
		FVCell c1 = get_cell(i,j,k);
		FVInterface c2 = get_ifi(i+1,j,k-1);
		FVInterface c3 = get_ifj(i,j+1,k-1);
		FVInterface c4 = get_ifi(i+1,j,k);
		FVInterface c5 = get_ifj(i,j+1,k);
		vtx.cloud_pos = [c0.pos[gtl], c1.pos[gtl], c2.pos, c3.pos, c4.pos, c5.pos];
		vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs];
	    }
	    // North-West edge [3]-->[7]
	    i = imin; j = jmax;
	    for ( k = kmin+1; k <= kmax; ++k ) {
		FVVertex vtx = get_vtx(i,j+1,k);
		FVCell c0 = get_cell(i,j,k-1);
		FVCell c1 = get_cell(i,j,k);
		FVInterface c2 = get_ifi(i,j,k-1);
		FVInterface c3 = get_ifj(i,j+1,k-1);
		FVInterface c4 = get_ifi(i,j,k);
		FVInterface c5 = get_ifj(i,j+1,k);
		vtx.cloud_pos = [c0.pos[gtl], c1.pos[gtl], c2.pos, c3.pos, c4.pos, c5.pos];
		vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs];
	    }
	    // Finally, the 8 corners.
	    // South-West-Bottom corner [0]
	    i = imin; j = jmin; k = kmin;
	    {
		FVVertex vtx = get_vtx(i,j,k);
		FVCell c0 = get_cell(i,j,k);
		FVInterface c1 = get_ifi(i,j,k);
		FVInterface c2 = get_ifj(i,j,k);
		FVInterface c3 = get_ifk(i,j,k);
		vtx.cloud_pos = [c0.pos[gtl], c1.pos, c2.pos, c3.pos];
		vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs];
	    }
	    // South-East-Bottom corner [1]
	    i = imax; j = jmin; k = kmin;
	    {
		FVVertex vtx = get_vtx(i+1,j,k);
		FVCell c0 = get_cell(i,j,k);
		FVInterface c1 = get_ifi(i+1,j,k);
		FVInterface c2 = get_ifj(i,j,k);
		FVInterface c3 = get_ifk(i,j,k);
		vtx.cloud_pos = [c0.pos[gtl], c1.pos, c2.pos, c3.pos];
		vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs];
	    }
	    // North-East-Bottom corner [2]
	    i = imax; j = jmax; k = kmin;
	    {
		FVVertex vtx = get_vtx(i+1,j+1,k);
		FVCell c0 = get_cell(i,j,k);
		FVInterface c1 = get_ifi(i+1,j,k);
		FVInterface c2 = get_ifj(i,j+1,k);
		FVInterface c3 = get_ifk(i,j,k);
		vtx.cloud_pos = [c0.pos[gtl], c1.pos, c2.pos, c3.pos];
		vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs];
	    }
	    // North-West-Bottom corner [3]
	    i = imin; j = jmax; k = kmin;
	    {
		FVVertex vtx = get_vtx(i,j+1,k);
		FVCell c0 = get_cell(i,j,k);
		FVInterface c1 = get_ifi(i,j,k);
		FVInterface c2 = get_ifj(i,j+1,k);
		FVInterface c3 = get_ifk(i,j,k);
		vtx.cloud_pos = [c0.pos[gtl], c1.pos, c2.pos, c3.pos];
		vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs];
	    }
	    // South-West-Top corner [4]
	    i = imin; j = jmin; k = kmax;
	    {
		FVVertex vtx = get_vtx(i,j,k+1);
		FVCell c0 = get_cell(i,j,k);
		FVInterface c1 = get_ifi(i,j,k);
		FVInterface c2 = get_ifj(i,j,k);
		FVInterface c3 = get_ifk(i,j,k+1);
		vtx.cloud_pos = [c0.pos[gtl], c1.pos, c2.pos, c3.pos];
		vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs];
	    }
	    // South-East-Top corner [5]
	    i = imax; j = jmin; k = kmax;
	    {
		FVVertex vtx = get_vtx(i+1,j,k+1);
		FVCell c0 = get_cell(i,j,k);
		FVInterface c1 = get_ifi(i+1,j,k);
		FVInterface c2 = get_ifj(i,j,k);
		FVInterface c3 = get_ifk(i,j,k+1);
		vtx.cloud_pos = [c0.pos[gtl], c1.pos, c2.pos, c3.pos];
		vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs];
	    }
	    // North-East-Top corner [6]
	    i = imax; j = jmax; k = kmax;
	    {
		FVVertex vtx = get_vtx(i+1,j+1,k+1);
		FVCell c0 = get_cell(i,j,k);
		FVInterface c1 = get_ifi(i+1,j,k);
		FVInterface c2 = get_ifj(i,j+1,k);
		FVInterface c3 = get_ifk(i,j,k+1);
		vtx.cloud_pos = [c0.pos[gtl], c1.pos, c2.pos, c3.pos];
		vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs];
	    }
	    // North-West-Top corner [7]
	    i = imin; j = jmax; k = kmax;
	    {
		FVVertex vtx = get_vtx(i,j+1,k+1);
		FVCell c0 = get_cell(i,j,k);
		FVInterface c1 = get_ifi(i,j,k);
		FVInterface c2 = get_ifj(i,j+1,k);
		FVInterface c3 = get_ifk(i,j,k+1);
		vtx.cloud_pos = [c0.pos[gtl], c1.pos, c2.pos, c3.pos];
		vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs];
	    }
	} // end if (myConfig.dimensions
    } // end assign_flow_locations_for_derivative_calc()

    override void read_grid(string filename, size_t gtl=0)
    // Read the grid vertices from a gzip file.
    {
	size_t nivtx, njvtx, nkvtx;
	double x, y, z;
	if (myConfig.verbosity_level >= 1) {
	    writeln("read_grid(): Start block ", id);
	}
	auto byLine = new GzipByLine(filename);
	auto line = byLine.front; byLine.popFront();
	formattedRead(line, "%d %d %d", &nivtx, &njvtx, &nkvtx);
	if ( myConfig.dimensions == 3 ) {
	    if ( nivtx-1 != nicell || njvtx-1 != njcell || nkvtx-1 != nkcell ) {
		throw new Error(text("For block[", id, "] we have a mismatch in 3D grid size.",
                                     " Have read nivtx=", nivtx, " njvtx=", njvtx,
				     " nkvtx=", nkvtx));
	    }
	    for ( size_t k = kmin; k <= kmax+1; ++k ) {
		for ( size_t j = jmin; j <= jmax+1; ++j ) {
		    for ( size_t i = imin; i <= imax+1; ++i ) {
			line = byLine.front; byLine.popFront();
			// Note that the line starts with whitespace.
			formattedRead(line, " %g %g %g", &x, &y, &z);
			auto vtx = get_vtx(i,j,k);
			vtx.pos[gtl].refx = x;
			vtx.pos[gtl].refy = y;
			vtx.pos[gtl].refz = z;
		    } // for i
		} // for j
	    } // for k
	} else { // 2D case
	    if ( nivtx-1 != nicell || njvtx-1 != njcell || nkvtx != 1 ) {
		throw new Error(text("For block[", id, "] we have a mismatch in 2D grid size.",
				     " Have read nivtx=", nivtx, " njvtx=", njvtx,
				     " nkvtx=", nkvtx));
	    }
	    for ( size_t j = jmin; j <= jmax+1; ++j ) {
		for ( size_t i = imin; i <= imax+1; ++i ) {
		    line = byLine.front; byLine.popFront();
		    // Note that the line starts with whitespace.
		    formattedRead(line, " %g %g", &x, &y);
		    auto vtx = get_vtx(i,j);
		    vtx.pos[gtl].refx = x;
		    vtx.pos[gtl].refy = y;
		    vtx.pos[gtl].refz = 0.0;
		} // for i
	    } // for j
	}
    } // end read_grid()

    override void write_grid(string filename, double sim_time, size_t gtl=0)
    {
	if (myConfig.verbosity_level >= 1) {
	    writeln("write_grid(): Start block ", id);
	}
	size_t kmaxrange;
	auto outfile = new GzipOut(filename);
	auto writer = appender!string();
	if ( myConfig.dimensions == 3 ) {
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
		    auto vtx = get_vtx(i,j,k);
		    writer = appender!string();
		    formattedWrite(writer, "%20.12e %20.12e %20.12e\n", vtx.pos[gtl].x,
				   vtx.pos[gtl].y, vtx.pos[gtl].z);
		    outfile.compress(writer.data);
		} // for i
	    } // for j
	} // for k
	outfile.finish();
    } // end write_grid()

    override double read_solution(string filename)
    // Note that the position data is read into grid-time-level 0
    // by scan_values_from_string(). 
    // Returns sim_time from file.
    {
	size_t ni, nj, nk;
	double sim_time;
	if (myConfig.verbosity_level >= 1) {
	    writeln("read_solution(): Start block ", id);
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
	     nk != ((myConfig.dimensions == 3) ? nkcell : 1) ) {
	    throw new Error(text("For block[", id, "] we have a mismatch in solution size.",
				 " Have read ni=", ni, " nj=", nj, " nk=", nk));
	}	
	for ( size_t k = kmin; k <= kmax; ++k ) {
	    for ( size_t j = jmin; j <= jmax; ++j ) {
		for ( size_t i = imin; i <= imax; ++i ) {
		    line = byLine.front; byLine.popFront();
		    get_cell(i,j,k).scan_values_from_string(line);
		} // for i
	    } // for j
	} // for k
	return sim_time;
    } // end read_solution()

    override void write_solution(string filename, double sim_time)
    // Write the flow solution (i.e. the primary variables at the cell centers)
    // for a single block.
    // This is almost Tecplot POINT format.
    {
	if (myConfig.verbosity_level >= 1) {
	    writeln("write_solution(): Start block ", id);
	}
	auto outfile = new GzipOut(filename);
	auto writer = appender!string();
	formattedWrite(writer, "%20.12e\n", sim_time);
	outfile.compress(writer.data);
	writer = appender!string();
	foreach(varname; variable_list_for_cell(myConfig.gmodel)) {
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
		    outfile.compress(" " ~ get_cell(i,j,k).write_values_to_string() ~ "\n");
		} // for i
	    } // for j
	} // for k
	outfile.finish();
    } // end write_solution()

    override void convective_flux()
    {
	FVCell cL1, cL0, cR0, cR1;
	FVInterface IFace;
	// ifi interfaces are East-facing interfaces.
	for ( size_t k = kmin; k <= kmax; ++k ) {
	    for ( size_t j = jmin; j <= jmax; ++j ) {
		for ( size_t i = imin; i <= imax+1; ++i ) {
		    IFace = get_ifi(i,j,k);
		    cL1 = get_cell(i-2,j,k);
		    cL0 = get_cell(i-1,j,k);
		    cR0 = get_cell(i,j,k);
		    cR1 = get_cell(i+1,j,k);
		    // Compute the flux from data on either-side of the interface.
		    // First, interpolate LEFT and RIGHT interface states from cell-center properties.
		    if ( (i == imin+1) && (bc[Face.west].ghost_cell_data_available == false) ) {
			one_d.interp_right(IFace, cL0, cR0, cR1, 
					   cL0.iLength, cR0.iLength, cR1.iLength,
					   Lft, Rght);
		    } else if ( (i == imax) && (bc[Face.east].ghost_cell_data_available == false) ) {
			one_d.interp_left(IFace, cL1, cL0, cR0, 
					  cL1.iLength, cL0.iLength, cR0.iLength,
					  Lft, Rght);
		    } else { // General symmetric reconstruction.
			one_d.interp_both(IFace, cL1, cL0, cR0, cR1,
					  cL1.iLength, cL0.iLength, cR0.iLength, cR1.iLength,
					  Lft, Rght);
		    }
		    // Second, save u, v, w, T for the viscous flux calculation by making a local average.
		    // The values for u, v and T may be updated subsequently by the interface-flux function.
		    if ( (i == imin) && (bc[Face.west].ghost_cell_data_available == false) ) {
			IFace.fs.copy_average_values_from(Rght, Rght);
		    } else if ( (i == imax+1) && (bc[Face.east].ghost_cell_data_available == false) ) {
			IFace.fs.copy_average_values_from(Lft, Lft);
		    } else {
			IFace.fs.copy_average_values_from(Lft, Rght);
		    }
		    compute_interface_flux(Lft, Rght, IFace, myConfig.gmodel, omegaz);
		} // i loop
	    } // j loop
	} // for k
	// ifj interfaces are North-facing interfaces.
	for ( size_t k = kmin; k <= kmax; ++k ) {
	    for ( size_t i = imin; i <= imax; ++i ) {
		for ( size_t j = jmin; j <= jmax+1; ++j ) {
		    IFace = get_ifj(i,j,k);
		    cL1 = get_cell(i,j-2,k);
		    cL0 = get_cell(i,j-1,k);
		    cR0 = get_cell(i,j,k);
		    cR1 = get_cell(i,j+1,k);
		    // Interpolate LEFT and RIGHT interface states from the cell-center properties.
		    if ( (j == jmin+1) && (bc[Face.south].ghost_cell_data_available == false) ) {
			one_d.interp_right(IFace, cL0, cR0, cR1, 
					   cL0.jLength, cR0.jLength, cR1.jLength,
					   Lft, Rght);
		    } else if ( (j == jmax) && (bc[Face.north].ghost_cell_data_available == false) ) {
			one_d.interp_left(IFace, cL1, cL0, cR0, 
					  cL1.jLength, cL0.jLength, cR0.jLength,
					  Lft, Rght);
		    } else { // General symmetric reconstruction.
			one_d.interp_both(IFace, cL1, cL0, cR0, cR1,
					  cL1.jLength, cL0.jLength, cR0.jLength, cR1.jLength,
					  Lft, Rght);
		    }
		    // Second, save u, v, w, T for the viscous flux calculation by making a local average.
		    // The values for u, v and T may be updated subsequently by the interface-flux function.
		    if ( (j == jmin) && (bc[Face.south].ghost_cell_data_available == false) ) {
			IFace.fs.copy_average_values_from(Rght, Rght);
		    } else if ( (j == jmax+1) && (bc[Face.north].ghost_cell_data_available == false) ) {
			IFace.fs.copy_average_values_from(Lft, Lft);
		    } else {
			IFace.fs.copy_average_values_from(Lft, Rght);
		    }
		    compute_interface_flux(Lft, Rght, IFace, myConfig.gmodel, omegaz);
		} // j loop
	    } // i loop
	} // for k
    
	if (myConfig.dimensions == 2) return;
    
	// ifk interfaces are Top-facing interfaces.
	for ( size_t i = imin; i <= imax; ++i ) {
	    for ( size_t j = jmin; j <= jmax; ++j ) {
		for ( size_t k = kmin; k <= kmax+1; ++k ) {
		    IFace = get_ifk(i,j,k);
		    cL1 = get_cell(i,j,k-2);
		    cL0 = get_cell(i,j,k-1);
		    cR0 = get_cell(i,j,k);
		    cR1 = get_cell(i,j,k+1);
		    // Interpolate LEFT and RIGHT interface states from the cell-center properties.
		    if ( (k == kmin+1) && (bc[Face.bottom].ghost_cell_data_available == false) ) {
			one_d.interp_right(IFace, cL0, cR0, cR1, 
					   cL0.kLength, cR0.kLength, cR1.kLength,
					   Lft, Rght);
		    } else if ( (k == kmax) && (bc[Face.top].ghost_cell_data_available == false) ) {
			one_d.interp_left(IFace, cL1, cL0, cR0, 
					  cL1.kLength, cL0.kLength, cR0.kLength,
					  Lft, Rght);
		    } else { // General symmetric reconstruction.
			one_d.interp_both(IFace, cL1, cL0, cR0, cR1,
					  cL1.kLength, cL0.kLength, cR0.kLength, cR1.kLength,
					  Lft, Rght);
		    }
		    // Second, save u, v, w, T for the viscous flux calculation by making a local average.
		    // The values for u, v and T may be updated subsequently by the interface-flux function.
		    if ( (k == kmin) && (bc[Face.bottom].ghost_cell_data_available == false) ) {
			IFace.fs.copy_average_values_from(Rght, Rght);
		    } else if ( (k == kmax+1) && (bc[Face.top].ghost_cell_data_available == false) ) {
			IFace.fs.copy_average_values_from(Lft, Lft);
		    } else {
			IFace.fs.copy_average_values_from(Lft, Rght);
		    }
		    compute_interface_flux(Lft, Rght, IFace, myConfig.gmodel, omegaz);
		} // for k 
	    } // j loop
	} // i loop
	return;
    } // end convective_flux()

    @nogc
    override void flow_property_derivatives(int gtl)
    {
	if (myConfig.dimensions == 2) {
	    size_t k = 0;
	    for ( size_t i = imin; i <= imax+1; ++i ) {
		for ( size_t j = jmin; j <= jmax+1; ++j ) {
		    FVVertex vtx = get_vtx(i,j);
		    // gradients_xy_leastsq(vtx, myConfig.diffusion);
		    gradients_xy_div(vtx, myConfig.diffusion);
		}
	    }
	} else {
	    // Flow quantity derivatives for 3D.
	    for ( size_t i = imin; i <= imax+1; ++i ) {
		for ( size_t j = jmin; j <= jmax+1; ++j ) {
		    for ( size_t k = kmin; k <= kmax+1; ++k ) {
			FVVertex vtx = get_vtx(i,j,k);
			gradients_xyz_leastsq(vtx, myConfig.diffusion);
		    }
		}
	    }
	} // end if (myConfig.dimensions
    } // end flow_property_derivatives()

    override void applyPreReconAction(double t, int gtl, int ftl)
    {
	bc[Face.north].applyPreReconAction(t, gtl, ftl);
	bc[Face.east].applyPreReconAction(t, gtl, ftl);
	bc[Face.south].applyPreReconAction(t, gtl, ftl);
	bc[Face.west].applyPreReconAction(t, gtl, ftl);
	if ( myConfig.dimensions == 3 ) {
	    bc[Face.top].applyPreReconAction(t, gtl, ftl);
	    bc[Face.bottom].applyPreReconAction(t, gtl, ftl);
	}
    } // end applyPreReconAction()

    override void applyPreSpatialDerivAction(double t, int gtl, int ftl)
    {
	bc[Face.north].applyPreSpatialDerivAction(t, gtl, ftl);
	bc[Face.east].applyPreSpatialDerivAction(t, gtl, ftl);
	bc[Face.south].applyPreSpatialDerivAction(t, gtl, ftl);
	bc[Face.west].applyPreSpatialDerivAction(t, gtl, ftl);
	if ( myConfig.dimensions == 3 ) {
	    bc[Face.top].applyPreSpatialDerivAction(t, gtl, ftl);
	    bc[Face.bottom].applyPreSpatialDerivAction(t, gtl, ftl);
	}
    } // end applyPreSpatialDerivAction()

    override void applyPostDiffFluxAction(double t, int gtl, int ftl)
    {
	bc[Face.north].applyPostDiffFluxAction(t, gtl, ftl);
	bc[Face.east].applyPostDiffFluxAction(t, gtl, ftl);
	bc[Face.south].applyPostDiffFluxAction(t, gtl, ftl);
	bc[Face.west].applyPostDiffFluxAction(t, gtl, ftl);
	if ( myConfig.dimensions == 3 ) {
	    bc[Face.top].applyPostDiffFluxAction(t, gtl, ftl);
	    bc[Face.bottom].applyPostDiffFluxAction(t, gtl, ftl);
	}
    } // end applyPreSpatialDerivAction()

    @nogc
    void copy_into_ghost_cells(int destination_face,
			       ref SBlock src_blk, int src_face, int src_orientation,
			       int type_of_copy, bool with_encode)
    {
	size_t i_dest, i_src, j_dest, j_src, k_dest, k_src, i, j, k;
	FVCell src0, dest0, src1, dest1;
	int gtl = 0; // Do the encode only for grid-time-level zero.
	//
	@nogc
	void copy_pair_of_cells(FVCell src0, FVCell dest0, 
				FVCell src1, FVCell dest1,
				bool with_encode)
	{
	    dest0.copy_values_from(src0, type_of_copy);
	    if (with_encode) dest0.encode_conserved(gtl, 0, omegaz);
	    dest1.copy_values_from(src1, type_of_copy);
	    if (with_encode) dest1.encode_conserved(gtl, 0, omegaz);
	}
	//
	if (myConfig.dimensions == 2) {
	    // Handle the 2D case separately.
	    switch (destination_face) {
	    case Face.north:
		j_dest = jmax;  // index of the north-most plane of active cells
		for (i = 0; i < nicell; ++i) {
		    i_dest = i + imin;
		    switch (src_face) {
		    case Face.north:
			i_src = (src_blk.nicell - i - 1) + src_blk.imin;
			j_src = src_blk.jmax; 
			src0 = src_blk.get_cell(i_src,j_src);
			src1 = src_blk.get_cell(i_src,j_src-1);
			break;
		    case Face.east:
			j_src = i + src_blk.jmin;
			i_src = src_blk.imax; 
			src0 = src_blk.get_cell(i_src,j_src);
			src1 = src_blk.get_cell(i_src-1,j_src);
			break;
		    case Face.south:
			i_src = i + src_blk.imin;
			j_src = src_blk.jmin; 
			src0 = src_blk.get_cell(i_src,j_src);
			src1 = src_blk.get_cell(i_src,j_src+1);
			break;
		    case Face.west:
			j_src = (src_blk.njcell - i - 1) + src_blk.jmin;
			i_src = src_blk.imin; 
			src0 = src_blk.get_cell(i_src,j_src,k_src);
			src1 = src_blk.get_cell(i_src+1,j_src,k_src);
			break;
		    default:
			assert(false, "Incorrect boundary connection, source face.");
		    } // end switch src_face
		    dest0 = get_cell(i_dest,j_dest+1);
		    dest1 = get_cell(i_dest,j_dest+2);
		    copy_pair_of_cells(src0, dest0, src1, dest1, with_encode);
		} // i loop
		break;
	    case Face.east:
		i_dest = imax;  // index of the east-most plane of active cells
		for (j = 0; j < njcell; ++j) {
		    j_dest = j + jmin;
		    switch (src_face) {
		    case Face.north:
			i_src = j + src_blk.imin;
			j_src = src_blk.jmax; 
			src0 = src_blk.get_cell(i_src,j_src);
			src1 = src_blk.get_cell(i_src,j_src-1);
			break;
		    case Face.east:
			j_src = (src_blk.njcell - j - 1) + src_blk.jmin;
			i_src = src_blk.imax; 
			src0 = src_blk.get_cell(i_src,j_src);
			src1 = src_blk.get_cell(i_src-1,j_src);
			break;
		    case Face.south:
			i_src = (src_blk.nicell - j - 1) + src_blk.imin;
			j_src = src_blk.jmin; 
			src0 = src_blk.get_cell(i_src,j_src);
			src1 = src_blk.get_cell(i_src,j_src+1);
			break;
		    case Face.west:
			j_src = j + src_blk.jmin;
			i_src = src_blk.imin; 
			src0 = src_blk.get_cell(i_src,j_src);
			src1 = src_blk.get_cell(i_src+1,j_src);
			break;
		    default:
			assert(false, "Incorrect boundary connection, source face.");
		    } // end switch src_face
		    dest0 = get_cell(i_dest+1,j_dest);
		    dest1 = get_cell(i_dest+2,j_dest);
		    copy_pair_of_cells(src0, dest0, src1, dest1, with_encode);
		} // j loop
		break;
	    case Face.south:
		j_dest = jmin;  // index of the south-most plane of active cells
		for (i = 0; i < nicell; ++i) {
		    i_dest = i + imin;
		    switch (src_face) {
		    case Face.north:
			i_src = i + src_blk.imin;
			j_src = src_blk.jmax; 
			src0 = src_blk.get_cell(i_src,j_src);
			src1 = src_blk.get_cell(i_src,j_src-1);
			break;
		    case Face.east:
			j_src = (src_blk.njcell - i - 1) + src_blk.jmin;
			i_src = src_blk.imax; 
			src0 = src_blk.get_cell(i_src,j_src);
			src1 = src_blk.get_cell(i_src-1,j_src);
			break;
		    case Face.south:
			i_src = (src_blk.nicell - i - 1) + src_blk.imin;
			j_src = src_blk.jmin; 
			src0 = src_blk.get_cell(i_src,j_src);
			src1 = src_blk.get_cell(i_src,j_src+1);
			break;
		    case Face.west:
			j_src = i + src_blk.jmin;
			i_src = src_blk.imin; 
			src0 = src_blk.get_cell(i_src,j_src);
			src1 = src_blk.get_cell(i_src+1,j_src);
			break;
		    default:
			assert(false, "Incorrect boundary connection, source face.");
		    } // end switch src_face
		    dest0 = get_cell(i_dest,j_dest-1);
		    dest1 = get_cell(i_dest,j_dest-2);
		    copy_pair_of_cells(src0, dest0, src1, dest1, with_encode);
		} // i loop
		break;
	    case Face.west:
		i_dest = imin;  // index of the west-most plane of active cells
		for (j = 0; j < njcell; ++j) {
		    j_dest = j + jmin;
		    switch (src_face) {
		    case Face.north:
			i_src = (src_blk.nicell - j - 1) + src_blk.imin;
			j_src = src_blk.jmax; 
			src0 = src_blk.get_cell(i_src,j_src);
			src1 = src_blk.get_cell(i_src,j_src-1);
			break;
		    case Face.east:
			j_src = j + src_blk.jmin;
			i_src = src_blk.imax; 
			src0 = src_blk.get_cell(i_src,j_src);
			src1 = src_blk.get_cell(i_src-1,j_src);
			break;
		    case Face.south:
			i_src = j + src_blk.imin;
			j_src = src_blk.jmin; 
			src0 = src_blk.get_cell(i_src,j_src);
			src1 = src_blk.get_cell(i_src,j_src+1);
			break;
		    case Face.west:
			j_src = (src_blk.njcell - j - 1) + src_blk.jmin;
			i_src = src_blk.imin; 
			src0 = src_blk.get_cell(i_src,j_src);
			src1 = src_blk.get_cell(i_src+1,j_src);
			break;
		    default:
			assert(false, "Incorrect boundary connection, source face.");
		    } // end switch src_face
		    dest0 = get_cell(i_dest-1,j_dest);
		    dest1 = get_cell(i_dest-2,j_dest);
		    copy_pair_of_cells(src0, dest0, src1, dest1, with_encode);
		} // j loop
		break;
	    default:
		assert(false, "Incorrect boundary connection, destination face.");
	    } // end switch destination_face
	    return;
	} // end if dimensions == 2

	// Continue on with 3D work...
	final switch (destination_face) {
	case Face.north:
	    j_dest = jmax;  // index of the north-most plane of active cells
	    for (i = 0; i < nicell; ++i) {
		i_dest = i + imin;
		for (k = 0; k < nkcell; ++k) {
		    k_dest = k + kmin;
		    final switch (src_face) {
		    case Face.north:
			final switch (src_orientation) {
			case 0: i_src = src_blk.nicell - i - 1; k_src = k; break;
			case 1: i_src = k; k_src = i; break;
			case 2: i_src = i; k_src = src_blk.nkcell - k - 1; break;
			case 3: i_src = src_blk.nicell - k - 1; k_src = src_blk.nkcell - i - 1;
			} // end switch (src_orientation)
			j_src = src_blk.jmax; 
			i_src += src_blk.imin; k_src += src_blk.kmin;
			src0 = src_blk.get_cell(i_src,j_src,k_src);
			src1 = src_blk.get_cell(i_src,j_src-1,k_src);
			break;
		    case Face.east:
			final switch (src_orientation) {
			case 0: j_src = i; k_src = k; break;
			case 1: j_src = src_blk.njcell - k - 1; k_src = i; break;
			case 2: j_src = src_blk.njcell - i - 1; k_src = src_blk.nkcell - k - 1; break;
			case 3: j_src = k; k_src = src_blk.nkcell - i - 1;
			} // end switch (src_orientation)
			i_src = src_blk.imax; 
			j_src += src_blk.jmin; k_src += src_blk.kmin;
			src0 = src_blk.get_cell(i_src,j_src,k_src);
			src1 = src_blk.get_cell(i_src-1,j_src,k_src);
			break;
		    case Face.south:
			final switch (src_orientation) {
			case 0: i_src = i; k_src = k; break;
			case 1: i_src = k; k_src = src_blk.nkcell - i - 1; break;
			case 2: i_src = src_blk.nicell - i - 1; k_src = src_blk.nkcell - k - 1; break;
			case 3: i_src = src_blk.nicell - k - 1; k_src = i;
			} // end switch (src_orientation)
			j_src = src_blk.jmin; 
			i_src += src_blk.imin; k_src += src_blk.kmin;
			src0 = src_blk.get_cell(i_src,j_src,k_src);
			src1 = src_blk.get_cell(i_src,j_src+1,k_src);
			break;
		    case Face.west:
			final switch (src_orientation) {
			case 0: j_src = src_blk.njcell - i - 1; k_src = k; break;
			case 1: j_src = k; k_src = i; break;
			case 2: j_src = i; k_src = src_blk.nkcell - k - 1; break;
			case 3: j_src = src_blk.njcell - k - 1; k_src = src_blk.nkcell - i - 1;
			} // end switch (src_orientation)
			i_src = src_blk.imin; 
			j_src += src_blk.jmin; k_src += src_blk.kmin;
			src0 = src_blk.get_cell(i_src,j_src,k_src);
			src1 = src_blk.get_cell(i_src+1,j_src,k_src);
			break;
		    case Face.top:
			final switch (src_orientation) {
			case 0: i_src = i; j_src = k; break;
			case 1: i_src = src_blk.nicell - k - 1; j_src = i; break;
			case 2: i_src = src_blk.nicell - i - 1; j_src = src_blk.njcell - k - 1; break;
			case 3: i_src = k; j_src = src_blk.njcell - i - 1;
			} // end switch (src_orientation)
			k_src = src_blk.kmax; 
			i_src += src_blk.imin; j_src += src_blk.jmin;
			src0 = src_blk.get_cell(i_src,j_src,k_src);
			src1 = src_blk.get_cell(i_src,j_src,k_src-1);
			break;
		    case Face.bottom:
			final switch (src_orientation) {
			case 0: i_src = src_blk.nicell - i - 1; j_src = k; break;
			case 1: i_src = k; j_src = i; break;
			case 2: i_src = src_blk.nicell - i - 1; j_src = src_blk.njcell - k - 1; break;
			case 3: i_src = src_blk.nicell - k - 1; j_src = src_blk.njcell - i - 1;
			} // end switch (src_orientation)
			k_src = src_blk.kmin; 
			i_src += src_blk.imin; j_src += src_blk.jmin;
			src0 = src_blk.get_cell(i_src,j_src,k_src);
			src1 = src_blk.get_cell(i_src,j_src,k_src+1);
		    } // end switch (src_face)
		    dest0 = get_cell(i_dest,j_dest+1,k_dest);
		    dest1 = get_cell(i_dest,j_dest+2,k_dest);
		    copy_pair_of_cells(src0, dest0, src1, dest1, with_encode);
		} // k loop
	    } // i loop
	    break;
	case Face.east:
	    i_dest = imax;  // index of the east-most plane of active cells
	    for (j = 0; j < njcell; ++j) {
		j_dest = j + jmin;
		for (k = 0; k < nkcell; ++k) {
		    k_dest = k + kmin;
		    final switch (src_face) {
		    case Face.north:
			final switch (src_orientation) {
			case 0: i_src = j; k_src = k; break;
			case 1: i_src = k; k_src = src_blk.nkcell - j - 1; break;
			case 2: i_src = src_blk.nicell - j - 1; k_src = src_blk.nkcell - k - 1; break;
			case 3: i_src = src_blk.nicell - k - 1; k_src = j;
			} // end switch (src_orientation)
			j_src = src_blk.jmax; 
			i_src += src_blk.imin; k_src += src_blk.kmin;
			src0 = src_blk.get_cell(i_src,j_src,k_src);
			src1 = src_blk.get_cell(i_src,j_src-1,k_src);
			break;
		    case Face.east:
			final switch (src_orientation) {
			case 0: j_src = src_blk.njcell - j - 1; k_src = k; break;
			case 1: j_src = src_blk.njcell - k - 1; k_src = src_blk.nkcell - j - 1; break;
			case 2: j_src = src_blk.njcell - j - 1; k_src = src_blk.nkcell - k - 1; break;
			case 3: j_src = j; k_src = src_blk.nkcell - k - 1;
			}
			i_src = src_blk.imax; 
			j_src += src_blk.jmin; k_src += src_blk.kmin;
			src0 = src_blk.get_cell(i_src,j_src,k_src);
			src1 = src_blk.get_cell(i_src-1,j_src,k_src);
			break;
		    case Face.south:
			final switch (src_orientation) {
			case 0: i_src = src_blk.nicell - j - 1; k_src = k; break;
			case 1: i_src = src_blk.nicell - k - 1; k_src = src_blk.nkcell - j - 1; break;
			case 2: i_src = j; k_src = src_blk.nkcell - k - 1; break;
			case 3: i_src = k; k_src = j;
			}
			j_src = src_blk.jmin; 
			i_src += src_blk.imin; k_src += src_blk.kmin;
			src0 = src_blk.get_cell(i_src,j_src,k_src);
			src1 = src_blk.get_cell(i_src,j_src+1,k_src);
			break;
		    case Face.west:
			final switch (src_orientation) {
			case 0: j_src = j; k_src = k; break;
			case 1: j_src = k; k_src = src_blk.nkcell - j - 1; break;
			case 2: j_src = src_blk.njcell - j - 1; k_src = src_blk.nkcell - k - 1; break;
			case 3: j_src = src_blk.njcell - k - 1; k_src = j;
			}
			i_src = src_blk.imin; 
			j_src += src_blk.jmin; k_src += src_blk.kmin;
			src0 = src_blk.get_cell(i_src,j_src,k_src);
			src1 = src_blk.get_cell(i_src+1,j_src,k_src);
			break;
		    case Face.top:
			final switch (src_orientation) {
			case 0: i_src = src_blk.nicell - j - 1; j_src = k; break;
			case 1: i_src = src_blk.nicell - k - 1; j_src = src_blk.njcell - j - 1; break;
			case 2: i_src = j; j_src = src_blk.njcell - k - 1; break;
			case 3: i_src = k; j_src = j;
			}
			k_src = src_blk.kmax; 
			i_src += src_blk.imin; j_src += src_blk.jmin;
			src0 = src_blk.get_cell(i_src,j_src,k_src);
			src1 = src_blk.get_cell(i_src,j_src+1,k_src-1);
			break;
		    case Face.bottom:
			final switch (src_orientation) {
			case 0: i_src = j; j_src = k; break;
			case 1: i_src = k; j_src = src_blk.njcell - j - 1; break;
			case 2: i_src = src_blk.nicell - j - 1; j_src = src_blk.njcell - k - 1; break;
			case 3: i_src = k; j_src = j;
			}
			k_src = src_blk.kmin; 
			i_src += src_blk.imin; j_src += src_blk.jmin;
			src0 = src_blk.get_cell(i_src,j_src,k_src);
			src1 = src_blk.get_cell(i_src,j_src+1,k_src+1);
		    } // end switch (src_face)
		    dest0 = get_cell(i_dest+1,j_dest,k_dest);
		    dest1 = get_cell(i_dest+2,j_dest,k_dest);
		    copy_pair_of_cells(src0, dest0, src1, dest1, with_encode);
		} // k loop
	    } // j loop
	    break;
	case Face.south:
	    j_dest = jmin;  // index of the south-most plane of active cells
	    for (i = 0; i < nicell; ++i) {
		i_dest = i + imin;
		for (k = 0; k < nkcell; ++k) {
		    k_dest = k + kmin;
		    final switch (src_face) {
		    case Face.north:
			final switch (src_orientation) {
			case 0: i_src = i; k_src = k; break;
			case 1: i_src = src_blk.nicell - k - 1; k_src = i; break;
			case 2: i_src = src_blk.nicell - i - 1; k_src = src_blk.nkcell - k - 1; break;
			case 3: i_src = k; k_src = src_blk.nkcell - i - 1;
			} // end switch (src_orientation)
			j_src = src_blk.jmax; 
			i_src += src_blk.imin; k_src += src_blk.kmin;
			src0 = src_blk.get_cell(i_src,j_src,k_src);
			src1 = src_blk.get_cell(i_src,j_src-1,k_src);
			break;
		    case Face.east:
			final switch (src_orientation) {
			case 0: j_src = src_blk.njcell - i - 1; k_src = k; break;
			case 1: j_src = src_blk.njcell - k - 1; k_src = src_blk.nkcell - i - 1; break;
			case 2: j_src = i; k_src = src_blk.nkcell - k - 1; break;
			case 3: j_src = k; k_src = i;
			} // end switch (src_orientation)
			i_src = src_blk.imax; 
			j_src += src_blk.jmin; k_src += src_blk.kmin;
			src0 = src_blk.get_cell(i_src,j_src,k_src);
			src1 = src_blk.get_cell(i_src-1,j_src,k_src);
			break;
		    case Face.south:
			final switch (src_orientation) {
			case 0: i_src = src_blk.nicell - i - 1; k_src = k; break;
			case 1: i_src = src_blk.nicell - k - 1; k_src = src_blk.nkcell - i - 1; break;
			case 2: i_src = i; k_src = src_blk.nkcell - k - 1; break;
			case 3: i_src = k; k_src = i;
			} // end switch (src_orientation)
			j_src = src_blk.jmin; 
			i_src += src_blk.imin; k_src += src_blk.kmin;
			src0 = src_blk.get_cell(i_src,j_src,k_src);
			src1 = src_blk.get_cell(i_src,j_src+1,k_src);
			break;
		    case Face.west:
			final switch (src_orientation) {
			case 0: j_src = i; k_src = k; break;
			case 1: j_src = k; k_src = src_blk.nkcell - i - 1; break;
			case 2: j_src = src_blk.njcell - i - 1; k_src = src_blk.nkcell - k - 1; break;
			case 3: j_src = src_blk.njcell - k - 1; k_src = i;
			} // end switch (src_orientation)
			i_src = src_blk.imin; 
			j_src += src_blk.jmin; k_src += src_blk.kmin;
			src0 = src_blk.get_cell(i_src,j_src,k_src);
			src1 = src_blk.get_cell(i_src+1,j_src,k_src);
			break;
		    case Face.top:
			final switch (src_orientation) {
			case 0: i_src = src_blk.nicell - i - 1; j_src = k; break;
			case 1: i_src = src_blk.nicell - k - 1; j_src = src_blk.njcell - i - 1; break;
			case 2: i_src = src_blk.nicell - i - 1; j_src = src_blk.njcell - k - 1; break;
			case 3: i_src = k; j_src = i;
			} // end switch (src_orientation)
			k_src = src_blk.kmax; 
			i_src += src_blk.imin; j_src += src_blk.jmin;
			src0 = src_blk.get_cell(i_src,j_src,k_src);
			src1 = src_blk.get_cell(i_src,j_src,k_src-1);
			break;
		    case Face.bottom:
			final switch (src_orientation) {
			case 0: i_src = i; j_src = k; break;
			case 1: i_src = k; j_src = src_blk.njcell - i - 1; break;
			case 2: i_src = src_blk.nicell - i - 1; j_src = src_blk.njcell - k - 1; break;
			case 3: i_src = src_blk.nicell - k - 1; j_src = i;
			} // end switch (src_orientation)
			k_src = src_blk.kmin; 
			i_src += src_blk.imin; j_src += src_blk.jmin;
			src0 = src_blk.get_cell(i_src,j_src,k_src);
			src1 = src_blk.get_cell(i_src,j_src,k_src+1);
		    } // end switch (src_face)
		    dest0 = get_cell(i_dest,j_dest-1,k_dest);
		    dest1 = get_cell(i_dest,j_dest-2,k_dest);
		    copy_pair_of_cells(src0, dest0, src1, dest1, with_encode);
		} // k loop
	    } // i loop
	    break;
	case Face.west:
	    i_dest = imin;  // index of the west-most plane of active cells
	    for (j = 0; j < njcell; ++j) {
		j_dest = j + jmin;
		for (k = 0; k < nkcell; ++k) {
		    k_dest = k + kmin;
		    final switch (src_face) {
		    case Face.north:
			final switch (src_orientation) {
			case 0: i_src = src_blk.nicell - j - 1; k_src = k; break;
			case 1: i_src = k; k_src = j; break;
			case 2: i_src = j; k_src = src_blk.nkcell - k - 1; break;
			case 3: i_src = src_blk.nicell - k - 1; k_src = src_blk.nkcell - j - 1;
			} // end switch (src_orientation)
			j_src = src_blk.jmax; 
			i_src += src_blk.imin; k_src += src_blk.kmin;
			src0 = src_blk.get_cell(i_src,j_src,k_src);
			src1 = src_blk.get_cell(i_src,j_src-1,k_src);
			break;
		    case Face.east:
			final switch (src_orientation) {
			case 0: j_src = j; k_src = k; break;
			case 1: j_src = src_blk.njcell - k - 1; k_src = j; break;
			case 2: j_src = src_blk.njcell - j - 1; k_src = src_blk.nkcell - k - 1; break;
			case 3: j_src = k; k_src = src_blk.nkcell - j - 1;
			} // end switch (src_orientation)
			i_src = src_blk.imax; 
			j_src += src_blk.jmin; k_src += src_blk.kmin;
			src0 = src_blk.get_cell(i_src,j_src,k_src);
			src1 = src_blk.get_cell(i_src-1,j_src,k_src);
			break;
		    case Face.south:
			final switch (src_orientation) {
			case 0: i_src = j; k_src = k; break;
			case 1: i_src = src_blk.nicell - k - 1; k_src = j; break;
			case 2: i_src = src_blk.nicell - j - 1; k_src = src_blk.nkcell - k - 1; break;
			case 3: i_src = k; k_src = src_blk.nkcell - j - 1;
			} // end switch (src_orientation)
			j_src = src_blk.jmin; 
			i_src += src_blk.imin; k_src += src_blk.kmin;
			src0 = src_blk.get_cell(i_src,j_src,k_src);
			src1 = src_blk.get_cell(i_src,j_src+1,k_src);
			break;
		    case Face.west:
			final switch (src_orientation) {
			case 0: j_src = src_blk.njcell - j - 1; k_src = src_blk.nkcell - k - 1; break;
			case 1: j_src = k; k_src = j; break;
			case 2: j_src = j; k_src = src_blk.nkcell - k - 1; break;
			case 3: j_src = src_blk.njcell - k - 1; k_src = src_blk.nkcell - j - 1;
			} // end switch (src_orientation)
			i_src = src_blk.imin; 
			j_src += src_blk.jmin; k_src += src_blk.kmin;
			src0 = src_blk.get_cell(i_src,j_src,k_src);
			src1 = src_blk.get_cell(i_src+1,j_src,k_src);
			break;
		    case Face.top:
			final switch (src_orientation) {
			case 0: i_src = j; j_src = k; break;
			case 1: i_src = src_blk.nicell - k - 1; j_src = j; break;
			case 2: i_src = src_blk.nicell - j - 1; j_src = src_blk.njcell - k - 1; break;
			case 3: i_src = k; j_src = src_blk.njcell - j - 1;
			} // end switch (src_orientation)
			k_src = src_blk.kmax; 
			i_src += src_blk.imin; j_src += src_blk.jmin;
			src0 = src_blk.get_cell(i_src,j_src,k_src);
			src1 = src_blk.get_cell(i_src,j_src,k_src-1);
			break;
		    case Face.bottom:
			final switch (src_orientation) {
			case 0: i_src = src_blk.nicell - j - 1; j_src = k; break;
			case 1: i_src = k; j_src = j; break;
			case 2: i_src = j; j_src = src_blk.njcell - k - 1; break;
			case 3: i_src = src_blk.nicell - k - 1; j_src = src_blk.njcell - j - 1;
			} // end switch (src_orientation)
			k_src = src_blk.kmin; 
			i_src += src_blk.imin; j_src += src_blk.jmin;
			src0 = src_blk.get_cell(i_src,j_src,k_src);
			src1 = src_blk.get_cell(i_src,j_src,k_src+1);
		    } // end switch (src_face)
		    dest0 = get_cell(i_dest-1,j_dest,k_dest);
		    dest1 = get_cell(i_dest-2,j_dest,k_dest);
		    copy_pair_of_cells(src0, dest0, src1, dest1, with_encode);
		} // k loop
	    } // j loop
	    break;
	case Face.top:
	    k_dest = kmax;  // index of the top-most plane of active cells
	    for (j = 0; j < njcell; ++j) {
		j_dest = j + jmin;
		for (i = 0; i < nicell; ++i) {
		    i_dest = i + imin;
		    final switch (src_face) {
		    case Face.north:
			final switch (src_orientation) {
			case 0: i_src = i; k_src = j; break;
			case 1: i_src = j; k_src = src_blk.nkcell - i - 1; break;
			case 2: i_src = src_blk.nicell - i - 1; k_src = src_blk.nkcell - j - 1; break;
			case 3: i_src = src_blk.nicell - j - 1; k_src = i;
			} // end switch (src_orientation)
			j_src = src_blk.jmax; 
			i_src += src_blk.imin; k_src += src_blk.kmin;
			src0 = src_blk.get_cell(i_src,j_src,k_src);
			src1 = src_blk.get_cell(i_src,j_src-1,k_src);
			break;
		    case Face.east:
			final switch (src_orientation) {
			case 0: j_src = src_blk.njcell - i - 1; k_src = j; break;
			case 1: j_src = src_blk.njcell - j - 1; k_src = src_blk.nkcell - i - 1; break;
			case 2: j_src = i; k_src = src_blk.nkcell - j - 1; break;
			case 3: j_src = j; k_src = i;
			} // end switch (src_orientation)
			i_src = src_blk.imax; 
			j_src += src_blk.jmin; k_src += src_blk.kmin;
			src0 = src_blk.get_cell(i_src,j_src,k_src);
			src1 = src_blk.get_cell(i_src-1,j_src,k_src);
			break;
		    case Face.south:
			final switch (src_orientation) {
			case 0: i_src = src_blk.nicell - i - 1; k_src = j; break;
			case 1: i_src = src_blk.nicell - j - 1; k_src = src_blk.nkcell - i - 1; break;
			case 2: i_src = i; k_src = src_blk.nkcell - j - 1; break;
			case 3: i_src = j; k_src = i;
			} // end switch (src_orientation)
			j_src = src_blk.jmin; 
			i_src += src_blk.imin; k_src += src_blk.kmin;
			src0 = src_blk.get_cell(i_src,j_src,k_src);
			src1 = src_blk.get_cell(i_src,j_src+1,k_src);
			break;
		    case Face.west:
			final switch (src_orientation) {
			case 0: j_src = i; k_src = j; break;
			case 1: j_src = j; k_src = src_blk.nkcell - i - 1; break;
			case 2: j_src = src_blk.njcell - i - 1; k_src = src_blk.nkcell - j - 1; break;
			case 3: j_src = src_blk.njcell - j - 1; k_src = i;
			} // end switch (src_orientation)
			i_src = src_blk.imin; 
			j_src += src_blk.jmin; k_src += src_blk.kmin;
			src0 = src_blk.get_cell(i_src,j_src,k_src);
			src1 = src_blk.get_cell(i_src+1,j_src,k_src);
			break;
		    case Face.top:
			final switch (src_orientation) {
			case 0: i_src = src_blk.nicell - i - 1; j_src = j; break;
			case 1: i_src = j; j_src = src_blk.njcell - i - 1; break;
			case 2: i_src = i; j_src = src_blk.njcell - j - 1; break;
			case 3: i_src = j; j_src = i;
			} // end switch (src_orientation)
			k_src = src_blk.kmax; 
			i_src += src_blk.imin; j_src += src_blk.jmin;
			src0 = src_blk.get_cell(i_src,j_src,k_src);
			src1 = src_blk.get_cell(i_src,j_src,k_src-1);
			break;
		    case Face.bottom:
			final switch (src_orientation) {
			case 0: i_src = i; j_src = j; break;
			case 1: i_src = j; j_src = src_blk.njcell - i - 1; break;
			case 2: i_src = src_blk.nicell - i - 1; j_src = src_blk.njcell - j - 1; break;
			case 3: i_src = src_blk.nicell - j - 1; j_src = i;
			} // end switch (src_orientation)
			k_src = src_blk.kmin; 
			i_src += src_blk.imin; j_src += src_blk.jmin;
			src0 = src_blk.get_cell(i_src,j_src,k_src);
			src1 = src_blk.get_cell(i_src,j_src,k_src+1);
		    } // end switch (src_face)
		    dest0 = get_cell(i_dest,j_dest,k_dest+1);
		    dest1 = get_cell(i_dest,j_dest,k_dest+2);
		    copy_pair_of_cells(src0, dest0, src1, dest1, with_encode);
		} // i loop
	    } // j loop
	    break;
	case Face.bottom:
	    k_dest = kmin;  // index of the bottom-most plane of active cells
	    for (j = 0; j < njcell; ++j) {
		j_dest = j + jmin;
		for (i = 0; i < nicell; ++i) {
		    i_dest = i + imin;
		    final switch (src_face) {
		    case Face.north:
			final switch (src_orientation) {
			case 0: i_src = src_blk.nicell - i - 1; k_src = j; break;
			case 1: i_src = j; k_src = i; break;
			case 2: i_src = i; k_src = src_blk.nkcell - j - 1; break;
			case 3: i_src = src_blk.nicell - j - 1; k_src = src_blk.nkcell - i - 1;
			} // end switch (src_orientation)
			j_src = src_blk.jmax; 
			i_src += src_blk.imin; k_src += src_blk.kmin;
			src0 = src_blk.get_cell(i_src,j_src,k_src);
			src1 = src_blk.get_cell(i_src,j_src-1,k_src);
			break;
		    case Face.east:
			final switch (src_orientation) {
			case 0: j_src = i; k_src = j; break;
			case 1: j_src = src_blk.njcell - j - 1; k_src = i; break;
			case 2: j_src = src_blk.njcell - i - 1; k_src = src_blk.nkcell - j - 1; break;
			case 3: j_src = j; k_src = src_blk.nkcell - i - 1;
			} // end switch (src_orientation)
			i_src = src_blk.imax; 
			j_src += src_blk.jmin; k_src += src_blk.kmin;
			src0 = src_blk.get_cell(i_src,j_src,k_src);
			src1 = src_blk.get_cell(i_src-1,j_src,k_src);
			break;
		    case Face.south:
			final switch (src_orientation) {
			case 0: i_src = src_blk.nicell - i - 1; k_src = j; break;
			case 1: i_src = src_blk.nicell - j - 1; k_src = i; break;
			case 2: i_src = src_blk.nicell - i - 1; k_src = src_blk.nkcell - j - 1; break;
			case 3: i_src = j; k_src = src_blk.nkcell - i - 1;
			} // end switch (src_orientation)
			j_src = src_blk.jmin; 
			i_src += src_blk.imin; k_src += src_blk.kmin;
			src0 = src_blk.get_cell(i_src,j_src,k_src);
			src1 = src_blk.get_cell(i_src,j_src+1,k_src);
			break;
		    case Face.west:
			final switch (src_orientation) {
			case 0: j_src = src_blk.njcell - i - 1; k_src = j; break;
			case 1: j_src = j; k_src = i; break;
			case 2: j_src = i; k_src = src_blk.nkcell - j - 1; break;
			case 3: j_src = src_blk.njcell - j - 1; k_src = src_blk.nkcell - i - 1;
			} // end switch (src_orientation)
			i_src = src_blk.imin; 
			j_src += src_blk.jmin; k_src += src_blk.kmin;
			src0 = src_blk.get_cell(i_src,j_src,k_src);
			src1 = src_blk.get_cell(i_src+1,j_src,k_src);
			break;
		    case Face.top:
			final switch (src_orientation) {
			case 0: i_src = i; j_src = j; break;
			case 1: i_src = src_blk.nicell - j - 1; j_src = i; break;
			case 2: i_src = src_blk.nicell - i - 1; j_src = src_blk.njcell - j - 1; break;
			case 3: i_src = j; j_src = src_blk.njcell - i - 1;
			} // end switch (src_orientation)
			k_src = src_blk.kmax; 
			i_src += src_blk.imin; j_src += src_blk.jmin;
			src0 = src_blk.get_cell(i_src,j_src,k_src);
			src1 = src_blk.get_cell(i_src,j_src,k_src-1);
			break;
		    case Face.bottom:
			final switch (src_orientation) {
			case 0: i_src = src_blk.nicell - i - 1; j_src = j; break;
			case 1: i_src = j; j_src = i; break;
			case 2: i_src = i; j_src = src_blk.njcell - j - 1; break;
			case 3: i_src = src_blk.nicell - j - 1; j_src = src_blk.njcell - i - 1;
			} // end switch (src_orientation)
			k_src = src_blk.kmin; 
			i_src += src_blk.imin; j_src += src_blk.jmin;
			src0 = src_blk.get_cell(i_src,j_src,k_src);
			src1 = src_blk.get_cell(i_src,j_src,k_src+1);
		    } // end switch src_face
		    dest0 = get_cell(i_dest,j_dest,k_dest-1);
		    dest1 = get_cell(i_dest,j_dest,k_dest-2);
		    copy_pair_of_cells(src0, dest0, src1, dest1, with_encode);
		} // i loop
	    } // j loop
	} // end switch destination_face
    } // end copy_to_ghost_cells()

    void propagate_inflow_data_west_to_east()
    {
	// Assume that the west-face ghost cells have appropriate data.
	for ( size_t k = kmin; k <= kmax; ++k ) {
	    for ( size_t j = jmin; j <= jmax; ++j ) {
		auto src_cell = get_cell(imin-1,j,k);
		for ( size_t i = imin; i <= imax; ++i ) {
		    auto dest_cell = get_cell(i,j,k);
		    dest_cell.copy_values_from(src_cell, CopyDataOption.all_flow);
		}
	    }
	}
    } // end propagate_inflow_data_west_to_east()

} // end class SBlock


/** Indexing of the data in 2D.
 *
 * \verbatim
 * The following figure shows cell [i,j] and its associated
 * vertices and faces. 
 * (New arrangement, planned August 2006, implemented Nov 2006)
 *
 *
 *
 *     Vertex 3         North face           Vertex 2 
 *   vtx[i,j+1]         ifj[i,j+1]           vtx[i+1,j+1]
 *             +--------------x--------------+
 *             |                             |
 *             |                             |
 *             |                             |
 *             |                             |
 *             |                             |
 *   West      |         cell center         |  East 
 *   face      |          ctr[i,j]           |  face
 *   ifi[i,j]  x              o              x  ifi[i+1,j]
 *             |                             |
 *             |                             |
 *             |                             |
 *             |                             |
 *             |                             |
 *             |                             |
 *             |                             |
 *             +--------------x--------------+
 *     Vertex 0           South face         Vertex 1
 *     vtx[i,j]           ifj[i,j]           vtx[i+1,j]
 *
 *
 * Thus...
 * ----
 * Active cells are indexed as ctr[i][i], where
 * imin <= i <= imax, jmin <= j <= jmax.
 *
 * Active east-facing interfaces are indexed as ifi[i][j], where
 * imin <= i <= imax+1, jmin <= j <= jmax.
 *
 * Active north-facing interfaces are indexed as ifj[i][j], where
 * imin <= i <= imax, jmin <= j <= jmax+1.
 *
 * Active vertices are indexed as vtx[i][j], where
 * imin <= i <= imax+1, jmin <= j <= jmax+1.
 *
 * Space for ghost cells is available outside these ranges.
 *
 * Indexing for the 3D data -- see page 8 in 3D CFD workbook
 * \endverbatim
 */
