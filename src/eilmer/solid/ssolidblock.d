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
import nm.complex;
import nm.number;
import std.json;
import util.lua;
import json_helper;
import gzip;
import geom;
import globalconfig;
import fluidblock;
import sfluidblock;
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
    SolidBoundaryCondition[] bc;

private:
    StructuredGrid grid; // for reading and writing

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
        foreach (boundary; 0 .. (myConfig.dimensions == 3 ? 6 : 4)) {
            string jsonKey = "face_" ~ face_name[boundary];
            auto bcJsonData = jsonData[jsonKey];
            bc ~= makeSolidBCFromJson(bcJsonData, id, boundary,
                                      nicell, njcell, nkcell);
        }
        foreach (bci; bc) bci.postBCconstruction();
    }

    void fillInOtherSizeData()
    // Helper function for the constructor
    {
        _nidim = nicell + 2 * n_ghost_cell_layers;
        _njdim = njcell + 2 * n_ghost_cell_layers;
        // Indices, in each grid direction for the active cells.
        // These limits are inclusive. The mincell and max cell
        // are both within the active set of cells.
        imin = n_ghost_cell_layers; imax = imin + nicell - 1;
        jmin = n_ghost_cell_layers; jmax = jmin + njcell - 1;
        if ( myConfig.dimensions == 2 ) {
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
            _nkdim = nkcell + 2 * n_ghost_cell_layers;
            kmin = n_ghost_cell_layers; kmax = kmin + nkcell - 1;
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
        if ( myConfig.verbosity_level >= 2 ) 
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
                _ctr ~= new SolidFVCell(myConfig); _ctr[gid].id = to!int(gid);
                auto ijk = toIJKIndices(gid);
                if ( ijk[0] >= imin && ijk[0] <= imax && 
                     ijk[1] >= jmin && ijk[1] <= jmax && 
                     ijk[2] >= kmin && ijk[2] <= kmax ) {
                    activeCells ~= _ctr[gid];
                }
                _ifi ~= new SolidFVInterface(); _ifi[gid].id = gid;
                _ifj ~= new SolidFVInterface(); _ifj[gid].id = gid;
                if ( myConfig.dimensions == 3 ) {
                    _ifk ~= new SolidFVInterface(); _ifk[gid].id = gid;
                }
                _vtx ~= new SolidFVVertex(); _vtx[gid].id = gid;
                _sifi ~= new SolidFVInterface(); _sifi[gid].id = gid;
                _sifj ~= new SolidFVInterface(); _sifj[gid].id = gid;
                if ( myConfig.dimensions == 3 ) {
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
        if ( myConfig.verbosity_level >= 2 ) {
            writefln("Done assembling arrays for %d solid cells.", ntot);
        }
    }
    override void bindFacesAndVerticesToCells()
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
                    SolidFVCell cell = getCell(i,j,k);
                    cell.iface ~= getIfj(i,j+1,k); // north
                    cell.iface ~= getIfi(i+1,j,k); // east
                    cell.iface ~= getIfj(i,j,k); // south
                    cell.iface ~= getIfi(i,j,k); // west
                    cell.vtx ~= getVtx(i,j,k);
                    cell.vtx ~= getVtx(i+1,j,k);
                    cell.vtx ~= getVtx(i+1,j+1,k);
                    cell.vtx ~= getVtx(i,j+1,k);
                    if ( myConfig.dimensions == 3 ) {
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
        if ( myConfig.verbosity_level >= 1 && id == 0 ) {
            writeln("SSolidBlock.readGrid(): Start block ", id);
        }
        grid = new StructuredGrid(filename, "gziptext");
        nivtx = grid.niv; njvtx = grid.njv; nkvtx = grid.nkv;
        if ( myConfig.dimensions == 3 ) {
            if ( nivtx-1 != nicell || njvtx-1 != njcell || nkvtx-1 != nkcell ) {
                throw new Error(text("For solid_block[", id, "] we have a mismatch in 3D grid size.",
                                     " Have read nivtx=", nivtx, " njvtx=", njvtx,
                                     " nkvtx=", nkvtx));
            }
            for ( size_t k = kmin; k <= kmax+1; ++k ) {
                for ( size_t j = jmin; j <= jmax+1; ++j ) {
                    for ( size_t i = imin; i <= imax+1; ++i ) {
                        auto src_vtx = grid[i-imin,j-jmin,k-kmin];
                        auto vtx = getVtx(i,j,k);
                        vtx.pos.refx = src_vtx.x;
                        vtx.pos.refy = src_vtx.y;
                        vtx.pos.refz = src_vtx.z;
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
                    auto vtx = getVtx(i,j);
                    auto src_vtx = grid[i-imin,j-jmin];
                    vtx.pos.refx = src_vtx.x;
                    vtx.pos.refy = src_vtx.y;
                    vtx.pos.refz = 0.0;
                } // for i
            } // for j
        }
    } // end read_grid()

    override void writeGrid(string filename, double sim_time)
    {
        if ( myConfig.verbosity_level >= 1 && id == 0 ) {
            writeln("SSolidBlock.writeGrid(): Start block ", id);
        }
        size_t kmaxrange;
        if ( myConfig.dimensions == 3 ) {
            kmaxrange = kmax + 1;
        } else { // 2D case
            kmaxrange = kmax;
        }
        for ( size_t k = kmin; k <= kmaxrange; ++k ) {
            for ( size_t j = jmin; j <= jmax+1; ++j ) {
                for ( size_t i = imin; i <= imax+1; ++i ) {
                    auto vtx = getVtx(i,j,k);
                    auto dest_vtx = grid[i-imin,j-jmin,k-kmin];
                    dest_vtx.refx = vtx.pos.x;
                    dest_vtx.refy = vtx.pos.y;
                    dest_vtx.refz = vtx.pos.z;
                } // for i
            } // for j
        } // for k
        grid.write_to_gzip_file(filename);
    } // end write_grid()

    override void readSolution(string filename)
    {
        size_t ni, nj, nk;
        double sim_time;
        if ( myConfig.verbosity_level >= 1 && id == 0 ) {
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
             nk != ((myConfig.dimensions == 3) ? nkcell : 1) ) {
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
        if ( myConfig.verbosity_level >= 1 && id == 0 ) {
            writeln("write_solution(): Start solid_block ", id);
        }
        auto outfile = new GzipOut(fileName);
        auto writer = appender!string();
        formattedWrite(writer, "%.18e\n", simTime);
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
        if ( myConfig.dimensions == 2 ) {
            calcVolumes2D();
            calcFaces2D();
            return;
        }
        throw new Error("SSolidBlock.computePrimaryCellGeometryData() not implemented for 3D.");

    } // end compute_primary_cell_geometric_data()

    void calcVolumes2D()
    {
        size_t i, j;
        number xA, yA, xB, yB, xC, yC, xD, yD;
        number xN, yN, xS, yS, xE, yE, xW, yW;
        number vol, max_vol, min_vol, xyarea;
        number dx, dy, dxN, dyN, dxE, dyE;
        number lengthN, lengthE, length_max, length_min, length_cross;
        number max_aspect, aspect_ratio;
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
                if ( myConfig.axisymmetric ) {
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
        number xA, xB, yA, yB, xC, yC;
        number LAB, LBC;

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
                if ( myConfig.axisymmetric ) {
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
                if ( myConfig.axisymmetric ) {
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

    override void assignVtxLocationsForDerivCalc()
    {
        // This code should follow very closesly the equivalent
        // code in sblock.d.
        size_t i, j, k;
        if (myConfig.dimensions == 2) {
            // First, do all of the internal secondary cells
            for ( i = imin+1; i <= imax; ++i ) {
                for ( j = jmin+1; j <= jmax; ++j ) {
                    // Secondary-cell centre IS a primary-cell vertex
                    // We are going to take a reference to this vertex
                    // and store some information there in it's role
                    // as a secondary-cell centre
                    SolidFVVertex vtx = getVtx(i, j);
                    // These are the corners of the secondary cell
                    SolidFVCell A = getCell(i, j-1);
                    SolidFVCell B = getCell(i, j);
                    SolidFVCell C = getCell(i-1, j);
                    SolidFVCell D = getCell(i-1, j-1);
                    // Retain locations and pointers to cell temperature for later
                    vtx.cloud_pos = [A.pos, B.pos, C.pos, D.pos];
                    vtx.cloud_T = [&(A.T), &(B.T), &(C.T), &(D.T)];
                } // j loop
            } // i loop
            // Half-cells along the edges of the block.
            // East boundary
            i = imax+1;
            for (j = jmin+1; j <= jmax; ++j) {
                SolidFVVertex vtx = getVtx(i, j);
                SolidFVInterface A = getIfi(i, j-1);
                SolidFVInterface B = getIfi(i, j);
                SolidFVCell C = getCell(i-1, j);
                SolidFVCell D = getCell(i-1, j-1);
                vtx.cloud_pos = [A.pos, B.pos, C.pos, D.pos];
                vtx.cloud_T = [&(A.T), &(B.T), &(C.T), &(D.T)];
            } // j loop
            // West boundary
            i = imin;
            for (j = jmin+1; j <= jmax; ++j ) {
                SolidFVVertex vtx = getVtx(i, j);
                SolidFVCell A = getCell(i, j-1);
                SolidFVCell B = getCell(i, j);
                SolidFVInterface C = getIfi(i, j);
                SolidFVInterface D = getIfi(i, j-1);
                vtx.cloud_pos = [A.pos, B.pos, C.pos, D.pos];
                vtx.cloud_T = [&(A.T), &(B.T), &(C.T), &(D.T)];
            } // j loop
            // North boundary
            j = jmax + 1;
            for (i = imin+1; i <= imax; ++i) {
                SolidFVVertex vtx = getVtx(i, j);
                SolidFVCell A = getCell(i, j-1);
                SolidFVInterface B = getIfj(i, j);
                SolidFVInterface C = getIfj(i-1, j);
                SolidFVCell D = getCell(i-1, j-1);
                vtx.cloud_pos = [A.pos, B.pos, C.pos, D.pos];
                vtx.cloud_T = [&(A.T), &(B.T), &(C.T), &(D.T)];
            } // i loop
            // South boundary
            j = jmin;
            for (i = imin+1; i <= imax; ++i) {
                SolidFVVertex vtx = getVtx(i, j);
                SolidFVInterface A = getIfj(i, j);
                SolidFVCell B = getCell(i, j);
                SolidFVCell C = getCell(i-1, j);
                SolidFVInterface D = getIfj(i-1, j);
                vtx.cloud_pos = [A.pos, B.pos, C.pos, D.pos];
                vtx.cloud_T = [&(A.T), &(B.T), &(C.T), &(D.T)];
            } // i loop
            // For the corners, we are going to use the same divergence-theorem-based
            // gradient calculator and let one edge collapse to a point, thus giving
            // it a triangle to compute over.  This should be fine. 
            // North-east corner
            {
                i = imax+1; j = jmax+1;
                SolidFVVertex vtx = getVtx(i, j);
                SolidFVInterface A = getIfi(i, j-1);
                SolidFVInterface B = getIfj(i-1, j);
                SolidFVCell C = getCell(i-1, j-1);
                vtx.cloud_pos = [A.pos, B.pos, C.pos];
                vtx.cloud_T = [&(A.T), &(B.T), &(C.T)];
            }
            // South-east corner
            {
                i = imax+1; j = jmin;
                SolidFVVertex vtx = getVtx(i, j);
                SolidFVInterface A = getIfi(i, j);
                SolidFVCell B = getCell(i-1, j);
                SolidFVInterface C = getIfj(i-1, j);
                vtx.cloud_pos = [A.pos, B.pos, C.pos];
                vtx.cloud_T = [&(A.T), &(B.T), &(C.T)];
            }
            // South-west corner
            {
                i = imin; j = jmin;
                SolidFVVertex vtx = getVtx(i, j);
                SolidFVInterface A = getIfj(i, j);
                SolidFVCell B = getCell(i, j);
                SolidFVInterface C = getIfi(i, j);
                vtx.cloud_pos = [A.pos, B.pos, C.pos];
                vtx.cloud_T = [&(A.T), &(B.T), &(C.T)];
            }
            // North-west corner
            {
                i = imin; j = jmax+1;
                SolidFVVertex vtx = getVtx(i, j);
                SolidFVCell A = getCell(i, j-1);
                SolidFVInterface B = getIfj(i, j);
                SolidFVInterface C = getIfi(i, j-1);
                vtx.cloud_pos = [A.pos, B.pos, C.pos];
                vtx.cloud_T = [&(A.T), &(B.T), &(C.T)];
            }
        } // end if 2D
    }

    override void applyPreSpatialDerivAction(double t, int tLevel)
    {
        bc[Face.north].applyPreSpatialDerivAction(t, tLevel);
        bc[Face.east].applyPreSpatialDerivAction(t, tLevel);
        bc[Face.south].applyPreSpatialDerivAction(t, tLevel);
        bc[Face.west].applyPreSpatialDerivAction(t, tLevel);
        if ( myConfig.dimensions == 3 ) {
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
        if ( myConfig.dimensions == 3 ) {
            bc[Face.top].applyPostFluxAction(t, tLevel);
            bc[Face.bottom].applyPostFluxAction(t, tLevel);
        }
    }

    override void computeSpatialDerivatives(int ftl)
    {
        if ( myConfig.dimensions == 3 ) {
            throw new Error("computeSpatialDerivatives() not implemented for 3D yet.");
        }
        
        size_t k = 0;
        for ( size_t i = imin; i <= imax+1; ++i ) {
            for ( size_t j = jmin; j <= jmax+1; ++j ) {
                SolidFVVertex vtx = getVtx(i, j);
                gradients_T_div(vtx);
            } // j loop
        } // i loop
    }

    override void computeFluxes()
    {
        size_t i, j, k;
        SolidFVInterface IFace;
        SolidFVVertex vtx1, vtx2;
        number dTdx, dTdy;
        number qx, qy;
        // East-facing interfaces
        for ( j = jmin; j <= jmax; ++j ) {
            for ( i = imin; i <= imax + 1; ++i ) {
                if ( i == imin && bc[Face.west].setsFluxDirectly )
                    continue;
                if ( i == imax+1 && bc[Face.east].setsFluxDirectly )
                    continue;
                IFace = getIfi(i, j);
                vtx1 = getVtx(i, j+1);
                vtx2 = getVtx(i, j);
                dTdx = 0.5*(vtx1.dTdx + vtx2.dTdx);
                dTdy = 0.5*(vtx1.dTdy + vtx2.dTdy);
                if (myConfig.solid_has_isotropic_properties) {
                    qx = -IFace.sp.k * dTdx;
                    qy = -IFace.sp.k * dTdy;
                }
                IFace.flux = qx * IFace.n.x + qy * IFace.n.y;
            }
        }
        // North-facing interfaces
        for ( j = jmin; j <= jmax + 1; ++j ) {
            for ( i = imin; i <= imax; ++i ) {
                if ( j == jmin && bc[Face.south].setsFluxDirectly ) 
                    continue;
                if ( j == jmax+1 && bc[Face.north].setsFluxDirectly )
                    continue;
                IFace = getIfj(i, j);
                vtx1 = getVtx(i, j);
                vtx2 = getVtx(i+1, j);
                dTdx = 0.5*(vtx1.dTdx + vtx2.dTdx);
                dTdy = 0.5*(vtx1.dTdy + vtx2.dTdy);
                if (myConfig.solid_has_isotropic_properties) {
                    qx = -IFace.sp.k * dTdx;
                    qy = -IFace.sp.k * dTdy;
                }
                IFace.flux = qx * IFace.n.x + qy * IFace.n.y;
            }
        }
    }

    override void clearSources()
    {
        foreach ( cell; activeCells ) cell.Q = 0.0;
    }
}

@nogc
void gradients_T_div(SolidFVVertex vtx)
{
    // Number of corners in our polygon.
    size_t n = vtx.cloud_pos.length;
    // Compute our own estimate of *twice* the area in xy plane here.
    // We can work with *twice* the area since it will cancel with
    // factor of 1/2 that appears in the contour integral.
    // Start with the contribution from the final segment of the bounding contour.
    number areaxy = (vtx.cloud_pos[0].x + vtx.cloud_pos[n-1].x) *
        (vtx.cloud_pos[0].y - vtx.cloud_pos[n-1].y);
    // Accumulate the contributions from the other segments.
    foreach (i; 0 .. n-1) {
        areaxy += (vtx.cloud_pos[i+1].x + vtx.cloud_pos[i].x) *
            (vtx.cloud_pos[i+1].y - vtx.cloud_pos[i].y);
    }
    number areaInv = 1.0 / areaxy;

    // Apply the divergence theorem to flow properties.
    //
    // Start with the contribution from the final segment of the bounding contour.
    number gradient_x = (*(vtx.cloud_T[0]) + *(vtx.cloud_T[n-1])) *
        (vtx.cloud_pos[0].y - vtx.cloud_pos[n-1].y);
    number gradient_y = (*(vtx.cloud_T[0]) + *(vtx.cloud_T[n-1])) *
        (vtx.cloud_pos[0].x - vtx.cloud_pos[n-1].x);
    // Accumulate the contributions from the other segments.
    foreach (i; 0 .. n-1) {
        gradient_x += (*(vtx.cloud_T[i+1]) + *(vtx.cloud_T[i])) *
            (vtx.cloud_pos[i+1].y - vtx.cloud_pos[i].y);
        gradient_y += (*(vtx.cloud_T[i+1]) + *(vtx.cloud_T[i])) *
            (vtx.cloud_pos[i+1].x - vtx.cloud_pos[i].x);
    }

    vtx.dTdx = gradient_x * areaInv;
    vtx.dTdy = -gradient_y * areaInv;
}
