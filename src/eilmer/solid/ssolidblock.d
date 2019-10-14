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
import nm.rsla;
import fvcore;
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
    size_t n_ghost_cell_layers;
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
        this.n_ghost_cell_layers = GlobalConfig.n_ghost_cell_layers;
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
        lua_pushinteger(myL, n_ghost_cell_layers); lua_setglobal(myL, "n_ghost_cell_layers");
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
        } else { // 3D
            calcVolumes3D();
            calcFaces3D();
            return;
            //throw new Error("SSolidBlock.computePrimaryCellGeometryData() not implemented for 3D.");
        }
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

    void calcVolumes3D()
    {
        size_t i, j, k;
        SolidFVCell cell;
        number volume;
        Vector3 centroid, p0, p1, p2, p3, p4, p5, p6, p7;
        //number iLength, jLength, kLength;
        for ( i = imin; i <= imax; ++i ) {
            for ( j = jmin; j <= jmax; ++j ) {
                for ( k = kmin; k <= kmax; ++k ) {
                    cell = getCell(i,j,k);
                    p0.set(cell.vtx[0].pos.x, cell.vtx[0].pos.y, cell.vtx[0].pos.z);
                    p1.set(cell.vtx[1].pos.x, cell.vtx[1].pos.y, cell.vtx[1].pos.z);
                    p2.set(cell.vtx[2].pos.x, cell.vtx[2].pos.y, cell.vtx[2].pos.z);
                    p3.set(cell.vtx[3].pos.x, cell.vtx[3].pos.y, cell.vtx[3].pos.z);
                    p4.set(cell.vtx[4].pos.x, cell.vtx[4].pos.y, cell.vtx[4].pos.z);
                    p5.set(cell.vtx[5].pos.x, cell.vtx[5].pos.y, cell.vtx[5].pos.z);
                    p6.set(cell.vtx[6].pos.x, cell.vtx[6].pos.y, cell.vtx[6].pos.z);
                    p7.set(cell.vtx[7].pos.x, cell.vtx[7].pos.y, cell.vtx[7].pos.z);
                    // Estimate the centroid so that we can use it as the peak
                    // of each of the pyramid sub-volumes.
                    centroid.set(0.125*(p0.x+p1.x+p2.x+p3.x+p4.x+p5.x+p6.x+p7.x),
                                 0.125*(p0.y+p1.y+p2.y+p3.y+p4.y+p5.y+p6.y+p7.y),
                                 0.125*(p0.z+p1.z+p2.z+p3.z+p4.z+p5.z+p6.z+p7.z));
                    // Mid-points of faces.
                    Vector3 pmN;
                    pmN.set(0.25*(p3.x+p2.x+p6.x+p7.x),
                            0.25*(p3.y+p2.y+p6.y+p7.y),
                            0.25*(p3.z+p2.z+p6.z+p7.z));
                    Vector3 pmE;
                    pmE.set(0.25*(p1.x+p2.x+p6.x+p5.x),
                            0.25*(p1.y+p2.y+p6.y+p5.y),
                            0.25*(p1.z+p2.z+p6.z+p5.z));
                    Vector3 pmS;
                    pmS.set(0.25*(p0.x+p1.x+p5.x+p4.x),
                            0.25*(p0.y+p1.y+p5.y+p4.y),
                            0.25*(p0.z+p1.z+p5.z+p4.z));
                    Vector3 pmW;
                    pmW.set(0.25*(p0.x+p3.x+p7.x+p4.x),
                            0.25*(p0.y+p3.y+p7.y+p4.y),
                            0.25*(p0.z+p3.z+p7.z+p4.z));
                    Vector3 pmT;
                    pmT.set(0.25*(p4.x+p5.x+p6.x+p7.x),
                            0.25*(p4.y+p5.y+p6.y+p7.y),
                            0.25*(p4.z+p5.z+p6.z+p7.z));
                    Vector3 pmB;
                    pmB.set(0.25*(p0.x+p1.x+p2.x+p3.x),
                            0.25*(p0.y+p1.y+p2.y+p3.y),
                            0.25*(p0.z+p1.z+p2.z+p3.z));
                    // Lengths between mid-points of faces.
                    // Note that we are assuming that the hexahedron is not very skewed
                    // when we later use these values as the widths of the hex cell.
                    //number dx, dy, dz;
                    //dx = pmE.x - pmW.x; dy = pmE.y - pmW.y; dz = pmE.z - pmW.z;
                    //iLen = sqrt(dx^^2 + dy^^2 + dz^^2);
                    //dx = pmN.x - pmS.x; dy = pmN.y - pmS.y; dz = pmN.z - pmS.z;
                    //jLen = sqrt(dx^^2 + dy^^2 + dz^^2);
                    //dx = pmT.x - pmB.x; dy = pmT.y - pmB.y; dz = pmT.z - pmB.z;
                    //kLen = sqrt(dx^^2 + dy^^2 + dz^^2);
                    // writeln("Single hexahedron divided into six tetragonal dipyramids.");
                    // J. Grandy (1997) Efficient Computation of Volume of Hexahedral Cells UCRL-ID-128886.
                    // Base of each dipyramid is specified clockwise from the outside.
                    number sub_volume; Vector3 sub_centroid;
                    volume = 0.0; Vector3 moment = Vector3(0.0, 0.0, 0.0);
                    pyramid_properties(p6, p7, p3, p2, centroid, sub_centroid, sub_volume);
                    volume += sub_volume; sub_centroid *= sub_volume; moment.add(sub_centroid);
                    pyramid_properties(p5, p6, p2, p1, centroid, sub_centroid, sub_volume);
                    volume += sub_volume; sub_centroid *= sub_volume; moment.add(sub_centroid);
                    pyramid_properties(p4, p5, p1, p0, centroid, sub_centroid, sub_volume);
                    volume += sub_volume; sub_centroid *= sub_volume; moment.add(sub_centroid);
                    pyramid_properties(p7, p4, p0, p3, centroid, sub_centroid, sub_volume);
                    volume += sub_volume; sub_centroid *= sub_volume; moment.add(sub_centroid);
                    pyramid_properties(p7, p6, p5, p4, centroid, sub_centroid, sub_volume);
                    volume += sub_volume; sub_centroid *= sub_volume; moment.add(sub_centroid);
                    pyramid_properties(p0, p1, p2, p3, centroid, sub_centroid, sub_volume);
                    volume += sub_volume; sub_centroid *= sub_volume; moment.add(sub_centroid);
                    //
                    if (volume < 0.0) {
                        throw new Error(text("Negative cell volume: solid_block ", id,
                                             " vol[", i, " ,", j, "]= ", volume));
                    }
                    //if (vol > max_vol) max_vol = vol;
                    //if (vol < min_vol) min_vol = vol;
                    cell.volume = volume;
                    //cell.areaxy = xyarea;
                } // j loop
            } // i loop
        } // k loop
    } // end calc_volumes_3D()
    
    void calcFaces3D()
    {
        SolidFVInterface iface;
        size_t i, j, k;
        Vector3 centroid, n, t1, t2, p0, p1, p2, p3;
        number area;
        // ifi interfaces are west interfaces, with their unit normal pointing east.
        for (k = kmin; k <= kmax; ++k) {
            for (j = jmin; j <= jmax; ++j) {
                for (i = imin; i <= imax+1; ++i) {
                    iface = getIfi(i,j,k);
                    // These are the corners.
                    p0.set(getVtx(i,j,k).pos.x, getVtx(i,j,k).pos.y, getVtx(i,j,k).pos.z);
                    p1.set(getVtx(i,j+1,k).pos.x, getVtx(i,j+1,k).pos.y, getVtx(i,j+1,k).pos.z);
                    p2.set(getVtx(i,j+1,k+1).pos.x, getVtx(i,j+1,k+1).pos.y, getVtx(i,j+1,k+1).pos.z);
                    p3.set(getVtx(i,j,k+1).pos.x, getVtx(i,j,k+1).pos.y, getVtx(i,j,k+1).pos.z);
                    quad_properties(p0, p1, p2, p3, centroid, n, t1, t2, area);

                    // set properties.
                    iface.n.refx = n.x;
                    iface.n.refy = n.y;
                    iface.n.refz = n.z;
                    iface.t2.refx = t2.x;
                    iface.t2.refy = t2.y;
                    iface.t2.refz = t2.z;
                    iface.t1.refx = t1.x;
                    iface.t1.refy = t1.y;
                    iface.t1.refz = t1.z;
                    iface.area = area;
                    iface.pos = (p0 + p1 + p2 + p3)/4.0;
                } // j loop
            } // j loop
        } // k loop
        // ifj interfaces are south interfaces, with their unit normal pointing north.
        for (k = kmin; k <= kmax; ++k) {
            for (i = imin; i <= imax; ++i) {
                for (j = jmin; j <= jmax+1; ++j) {
                    iface = getIfj(i,j,k);
                    // These are the corners.
                    p0.set(getVtx(i,j,k).pos.x, getVtx(i,j,k).pos.y, getVtx(i,j,k).pos.z);
                    p1.set(getVtx(i,j,k+1).pos.x, getVtx(i,j,k+1).pos.y, getVtx(i,j,k+1).pos.z);
                    p2.set(getVtx(i+1,j,k+1).pos.x, getVtx(i+1,j,k+1).pos.y, getVtx(i+1,j,k+1).pos.z);
                    p3.set(getVtx(i+1,j,k).pos.x, getVtx(i+1,j,k).pos.y, getVtx(i+1,j,k).pos.z);
                    quad_properties(p0, p1, p2, p3, centroid, n, t1, t2, area);

                    // set properties.
                    iface.n.refx = n.x;
                    iface.n.refy = n.y;
                    iface.n.refz = n.z;
                    iface.t2.refx = t2.x;
                    iface.t2.refy = t2.y;
                    iface.t2.refz = t2.z;
                    iface.t1.refx = t1.x;
                    iface.t1.refy = t1.y;
                    iface.t1.refz = t1.z;
                    iface.area = area;
                    iface.pos = (p0 + p1 + p2 + p3)/4.0;
                } // j loop
            } // i loop
        } // k loop
        // ifk interfaces are bottom interfaces, with unit normal pointing to top.
        for (i = imin; i <= imax; ++i) {
            for (j = jmin; j <= jmax; ++j) {
                for (k = kmin; k <= kmax+1; ++k) {
                    iface = getIfk(i,j,k);
                    // These are the corners.
                    p0.set(getVtx(i,j,k).pos.x, getVtx(i,j,k).pos.y, getVtx(i,j,k).pos.z);
                    p1.set(getVtx(i+1,j,k).pos.x, getVtx(i+1,j,k).pos.y, getVtx(i+1,j,k).pos.z);
                    p2.set(getVtx(i+1,j+1,k).pos.x, getVtx(i+1,j+1,k).pos.y, getVtx(i+1,j+1,k).pos.z);
                    p3.set(getVtx(i,j+1,k).pos.x, getVtx(i,j+1,k).pos.y, getVtx(i,j+1,k).pos.z);
                    quad_properties(p0, p1, p2, p3, centroid, n, t1, t2, area);

                    // set properties.
                    iface.n.refx = n.x;
                    iface.n.refy = n.y;
                    iface.n.refz = n.z;
                    iface.t2.refx = t2.x;
                    iface.t2.refy = t2.y;
                    iface.t2.refz = t2.z;
                    iface.t1.refx = t1.x;
                    iface.t1.refy = t1.y;
                    iface.t1.refz = t1.z;
                    iface.area = area;
                    iface.pos = (p0 + p1 + p2 + p3)/4.0;
                } // j loop
            } // i loop
        } // k loop
    } // end calc_faces_2D()

    void setupSpatialDerivativeCalc()
    {
        //if ( myConfig.dimensions == 3 ) {
        //    throw new Error("computeSpatialDerivatives() not implemented for 3D yet.");
        //}
        
        if ( myConfig.dimensions == 2 ) {
            size_t k = 0;
            for ( size_t i = imin; i <= imax+1; ++i ) {
                for ( size_t j = jmin; j <= jmax+1; ++j ) {
                    SolidFVVertex vtx = getVtx(i, j);
                    if (myConfig.spatial_deriv_calc == SpatialDerivCalc.divergence)
                        { }// nothing to do here
                    else // myConfig.spatial_deriv_calc == SpatialDerivCalc.least_squares
                        { gradients_T_lsq_setup(vtx, 2); }
                } // j loop
            } // i loop
        } else { // 3D
            for ( size_t i = imin; i <= imax+1; ++i ) {
                for ( size_t j = jmin; j <= jmax+1; ++j ) {
                    for ( size_t k = kmin; k <= kmax+1; ++k ) {
                        SolidFVVertex vtx = getVtx(i, j, k);
                        gradients_T_lsq_setup(vtx, 3);
                    } // k loop
                } // j loop
            } // i loop
        }
    }

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
        } else { // Flow quantity derivatives for 3D.
            // Internal secondary cell geometry information
            for ( i = imin; i <= imax-1; ++i ) {
                for ( j = jmin; j <= jmax-1; ++j ) {
                    for ( k = kmin; k <= kmax-1; ++k ) {
                        SolidFVVertex vtx = getVtx(i+1,j+1,k+1);
                        SolidFVCell c0 = getCell(i,j,k);
                        SolidFVCell c1 = getCell(i+1,j,k);
                        SolidFVCell c2 = getCell(i+1,j+1,k);
                        SolidFVCell c3 = getCell(i,j+1,k);
                        SolidFVCell c4 = getCell(i,j,k+1);
                        SolidFVCell c5 = getCell(i+1,j,k+1);
                        SolidFVCell c6 = getCell(i+1,j+1,k+1);
                        SolidFVCell c7 = getCell(i,j+1,k+1);
                        vtx.cloud_pos = [c0.pos, c1.pos, c2.pos, c3.pos,
                                         c4.pos, c5.pos, c6.pos, c7.pos];
                        vtx.cloud_T = [&(c0.T), &(c1.T), &(c2.T), &(c3.T), &(c4.T), &(c5.T), &(c6.T), &(c7.T)];
                    }
                }
            }
            // East boundary secondary cell geometry information
            i = imax;
            for ( j = jmin; j <= jmax-1; ++j ) {
                for ( k = kmin; k <= kmax-1; ++k ) {
                    SolidFVVertex vtx = getVtx(i+1,j+1,k+1);
                    SolidFVCell c0 = getCell(i,j,k);
                    SolidFVInterface c1 = getIfi(i+1,j,k);
                    SolidFVInterface c2 = getIfi(i+1,j+1,k);
                    SolidFVCell c3 = getCell(i,j+1,k);
                    SolidFVCell c4 = getCell(i,j,k+1);
                    SolidFVInterface c5 = getIfi(i+1,j,k+1);
                    SolidFVInterface c6 = getIfi(i+1,j+1,k+1);
                    SolidFVCell c7 = getCell(i,j+1,k+1);
                    vtx.cloud_pos = [c0.pos, c1.pos, c2.pos, c3.pos,
                                     c4.pos, c5.pos, c6.pos, c7.pos];
                    vtx.cloud_T = [&(c0.T), &(c1.T), &(c2.T), &(c3.T), &(c4.T), &(c5.T), &(c6.T), &(c7.T)];
                }
            }
            // West boundary secondary cell geometry information
            i = imin - 1;
            for ( j = jmin; j <= jmax-1; ++j ) {
                for ( k = kmin; k <= kmax-1; ++k ) {
                    SolidFVVertex vtx = getVtx(i+1,j+1,k+1);
                    SolidFVInterface c0 = getIfi(i+1,j,k);
                    SolidFVCell c1 = getCell(i+1,j,k);
                    SolidFVCell c2 = getCell(i+1,j+1,k);
                    SolidFVInterface c3 = getIfi(i+1,j+1,k);
                    SolidFVInterface c4 = getIfi(i+1,j,k+1);
                    SolidFVCell c5 = getCell(i+1,j,k+1);
                    SolidFVCell c6 = getCell(i+1,j+1,k+1);
                    SolidFVInterface c7 = getIfi(i+1,j+1,k+1);
                    vtx.cloud_pos = [c0.pos, c1.pos, c2.pos, c3.pos,
                                     c4.pos, c5.pos, c6.pos, c7.pos];
                    vtx.cloud_T = [&(c0.T), &(c1.T), &(c2.T), &(c3.T), &(c4.T), &(c5.T), &(c6.T), &(c7.T)];
                }
            }
            // North boundary secondary cell geometry information
            j = jmax;
            for ( i = imin; i <= imax-1; ++i ) {
                for ( k = kmin; k <= kmax-1; ++k ) {
                    SolidFVVertex vtx = getVtx(i+1,j+1,k+1);
                    SolidFVCell c0 = getCell(i,j,k);
                    SolidFVCell c1 = getCell(i+1,j,k);
                    SolidFVInterface c2 = getIfj(i+1,j+1,k);
                    SolidFVInterface c3 = getIfj(i,j+1,k);
                    SolidFVCell c4 = getCell(i,j,k+1);
                    SolidFVCell c5 = getCell(i+1,j,k+1);
                    SolidFVInterface c6 = getIfj(i+1,j+1,k+1);
                    SolidFVInterface c7 = getIfj(i,j+1,k+1);
                    vtx.cloud_pos = [c0.pos, c1.pos, c2.pos, c3.pos,
                                     c4.pos, c5.pos, c6.pos, c7.pos];
                    vtx.cloud_T = [&(c0.T), &(c1.T), &(c2.T), &(c3.T), &(c4.T), &(c5.T), &(c6.T), &(c7.T)];
                }
            }
            // South boundary secondary cell geometry information
            j = jmin - 1;
            for ( i = imin; i <= imax-1; ++i ) {
                for ( k = kmin; k <= kmax-1; ++k ) {
                    SolidFVVertex vtx = getVtx(i+1,j+1,k+1);
                    SolidFVInterface c0 = getIfj(i,j+1,k);
                    SolidFVInterface c1 = getIfj(i+1,j+1,k);
                    SolidFVCell c2 = getCell(i+1,j+1,k);
                    SolidFVCell c3 = getCell(i,j+1,k);
                    SolidFVInterface c4 = getIfj(i,j+1,k+1);
                    SolidFVInterface c5 = getIfj(i+1,j+1,k+1);
                    SolidFVCell c6 = getCell(i+1,j+1,k+1);
                    SolidFVCell c7 = getCell(i,j+1,k+1);
                    vtx.cloud_pos = [c0.pos, c1.pos, c2.pos, c3.pos,
                                     c4.pos, c5.pos, c6.pos, c7.pos];
                    vtx.cloud_T = [&(c0.T), &(c1.T), &(c2.T), &(c3.T), &(c4.T), &(c5.T), &(c6.T), &(c7.T)];
                }
            }
            // Top boundary secondary cell geometry information
            k = kmax;
            for ( i = imin; i <= imax-1; ++i ) {
                for ( j = jmin; j <= jmax-1; ++j ) {
                    SolidFVVertex vtx = getVtx(i+1,j+1,k+1);
                    SolidFVCell c0 = getCell(i,j,k);
                    SolidFVCell c1 = getCell(i+1,j,k);
                    SolidFVCell c2 = getCell(i+1,j+1,k);
                    SolidFVCell c3 = getCell(i,j+1,k);
                    SolidFVInterface c4 = getIfk(i,j,k+1);
                    SolidFVInterface c5 = getIfk(i+1,j,k+1);
                    SolidFVInterface c6 = getIfk(i+1,j+1,k+1);
                    SolidFVInterface c7 = getIfk(i,j+1,k+1);
                    vtx.cloud_pos = [c0.pos, c1.pos, c2.pos, c3.pos,
                                     c4.pos, c5.pos, c6.pos, c7.pos];
                    vtx.cloud_T = [&(c0.T), &(c1.T), &(c2.T), &(c3.T), &(c4.T), &(c5.T), &(c6.T), &(c7.T)];
                }
            }
            // Bottom boundary secondary cell geometry information
            k = kmin - 1;
            for ( i = imin; i <= imax-1; ++i ) {
                for ( j = jmin; j <= jmax-1; ++j ) {
                    SolidFVVertex vtx = getVtx(i+1,j+1,k+1);
                    SolidFVInterface c0 = getIfk(i,j,k+1);
                    SolidFVInterface c1 = getIfk(i+1,j,k+1);
                    SolidFVInterface c2 = getIfk(i+1,j+1,k+1);
                    SolidFVInterface c3 = getIfk(i,j+1,k+1);
                    SolidFVCell c4 = getCell(i,j,k+1);
                    SolidFVCell c5 = getCell(i+1,j,k+1);
                    SolidFVCell c6 = getCell(i+1,j+1,k+1);
                    SolidFVCell c7 = getCell(i,j+1,k+1);
                    vtx.cloud_pos = [c0.pos, c1.pos, c2.pos, c3.pos,
                                     c4.pos, c5.pos, c6.pos, c7.pos];
                    vtx.cloud_T = [&(c0.T), &(c1.T), &(c2.T), &(c3.T), &(c4.T), &(c5.T), &(c6.T), &(c7.T)];
                }
            }
            // Now, do the 4 edges around the bottom face.
            // Bottom-South edge [0]-->[1]
            j = jmin; k = kmin;         
            for ( i = imin+1; i <= imax; ++i ) {
                SolidFVVertex vtx = getVtx(i,j,k);
                SolidFVCell c0 = getCell(i-1,j,k);
                SolidFVCell c1 = getCell(i,j,k);
                SolidFVInterface c2 = getIfj(i-1,j,k);
                SolidFVInterface c3 = getIfk(i-1,j,k);
                SolidFVInterface c4 = getIfj(i,j,k);
                SolidFVInterface c5 = getIfk(i,j,k);
                vtx.cloud_pos = [c0.pos, c1.pos, c2.pos,
                                 c3.pos, c4.pos, c5.pos];
                vtx.cloud_T = [&(c0.T), &(c1.T), &(c2.T), &(c3.T), &(c4.T), &(c5.T)];
            }
            // Bottom-North edge [3]-->[2]
            j = jmax; k = kmin;
            for ( i = imin+1; i <= imax; ++i ) {
                SolidFVVertex vtx = getVtx(i,j+1,k);
                SolidFVCell c0 = getCell(i-1,j,k);
                SolidFVCell c1 = getCell(i,j,k);
                SolidFVInterface c2 = getIfj(i-1,j+1,k);
                SolidFVInterface c3 = getIfk(i-1,j,k);
                SolidFVInterface c4 = getIfj(i,j+1,k);
                SolidFVInterface c5 = getIfk(i,j,k);
                vtx.cloud_pos = [c0.pos, c1.pos, c2.pos,
                                 c3.pos, c4.pos, c5.pos];
                vtx.cloud_T = [&(c0.T), &(c1.T), &(c2.T), &(c3.T), &(c4.T), &(c5.T)];
            }
            // Bottom-West edge [0]-->[3]
            i = imin; k = kmin;
            for ( j = jmin+1; j <= jmax; ++j ) {
                SolidFVVertex vtx = getVtx(i,j,k);
                SolidFVCell c0 = getCell(i,j-1,k);
                SolidFVCell c1 = getCell(i,j,k);
                SolidFVInterface c2 = getIfi(i,j-1,k);
                SolidFVInterface c3 = getIfk(i,j-1,k);
                SolidFVInterface c4 = getIfi(i,j,k);
                SolidFVInterface c5 = getIfk(i,j,k);
                vtx.cloud_pos = [c0.pos, c1.pos, c2.pos,
                                 c3.pos, c4.pos, c5.pos];
                vtx.cloud_T = [&(c0.T), &(c1.T), &(c2.T), &(c3.T), &(c4.T), &(c5.T)];
            }
            // Bottom-East edge [1]-->[2]
            i = imax; k = kmin;
            for ( j = jmin+1; j <= jmax; ++j ) {
                SolidFVVertex vtx = getVtx(i+1,j,k);
                SolidFVCell c0 = getCell(i,j-1,k);
                SolidFVCell c1 = getCell(i,j,k);
                SolidFVInterface c2 = getIfi(i+1,j-1,k);
                SolidFVInterface c3 = getIfk(i,j-1,k);
                SolidFVInterface c4 = getIfi(i+1,j,k);
                SolidFVInterface c5 = getIfk(i,j,k);
                vtx.cloud_pos = [c0.pos, c1.pos, c2.pos,
                                 c3.pos, c4.pos, c5.pos];
                vtx.cloud_T = [&(c0.T), &(c1.T), &(c2.T), &(c3.T), &(c4.T), &(c5.T)];
            }
            // 4 edges around the top face.
            // Top-South edge [4]-->[5]
            j = jmin; k = kmax;
            for ( i = imin+1; i <= imax; ++i ) {
                SolidFVVertex vtx = getVtx(i,j,k+1);
                SolidFVCell c0 = getCell(i-1,j,k);
                SolidFVCell c1 = getCell(i,j,k);
                SolidFVInterface c2 = getIfj(i-1,j,k);
                SolidFVInterface c3 = getIfk(i-1,j,k+1);
                SolidFVInterface c4 = getIfj(i,j,k);
                SolidFVInterface c5 = getIfk(i,j,k+1);
                vtx.cloud_pos = [c0.pos, c1.pos, c2.pos,
                                 c3.pos, c4.pos, c5.pos];
                vtx.cloud_T = [&(c0.T), &(c1.T), &(c2.T), &(c3.T), &(c4.T), &(c5.T)];
            }
            // Top-North edge [7]-->[6]
            j = jmax; k = kmax;
            for ( i = imin+1; i <= imax; ++i ) {
                SolidFVVertex vtx = getVtx(i,j+1,k+1);
                SolidFVCell c0 = getCell(i-1,j,k);
                SolidFVCell c1 = getCell(i,j,k);
                SolidFVInterface c2 = getIfj(i-1,j+1,k);
                SolidFVInterface c3 = getIfk(i-1,j,k+1);
                SolidFVInterface c4 = getIfj(i,j+1,k);
                SolidFVInterface c5 = getIfk(i,j,k+1);
                vtx.cloud_pos = [c0.pos, c1.pos, c2.pos,
                                 c3.pos, c4.pos, c5.pos];
                vtx.cloud_T = [&(c0.T), &(c1.T), &(c2.T), &(c3.T), &(c4.T), &(c5.T)];
            }
            // Top-West edge [4]-->[7]
            i = imin; k = kmax;
            for ( j = jmin+1; j <= jmax; ++j ) {
                SolidFVVertex vtx = getVtx(i,j,k+1);
                SolidFVCell c0 = getCell(i,j-1,k);
                SolidFVCell c1 = getCell(i,j,k);
                SolidFVInterface c2 = getIfi(i,j-1,k);
                SolidFVInterface c3 = getIfk(i,j-1,k+1);
                SolidFVInterface c4 = getIfi(i,j,k);
                SolidFVInterface c5 = getIfk(i,j,k+1);
                vtx.cloud_pos = [c0.pos, c1.pos, c2.pos,
                                 c3.pos, c4.pos, c5.pos];
                vtx.cloud_T = [&(c0.T), &(c1.T), &(c2.T), &(c3.T), &(c4.T), &(c5.T)];
            }
            // Top-East edge [5]-->[6]
            i = imax; k = kmax;
            for ( j = jmin+1; j <= jmax; ++j ) {
                SolidFVVertex vtx = getVtx(i+1,j,k+1);
                SolidFVCell c0 = getCell(i,j-1,k);
                SolidFVCell c1 = getCell(i,j,k);
                SolidFVInterface c2 = getIfi(i+1,j-1,k);
                SolidFVInterface c3 = getIfk(i,j-1,k+1);
                SolidFVInterface c4 = getIfi(i+1,j,k);
                SolidFVInterface c5 = getIfk(i,j,k+1);
                vtx.cloud_pos = [c0.pos, c1.pos, c2.pos,
                                 c3.pos, c4.pos, c5.pos];
                vtx.cloud_T = [&(c0.T), &(c1.T), &(c2.T), &(c3.T), &(c4.T), &(c5.T)];
            }
            // 4 edges running from bottom to top.
            // South-West edge [0]-->[4]
            i = imin; j = jmin;
            for ( k = kmin+1; k <= kmax; ++k ) {
                SolidFVVertex vtx = getVtx(i,j,k);
                SolidFVCell c0 = getCell(i,j,k-1);
                SolidFVCell c1 = getCell(i,j,k);
                SolidFVInterface c2 = getIfi(i,j,k-1);
                SolidFVInterface c3 = getIfj(i,j,k-1);
                SolidFVInterface c4 = getIfi(i,j,k);
                SolidFVInterface c5 = getIfj(i,j,k);
                vtx.cloud_pos = [c0.pos, c1.pos, c2.pos,
                                 c3.pos, c4.pos, c5.pos];
                vtx.cloud_T = [&(c0.T), &(c1.T), &(c2.T), &(c3.T), &(c4.T), &(c5.T)];
            }
            // South-East edge [1]-->[5]
            i = imax; j = jmin;
            for ( k = kmin+1; k <= kmax; ++k ) {
                SolidFVVertex vtx = getVtx(i+1,j,k);
                SolidFVCell c0 = getCell(i,j,k-1);
                SolidFVCell c1 = getCell(i,j,k);
                SolidFVInterface c2 = getIfi(i+1,j,k-1);
                SolidFVInterface c3 = getIfj(i,j,k-1);
                SolidFVInterface c4 = getIfi(i+1,j,k);
                SolidFVInterface c5 = getIfj(i,j,k);
                vtx.cloud_pos = [c0.pos, c1.pos, c2.pos,
                                 c3.pos, c4.pos, c5.pos];
                vtx.cloud_T = [&(c0.T), &(c1.T), &(c2.T), &(c3.T), &(c4.T), &(c5.T)];
            }
            // North-East edge [2]-->[6]
            i = imax; j = jmax;
            for ( k = kmin+1; k <= kmax; ++k ) {
                SolidFVVertex vtx = getVtx(i+1,j+1,k);
                SolidFVCell c0 = getCell(i,j,k-1);
                SolidFVCell c1 = getCell(i,j,k);
                SolidFVInterface c2 = getIfi(i+1,j,k-1);
                SolidFVInterface c3 = getIfj(i,j+1,k-1);
                SolidFVInterface c4 = getIfi(i+1,j,k);
                SolidFVInterface c5 = getIfj(i,j+1,k);
                vtx.cloud_pos = [c0.pos, c1.pos, c2.pos,
                                 c3.pos, c4.pos, c5.pos];
                vtx.cloud_T = [&(c0.T), &(c1.T), &(c2.T), &(c3.T), &(c4.T), &(c5.T)];
            }
            // North-West edge [3]-->[7]
            i = imin; j = jmax;
            for ( k = kmin+1; k <= kmax; ++k ) {
                SolidFVVertex vtx = getVtx(i,j+1,k);
                SolidFVCell c0 = getCell(i,j,k-1);
                SolidFVCell c1 = getCell(i,j,k);
                SolidFVInterface c2 = getIfi(i,j,k-1);
                SolidFVInterface c3 = getIfj(i,j+1,k-1);
                SolidFVInterface c4 = getIfi(i,j,k);
                SolidFVInterface c5 = getIfj(i,j+1,k);
                vtx.cloud_pos = [c0.pos, c1.pos, c2.pos,
                                 c3.pos, c4.pos, c5.pos];
                vtx.cloud_T = [&(c0.T), &(c1.T), &(c2.T), &(c3.T), &(c4.T), &(c5.T)];
            }
            // Finally, the 8 corners.
            // South-West-Bottom corner [0]
            i = imin; j = jmin; k = kmin;
            {
                SolidFVVertex vtx = getVtx(i,j,k);
                SolidFVCell c0 = getCell(i,j,k);
                SolidFVInterface c1 = getIfi(i,j,k);
                SolidFVInterface c2 = getIfj(i,j,k);
                SolidFVInterface c3 = getIfk(i,j,k);
                vtx.cloud_pos = [c0.pos, c1.pos, c2.pos, c3.pos];
                vtx.cloud_T = [&(c0.T), &(c1.T), &(c2.T), &(c3.T)];
            }
            // South-East-Bottom corner [1]
            i = imax; j = jmin; k = kmin;
            {
                SolidFVVertex vtx = getVtx(i+1,j,k);
                SolidFVCell c0 = getCell(i,j,k);
                SolidFVInterface c1 = getIfi(i+1,j,k);
                SolidFVInterface c2 = getIfj(i,j,k);
                SolidFVInterface c3 = getIfk(i,j,k);
                vtx.cloud_pos = [c0.pos, c1.pos, c2.pos, c3.pos];
                vtx.cloud_T = [&(c0.T), &(c1.T), &(c2.T), &(c3.T)];
            }
            // North-East-Bottom corner [2]
            i = imax; j = jmax; k = kmin;
            {
                SolidFVVertex vtx = getVtx(i+1,j+1,k);
                SolidFVCell c0 = getCell(i,j,k);
                SolidFVInterface c1 = getIfi(i+1,j,k);
                SolidFVInterface c2 = getIfj(i,j+1,k);
                SolidFVInterface c3 = getIfk(i,j,k);
                vtx.cloud_pos = [c0.pos, c1.pos, c2.pos, c3.pos];
                vtx.cloud_T = [&(c0.T), &(c1.T), &(c2.T), &(c3.T)];
            }
            // North-West-Bottom corner [3]
            i = imin; j = jmax; k = kmin;
            {
                SolidFVVertex vtx = getVtx(i,j+1,k);
                SolidFVCell c0 = getCell(i,j,k);
                SolidFVInterface c1 = getIfi(i,j,k);
                SolidFVInterface c2 = getIfj(i,j+1,k);
                SolidFVInterface c3 = getIfk(i,j,k);
                vtx.cloud_pos = [c0.pos, c1.pos, c2.pos, c3.pos];
                vtx.cloud_T = [&(c0.T), &(c1.T), &(c2.T), &(c3.T)];
            }
            // South-West-Top corner [4]
            i = imin; j = jmin; k = kmax;
            {
                SolidFVVertex vtx = getVtx(i,j,k+1);
                SolidFVCell c0 = getCell(i,j,k);
                SolidFVInterface c1 = getIfi(i,j,k);
                SolidFVInterface c2 = getIfj(i,j,k);
                SolidFVInterface c3 = getIfk(i,j,k+1);
                vtx.cloud_pos = [c0.pos, c1.pos, c2.pos, c3.pos];
                vtx.cloud_T = [&(c0.T), &(c1.T), &(c2.T), &(c3.T)];
            }
            // South-East-Top corner [5]
            i = imax; j = jmin; k = kmax;
            {
                SolidFVVertex vtx = getVtx(i+1,j,k+1);
                SolidFVCell c0 = getCell(i,j,k);
                SolidFVInterface c1 = getIfi(i+1,j,k);
                SolidFVInterface c2 = getIfj(i,j,k);
                SolidFVInterface c3 = getIfk(i,j,k+1);
                vtx.cloud_pos = [c0.pos, c1.pos, c2.pos, c3.pos];
                vtx.cloud_T = [&(c0.T), &(c1.T), &(c2.T), &(c3.T)];
            }
            // North-East-Top corner [6]
            i = imax; j = jmax; k = kmax;
            {
                SolidFVVertex vtx = getVtx(i+1,j+1,k+1);
                SolidFVCell c0 = getCell(i,j,k);
                SolidFVInterface c1 = getIfi(i+1,j,k);
                SolidFVInterface c2 = getIfj(i,j+1,k);
                SolidFVInterface c3 = getIfk(i,j,k+1);
                vtx.cloud_pos = [c0.pos, c1.pos, c2.pos, c3.pos];
                vtx.cloud_T = [&(c0.T), &(c1.T), &(c2.T), &(c3.T)];
            }
            // North-West-Top corner [7]
            i = imin; j = jmax; k = kmax;
            {
                SolidFVVertex vtx = getVtx(i,j+1,k+1);
                SolidFVCell c0 = getCell(i,j,k);
                SolidFVInterface c1 = getIfi(i,j,k);
                SolidFVInterface c2 = getIfj(i,j+1,k);
                SolidFVInterface c3 = getIfk(i,j,k+1);
                vtx.cloud_pos = [c0.pos, c1.pos, c2.pos, c3.pos];
                vtx.cloud_T = [&(c0.T), &(c1.T), &(c2.T), &(c3.T)];
            }
        } // end if (myConfig.dimensions
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
        //if ( myConfig.dimensions == 3 ) {
        //    throw new Error("computeSpatialDerivatives() not implemented for 3D yet.");
        //}

        if ( myConfig.dimensions == 2 ) {
            size_t k = 0;
            for ( size_t i = imin; i <= imax+1; ++i ) {
                for ( size_t j = jmin; j <= jmax+1; ++j ) {
                    SolidFVVertex vtx = getVtx(i, j);
                    if (myConfig.spatial_deriv_calc == SpatialDerivCalc.divergence)
                        { gradients_T_div(vtx); }
                    else // myConfig.spatial_deriv_calc == SpatialDerivCalc.least_squares
                        { gradients_T_lsq(vtx, 2); }
                } // j loop
            } // i loop
        } else { // 3D
            for ( size_t i = imin; i <= imax+1; ++i ) {
                for ( size_t j = jmin; j <= jmax+1; ++j ) {
                    for ( size_t k = kmin; k <= kmax+1; ++k ) {
                        SolidFVVertex vtx = getVtx(i, j, k);
                        gradients_T_lsq(vtx, 3);
                    } // k loop
                } // j loop
            } // i loop
        }
    }

    override void computeFluxes()
    {
        size_t i, j, k;
        SolidFVInterface IFace;
        SolidFVVertex vtx1, vtx2, vtx3, vtx4;
        number dTdx, dTdy, dTdz;
        number qx, qy, qz;
        if ( myConfig.dimensions == 2 ) {
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
                    else if (myConfig.solid_has_homogeneous_properties) { 
                        qx = -IFace.sp.k11 * dTdx - IFace.sp.k12 * dTdy; 
                        qy = -IFace.sp.k21 * dTdx - IFace.sp.k22 * dTdy;  
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
                    else if (myConfig.solid_has_homogeneous_properties) { 
                        qx = -IFace.sp.k11 * dTdx - IFace.sp.k12 * dTdy; 
                        qy = -IFace.sp.k21 * dTdx - IFace.sp.k22 * dTdy;  
                    }
                    IFace.flux = qx * IFace.n.x + qy * IFace.n.y;
                }
            }
        } else { // 3D
            // ifi interfaces are west interfaces, with their unit normal pointing east.
            for (k = kmin; k <= kmax; ++k) {
                for (j = jmin; j <= jmax; ++j) {
                    for (i = imin; i <= imax+1; ++i) {
                    if ( i == imin && bc[Face.west].setsFluxDirectly )
                        continue;
                    if ( i == imax+1 && bc[Face.east].setsFluxDirectly )
                        continue;
                    IFace = getIfi(i, j, k);
                    vtx1 = getVtx(i,j,k);
                    vtx2 = getVtx(i,j+1,k);
                    vtx3 = getVtx(i,j+1,k+1);
                    vtx4 = getVtx(i,j,k+1);
                    dTdx = 0.25*(vtx1.dTdx + vtx2.dTdx + vtx3.dTdx + vtx4.dTdx);
                    dTdy = 0.25*(vtx1.dTdy + vtx2.dTdy + vtx3.dTdy + vtx4.dTdy);
                    dTdz = 0.25*(vtx1.dTdz + vtx2.dTdz + vtx3.dTdz + vtx4.dTdz);
                    if (myConfig.solid_has_isotropic_properties) {
                        qx = -IFace.sp.k * dTdx;
                        qy = -IFace.sp.k * dTdy;
                        qz = -IFace.sp.k * dTdz;
                    }
                    else if (myConfig.solid_has_homogeneous_properties) { 
                        throw new Error("solid_has_homogeneous_properties not implemented for 3D yet.");
                        //qx = -IFace.sp.k11 * dTdx - IFace.sp.k12 * dTdy; 
                        //qy = -IFace.sp.k21 * dTdx - IFace.sp.k22 * dTdy;
                    }
                    IFace.flux = qx * IFace.n.x + qy * IFace.n.y + qz * IFace.n.z;
                    }  // i loop
                } // j loop
            } // k loop
            
            // ifj interfaces are south interfaces, with their unit normal pointing north.
            for (k = kmin; k <= kmax; ++k) {
                for (i = imin; i <= imax; ++i) {
                    for (j = jmin; j <= jmax+1; ++j) {
                        if ( j == jmin && bc[Face.south].setsFluxDirectly ) 
                            continue;
                        if ( j == jmax+1 && bc[Face.north].setsFluxDirectly )
                            continue;
                        IFace = getIfj(i, j, k);
                        vtx1 = getVtx(i,j,k);
                        vtx2 = getVtx(i,j,k+1);
                        vtx3 = getVtx(i+1,j,k+1);
                        vtx4 = getVtx(i+1,j,k);
                        dTdx = 0.25*(vtx1.dTdx + vtx2.dTdx + vtx3.dTdx + vtx4.dTdx);
                        dTdy = 0.25*(vtx1.dTdy + vtx2.dTdy + vtx3.dTdy + vtx4.dTdy);
                        dTdz = 0.25*(vtx1.dTdz + vtx2.dTdz + vtx3.dTdz + vtx4.dTdz);
                        if (myConfig.solid_has_isotropic_properties) {
                            qx = -IFace.sp.k * dTdx;
                            qy = -IFace.sp.k * dTdy;
                            qz = -IFace.sp.k * dTdz;
                        }
                        else if (myConfig.solid_has_homogeneous_properties) { 
                            throw new Error("solid_has_homogeneous_properties not implemented for 3D yet.");
                            //qx = -IFace.sp.k11 * dTdx - IFace.sp.k12 * dTdy; 
                            //qy = -IFace.sp.k21 * dTdx - IFace.sp.k22 * dTdy;
                        }
                        IFace.flux = qx * IFace.n.x + qy * IFace.n.y + qz * IFace.n.z;
                    } // j loop
                } // i loop
            } // k loop
            // ifk interfaces are bottom interfaces, with unit normal pointing to top.
            for (i = imin; i <= imax; ++i) {
                for (j = jmin; j <= jmax; ++j) {
                    for (k = kmin; k <= kmax+1; ++k) {
                        if ( k == kmin && bc[Face.bottom].setsFluxDirectly ) 
                            continue;
                        if ( k == kmax+1 && bc[Face.top].setsFluxDirectly )
                            continue;
                        IFace = getIfk(i, j, k);
                        vtx1 = getVtx(i,j,k);
                        vtx2 = getVtx(i+1,j,k);
                        vtx3 = getVtx(i+1,j+1,k);
                        vtx4 = getVtx(i,j+1,k);
                        dTdx = 0.25*(vtx1.dTdx + vtx2.dTdx + vtx3.dTdx + vtx4.dTdx);
                        dTdy = 0.25*(vtx1.dTdy + vtx2.dTdy + vtx3.dTdy + vtx4.dTdy);
                        dTdz = 0.25*(vtx1.dTdz + vtx2.dTdz + vtx3.dTdz + vtx4.dTdz);
                        if (myConfig.solid_has_isotropic_properties) {
                            qx = -IFace.sp.k * dTdx;
                            qy = -IFace.sp.k * dTdy;
                            qz = -IFace.sp.k * dTdz;
                        }
                        else if (myConfig.solid_has_homogeneous_properties) { 
                            throw new Error("solid_has_homogeneous_properties not implemented for 3D yet.");
                            //qx = -IFace.sp.k11 * dTdx - IFace.sp.k12 * dTdy; 
                            //qy = -IFace.sp.k21 * dTdx - IFace.sp.k22 * dTdy;
                        }
                        IFace.flux = qx * IFace.n.x + qy * IFace.n.y + qz * IFace.n.z;
                    } // j loop
                } // i loop
            } // k loop
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

@nogc
void gradients_T_lsq_setup(SolidFVVertex vtx, int dimensions)
{
    size_t n = vtx.cloud_pos.length;
    number[12] weights2;
    size_t loop_init;
    loop_init = 0; // All points count.
    //
    // Calculate weights used in the least-squares gradient calculation.
    // These are the square of the weights on the original linear constraint eqns.
    // These weights are calculated with the current interface/vertex
    // as the reference point (pos).
    // For the "faces" spatial location we are expecting the primary point
    // (i.e. the face at which we are calculating the gradients) to be in
    // the first cloud position. 
    number x0 = vtx.pos.x; number y0 = vtx.pos.y; number z0 = vtx.pos.z;
    if (dimensions == 2) {
        foreach (i; loop_init .. n) {
            number dx = vtx.cloud_pos[i].x - x0;
            number dy = vtx.cloud_pos[i].y - y0;
            weights2[i] = 1.0/(dx*dx+dy*dy);
        }
    } else { //3D
        foreach (i; loop_init .. n) {
            number dx = vtx.cloud_pos[i].x - x0;
            number dy = vtx.cloud_pos[i].y - y0;
            number dz = vtx.cloud_pos[i].z - z0;
            weights2[i] = 1.0/(dx*dx+dy*dy+dz*dz);
        }
    }
    x0 = 0.0; y0 = 0.0; z0 = 0.0;
    foreach (i; 0 .. n) {
        x0 += vtx.cloud_pos[i].x; y0 += vtx.cloud_pos[i].y; z0 += vtx.cloud_pos[i].z;
    }
    x0 /= n; y0 /= n; z0 /= n; // midpoint
    
    number[12] dx, dy, dz;
    //
    // Assemble and invert the normal matrix.
    // We'll reuse the resulting inverse for each flow-field quantity.
    number[12] wx, wy, wz; 
    if (dimensions == 3) {
        number[6][3] xTx; // normal matrix, augmented to give 6 entries per row
        number xx = 0.0; number xy = 0.0; number xz = 0.0;
        number yy = 0.0; number yz = 0.0; number zz = 0.0;
        foreach (i; loop_init .. n) {
            dx[i] = vtx.cloud_pos[i].x - x0;
            dy[i] = vtx.cloud_pos[i].y - y0;
            dz[i] = vtx.cloud_pos[i].z - z0;
            xx += weights2[i]*dx[i]*dx[i];
            xy += weights2[i]*dx[i]*dy[i];
            xz += weights2[i]*dx[i]*dz[i];
            yy += weights2[i]*dy[i]*dy[i];
            yz += weights2[i]*dy[i]*dz[i];
            zz += weights2[i]*dz[i]*dz[i];
        }
        xTx[0][0] = xx; xTx[0][1] = xy; xTx[0][2] = xz;
        xTx[1][0] = xy; xTx[1][1] = yy; xTx[1][2] = yz;
        xTx[2][0] = xz; xTx[2][1] = yz; xTx[2][2] = zz;
        xTx[0][3] = 1.0; xTx[0][4] = 0.0; xTx[0][5] = 0.0;
        xTx[1][3] = 0.0; xTx[1][4] = 1.0; xTx[1][5] = 0.0;
        xTx[2][3] = 0.0; xTx[2][4] = 0.0; xTx[2][5] = 1.0;
        double very_small_value = 1.0e-16*(normInf!(3,3,6,number)(xTx).re)^^3;
        if (0 != computeInverse!(3,3,6,number)(xTx, very_small_value)) {
            throw new FlowSolverException("Failed to invert LSQ normal matrix");
            // Assume that the rows are linearly dependent 
            // because the sample points are coplanar or colinear.
        }
        // Prepare final weights for later use in the reconstruction phase.
        foreach (i; loop_init .. n) {
            vtx.wx[i] = xTx[0][3]*dx[i] + xTx[0][4]*dy[i] + xTx[0][5]*dz[i];
            vtx.wx[i] *= weights2[i];
            vtx.wy[i] = xTx[1][3]*dx[i] + xTx[1][4]*dy[i] + xTx[1][5]*dz[i];
            vtx.wy[i] *= weights2[i];
            vtx.wz[i] = xTx[2][3]*dx[i] + xTx[2][4]*dy[i] + xTx[2][5]*dz[i];
            vtx.wz[i] *= weights2[i];
        }
    } else {
        // dimensions == 2
        number[4][2] xTx; // normal matrix, augmented to give 4 entries per row
        number xx = 0.0; number xy = 0.0; number yy = 0.0;
        foreach (i; loop_init .. n) {
            dx[i] = vtx.cloud_pos[i].x - x0;
            dy[i] = vtx.cloud_pos[i].y - y0;
            xx += weights2[i]*dx[i]*dx[i];
            xy += weights2[i]*dx[i]*dy[i];
            yy += weights2[i]*dy[i]*dy[i];
        }
        xTx[0][0] = xx; xTx[0][1] = xy;
        xTx[1][0] = xy; xTx[1][1] = yy;
        xTx[0][2] = 1.0; xTx[0][3] = 0.0;
        xTx[1][2] = 0.0; xTx[1][3] = 1.0;
        double very_small_value = 1.0e-16*(normInf!(2,2,4,number)(xTx).re)^^2;
        if (0 != computeInverse!(2,2,4,number)(xTx, very_small_value)) {
            throw new FlowSolverException("Failed to invert LSQ normal matrix");
            // Assume that the rows are linearly dependent 
            // because the sample points are colinear.
        }
        //number[12] wx, wy, wz; 
        // Prepare final weights for later use in the reconstruction phase.
        foreach (i; loop_init .. n) {
            vtx.wx[i] = xTx[0][2]*dx[i] + xTx[0][3]*dy[i];
            vtx.wx[i] *= weights2[i];
            vtx.wy[i] = xTx[1][2]*dx[i] + xTx[1][3]*dy[i];
            vtx.wy[i] *= weights2[i];
            vtx.wz[i] = 0.0;
        }
    }
}

@nogc
void gradients_T_lsq(SolidFVVertex vtx, int dimensions)
{
    size_t n = vtx.cloud_pos.length;
    size_t loop_init;
    loop_init = 0; // All points count.
 
    number T0;
    number[3] gradT;
    T0 = 0.0;
    foreach (i; loop_init .. n) { T0 += *(vtx.cloud_T[i]); }
    T0 /= n;
    gradT[0] = 0.0; gradT[1] = 0.0; gradT[2] = 0.0;
    foreach (i; loop_init .. n) {
        number dT = *(vtx.cloud_T[i]) - T0;
        gradT[0] += vtx.wx[i] * dT;
        gradT[1] += vtx.wy[i] * dT;
        if (dimensions == 3) { gradT[2] += vtx.wz[i] * dT; }
    }
     
    vtx.dTdx = gradT[0];
    vtx.dTdy = gradT[1];
    if (dimensions == 3) vtx.dTdz = gradT[2];
}

