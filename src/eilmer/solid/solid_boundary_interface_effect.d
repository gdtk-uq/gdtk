module solid_boundary_interface_effect;

import std.stdio;
import std.json;
import std.string;
import std.format;
import util.lua;
import std.math;
import nm.complex;
import nm.number;

import simcore;
import json_helper;
import geom;
import globaldata;
import globalconfig;
import solidfvinterface;
import ssolidblock;
import solidfvcell;

SolidBoundaryInterfaceEffect makeSolidBIEfromJson(JSONValue jsonData, int blk_id, int boundary)
{
    string bieType = jsonData["type"].str;
    SolidBoundaryInterfaceEffect newBIE;
    switch (bieType) {
    case "fixed_temperature":
        double Twall = getJSONdouble(jsonData, "Twall", 300.0);
        newBIE = new SolidBIE_FixedT(blk_id, boundary, Twall);
        break;
    case "copy_adjacent_cell_temperature":
        newBIE = new SolidBIE_CopyAdjacentCellT(blk_id, boundary);
        break;
    case "user_defined":
        string fname = getJSONstring(jsonData, "filename", "none");
        newBIE = new SolidBIE_UserDefined(blk_id, boundary, fname);
        break;
    case "connection_boundary":
        int otherBlk = getJSONint(jsonData, "otherBlock", -1);
        string otherFace = getJSONstring(jsonData, "otherFace", "none");
        int orientation = getJSONint(jsonData, "orientation", -1);
        newBIE = new SolidBIE_ConnectionBoundary(blk_id, boundary,
                                                 otherBlk, face_index(otherFace), orientation);
        break;
    case "temperature_and_flux_from_gas_solid_interface":
        newBIE = new SolidBIE_TemperatureAndFluxFromSolidGasInterface(blk_id, boundary);
        break;
    default:
        string errMsg = format("ERROR: The SolidBoundaryInterfaceEffect type: '%s' is unknown.", bieType);
        throw new Exception(errMsg);
    }
    return newBIE;
}

class SolidBoundaryInterfaceEffect {
public:
    SSolidBlock blk;
    int whichBoundary;
    string desc;
    
    this(int id, int boundary, string description) {
        blk = cast(SSolidBlock) globalBlocks[id];
        assert(blk !is null, "Oops, this should be a SSolidBlock object.");
        whichBoundary = boundary;
        desc = description;
    }
    void postBCconstruction() {}
    abstract void apply(double t, int tLevel);
}

class SolidBIE_FixedT : SolidBoundaryInterfaceEffect {
public:
    this(int id, int boundary, double Twall)
    {
        super(id, boundary, "FixedT");
        _Twall = Twall;
    }

    override void apply(double t, int tLevel)
    {
        size_t i, j, k;
        SolidFVInterface IFace;
        
        final switch (whichBoundary) {
        case Face.north:
            j = blk.jmax + 1;
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (i = blk.imin; i <= blk.imax; ++i) {
                    IFace = blk.getIfj(i, j, k);
                    IFace.T = _Twall;
                }
            }
            break;
        case Face.east:
            i = blk.imax + 1;
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (j = blk.jmin; j <= blk.jmax; ++j) {
                    IFace = blk.getIfi(i, j, k);
                    IFace.T = _Twall;
                }
            }
            break;
        case Face.south:
            j = blk.jmin;
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (i = blk.imin; i <= blk.imax; ++i) {
                    IFace = blk.getIfj(i, j, k);
                    IFace.T = _Twall;
                }
            }
            break;
        case Face.west:
            i = blk.imin;
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (j = blk.jmin; j <= blk.jmax; ++j) {
                    IFace = blk.getIfi(i, j, k);
                    IFace.T = _Twall;
                }
            }
            break;
        case Face.top:
            k = blk.kmax + 1;
            for (i = blk.imin; i <= blk.imax; ++i) {
                for (j = blk.jmin; j <= blk.jmax; ++j) {
                    IFace = blk.getIfk(i, j, k);
                    IFace.T = _Twall;
                } // end j loop
            } // end for i
            break;
        case Face.bottom:
            k = blk.kmin;
            for (i = blk.imin; i <= blk.imax; ++i) {
                for (j = blk.jmin; j <= blk.jmax; ++j) {
                    IFace = blk.getIfk(i, j, k);
                    IFace.T = _Twall;
                } // end j loop
            } // end for i
            break;
        }   
    }
    
private:
    double _Twall;
}

class SolidBIE_CopyAdjacentCellT : SolidBoundaryInterfaceEffect {
public:
    this(int id, int boundary)
    {
        super(id, boundary, "CopyAdjacentCellT");
    }

    override void apply(double t, int tLevel)
    {
        size_t i, j, k;
        SolidFVCell cell;
        SolidFVInterface IFace;
        
        final switch (whichBoundary) {
        case Face.north:
            j = blk.jmax + 1;
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (i = blk.imin; i <= blk.imax; ++i) {
                    cell = blk.getCell(i, j-1, k);
                    IFace = blk.getIfj(i, j, k);
                    IFace.T = cell.T;
                }
            }
            break;
        case Face.east:
            i = blk.imax + 1;
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (j = blk.jmin; j <= blk.jmax; ++j) {
                    cell = blk.getCell(i-1, j, k);
                    IFace = blk.getIfi(i, j, k);
                    IFace.T = cell.T;
                }
            }
            break;
        case Face.south:
            j = blk.jmin;
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (i = blk.imin; i <= blk.imax; ++i) {
                    cell = blk.getCell(i, j, k);
                    IFace = blk.getIfj(i, j, k);
                    IFace.T = cell.T;
                }
            }
            break;
        case Face.west:
            i = blk.imin;
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (j = blk.jmin; j <= blk.jmax; ++j) {
                    cell = blk.getCell(i, j, k);
                    IFace = blk.getIfi(i, j, k);
                    IFace.T = cell.T;
                }
            }
            break;
        case Face.top:
            k = blk.kmax + 1;
            for (i = blk.imin; i <= blk.imax; ++i) {
                for (j = blk.jmin; j <= blk.jmax; ++j) {
                    cell = blk.getCell(i, j, k-1);
                    IFace = blk.getIfk(i, j, k);
                    IFace.T = cell.T;
                } // end j loop
            } // end for i
            break;
        case Face.bottom:
            k = blk.kmin;
            for (i = blk.imin; i <= blk.imax; ++i) {
                for (j = blk.jmin; j <= blk.jmax; ++j) {
                    cell = blk.getCell(i, j, k);
                    IFace = blk.getIfk(i, j, k);
                    IFace.T = cell.T;
                } // end j loop
            } // end for i
            break;
        }
    }
}


class SolidBIE_UserDefined : SolidBoundaryInterfaceEffect {
public:
    string luafname;
    this(int id, int boundary, string fname)
    {
        super(id, boundary, "UserDefined");
        luafname = fname;
    }
    override void postBCconstruction()
    {
        if (blk.bc[whichBoundary].myL == null) {
            blk.bc[whichBoundary].myL = luaL_newstate();
            auto L = blk.bc[whichBoundary].myL;
            luaL_openlibs(L);
            lua_pushinteger(L, blk.id); lua_setglobal(L, "blkId");
            lua_pushinteger(L, blk.nicell); lua_setglobal(L, "nicell");
            lua_pushinteger(L, blk.njcell); lua_setglobal(L, "njcell");
            lua_pushinteger(L, blk.nkcell); lua_setglobal(L, "nkcell");
            lua_pushinteger(L, Face.north); lua_setglobal(L, "north");
            lua_pushinteger(L, Face.east); lua_setglobal(L, "east");
            lua_pushinteger(L, Face.south); lua_setglobal(L, "south");
            lua_pushinteger(L, Face.west); lua_setglobal(L, "west");
            lua_pushinteger(L, Face.top); lua_setglobal(L, "top");
            lua_pushinteger(L, Face.bottom); lua_setglobal(L, "bottom");
        }
        if ( luaL_dofile(blk.bc[whichBoundary].myL, luafname.toStringz) != 0 ) {
            luaL_error(blk.bc[whichBoundary].myL, "error while loading user-defined b.c. file '%s':\n %s\n",
                       luafname.toStringz, lua_tostring(blk.bc[whichBoundary].myL, -1));
        }
    }
    override void apply(double t, int tLevel)
    {
        size_t i, j, k;
        SolidFVInterface IFace;

        final switch (whichBoundary) {
        case Face.north:
            j = blk.jmax + 1;
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (i = blk.imin; i <= blk.imax; ++i) {
                    IFace = blk.getIfj(i, j, k);
                    callSolidIfaceUDF(t, tLevel, i, j, k, IFace, "north");
                }
            }
            break;
        case Face.east:
            i = blk.imax + 1;
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (j = blk.jmin; j <= blk.jmax; ++j) {
                    IFace = blk.getIfi(i, j, k);
                    callSolidIfaceUDF(t, tLevel, i, j, k, IFace, "east");
                }
            }
            break;
        case Face.south:
            j = blk.jmin;
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (i = blk.imin; i <= blk.imax; ++i) {
                    IFace = blk.getIfj(i, j, k);
                    callSolidIfaceUDF(t, tLevel, i, j, k, IFace, "south");
                }
            }
            break;
        case Face.west:
            i = blk.imin;
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (j = blk.jmin; j <= blk.jmax; ++j) {
                    IFace = blk.getIfi(i, j, k);
                    callSolidIfaceUDF(t, tLevel, i, j, k, IFace, "west");
                }
            }
            break;
        case Face.top:
            k = blk.kmax + 1;
            for (i = blk.imin; i <= blk.imax; ++i) {
                for (j = blk.jmin; j <= blk.jmax; ++j) {
                    IFace = blk.getIfk(i, j, k);
                    callSolidIfaceUDF(t, tLevel, i, j, k, IFace, "top");
                } // end j loop
            } // end for i
            break;
        case Face.bottom:
            k = blk.kmin;
            for (i = blk.imin; i <= blk.imax; ++i) {
                for (j = blk.jmin; j <= blk.jmax; ++j) {
                    IFace = blk.getIfk(i, j, k);
                    callSolidIfaceUDF(t, tLevel, i, j, k, IFace, "bottom");
                } // end j loop
            } // end for i
            break;
        }
    }
    
    void callSolidIfaceUDF(double t, int tLevel, size_t i, size_t j, size_t k,
                           SolidFVInterface IFace, string boundaryName)
    {
        auto L = blk.bc[whichBoundary].myL;
        lua_getglobal(L, toStringz("solidInterface"));
        // Set some userful values for the caller in table
        lua_newtable(L);
        lua_pushnumber(L, t); lua_setfield(L, -2, "t");
        lua_pushnumber(L, SimState.dt_global); lua_setfield(L, -2, "dt");
        lua_pushinteger(L, SimState.step); lua_setfield(L, -2, "timeStep");
        lua_pushinteger(L, tLevel); lua_setfield(L, -2, "timeLevel");
        lua_pushinteger(L, whichBoundary); lua_setfield(L, -2, "boundaryId");
        lua_pushstring(L, boundaryName.toStringz); lua_setfield(L, -2, "boundaryName");
        lua_pushnumber(L, IFace.pos.x); lua_setfield(L, -2, "x");
        lua_pushnumber(L, IFace.pos.y); lua_setfield(L, -2, "y");
        lua_pushnumber(L, IFace.pos.z); lua_setfield(L, -2, "z");
        lua_pushnumber(L, IFace.n.x); lua_setfield(L, -2, "csX");
        lua_pushnumber(L, IFace.n.y); lua_setfield(L, -2, "csY");
        lua_pushnumber(L, IFace.n.z); lua_setfield(L, -2, "csZ");
        lua_pushnumber(L, IFace.t1.x); lua_setfield(L, -2, "csX1");
        lua_pushnumber(L, IFace.t1.y); lua_setfield(L, -2, "csY1");
        lua_pushnumber(L, IFace.t1.z); lua_setfield(L, -2, "csZ1");
        lua_pushnumber(L, IFace.t2.x); lua_setfield(L, -2, "csX2");
        lua_pushnumber(L, IFace.t2.y); lua_setfield(L, -2, "csY2");
        lua_pushnumber(L, IFace.t2.z); lua_setfield(L, -2, "csZ2");
        lua_pushinteger(L, i); lua_setfield(L, -2, "i");
        lua_pushinteger(L, j); lua_setfield(L, -2, "j");
        lua_pushinteger(L, k); lua_setfield(L, -2, "k");

        // Call function and expect back a temperature value.
        int number_args = 1;
        int number_results = 1;
        if ( lua_pcall(L, number_args, number_results, 0) != 0 ) {
            luaL_error(L, "error running user user-defined b.c. solidInterface function: %s\n",
                       lua_tostring(L, -1));
        }
        
        IFace.T = luaL_checknumber(L, -1);
        lua_pop(L, 1);
    }
}

class SolidBIE_ConnectionBoundary : SolidBoundaryInterfaceEffect {
public:
    SSolidBlock neighbourBlk;
    int neighbourFace;
    int neighbourOrientation;

    this(int id, int boundary,
         int otherBlock, int otherFace, int orient)
    {
        super(id, boundary, "ConnectionBoundary");
        neighbourBlk = cast(SSolidBlock) globalBlocks[otherBlock];
        assert(neighbourBlk !is null, "Oops, this should be a SSolidBlock object.");
        neighbourFace = otherFace;
        neighbourOrientation = orient;
    }

    override void apply(double t, int tLevel) {
        size_t i, j, k;
        SolidFVCell Lft, Rght;
        SolidFVInterface IFace;

        if (GlobalConfig.dimensions == 2) {
            switch (whichBoundary) {
            case Face.north:
                switch (neighbourFace) {
                case Face.north:
                    throw new Error("SolidBFE_ConnectionBoundary not implemented for NORTH-NORTH connections.");
                case Face.east:
                    throw new Error("SolidBFE_ConnectionBoundary not implemented for NORTH-EAST connections.");
                case Face.south:
                    for (i = blk.imin; i <= blk.imax; ++i) {
                        Lft = blk.getCell(i, blk.jmax);
                        Rght = neighbourBlk.getCell(i, neighbourBlk.jmin);
                        IFace = blk.getIfj(i, blk.jmax+1);
                        computeBoundaryFlux(Lft, Rght, IFace,
                                                Lft.sp.k, Rght.sp.k);
                    }
                    break;
                case Face.top:
                    throw new Error("SolidBFE_ConnectionBoundary not implemented for NORTH-SOUTH connections.");
                case Face.bottom:
                    throw new Error("SolidBFE_ConnectionBoundary not implemented for NORTH-SOUTH connections.");
                case Face.west:
                    throw new Error("SolidBFE_ConnectionBoundary not implemented for NORTH-WEST connections.");
                default:
                    throw new Error("SolidBFE_ConnectionBoundary: connection type not available in 2D.");
                }
                break;
            case Face.east:
                switch (neighbourFace) {
                case Face.north:
                    throw new Error("SolidBFE_ConnectionBoundary not implemented for EAST-NORTH connections.");
                case Face.east:
                    throw new Error("SolidBFE_ConnectionBoundary not implemented for EAST-EAST connections.");
                case Face.south:
                    throw new Error("SolidBFE_ConnectionBoundary not implemented for EAST-SOUTH connections.");
                case Face.west:
                    for (j = blk.jmin; j <= blk.jmax; ++j) {
                        Lft = blk.getCell(blk.imax, j);
                        Rght = neighbourBlk.getCell(neighbourBlk.imin, j);
                        IFace = blk.getIfi(blk.imax+1, j);
                        computeBoundaryFlux(Lft, Rght, IFace,
                                            Lft.sp.k, Rght.sp.k);
                    }
                    break;
                default:
                    throw new Error("SolidBFE_ConnectionBoundary: connection type not available in 2D.");
                }
                break;
            case Face.south:
                switch (neighbourFace) {
                case Face.north:
                    for (i = blk.imin; i <= blk.imax; ++i) {
                        Lft = neighbourBlk.getCell(i, neighbourBlk.jmax);
                        Rght = blk.getCell(i, blk.jmin);
                        IFace = blk.getIfj(i, blk.jmin);
                        computeBoundaryFlux(Lft, Rght, IFace,
                                            Lft.sp.k, Rght.sp.k);
                    }
                    break;
                case Face.east:
                    throw new Error("SolidBFE_ConnectionBoundary not implemented for SOUTH-EAST connections.");
                case Face.south:                    
                    throw new Error("SolidBFE_ConnectionBoundary not implemented for SOUTH-SOUTH connections.");
                case Face.top:
                    throw new Error("SolidBFE_ConnectionBoundary not implemented for SOUTH-TOP connections.");
                case Face.bottom:
                    throw new Error("SolidBFE_ConnectionBoundary not implemented for SOUTH-BOTTOM connections.");
                case Face.west:
                    throw new Error("SolidBFE_ConnectionBoundary not implemented for SOUTH-WEST connections.");
                default:
                    throw new Error("SolidBFE_ConnectionBoundary: connection type not available in 2D.");
                }
                break;
            case Face.west:
                switch (neighbourFace) {
                case Face.north:
                    throw new Error("SolidBFE_ConnectionBoundary not implemented for WEST-NORTH connections.");
                case Face.east:
                    for (j = blk.jmin; j <= blk.jmax; ++j) {
                        Lft = neighbourBlk.getCell(neighbourBlk.imax, j);
                        Rght = blk.getCell(blk.imin, j);
                        IFace = blk.getIfi(blk.imin, j);
                        computeBoundaryFlux(Lft, Rght, IFace,
                                            Lft.sp.k, Rght.sp.k);
                    }
                    break;
                case Face.south:
                    throw new Error("SolidBFE_ConnectionBoundary not implemented for WEST-SOUTH connections.");
                case Face.west:
                    throw new Error("SolidBFE_ConnectionBoundary not implemented for WEST-WEST connections.");
                default:
                    throw new Error("SolidBFE_ConnectionBoundary: boundary connection not available in 2D.");
                }
                break;
            default:
                throw new Error("SolidBFE_ConnectionBoundary: boundary not available in 2D.");
            }
        }
        else { // 3D
            switch (whichBoundary) {
            case Face.north:
                switch (neighbourFace) {
                case Face.north:
                    throw new Error("SolidBFE_ConnectionBoundary not implemented for NORTH-NORTH connections.");
                case Face.east:
                    throw new Error("SolidBFE_ConnectionBoundary not implemented for NORTH-EAST connections.");
                case Face.south:
                    for (k = blk.kmin; k <= blk.kmax; ++k) {
                        for (i = blk.imin; i <= blk.imax; ++i) {
                            Lft = blk.getCell(i, blk.jmax, k);
                            Rght = neighbourBlk.getCell(i, neighbourBlk.jmin, k);
                            IFace = blk.getIfj(i, blk.jmax+1, k);
                            computeBoundaryFlux(Lft, Rght, IFace,
                                                Lft.sp.k, Rght.sp.k);
                        }
                    }
                    break;
                case Face.west:
                    throw new Error("SolidBFE_ConnectionBoundary not implemented for NORTH-WEST connections.");
                case Face.top:
                    throw new Error("SolidBFE_ConnectionBoundary not implemented for NORTH-TOP connections.");
                case Face.bottom:
                    throw new Error("SolidBFE_ConnectionBoundary not implemented for NORTH-BOTTOM connections.");
                default:
                    throw new Error("SolidBFE_ConnectionBoundary: connection type not available in 2D.");
                }
                break;
            case Face.east:
                switch (neighbourFace) {
                case Face.north:
                    throw new Error("SolidBFE_ConnectionBoundary not implemented for EAST-NORTH connections.");
                case Face.east:
                    throw new Error("SolidBFE_ConnectionBoundary not implemented for EAST-EAST connections.");
                case Face.south:
                    throw new Error("SolidBFE_ConnectionBoundary not implemented for EAST-SOUTH connections.");
                case Face.west:
                    for (k = blk.kmin; k <= blk.kmax; ++k) {
                        for (j = blk.jmin; j <= blk.jmax; ++j) {
                            Lft = blk.getCell(blk.imax, j, k);
                            Rght = neighbourBlk.getCell(neighbourBlk.imin, j, k);
                            IFace = blk.getIfi(blk.imax+1, j, k);
                            computeBoundaryFlux(Lft, Rght, IFace,
                                                 Lft.sp.k, Rght.sp.k);
                        }
                    }
                    break;
                case Face.top:
                    throw new Error("SolidBFE_ConnectionBoundary not implemented for EAST-TOP connections.");
                case Face.bottom:
                    throw new Error("SolidBFE_ConnectionBoundary not implemented for EAST-BOTTOM connections.");
                default:
                    throw new Error("SolidBFE_ConnectionBoundary: connection type not available in 2D.");
                }
                break;
            case Face.south:
                switch (neighbourFace) {
                case Face.north:
                    for (k = blk.kmin; k <= blk.kmax; ++k) {
                        for (i = blk.imin; i <= blk.imax; ++i) {
                            Lft = neighbourBlk.getCell(i, neighbourBlk.jmax, k);
                            Rght = blk.getCell(i, blk.jmin, k);
                            IFace = blk.getIfj(i, blk.jmin, k);
                            computeBoundaryFlux(Lft, Rght, IFace,
                                                Lft.sp.k, Rght.sp.k);
                        }
                    }
                    break;
                case Face.east:
                    throw new Error("SolidBFE_ConnectionBoundary not implemented for SOUTH-EAST connections.");
                case Face.south:                    
                    throw new Error("SolidBFE_ConnectionBoundary not implemented for SOUTH-SOUTH connections.");
                case Face.west:
                    throw new Error("SolidBFE_ConnectionBoundary not implemented for SOUTH-WEST connections.");
                case Face.top:
                    throw new Error("SolidBFE_ConnectionBoundary not implemented for SOUTH-TOP connections.");
                case Face.bottom:
                    throw new Error("SolidBFE_ConnectionBoundary not implemented for SOUTH-BOTTOM connections.");
                default:
                    throw new Error("SolidBFE_ConnectionBoundary: connection type not available in 2D.");
                }
                break;
            case Face.west:
                switch (neighbourFace) {
                case Face.north:
                    throw new Error("SolidBFE_ConnectionBoundary not implemented for WEST-NORTH connections.");
                case Face.east:
                    for (k = blk.kmin; k <= blk.kmax; ++k) {
                        for (j = blk.jmin; j <= blk.jmax; ++j) {
                            Lft = neighbourBlk.getCell(neighbourBlk.imax, j, k);
                            Rght = blk.getCell(blk.imin, j, k);
                            IFace = blk.getIfi(blk.imin, j, k);
                            computeBoundaryFlux(Lft, Rght, IFace,
                                                Lft.sp.k, Rght.sp.k);
                        }
                    }
                    break;
                case Face.south:
                    throw new Error("SolidBFE_ConnectionBoundary not implemented for WEST-SOUTH connections.");
                case Face.west:
                    throw new Error("SolidBFE_ConnectionBoundary not implemented for WEST-WEST connections.");
                case Face.top:
                    throw new Error("SolidBFE_ConnectionBoundary not implemented for EAST-SOUTH connections.");
                case Face.bottom:
                    throw new Error("SolidBFE_ConnectionBoundary not implemented for EAST-SOUTH connections.");
                default:
                    throw new Error("SolidBFE_ConnectionBoundary: boundary connection not available in 2D.");
                }
                break;
            case Face.top:
                switch (neighbourFace) {
                case Face.north:
                    throw new Error("SolidBFE_ConnectionBoundary not implemented for TOP-NORTH connections.");
                case Face.east:
                    throw new Error("SolidBFE_ConnectionBoundary not implemented for TOP-EAST connections.");
                case Face.south:
                    throw new Error("SolidBFE_ConnectionBoundary not implemented for TOP-SOUTH connections.");
                case Face.west:
                    throw new Error("SolidBFE_ConnectionBoundary not implemented for TOP-WEST connections.");
                case Face.top:
                    throw new Error("SolidBFE_ConnectionBoundary not implemented for TOP-TOP connections.");
                case Face.bottom:
                    for (i = blk.imin; i <= blk.imax; ++i) {
                        for (j = blk.jmin; j <= blk.jmax; ++j) {
                            Lft = blk.getCell(i, j, blk.kmax);
                            Rght = neighbourBlk.getCell(i, j, neighbourBlk.kmin);
                            IFace = blk.getIfk(i, j, blk.kmax+1);
                            computeBoundaryFlux(Lft, Rght, IFace,
                                                Lft.sp.k, Rght.sp.k);
                        }
                    }
                    break;
                default:
                    throw new Error("SolidBFE_ConnectionBoundary: boundary connection not available in 2D.");
                }
                break;
            case Face.bottom:
                switch (neighbourFace) {
                case Face.north:
                    throw new Error("SolidBFE_ConnectionBoundary not implemented for BOTTOM-NORTH connections.");
                case Face.east:
                    throw new Error("SolidBFE_ConnectionBoundary not implemented for BOTTOM-EAST connections.");
                case Face.south:
                    throw new Error("SolidBFE_ConnectionBoundary not implemented for BOTTOM-SOUTH connections.");
                case Face.west:
                    throw new Error("SolidBFE_ConnectionBoundary not implemented for BOTTOM-WEST connections.");
                case Face.top:
                    for (i = blk.imin; i <= blk.imax; ++i) {
                        for (j = blk.jmin; j <= blk.jmax; ++j) {
                            Lft = neighbourBlk.getCell(i, j, neighbourBlk.kmax);
                            Rght = blk.getCell(i, j, blk.kmin);
                            IFace = blk.getIfk(i, j, blk.kmin);
                            computeBoundaryFlux(Lft, Rght, IFace,
                                                Lft.sp.k, Rght.sp.k);
                        }
                    }
                    break;
                case Face.bottom:
                    throw new Error("SolidBFE_ConnectionBoundary not implemented for BOTTOM-BOTTOM connections.");
                default:
                    throw new Error("SolidBFE_ConnectionBoundary: boundary not available in 3D.");
                }
                break;
            default:
                throw new Error("SolidBFE_ConnectionBoundary: boundary not available in 3D.");
            }
        }
    }
    
    void computeBoundaryFlux(SolidFVCell Lft, SolidFVCell Rght, SolidFVInterface IFace, double kL, double kR)
    {
        Vector3 LI, IR;
        number dL, dR;
        number kL_dL, kR_dR;
        number T, q;
        
        LI = IFace.pos - Lft.pos;
        IR = Rght.pos - IFace.pos;
        dL = fabs(dot(LI, IFace.n).re);
        dR = fabs(dot(IR, IFace.n).re);
        kL_dL = kL/dL;
        kR_dR = kR/dR;

        T = (Lft.T*kL_dL + Rght.T*kR_dR)/(kL_dL + kR_dR);
        q = -kL_dL*(T - Lft.T);
        
        IFace.T = T;
        IFace.flux = q;
    }

}

class SolidBIE_TemperatureAndFluxFromSolidGasInterface : SolidBoundaryInterfaceEffect {
public:
    this(int id, int boundary)
    {
        super(id, boundary, "SolidGasInterface");
    }

    override void apply(double t, int tLevel)
    {

        auto myBC = blk.bc[whichBoundary];
        number dxG, dyG, dzG, dnG, dxS, dyS, dzS, dnS;
        number kG_dnG, kS_dnS, cosA, cosB, cosC;
        number T, q;
        int outsign;
        
        switch(whichBoundary){
        case Face.north:
            outsign = 1;
            break;
        case Face.east:
            outsign = 1;
            break;
        case Face.south:
            outsign = -1;
            break;
        case Face.west:
            outsign = -1;
            break;
        case Face.top:
            outsign = 1;
            break;
        case Face.bottom:
            outsign = -1;
            break;
        default:
            throw new Error("oops, wrong boundary id");
        } // end switch

        foreach ( i; 0 .. myBC.ifaces.length ) {
            cosA = myBC.ifaces[i].n.x;
            cosB = myBC.ifaces[i].n.y;
            cosC = myBC.ifaces[i].n.z;
            
            dxG = myBC.ifaces[i].pos.x - myBC.gasCells[i].pos[0].x;
            dyG = myBC.ifaces[i].pos.y - myBC.gasCells[i].pos[0].y;
            dzG = myBC.ifaces[i].pos.z - myBC.gasCells[i].pos[0].z;
            dnG = fabs(cosA*dxG + cosB*dyG + cosC*dzG);
            
            dxS = myBC.ifaces[i].pos.x - myBC.solidCells[i].pos.x;
            dyS = myBC.ifaces[i].pos.y - myBC.solidCells[i].pos.y;
            dzS = myBC.ifaces[i].pos.z - myBC.solidCells[i].pos.z;
            dnS = fabs(cosA*dxS + cosB*dyS + cosC*dzS);
            
            kG_dnG = myBC.gasCells[i].fs.gas.k / dnG;
            kS_dnS = myBC.solidCells[i].sp.k / dnS;
            
            T = (myBC.gasCells[i].fs.gas.T*kG_dnG + myBC.solidCells[i].T*kS_dnS) / (kG_dnG + kS_dnS);
            q = -kS_dnS * (T - myBC.solidCells[i].T);

            // Finally update properties in interfaces
            myBC.ifaces[i].T = T;
            myBC.ifaces[i].flux = outsign*q;
        }
    }
}

