module solid_boundary_flux_effect;

import std.math;
import std.stdio;
import std.json;
import std.string;
import std.format;
import util.lua;

import simcore;
import util.json_helper;
import lua_helper;
import geom;
import globaldata;
import globalconfig;
import solidfvinterface;
import ssolidblock;
import solidfvcell;
import ntypes.complex;
import nm.number;
import solidbc;

SolidBoundaryFluxEffect makeSolidBFEfromJson(JSONValue jsonData, int blk_id, int boundary)
{
    string bfeType = jsonData["type"].str;
    SolidBoundaryFluxEffect newBFE;
    switch (bfeType) {
    case "zero_flux":
        newBFE = new SolidBFE_ZeroFlux(blk_id, boundary);
        break;
    case "constant_flux":
        double fluxValue = getJSONdouble(jsonData, "flux_value", 0.0);
        newBFE = new SolidBFE_ConstantFlux(blk_id, boundary, fluxValue);
        break;
    case "constant_flux_from_solid_gas_interface":
        newBFE = new SolidBFE_ConstantFluxFromSolidGasInterface(blk_id, boundary);
        break;
    case "user_defined":
        string fname = getJSONstring(jsonData, "filename", "none");
        newBFE = new SolidBFE_UserDefined(blk_id, boundary, fname);
        break;
    default:
        string errMsg = format("ERROR: The SolidBoundaryFluxEffect type: '%s' is unknown.", bfeType);
        throw new Exception(errMsg);
    }
    return newBFE;
}

class SolidBoundaryFluxEffect {
public:
    SSolidBlock blk;
    int whichBoundary;
    string desc;

    this(int id, int boundary, string description)
    {
        blk = cast(SSolidBlock) globalBlocks[id];
        assert(blk !is null, "Oops, this should be a SolidBlock object.");
        whichBoundary = boundary;
        desc = description;
    }
    void postBCconstruction() {}
    abstract void apply(double t, int tLevel);
    abstract void apply_for_interface(double t, int tLevel, SolidFVInterface f);
}

class SolidBFE_ZeroFlux : SolidBoundaryFluxEffect {
public:
    this(int id, int boundary)
    {
        super(id, boundary, "ZeroFlux");
    }

    override void apply_for_interface(double t, int tLevel, SolidFVInterface f)
    {
        f.flux = 0.0;
    }

    override void apply(double t, int tLevel)
    {
        SolidBoundaryCondition bc = blk.bc[whichBoundary];
        foreach (face; bc.faces) { apply_for_interface(t, tLevel, face); }
    }
}

class SolidBFE_ConstantFlux : SolidBoundaryFluxEffect {
public:
    this(int id, int boundary, double fluxValue)
    {
        super(id, boundary, "ConstantFlux");
        _fluxValue = fluxValue;
    }

    override void apply_for_interface(double t, int tLevel, SolidFVInterface f)
    {
        SolidBoundaryCondition bc = blk.bc[whichBoundary];
        f.flux = bc.outsigns[f.bc_idx]*_fluxValue;
    }

    override void apply(double t, int tLevel)
    {
        SolidBoundaryCondition bc = blk.bc[whichBoundary];
        foreach (face; bc.faces) { apply_for_interface(t, tLevel, face); }
    }

private:
    double _fluxValue;
}

class SolidBFE_ConstantFluxFromSolidGasInterface : SolidBoundaryFluxEffect {
public:
    this(int id, int boundary)
    {
        super(id, boundary, "SolidGasInterface");
    }

    override void apply_for_interface(double t, int tLevel, SolidFVInterface f)
    {
        auto myBC = blk.bc[whichBoundary];
        if (blk.myConfig.fluid_solid_bc_use_heat_transfer_coeff) {
            // Note that the heat_transfer_into_solid is the total heat transfer (e.g. q_conduction + q_diffusion),
            // so we are assuming here that both q_conduction and q_diffusion will diminish as the fluid and solid
            // temperatures equilibrate. TODO: we should think about whether this is a valid assumption. KAD 2023-05-09
            number htc = myBC.gasCells[f.bc_idx].heat_transfer_into_solid;
            number dT = myBC.gasCells[f.bc_idx].fs.gas.T - myBC.solidCells[f.bc_idx].T;
            myBC.ifaces[f.bc_idx].flux = -htc*dT;
        } else {
            number q = myBC.gasCells[f.bc_idx].heat_transfer_into_solid;
            myBC.ifaces[f.bc_idx].flux = -q;
        }
    }

    override void apply(double t, int tLevel)
    {
        auto myBC = blk.bc[whichBoundary];
        foreach ( idx; 0 .. myBC.ifaces.length ) {
            SolidFVInterface face = myBC.ifaces[idx];
            apply_for_interface(t, tLevel, face);
        }
    }

    /*
    override void apply(double t, int tLevel)
    {
        auto myBC = blk.bc[whichBoundary];
        if (blk.myConfig.fluid_solid_bc_use_heat_transfer_coeff) {
            foreach ( idx; 0 .. myBC.ifaces.length ) {
                // Note that the heat_transfer_into_solid is the total heat transfer (e.g. q_conduction + q_diffusion),
                // so we are assuming here that both q_conduction and q_diffusion will diminish as the fluid and solid
                // temperatures equilibrate. TODO: we should think about whether this is a valid assumption. KAD 2023-05-09
                number htc = myBC.gasCells[idx].heat_transfer_into_solid;
                number dT = myBC.gasCells[idx].fs.gas.T - myBC.solidCells[idx].T;
                myBC.ifaces[idx].flux = -htc*dT;
            }
        } else {
            foreach ( idx; 0 .. myBC.ifaces.length ) {
                number q = myBC.gasCells[idx].heat_transfer_into_solid;
                myBC.ifaces[idx].flux = -q;
            }
        }
    }
    */
}

class SolidBFE_UserDefined : SolidBoundaryFluxEffect {
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
            lua_pushcfunction(L, &luafn_sampleSolidCell);
            lua_setglobal(L, "sampleSolidCell");
        }
        if ( luaL_dofile(blk.bc[whichBoundary].myL, luafname.toStringz) != 0 ) {
            luaL_error(blk.bc[whichBoundary].myL, "error while loading user-defined b.c. file '%s':\n %s\n",
                       luafname.toStringz, lua_tostring(blk.bc[whichBoundary].myL, -1));
        }
    }

    override void apply_for_interface(double t, int tLevel, SolidFVInterface f)
    {
        // This is an approximation to allow for steady-state simulations to operate with a user-defined BIE,
        // it will not give an analytically accurate contribution to the Jacobian.
        size_t j = 0, k = 0;
        callSolidFluxUDF(t, tLevel, f.bc_idx, j, k, f, "null");
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
                    callSolidFluxUDF(t, tLevel, i, j-1, k, IFace, "north");
                }
            }
            break;
        case Face.east:
            i = blk.imax + 1;
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (j = blk.jmin; j <= blk.jmax; ++j) {
                    IFace = blk.getIfi(i, j, k);
                    callSolidFluxUDF(t, tLevel, i-1, j, k, IFace, "east");
                }
            }
            break;
        case Face.south:
            j = blk.jmin;
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (i = blk.imin; i <= blk.imax; ++i) {
                    IFace = blk.getIfj(i, j, k);
                    callSolidFluxUDF(t, tLevel, i, j, k, IFace, "south");
                }
            }
            break;
        case Face.west:
            i = blk.imin;
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (j = blk.jmin; j <= blk.jmax; ++j) {
                    IFace = blk.getIfi(i, j, k);
                    callSolidFluxUDF(t, tLevel, i, j, k, IFace, "west");
                }
            }
            break;
        case Face.top:
            k = blk.kmax + 1;
            for (i = blk.imin; i <= blk.imax; ++i) {
                for (j = blk.jmin; j <= blk.jmax; ++j) {
                    IFace = blk.getIfk(i, j, k);
                    callSolidFluxUDF(t, tLevel, i, j, k-1, IFace, "top");
                } // end j loop
            } // end for i
            break;
        case Face.bottom:
            k = blk.kmin;
            for (i = blk.imin; i <= blk.imax; ++i) {
                for (j = blk.jmin; j <= blk.jmax; ++j) {
                    IFace = blk.getIfk(i, j, k);
                    callSolidFluxUDF(t, tLevel, i, j, k, IFace, "bottom");
                } // end j loop
            } // end for i
            break;       
        }
    }

    void callSolidFluxUDF(double t, int tLevel, size_t i, size_t j, size_t k,
                           SolidFVInterface IFace, string boundaryName)
    {
        auto L = blk.bc[whichBoundary].myL;
        lua_getglobal(L, toStringz("solidFlux"));
        // Set some userful values for the caller in table
        lua_newtable(L);
        lua_pushnumber(L, t); lua_setfield(L, -2, "t");
        lua_pushnumber(L, SimState.dt_global); lua_setfield(L, -2, "dt");
        lua_pushinteger(L, SimState.step); lua_setfield(L, -2, "timeStep");
        lua_pushinteger(L, tLevel); lua_setfield(L, -2, "timeLevel");
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
            luaL_error(L, "error running user user-defined b.c. solidFlux function: %s\n",
                       lua_tostring(L, -1));
        }
        
        IFace.flux = luaL_checknumber(L, -1);
        lua_pop(L, 1);
    }
}
