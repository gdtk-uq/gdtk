/**
 * boundary_flux_effect.d
 *
 * Authors: RG and PJ
 * Date: 2015-05-07
 * Author: KD added ConstFlux
 * Date: 2015-11-10
 **/

module bc.boundary_flux_effect;

import std.stdio;
import std.json;
import std.string;
import std.conv;
import std.math;
import nm.complex;
import nm.number;
import nm.bbla;
import nm.brent;
import nm.bracketing;

import geom;
import json_helper;
import globalconfig;
import globaldata;
import fluidblock;
import sfluidblock;
import ufluidblock;
import fvcore;
import fvcell;
import fvinterface;
import solidfvcell;
import solidfvinterface;
import gas_solid_interface;
import flowstate;
import gas;
import bc;
import flowgradients;
import mass_diffusion;
//import nm.ridder;
//import nm.bracketing;

BoundaryFluxEffect make_BFE_from_json(JSONValue jsonData, int blk_id, int boundary)
{
    string bfeType = jsonData["type"].str;
    BoundaryFluxEffect newBFE;
    auto gmodel = GlobalConfig.gmodel_master;

    switch ( bfeType ) {
    case "const_flux":
        auto flowstate = new FlowState(jsonData["flowstate"], gmodel);
        newBFE = new BFE_ConstFlux(blk_id, boundary, flowstate);
        break;
    case "simple_outflow_flux":
        newBFE = new BFE_SimpleOutflowFlux(blk_id, boundary);
        break;
    case "user_defined":
        string fname = getJSONstring(jsonData, "filename", "none");
        string funcName = getJSONstring(jsonData, "function_name", "none");
        newBFE = new BFE_UserDefined(blk_id, boundary, fname, funcName);
        break;
    case "energy_flux_from_adjacent_solid":
        int otherBlock = getJSONint(jsonData, "other_block", -1);
        string otherFaceName = getJSONstring(jsonData, "other_face", "none");
        int neighbourOrientation = getJSONint(jsonData, "neighbour_orientation", 0);
        newBFE = new BFE_EnergyFluxFromAdjacentSolid(blk_id, boundary,
                                                     otherBlock, face_index(otherFaceName),
                                                     neighbourOrientation);
        break;
    case "update_energy_wall_normal_velocity":
        newBFE = new BFE_UpdateEnergyWallNormalVelocity(blk_id, boundary);
        break;
    case "thermionic_electron_flux":
        double Ar = getJSONdouble(jsonData, "Ar", 0.0);
        double phi = getJSONdouble(jsonData, "phi", 0.0);
        newBFE = new BFE_ThermionicElectronFlux(blk_id, boundary, Ar, phi, GlobalConfig.gmodel_master);
        break;
    default:
        string errMsg = format("ERROR: The BoundaryFluxEffect type: '%s' is unknown.", bfeType);
        throw new Error(errMsg);
    }

    return newBFE;
}

class BoundaryFluxEffect {
public:
    FluidBlock blk;
    int which_boundary;
    string desc;

    this(int id, int boundary, string description)
    {
        blk = cast(FluidBlock) globalBlocks[id];
        assert(blk !is null, "Oops, this should be a FluidBlock object.");
        which_boundary = boundary;
        desc = description;
    }
    void post_bc_construction() {}
    override string toString() const
    {
        return "BoundaryFluxEffect()";
    }
    void apply(double t, int gtl, int ftl)
    {
        final switch (blk.grid_type) {
        case Grid_t.unstructured_grid:
            apply_unstructured_grid(t, gtl, ftl);
            break;
        case Grid_t.structured_grid:
            apply_structured_grid(t, gtl, ftl);
        }
    }
    void apply_for_interface(double t, int gtl, int ftl, FVInterface f)
    {
        final switch (blk.grid_type) {
        case Grid_t.unstructured_grid:
            apply_for_interface_unstructured_grid(t, gtl, ftl, f);
            break;
        case Grid_t.structured_grid:
            apply_for_interface_structured_grid(t, gtl, ftl, f);
            break;
        }
    }
    abstract void apply_for_interface_unstructured_grid(double t, int gtl, int ftl, FVInterface f);
    abstract void apply_unstructured_grid(double t, int gtl, int ftl);
    abstract void apply_for_interface_structured_grid(double t, int gtl, int ftl, FVInterface f);
    abstract void apply_structured_grid(double t, int gtl, int ftl);
} // end class BoundaryFluxEffect()

class BFE_EnergyFluxFromAdjacentSolid : BoundaryFluxEffect {
public:
    int neighbourSolidBlk;
    int neighbourSolidFace;
    int neighbourOrientation;

    this(int id, int boundary,
         int otherBlock, int otherFace, int orient)
    {
        super(id, boundary, "EnergyFluxFromAdjacentSolid");
        neighbourSolidBlk = otherBlock;
        neighbourSolidFace = otherFace;
        neighbourOrientation = orient;
    }

    override string toString() const
    {
        return "BFE_EnergyFluxFromAdjacentSolid()";
    }

    @nogc
    override void apply_for_interface_unstructured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        throw new Error("BFE_EnergyFluxFromAdjacentSolid.apply_for_interface_unstructured_grid() not yet implemented");
    }

    @nogc
    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
        throw new Error("BFE_EnergyFluxFromAdjacentSolid.apply_unstructured_grid() not yet implemented");
    }

    @nogc
    override void apply_for_interface_structured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        auto myBC = blk.bc[which_boundary];
        number dxG, dyG, dzG, dnG, dxS, dyS, dzS, dnS;
        number kG_dnG, kS_dnS, cosA, cosB, cosC;
        number T, q;
        int outsign;

        switch(which_boundary){
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

        cosA = myBC.ifaces[f.i_bndry].n.x;
        cosB = myBC.ifaces[f.i_bndry].n.y;
        cosC = myBC.ifaces[f.i_bndry].n.z;

        dxG = myBC.ifaces[f.i_bndry].pos.x - myBC.gasCells[f.i_bndry].pos[0].x;
        dyG = myBC.ifaces[f.i_bndry].pos.y - myBC.gasCells[f.i_bndry].pos[0].y;
        dzG = myBC.ifaces[f.i_bndry].pos.z - myBC.gasCells[f.i_bndry].pos[0].z;
        dnG = fabs(cosA*dxG + cosB*dyG + cosC*dzG);

        dxS = myBC.ifaces[f.i_bndry].pos.x - myBC.solidCells[f.i_bndry].pos.x;
        dyS = myBC.ifaces[f.i_bndry].pos.y - myBC.solidCells[f.i_bndry].pos.y;
        dzS = myBC.ifaces[f.i_bndry].pos.z - myBC.solidCells[f.i_bndry].pos.z;
        dnS = fabs(cosA*dxS + cosB*dyS + cosC*dzS);

        kG_dnG = myBC.gasCells[f.i_bndry].fs.gas.k / dnG;
        kS_dnS = myBC.solidCells[f.i_bndry].sp.k / dnS;

        T = (myBC.gasCells[f.i_bndry].fs.gas.T*kG_dnG + myBC.solidCells[f.i_bndry].T*kS_dnS) / (kG_dnG + kS_dnS);
        q = -kG_dnG * (T - myBC.gasCells[f.i_bndry].fs.gas.T);

        // Finally update properties in interfaces
        myBC.ifaces[f.i_bndry].fs.gas.T = T;
        myBC.ifaces[f.i_bndry].F.total_energy = outsign*q;
    }

    // not @nogc
    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        auto myBC = blk.bc[which_boundary];
        number dxG, dyG, dzG, dnG, dxS, dyS, dzS, dnS;
        number kG_dnG, kS_dnS, cosA, cosB, cosC;
        number T, q;
        int outsign;

        switch(which_boundary){
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
            q = -kG_dnG * (T - myBC.gasCells[i].fs.gas.T);

            // Finally update properties in interfaces
            myBC.ifaces[i].fs.gas.T = T;
            myBC.ifaces[i].F.total_energy = outsign*q;
        }

        //auto myBC = blk.bc[which_boundary];
        // computeFluxesAndTemperatures(ftl, myBC.gasCells, myBC.faces, myBC.solidCells, myBC.solidIFaces);
        // RJG 2019-09-10:
        // Disable non-isotropic code
        /*
          if (blk.myConfig.solid_has_isotropic_properties) {
          computeFluxesAndTemperatures(ftl, _gasCells, _gasIFaces, _solidCells, _solidIFaces);
          }
          else {
          computeFluxesAndTemperatures2(ftl, _gasCells, _gasIFaces, _solidCells, _solidIFaces,
          _T, _B, _A, _pivot);
          }
        */
    }

private:

    // RJG 2019-09-10:
    // Disable non-isotropic code
    /*
    // Some private working arrays.
    // We'll pack data into these can pass out
    // to a routine that can compute the flux and
    // temperatures that balance at the interface.
    FVCell[] _gasCells;
    FVInterface[] _gasIFaces;
    SolidFVCell[] _solidCells;
    SolidFVInterface[] _solidIFaces;
    number[] _T;
    number[] _B;
    Matrix!number _A;
    int[] _pivot;
    */

public:
    // RJG 2019-09-10:
    // Disable non-isotropic code
    /*
    void initSolidCellsAndIFaces()
    {
        size_t i, j, k;
        auto blk = solidBlocks[neighbourSolidBlk];
        switch ( neighbourSolidFace ) {
        case Face.south:
            j = blk.jmin;
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (i = blk.imin; i <= blk.imax; ++i) {
                    _solidCells ~= blk.getCell(i, j, k);
                    _solidIFaces ~= _solidCells[$-1].iface[Face.south];
                }
            }
            if (!blk.myConfig.solid_has_isotropic_properties) {
                // We'll need to initialise working space for the
                // linear system solve.
                auto n = _solidIFaces.length;
                _T.length = n;
                _B.length = n;
                _A = new Matrix!number(n);
                _pivot.length = n;
            }
            break;
        default:
            throw new Error("initSolidCellsAndIFaces() only implemented for SOUTH face.");
        }
    }

    void initGasCellsAndIFaces()
    {
        size_t i, j, k;
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");
        switch ( which_boundary ) {
        case Face.north:
            j = blk.jmax;
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (i = blk.imin; i <= blk.imax; ++i) {
                    _gasCells ~= blk.get_cell(i, j, k);
                    _gasIFaces ~= _gasCells[$-1].iface[Face.north];
                }
            }
            break;
        default:
            throw new Error("initGasCellsAndIFaces() only implemented for NORTH gas face.");
        }
    }
    */
}

class BFE_ConstFlux : BoundaryFluxEffect {
public:
    FlowState fstate;

private:
    number[] _massf;
    number _e, _rho, _p, _u, _v;
    FlowState _fstate;
    int _nsp;

public:
    this(int id, int boundary, in FlowState fstate)
    {
        /+ We only need to gather the freestream values once at
         + the start of simulation since we are interested in
         + applying a constant flux as the incoming boundary condition.
         + Note that, at this time, the gmodel held by the block is notavailable.
        +/
        auto gmodel = GlobalConfig.gmodel_master;
        super(id, boundary, "Const_Flux");
        _u = fstate.vel.x;
        _v = fstate.vel.y;
        // [TODO]: Kyle, think about z component.
        _p = fstate.gas.p;
        _rho = fstate.gas.rho;
        _e = gmodel.internal_energy(fstate.gas);
        _nsp = gmodel.n_species;
        version(multi_species_gas) {
            _massf.length = _nsp;
            for (int _isp=0; _isp < _nsp; _isp++) {
                _massf[_isp] = fstate.gas.massf[_isp];
            }
        }
        this.fstate = fstate.dup();
    }

    override string toString() const
    {
        return "BFE_ConstFlux";
    }

    @nogc
    override void apply_for_interface_unstructured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        throw new Error("BFE_ConstFlux.apply_for_interface_unstructured_grid() not yet implemented");
    }

    @nogc
    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
        throw new Error("BFE_ConstFlux.apply_unstructured_grid() not yet implemented");
    }

    @nogc
    override void apply_for_interface_structured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        throw new Error("BFE_ConstFlux.apply_for_interface_structured_grid() not yet implemented");
    }

    @nogc
    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");
        BoundaryCondition bc = blk.bc[which_boundary];
        //
        foreach (i, f; bc.faces) {
            // for a moving grid we need vel relative to the interface
            number _u_rel = _u - f.gvel.x;
            number _v_rel = _v - f.gvel.y;
            f.F.mass = _rho * ( _u_rel*f.n.x + _v_rel*f.n.y );
            /++ when the boundary is moving we use the relative velocity
             + between the fluid and the boundary interface to determine
             + the amount of mass flux across the cell face (above).
             + Alternatively momentum is a fluid property hence we use the
             + fluid velocity in determining the momentum flux -- this is
             + akin to saying we know how much mass flux is crossing
             + the cell face of which this mass has a momentum dependant
             + on its velocity. Since we we want this momentum flux in global
             + coordinates there is no need to rotate the velocity.
             ++/
            f.F.momentum.refx = _p * f.n.x + _u*f.F.mass;
            f.F.momentum.refy = _p * f.n.y + _v*f.F.mass;
            f.F.momentum.refz = 0.0;
            // [TODO]: Kyle, think about z component.
            f.F.total_energy = f.F.mass * (_e + 0.5*(_u*_u+_v*_v)) + _p*(_u*f.n.x+_v*f.n.y);
            version(multi_species_gas) {
                for ( int _isp = 0; _isp < _nsp; _isp++ ){
                    f.F.massf[_isp] = f.F.mass * _massf[_isp];
                }
            }
            version(multi_T_gas) {
                // [TODO]: Kyle, separate energy modes for multi-species simulations.
            }
        }
    }
}

class BFE_SimpleOutflowFlux : BoundaryFluxEffect {
public:
    this(int id, int boundary)
    {
        super(id, boundary, "Simple_Outflow_Flux");
    }

    override string toString() const
    {
        return "BFE_SimpleOutflowFlux";
    }

    @nogc
    override void apply_for_interface_unstructured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        auto blk = cast(UFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be a UFluidBlock object.");
        assert(!(blk.myConfig.MHD), "Oops, not implemented for MHD.");

        BoundaryCondition bc = blk.bc[which_boundary];
	int outsign = bc.outsigns[f.i_bndry];
	FVCell interior_cell = (outsign == 1) ? f.left_cell : f.right_cell;
	compute_outflow_flux(interior_cell.fs, outsign, blk.omegaz, f);
    }

    @nogc
    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
        auto blk = cast(UFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be a UFluidBlock object.");
        assert(!(blk.myConfig.MHD), "Oops, not implemented for MHD.");

        BoundaryCondition bc = blk.bc[which_boundary];
        foreach (i, face; bc.faces) {
            int outsign = bc.outsigns[i];
            FVCell interior_cell = (outsign == 1) ? face.left_cell : face.right_cell;
            compute_outflow_flux(interior_cell.fs, outsign, blk.omegaz, face);
        }
    }

    @nogc
    override void apply_for_interface_structured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");
        assert(!(blk.myConfig.MHD), "Oops, not implemented for MHD.");
        BoundaryCondition bc = blk.bc[which_boundary];
	int outsign = bc.outsigns[f.i_bndry];
	FVCell interior_cell = (outsign == 1) ? f.left_cells[0] : f.right_cells[0];
	compute_outflow_flux(interior_cell.fs, outsign, blk.omegaz, f);
    } // end apply_for_interface_structured_grid()

    @nogc
    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");
        assert(!(blk.myConfig.MHD), "Oops, not implemented for MHD.");
        BoundaryCondition bc = blk.bc[which_boundary];
        //
        foreach (i, f; bc.faces) {
            int outsign = bc.outsigns[i];
            FVCell interior_cell = (outsign == 1) ? f.left_cells[0] : f.right_cells[0];
            compute_outflow_flux(interior_cell.fs, outsign, blk.omegaz, f);
        }
    } // end apply_structured_grid()

private:
    @nogc
    void compute_outflow_flux(ref const(FlowState) fs, int outsign, double omegaz, ref FVInterface face)
    {
        // Flux equations use the local flow state information.
        // For a moving grid we need vel relative to the interface.
        Vector3 v_rel; v_rel.set(fs.vel); v_rel -= face.gvel;
        number mass_flux = fs.gas.rho * dot(v_rel, face.n);
        if ((outsign*mass_flux) > 0.0) {
            // We have a true outflow flux.
            face.F.mass = mass_flux;
            face.F.momentum.refx = fs.gas.p * face.n.x + fs.vel.x * mass_flux;
            face.F.momentum.refy = fs.gas.p * face.n.y + fs.vel.y * mass_flux;
            face.F.momentum.refz = fs.gas.p * face.n.z + fs.vel.z * mass_flux;
            number utot = fs.gas.u + 0.5*dot(fs.vel,fs.vel);
            version(multi_T_gas) {
                foreach (umode; fs.gas.u_modes) { utot += umode; }
            }
            version(turbulence) {
                utot += blk.myConfig.turb_model.turbulent_kinetic_energy(fs);
            }
            face.F.total_energy = mass_flux*utot + fs.gas.p*dot(fs.vel,face.n);
            if (omegaz != 0.0) {
                // Rotating frame.
                number x = face.pos.x;
                number y = face.pos.y;
                number rsq = x*x + y*y;
                // The conserved quantity is rotating-frame total energy,
                // so we need to take -(u**2)/2 off the total energy.
                // Note that rotating frame velocity u = omegaz * r.
                face.F.total_energy -= mass_flux * 0.5*omegaz*omegaz*rsq;
            }
            version(turbulence) {
                foreach (i; 0 .. blk.myConfig.turb_model.nturb) { face.F.rhoturb[i] = mass_flux * fs.turb[i]; }
            }
            version(MHD) {
                // [TODO] PJ 2018-10-25 MHD?
            }
            version(multi_species_gas) {
                foreach (i; 0 .. face.F.massf.length) { face.F.massf[i] = mass_flux * fs.gas.massf[i]; }
            }
            version(multi_T_gas) {
                foreach (i; 0 .. face.F.energies.length) { face.F.energies[i] = mass_flux * fs.gas.u_modes[i]; }
            }
        } else {
            // We have a situation where the nominal mass flux
            // indicates that flow should be coming into the domain.
            // Since we really do not want to have this happen,
            // we close off the face and think of it as a wall.
            face.F.mass = 0.0;
            face.F.momentum.set(face.n); face.F.momentum *= fs.gas.p;
            face.F.total_energy = fs.gas.p*dot(face.gvel,face.n);
            version(turbulence) {
                foreach (i; 0 .. blk.myConfig.turb_model.nturb) { face.F.rhoturb[i] = 0.0; }
            }
            version(MHD) {
                // [TODO] PJ 2018-10-25 MHD?
            }
            version(multi_species_gas) {
                foreach (i; 0 .. face.F.massf.length) { face.F.massf[i] = 0.0; }
            }
            version(multi_T_gas) {
                foreach (i; 0 .. face.F.energies.length) { face.F.energies[i] = 0.0; }
            }
        }
        return;
    }
} // end class BFE_SimpleOutflowFlux


/**
 * BFE_UpdateEnergyWallNormalVelocity is a boundary Flux Effect
 * that can be called for moving walls that have a wall normal
 * velocity component.
 * It operates by incrementing total_energy to correct for work
 * done on fluid by moving wall:
 * total_energy += pressure * Wall_normal_velocity
 *
*/
class BFE_UpdateEnergyWallNormalVelocity : BoundaryFluxEffect {
public:
    this(int id, int boundary)
    {
        // Don't need to do anything specific
        super(id, boundary, "UpdateEnergyWallNormalVelocity");
    }

    override string toString() const
    {
        return "BFE_UpdateEnergyWallNormalVelocity";
    }

    @nogc
    override void apply_for_interface_unstructured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        throw new Error("BFE_UpdateEnergyWallNormalVelocity.apply_for_interface_unstructured_grid() not yet implemented");
    }

    @nogc
    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
        throw new Error("BFE_UpdateEnergyWallNormalVelocity.apply_unstructured_grid() not yet implemented");
    }

    @nogc
    override void apply_for_interface_structured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        throw new Error("BFE_UpdateEnergyWallNormalVelocity.apply_for_interface_structured_grid() not yet implemented");
    }

    @nogc
    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        Vector3 nx,ny,nz;
        nx.set(1,0,0); ny.set(0,1,0); nz.set(0,0,1);
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");
        BoundaryCondition bc = blk.bc[which_boundary];
        //
        foreach (i, f; bc.faces) {
            f.F.total_energy = f.fs.gas.p * dot(f.n, f.gvel);
            f.F.momentum.refx = f.fs.gas.p * dot(f.n, nx);
            f.F.momentum.refy = f.fs.gas.p * dot(f.n, ny);
            f.F.momentum.refz = f.fs.gas.p * dot(f.n, nz);
            f.F.mass = 0.;
        }
    } // end apply_structured_grid()
} // end BFE_UpdateEnergyWallNormalVelocity

class BFE_ThermionicElectronFlux : BoundaryFluxEffect {
    this(int id, int boundary, double Ar, double phi, GasModel gmodel)
    {
        super(id, boundary, "ThermionicElectronFlux");
        this.Ar = Ar;
        this.phi = phi*Qe;  // Convert phi from input 'eV' to 'J'
        if (!gmodel.is_plasma)
            throw new Error("ThermionicElectronFlux Flux Effect requires a gas model with electrons");
        is_one_temperature = gmodel.n_modes==0;
        electron_index = gmodel.species_index("e-");
        Me = gmodel.mol_masses[electron_index];
    }

    override string toString() const
    {
        return "BFE_ThermionicElectronEmission(" ~
            "Work Function =" ~ to!string(phi/Qe) ~
            ", Richardson Constant=" ~ to!string(Ar) ~
            ")";
    }

    @nogc
    override void apply_for_interface_unstructured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        throw new Error("BIE_ThermionicElectronEmission.apply_for_interface_unstructured_grid() not yet implemented");
    }

    @nogc
    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
        throw new Error("BIE_ThermionicElectronEmission.apply_unstructured_grid() not yet implemented");
    }

    @nogc
    override void apply_for_interface_structured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        throw new Error("BIE_ThermionicElectronEmission.apply_for_interface_structured_grid() not yet implemented");
    }

    @nogc
    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");
        BoundaryCondition bc = blk.bc[which_boundary];

        foreach (i, f; bc.faces) {
            int sign = bc.outsigns[i];
            //version(multi_species_gas){
            //    number emf = electron_mass_flux(f);
            //    f.F.massf[electron_index] -= sign*emf;
            //}

            number eef = electron_energy_flux(f);
            f.F.total_energy -= sign*eef;
            version(multi_T_gas) { f.F.energies[0] -= sign*eef; }
        }
    }

protected:
    // Function inputs from Eilmer4 .lua simulation input
    double Ar;          // Richardson constant, material-dependent
    double phi;         // Work function, material dependent. Input units in eV,
                        // this gets converted to Joules by multiplying by Elementary charge, Qe
    size_t electron_index;
    number Me;
    bool is_one_temperature;
    // Constants used in analysis
    immutable double kb = 1.38064852e-23;     // Boltzmann constant.          Units: (m^2 kg)/(s^2 K^1)
    immutable double Qe = 1.60217662e-19;     // Elementary charge.           Units: C
    // Faraday's constant Qe*Na, but exact as of the 2019 redefinition of base SI units
    immutable double Faraday  = 96485.3321233100184;

    @nogc
    number electron_mass_flux(const FVInterface f)
    {
    /*
        Compute the electron flux per unit area due to thermionic emission
    */
        number Tw;
        version(multi_T_gas) {
            Tw = is_one_temperature ? f.fs.gas.T : f.fs.gas.T_modes[0];
        } else {
            Tw = f.fs.gas.T;
        }
        return Ar*Tw*Tw*exp(-phi/kb/Tw)/Faraday*Me;
    }

    @nogc
    number electron_energy_flux(const FVInterface f)
    {
    /*
        Compute the energy flux per unit area due to thermionic emission
    */
        number Tw;
        version(multi_T_gas) {
            Tw = is_one_temperature ? f.fs.gas.T : f.fs.gas.T_modes[0];
        } else {
            Tw = f.fs.gas.T;
        }
        return Ar*Tw*Tw/Qe*exp(-phi/kb/Tw)*(phi + 2*kb*Tw);
    }

} // end class BFE_ThermionicElectronFlux
