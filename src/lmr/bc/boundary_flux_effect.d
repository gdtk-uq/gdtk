/**
 * boundary_flux_effect.d
 *
 * Authors: RG and PJ
 * Date: 2015-05-07
 * Author: KD added ConstFlux
 * Date: 2015-11-10
 * 2021-07-29 PJ rework ConstFlux
 **/

module bc.boundary_flux_effect;

import std.stdio;
import std.json;
import std.string;
import std.conv;
import std.math;
import ntypes.complex;
import nm.number;
import nm.bbla;
import nm.brent;
import nm.bracketing;

import geom;
import util.json_helper;
import globalconfig;
import globaldata;
import fluidblock;
import sfluidblock;
import ufluidblock;
import lmr.fluidfvcell;
import fvinterface;
import solidfvcell;
import solidfvinterface;
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
        auto flowstate = FlowState(jsonData["flowstate"], gmodel);
        double x0 = getJSONdouble(jsonData, "x0", 0.0);
        double y0 = getJSONdouble(jsonData, "y0", 0.0);
        double z0 = getJSONdouble(jsonData, "z0", 0.0);
        double r = getJSONdouble(jsonData, "r", 0.0);
        newBFE = new BFE_ConstFlux(blk_id, boundary, flowstate, x0, y0, z0, r);
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
} // end class BoundaryFluxEffect


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
        auto cqi = blk.myConfig.cqi;

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
        kS_dnS = myBC.solidCells[f.i_bndry].ss.k / dnS;

        T = (myBC.gasCells[f.i_bndry].fs.gas.T*kG_dnG + myBC.solidCells[f.i_bndry].T*kS_dnS) / (kG_dnG + kS_dnS);
        q = -kG_dnG * (T - myBC.gasCells[f.i_bndry].fs.gas.T);

        // Finally update properties in interfaces
        myBC.ifaces[f.i_bndry].fs.gas.T = T;
        myBC.ifaces[f.i_bndry].F[cqi.totEnergy] = outsign*q;
    }

    // not @nogc
    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        auto myBC = blk.bc[which_boundary];
        number dxG, dyG, dzG, dnG, dxS, dyS, dzS, dnS;
        number kG_dnG, kS_dnS, cosA, cosB, cosC;
        number T, q;
        int outsign;
        auto cqi = blk.myConfig.cqi;

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
            kS_dnS = myBC.solidCells[i].ss.k / dnS;

            T = (myBC.gasCells[i].fs.gas.T*kG_dnG + myBC.solidCells[i].T*kS_dnS) / (kG_dnG + kS_dnS);
            q = -kG_dnG * (T - myBC.gasCells[i].fs.gas.T);

            // Finally update properties in interfaces
            myBC.ifaces[i].fs.gas.T = T;
            myBC.ifaces[i].F[cqi.totEnergy] = outsign*q;
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
    FluidFVCell[] _gasCells;
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
} // end class BFE_EnergyFluxFromAdjacentSolid


class BFE_ConstFlux : BoundaryFluxEffect {
public:
    FlowState fstate;
    SourceFlow sflow;
    // The shock-fitting functions need to see the following parameters.
    double x0, y0, z0, r; // conical-flow parameters

public:
    this(int id, int boundary, in FlowState fstate, double x0, double y0, double z0, double r)
    {
        super(id, boundary, "Const_Flux");
        /+ We only need to gather the freestream values once at
         + the start of simulation since we are interested in
         + applying a constant flux as the incoming boundary condition.
         + Note that, at this time, the gmodel held by the block is not available.
        +/
        auto gmodel = GlobalConfig.gmodel_master;
        this.fstate = fstate.dup();
        this.x0 = x0;
        this.y0 = y0;
        this.z0 = z0;
        this.r = r;
        sflow = new SourceFlow(gmodel, fstate, r);
        //
        auto myblk = cast(FluidBlock) globalBlocks[id];
        if (myblk.omegaz != 0.0) {
            throw new Error("BFE_ConstFlux not implemented for rotating frame.");
        }
    }

    override string toString() const
    {
        return "BFE_ConstFlux";
    }

    @nogc
    override void apply_for_interface_unstructured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        apply_to_single_face(f);
    }

    @nogc
    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
        BoundaryCondition bc = blk.bc[which_boundary];
        foreach (i, f; bc.faces) { apply_to_single_face(f); }
    }

    @nogc
    override void apply_for_interface_structured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        apply_to_single_face(f);
    }

    @nogc
    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        BoundaryCondition bc = blk.bc[which_boundary];
        foreach (i, f; bc.faces) { apply_to_single_face(f); }
    }

private:
    @nogc
    void apply_to_single_face(FVInterface f)
    {
        auto cqi = blk.myConfig.cqi;
        // Start by assuming uniform, parallel flow.
        number p = fstate.gas.p;
        number rho = fstate.gas.rho;
        auto gmodel = blk.myConfig.gmodel;
        number u = gmodel.internal_energy(fstate.gas);
        number velx = fstate.vel.x;
        number vely = fstate.vel.y;
        number velz = (cqi.threeD) ? fstate.vel.z : to!number(0.0);
        if (r > 0.0) {
            // (Approximate) conical inflow.
            double dx = f.pos.x.re - x0;
            double dy = f.pos.y.re - y0;
            double dz = (cqi.threeD) ? f.pos.z.re - z0 : 0.0;
            double hypot = sqrt(dx*dx + dy*dy + dz*dz);
            double[4] deltas = sflow.get_rho_v_p_u_increments(hypot-r);
            rho += deltas[0];
            double v = fstate.vel.x.re + deltas[1];
            velx = v * dx/hypot;
            vely = v * dy/hypot;
            velz = (cqi.threeD) ? v * dz/hypot : 0.0;
            p += deltas[2];
            u += deltas[3];
        }
        // For a moving grid, we need gas velocity relative to the interface.
        number velx_rel = velx - f.gvel.x;
        number vely_rel = vely - f.gvel.y;
        number velz_rel = velz - f.gvel.z;
        number massFlux = rho * (velx_rel*f.n.x + vely_rel*f.n.y + velz_rel*f.n.z);
        if (cqi.mass==0) f.F[cqi.mass] = massFlux;
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
        f.F[cqi.xMom] = p*f.n.x + velx*massFlux;
        f.F[cqi.yMom] = p*f.n.y + vely*massFlux;
        if (cqi.threeD) {
            f.F[cqi.zMom] = p*f.n.z + velz*massFlux;
        }
        f.F[cqi.totEnergy] = massFlux*(u + 0.5*(velx*velx + vely*vely + velz*velz)) +
            p*(velx*f.n.x + vely*f.n.y + velz*f.n.z);
        version(multi_species_gas) {
            if (cqi.n_species > 1) {
                foreach (i; 0 .. cqi.n_species) {
                    f.F[cqi.species+i] = massFlux*fstate.gas.massf[i];
                }
            }
        }
        version(multi_T_gas) {
            foreach (i; 0 .. cqi.n_modes){
                f.F[cqi.modes+i] = massFlux*fstate.gas.u_modes[i];
            }
        }
    } // end apply_to_single_face()

} // end class BFE_ConstFlux


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
        //
        BoundaryCondition bc = blk.bc[which_boundary];
	int outsign = bc.outsigns[f.i_bndry];
	FluidFVCell interior_cell = (outsign == 1) ? f.left_cell : f.right_cell;
	compute_outflow_flux(interior_cell.fs, outsign, blk.omegaz, f);
    }

    @nogc
    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
        auto blk = cast(UFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be a UFluidBlock object.");
        assert(!(blk.myConfig.MHD), "Oops, not implemented for MHD.");
        //
        BoundaryCondition bc = blk.bc[which_boundary];
        foreach (i, face; bc.faces) {
            int outsign = bc.outsigns[i];
            FluidFVCell interior_cell = (outsign == 1) ? face.left_cell : face.right_cell;
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
	FluidFVCell interior_cell = (outsign == 1) ? f.left_cells[0] : f.right_cells[0];
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
            FluidFVCell interior_cell = (outsign == 1) ? f.left_cells[0] : f.right_cells[0];
            compute_outflow_flux(interior_cell.fs, outsign, blk.omegaz, f);
        }
    } // end apply_structured_grid()

private:
    @nogc
    void compute_outflow_flux(FlowState* fs, int outsign, double omegaz, ref FVInterface face)
    {
        // Flux equations use the local flow state information.
        // For a moving grid we need vel relative to the interface.
        auto cqi = blk.myConfig.cqi;
        Vector3 v_rel; v_rel.set(fs.vel); v_rel -= face.gvel;
        number mass_flux = fs.gas.rho * dot(v_rel, face.n);
        if ((outsign*mass_flux) > 0.0) {
            // We have a true outflow flux.
            if (cqi.mass==0) face.F[cqi.mass] = mass_flux;
            face.F[cqi.xMom] = fs.gas.p * face.n.x + fs.vel.x * mass_flux;
            face.F[cqi.yMom] = fs.gas.p * face.n.y + fs.vel.y * mass_flux;
            if (cqi.threeD) { face.F[cqi.zMom] = fs.gas.p * face.n.z + fs.vel.z * mass_flux; }
            number utot = fs.gas.u + 0.5*dot(fs.vel,fs.vel);
            version(multi_T_gas) {
                foreach (umode; fs.gas.u_modes) { utot += umode; }
            }
            version(turbulence) {
                utot += blk.myConfig.turb_model.turbulent_kinetic_energy(*fs);
            }
            face.F[cqi.totEnergy] = mass_flux*utot + fs.gas.p*dot(fs.vel,face.n);
            if (omegaz != 0.0) {
                // Rotating frame.
                number x = face.pos.x;
                number y = face.pos.y;
                number rsq = x*x + y*y;
                // The conserved quantity is rotating-frame total energy,
                // so we need to take -(u**2)/2 off the total energy.
                // Note that rotating frame velocity u = omegaz * r.
                face.F[cqi.totEnergy] -= mass_flux * 0.5*omegaz*omegaz*rsq;
            }
            version(turbulence) {
                foreach (i; 0 .. blk.myConfig.turb_model.nturb) { face.F[cqi.rhoturb+i] = mass_flux * fs.turb[i]; }
            }
            version(MHD) {
                // [TODO] PJ 2018-10-25 MHD?
            }
            version(multi_species_gas) {
                if (cqi.n_species > 1) {
                    foreach (i; 0 .. cqi.n_species) { face.F[cqi.species+i] = mass_flux * fs.gas.massf[i]; }
                }
            }
            version(multi_T_gas) {
                foreach (i; 0 .. cqi.n_modes) { face.F[cqi.modes+i] = mass_flux * fs.gas.u_modes[i]; }
            }
        } else {
            // We have a situation where the nominal mass flux
            // indicates that flow should be coming into the domain.
            // Since we really do not want to have this happen,
            // we close off the face and think of it as a wall.
            if (cqi.mass==0) face.F[cqi.mass] = 0.0;
            face.F[cqi.xMom] = face.n.x * fs.gas.p;
            face.F[cqi.yMom] = face.n.y * fs.gas.p;
            if (cqi.threeD) { face.F[cqi.zMom] = face.n.z * fs.gas.p; }
            face.F[cqi.totEnergy] = fs.gas.p * dot(face.gvel,face.n);
            version(turbulence) {
                foreach (i; 0 .. blk.myConfig.turb_model.nturb) { face.F[cqi.rhoturb+i] = 0.0; }
            }
            version(MHD) {
                // [TODO] PJ 2018-10-25 MHD?
            }
            version(multi_species_gas) {
                if (cqi.n_species > 1) {
                    foreach (i; 0 .. cqi.n_species) { face.F[cqi.species+i] = 0.0; }
                }
            }
            version(multi_T_gas) {
                foreach (i; 0 .. cqi.n_modes) { face.F[cqi.modes+i] = 0.0; }
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
        auto cqi = blk.myConfig.cqi;
        //
        foreach (i, f; bc.faces) {
            f.F[cqi.totEnergy] = f.fs.gas.p * dot(f.n, f.gvel);
            f.F[cqi.xMom] = f.fs.gas.p * dot(f.n, nx);
            f.F[cqi.yMom] = f.fs.gas.p * dot(f.n, ny);
            if (cqi.threeD) { f.F[cqi.zMom] = f.fs.gas.p * dot(f.n, nz); }
            if (cqi.mass==0) f.F[cqi.mass] = 0.;
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
        auto cqi = blk.myConfig.cqi;
        //
        foreach (i, f; bc.faces) {
            int sign = bc.outsigns[i];
            //version(multi_species_gas){
            //    number emf = electron_mass_flux(f);
            //    f.F.massf[electron_index] -= sign*emf;
            //}
            number eef = electron_energy_flux(f);
            f.F[cqi.totEnergy] -= sign*eef;
            version(multi_T_gas) { f.F[cqi.modes+0] -= sign*eef; }
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
