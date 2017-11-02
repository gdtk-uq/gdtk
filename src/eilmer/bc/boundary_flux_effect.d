/**
 * boundary_flux_effct.d
 *
 * Authors: RG and PJ
 * Date: 2015-05-07
 * Author: KD added ConstFlux
 * Date: 2015-11-10
 **/

module boundary_flux_effect;

import std.stdio;
import std.json;
import std.string;
import std.conv;
import std.math;

import geom;
import grid;
import json_helper;
import globalconfig;
import globaldata;
import block;
import sblock;
import fvcore;
import fvcell;
import fvinterface;
import solidfvcell;
import solidfvinterface;
import gas_solid_interface;
import flowstate;
import gas;
import user_defined_effects;
import flowgradients;

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
	case "energy_balance_thermionic":
	double emissivity = getJSONdouble(jsonData, "emissivity", 0.0);
	double Ar = getJSONdouble(jsonData, "Ar", 0.0);
	double phi = getJSONdouble(jsonData, "phi", 0.0);
	int ThermionicEmissionActive = getJSONint(jsonData, "ThermionicEmissionActive", 1);
	int Twall_iterations = getJSONint(jsonData, "Twall_iterations", 200);
	int Twall_subiterations = getJSONint(jsonData, "Twall_subiterations", 20);
	newBFE = new BFE_EnergyBalanceThermionic(blk_id, boundary, emissivity, Ar, phi,
				 ThermionicEmissionActive, Twall_iterations, Twall_subiterations);
	break;
    default:
	string errMsg = format("ERROR: The BoundaryFluxEffect type: '%s' is unknown.", bfeType);
	throw new Error(errMsg);
    }
    
    return newBFE;
}

class BoundaryFluxEffect {
public:
    Block blk;
    int which_boundary;
    string type;
    
    this(int id, int boundary, string _type)
    {
	blk = gasBlocks[id];
	which_boundary = boundary;
	type = _type;
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
    abstract void apply_unstructured_grid(double t, int gtl, int ftl);
    abstract void apply_structured_grid(double t, int gtl, int ftl);
} // end class BoundaryFluxEffect()

// NOTE: This GAS DOMAIN boundary effect has a large
//       and important side-effect:
//       IT ALSO SETS THE FLUX IN THE ADJACENT SOLID DOMAIN
//       AT THE TIME IT IS CALLED.

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
    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
	throw new Error("BFE_EnergyFluxFromAdjacentSolid.apply_unstructured_grid() not yet implemented");
    }

    override void apply_structured_grid(double t, int gtl, int ftl)
    {
	computeFluxesAndTemperatures(ftl, _gasCells, _gasIFaces, _solidCells, _solidIFaces);
    }

private:
    // Some private working arrays.
    // We'll pack data into these can pass out
    // to a routine that can compute the flux and
    // temperatures that balance at the interface.
    FVCell[] _gasCells;
    FVInterface[] _gasIFaces;
    SolidFVCell[] _solidCells;
    SolidFVInterface[] _solidIFaces;

public:
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
	    break;
	default:
	    throw new Error("initSolidCellsAndIFaces() only implemented for SOUTH face.");
	}
    }

    void initGasCellsAndIFaces()
    {
	size_t i, j, k;
	SBlock blk = cast(SBlock) this.blk;
	assert(blk !is null, "Oops, this should be an SBlock object.");
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
}
class BFE_ConstFlux : BoundaryFluxEffect {
public:
    FlowState fstate;

private:
    double[] _massf;
    double _e, _rho, _p, _u, _v;
    FlowState _fstate;
    int _nsp;
   
public:  
    this(int id, int boundary, in FlowState fstate)
    {
	/+ We only need to gather the freestream values once at
	 + the start of simulation since we are interested in
	 + applying a constant flux as the incoming boundary
	 + condition.
	+/
	//auto gmodel = blk.myConfig.gmodel;
	auto gmodel = GlobalConfig.gmodel_master;
	super(id, boundary, "Const_Flux");
	_u = fstate.vel.x;
	_v = fstate.vel.y;
	// [TODO]: Kyle, think about z component.
	_p = fstate.gas.p;
	_rho = fstate.gas.rho;
	_e = gmodel.internal_energy(fstate.gas);
	_nsp = gmodel.n_species;
	_massf.length = _nsp;
	for (int _isp=0; _isp < _nsp; _isp++) {
	    _massf[_isp] = fstate.gas.massf[_isp];
	}
	this.fstate = fstate.dup();
    }

    override string toString() const 
    {
	return "BFE_ConstFlux";
    }

    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
	throw new Error("BFE_ConstFlux.apply_unstructured_grid() not yet implemented");
    }
    
    override void apply_structured_grid(double t, int gtl, int ftl)
    {
	FVInterface IFace;
	size_t i, j, k;
	double _u_rel, _v_rel;
	SBlock blk = cast(SBlock) this.blk;
	assert(blk !is null, "Oops, this should be an SBlock object.");

        switch(which_boundary){
	case Face.west:
	    i = blk.imin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    // Flux equations
		    IFace = blk.get_cell(i,j,k).iface[Face.west];
		    // for a moving grid we need vel relative to the interface
		    _u_rel = _u - IFace.gvel.x;
		    _v_rel = _v - IFace.gvel.y;
		    IFace.F.mass = _rho * ( _u_rel*IFace.n.x + _v_rel*IFace.n.y );
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
		    IFace.F.momentum.refx = _p * IFace.n.x + _u*IFace.F.mass;
		    IFace.F.momentum.refy = _p * IFace.n.y + _v*IFace.F.mass;
		    IFace.F.momentum.refz = 0.0;
		    // [TODO]: Kyle, think about z component.
		    IFace.F.total_energy = IFace.F.mass * (_e + 0.5*(_u*_u+_v*_v)) + _p*(_u*IFace.n.x+_v*IFace.n.y);
		    for ( int _isp = 0; _isp < _nsp; _isp++ ){
			IFace.F.massf[_isp] = IFace.F.mass * _massf[_isp];
		    }
		    // [TODO]: Kyle, separate energy modes for multi-species simulations.
		} // end j loop
	    } // end k loop
	    break;
	default:
	    throw new Error("Const_Flux only implemented for WEST gas face.");
	}
    }
}

class BFE_EnergyBalanceThermionic : BoundaryFluxEffect {
public:
	// Function inputs from Eilmer4 .lua simulation input
    double emissivity;					// Input emissivity, 0<e<=1.0. Assumed black body radiation out from wall
    double Ar;							// Richardson constant, material-dependent
    double phi;							// Work function, material dependent. Input units in eV, 
    									// this gets converted to Joules by multiplying by Elementary charge, Qe
    int ThermionicEmissionActive;       // Whether or not Thermionic Emission is active. Default is 'on'

    // Solver iteration counts
    int Twall_iterations;				// Iterations for primary Twall calculations. Default = 200
    int Twall_subiterations;			// Iterations for newton method when ThermionicEmissionActive==1. Default = 20

	// Constants used in analysis
    double SB_sigma = 5.670373e-8;   	// Stefan-Boltzmann constant. 	Units: W/(m^2 K^4)
    double kb = 1.38064852e-23;			// Boltzmann constant. 			Units: (m^2 kg)/(s^2 K^1)
    double Qe = 1.60217662e-19; 		// Elementary charge. 			Units: C
 
    this(int id, int boundary, double emissivity, double Ar, double phi, int ThermionicEmissionActive,
    	int Twall_iterations, int Twall_subiterations)
    {
	super(id, boundary, "EnergyBalanceThermionic");
	this.emissivity = emissivity;
	this.Ar = Ar;
	this.phi = phi*Qe;  // Convert phi from input 'eV' to 'J'
	this.ThermionicEmissionActive = ThermionicEmissionActive;
	this.Twall_iterations = Twall_iterations;
	this.Twall_subiterations = Twall_subiterations;
    }

    override string toString() const 
    {
		return "BFE_EnergyBalanceThermionic(ThermionicEmissionActive=" ~ to!string(ThermionicEmissionActive) ~ 
		", Work Function =" ~ to!string(phi/Qe) ~
		"eV , emissivity=" ~ to!string(emissivity) ~ 
		", Richardson Constant=" ~ to!string(Ar) ~
		")";
    }

    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
	throw new Error("BFE_EnergyBalanceThermionic.apply_unstructured_grid() not yet implemented");
    }
    
    override void apply_structured_grid(double t, int gtl, int ftl)
    {
    	SBlock blk = cast(SBlock) this.blk;
		assert(blk !is null, "Oops, this should be an SBlock object.");
		if ( emissivity <= 0.0 || emissivity > 1.0 ) {
			// Check if emissivity value is valid
			throw new Error("emissivity should be 0.0<e<=1.0\n");
		} else if ( Ar == 0.0){
			throw new Error("Ar should be set!\n");			
		} else if ( phi == 0.0){
			throw new Error("phi should be set!\n");
		} else if (blk.myConfig.turbulence_model != TurbulenceModel.none) {
			throw new Error("WallBC_ThermionicEmission only implemented for laminar flow\n");
		} else {

		FVInterface IFace;
		size_t i, j, k;
		FVCell cell;
	
    	switch(which_boundary){
			case Face.south:
			    j = blk.jmin;
			    for (k = blk.kmin; k <= blk.kmax; ++k) {
				for (i = blk.imin; i <= blk.imax; ++i) {
					// Set cell/face properties
		    		cell = blk.get_cell(i,j,k);
		    		IFace = cell.iface[Face.south];
		    		double dn = distance_between(cell.pos[0], IFace.pos);
                	solve_for_wall_temperature(cell, IFace, dn, false);
				} // end j loop
			    } // end k loop
			    break;
			default:
			    throw new Error("WallBC_ThermionicEmission only implemented for SOUTH face.");
		} // end switch which_boundary
		// Start by making bounds of SOUTH FACE work	
	    } // end apply_structured_grid()
	}

	void solve_for_wall_temperature(const FVCell cell, FVInterface IFace, double dn, bool outwardNormal)
    // Iteratively converge on wall temp
    {

	auto gmodel = blk.myConfig.gmodel; 
    double dT, dTdn, dTdnx, dTdny, dTdnz, qx, qy, qz, k_eff, k_lam_wall, h, q_cond, Twall;
	double f_relax = 0.05;
    double tolerance = 1.0e-3; 
	double Twall_prev = IFace.fs.gas.Ttr;
	double Twall_prev_backup = IFace.fs.gas.Ttr;
    double q_total, q_total_prev, q_total_prev_backup = 0.0;
    size_t Twall_iteration_count;
    int iteration_check, subiteration_check = 0;
    // IFace orientation
    double nx = IFace.n.x; double ny = IFace.n.y; double nz = IFace.n.z;
    // IFace properties.
	FlowGradients grad = IFace.grad;

	// Mass diffusion things
    double[] jx; // diffusive mass flux in x
    double[] jy; // diffusive mass flux in y
    double[] jz; // diffusive mass flux in z
	size_t n_species = gmodel.n_species;
	double viscous_factor = blk.myConfig.viscous_factor;

	// Mass diffusion
	jx = IFace.jx.dup();
	jy = IFace.jy.dup();
	jz = IFace.jz.dup();

	// Newton method to find Twall when thermionic active
	double f_rad, f_thermionic, f_drv, Twall_1, Twall_0;
	size_t Twall_subiteration_count;

	// Electron creation

	// Iteratively solve for the wall temperature
	for (Twall_iteration_count=0; Twall_iteration_count <= Twall_iterations; ++Twall_iteration_count) {

		// Update the thermodynamic and transport properties at IFace
	    IFace.fs.gas.Ttr = Twall_prev;
	    IFace.fs.gas.p = cell.fs.gas.p;
	    gmodel.update_thermo_from_pT(IFace.fs.gas);
	    gmodel.update_trans_coeffs(IFace.fs.gas);

	     //Determine convective heat flux at current iteration with Twall_prev    
	    dT = (cell.fs.gas.Ttr - IFace.fs.gas.Ttr); // Positive is heat into wall
	    if (dT < 0.0){
	    	// Catch in case the iteration goes negative (radiation would become imaginary)
			IFace.fs.gas.Ttr = 0.9*cell.fs.gas.Ttr;
			Twall_prev = IFace.fs.gas.Ttr;
			dT = (cell.fs.gas.Ttr - IFace.fs.gas.Ttr);
	    }

	    // Calculate thermal conductivity 
    	k_eff = viscous_factor * (IFace.fs.gas.k + IFace.fs.k_t);
    	// Temperature gradient of cell centred value, to the present wall temperature
    	dTdn = dT / dn;
    	// Break up into components to remove mass diffusion Cartesian components
    	dTdnx = dTdn*nx; dTdny = dTdn*ny; dTdnz = dTdn*nz; 
    	qx = k_eff*dTdnx; qy = k_eff*dTdny; qz = k_eff*dTdnz; 

		// Apply molecular diffusion for the present laminar flow
		foreach (isp; 0 .. n_species) {
		    jx[isp] *= viscous_factor;
		    jy[isp] *= viscous_factor;
		    jz[isp] *= viscous_factor;
		}
		// Subtract mass diffusion contribution to heat transfer
		foreach (isp; 0 .. n_species) {
			h = gmodel.enthalpy(IFace.fs.gas, to!int(isp));
			qx -= jx[isp] * h;
			qy -= jy[isp] * h;
			qz -= jz[isp] * h;
		}

		// Don't nest pow functions. It is very slow....
		q_total = pow(qx*qx + qy*qy + qz*qz, 0.5);

		if (outwardNormal) {
				q_total = -q_total;
		}

	    // If we reach our tolerance value, we break. Tolerance based off heat loads
	    // as it is very susceptible to minor changes in Twall
	    if ( fabs((q_total - q_total_prev)/q_total) < tolerance) {
			//printf("Convergence reached\n");
		 //   printf("Twall = %.4g\n", Twall);
		    //printf("Calculated q_rad out = %.4g\n", pow(Twall,4)*emissivity*SB_sigma);
		 //   printf("cell.fs.gas.Ttr = %.4g\n", cell.fs.gas.Ttr);
		 //   printf("q_total = %.4g\n", q_total);
		 //   printf("q_total_prev = %.4g\n", q_total_prev);
		 //   printf("dTdn = %.4g\n", dTdn);
		 //   printf("k = %.4g\n", k_eff);
		 //   printf("\n");
		    // Update your wall temp with your final value
    	    IFace.fs.gas.Ttr = Twall;
		    IFace.fs.gas.p = cell.fs.gas.p;
		    gmodel.update_thermo_from_pT(IFace.fs.gas);
		    gmodel.update_trans_coeffs(IFace.fs.gas);

		    // Update flow gradients file
		    // The WallBC_ThermionicEmission is a copy of the adiabatic wall BC
		    // Hence, unless we alter these parameters, the exported heat transfer 
		    // will equal zero. It is assumed that there is no heat transfer
		    // between boundary cells (ie no conjugate heat transfer)
		    grad.Ttr[0] = dTdnx;
		    grad.Ttr[1] = dTdny;
		    grad.Ttr[2] = dTdnz;
		    iteration_check = 1;
	    	break;
		}

    	// What wall temperature would reach radiative equilibrium at current q_total
    	if (ThermionicEmissionActive == 0){
    		Twall = pow( q_total / ( emissivity * SB_sigma ), 0.25 );
    	} else {
    		Twall_0 = Twall_prev;
    		for (Twall_subiteration_count=0; Twall_subiteration_count < Twall_subiterations; ++Twall_subiteration_count) {
	    		f_rad = emissivity*SB_sigma*Twall_0*Twall_0*Twall_0*Twall_0;
	    		f_thermionic = phi/Qe*Ar*Twall_0*Twall_0*exp(-phi/(kb*Twall_0));
	    		f_drv = f_rad*4/Twall_0 + Ar*exp(-phi/(kb*Twall_0))/(Qe*kb)*(phi*(phi + 2*kb*Twall_0) + 2*kb*Twall_0*(phi + 3*kb*Twall_0));
	    		Twall_1 = Twall_0 - (f_rad + f_thermionic - q_total)/f_drv;
	    		if (fabs((Twall_1 - Twall_0))/Twall_0 <= 0.001) {
	    			subiteration_check = 1;
	    			break;
	    		}
	    		Twall_0 = Twall_1;
    		}
    		Twall = Twall_1;
    		//printf("\n");
    	}
	    // Determine new guess as a weighted average so we don't take too much of a step.
        Twall = f_relax * Twall + ( 1.0 - f_relax ) * Twall_prev;
    	Twall_prev_backup = Twall_prev;
    	Twall_prev = Twall;
    	q_total_prev_backup = q_total_prev;
    	q_total_prev = q_total;
	}
	if (iteration_check == 0){
		printf("Iteration's didn't converge\n");
        printf("Increase Twall_iterations from default 200 value\n");
        printf("\n");

	} else if (subiteration_check == 0 && ThermionicEmissionActive == 1) {
		printf("Subiteration's didn't converge\n");
        printf("Increase Twall_subiterations from default 20 value\n");
        printf("\n");
	}
return;
    } // end solve_for_wall_temperature()
} // end class BIE_EnergyBalanceThermionic

