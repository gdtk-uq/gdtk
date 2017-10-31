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
	newBFE = new BFE_EnergyBalanceThermionic(blk_id, boundary, emissivity, Ar, phi, ThermionicEmissionActive);
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
    FlowState fstate;

	// Function inputs from Eilmer4 .lua simulation input
    double emissivity;					// Input emissivity, 0<e<=1.0. Assumed black body radiation out from wall
    double Ar;							// Richardson constant, material-dependent
    double phi;							// Work function, material dependent. Input units in eV, 
    									// this get's converted to Joules by multiplying by Elementary charge, Qe
    int ThermionicEmissionActive;       // Whether or not Thermionic Emmision is active. Default is 'on'

	// Constants used in analysis
    double SB_sigma_SI = 5.670373e-8;   // Stefan-Boltzmann constant. 	Units: W/(m^2 K^4)
    double kb = 1.38064852e-23;			// Boltzmann constant. 			Units: (m^2 kg)/(s^2 K^1)
    double Qe = 1.60217662e-19; 		// Elementary charge. 			Units: C
 
    this(int id, int boundary, double emissivity, double Ar, double phi, int ThermionicEmissionActive)
    {
	super(id, boundary, "EnergyBalanceThermionic");
	this.emissivity = emissivity;
	this.Ar = Ar;
	this.phi = phi*Qe;  // Convert phi from input 'eV' to 'J'
	this.ThermionicEmissionActive = ThermionicEmissionActive;
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
		if ( emissivity < 0.0 || emissivity > 1.0 ) {
			// Check if emissivity value is valid
			throw new Error("emissivity should be 0.0<=e<=1.0");
		} else if ( Ar == 0.0){
			throw new Error("Ar should be set!");			
		} else if ( phi == 0.0){
			throw new Error("Phi should be set!");			
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
		    		FlowState fs = IFace.fs;
		    		double dn = distance_between(cell.pos[0], IFace.pos);
                	solve_for_wall_temperature(cell, IFace, dn, false);
				    // [TODO]: Kyle, separate energy modes for multi-species simulations.
				} // end j loop
			    } // end k loop
			    break;
			default:
			    throw new Error("Const_Flux only implemented for WEST gas face.");
		} // end switch which_boundary
		// Start by making bounds of SOUTH FACE work	
	    } // end apply()
	}

	void solve_for_wall_temperature(const FVCell cell, FVInterface IFace, double dn, bool outwardNormal)
    // Modify the wall function method from BIE_WallFunction to 
    // iteratively converge on wall temp
    {

	auto gmodel = blk.myConfig.gmodel; 
	double f_relax = 0.01;

    double dT, dTdy, k_lam_wall;

    double tolerance = 1.0e-4; 
    size_t counter = 0;
    size_t Twall_counter = 0;
    size_t max_iterations = 100;
    size_t Twall_max_iterations = 10;
    double q_total, q_total_prev = 0.0;
    double Twall;

	double Twall_prev = IFace.fs.gas.Ttr;
    double q_cond; 

    // Newtons method stuff
	double Tw_1;
	double Tw_0;
	double funcTw, dTw_funcTw;

	// Iteratively solve for the wall temperature
	while ( counter < max_iterations) {

		// Update the thermodynamic and transport properties at IFace
	    IFace.fs.gas.Ttr = Twall_prev;
	    IFace.fs.gas.p = cell.fs.gas.p;
	    gmodel.update_thermo_from_pT(IFace.fs.gas);
	    gmodel.update_trans_coeffs(IFace.fs.gas);

	     //Determine convective heat flux at current iteration with Twall_prev    
	    dT = (cell.fs.gas.Ttr - IFace.fs.gas.Ttr); // Positive is heat into wall
	    if (dT < 0.0){
			IFace.fs.gas.Ttr = 0.9*cell.fs.gas.Ttr;
			Twall_prev = IFace.fs.gas.Ttr;
			dT = (cell.fs.gas.Ttr - IFace.fs.gas.Ttr);
	    }

    	dTdy = dT / dn;
    	k_lam_wall = IFace.fs.gas.k;
    	q_cond = k_lam_wall * dTdy;

		if (outwardNormal) {
				q_cond = -q_cond;
		}

	 	// Get total q into wall at current Twall
	 	// To do: Mass diffusion heat transfer
    	q_total = (q_cond);

	    // If we reach our tolerance value, we break. Tolerance based off heat loads
	    // as it is apparently very surceptable to minor changes in Twall
	    if ( fabs((q_total - q_total_prev)/q_total) < tolerance) {
			printf("Convergence reached?\n");
		    printf("Twall = %.4g\n", Twall);
		    printf("Calculated q_rad out = %.4g\n", pow(Twall,4)*emissivity*SB_sigma_SI);
		    printf("cell.fs.gas.Ttr = %.4g\n", cell.fs.gas.Ttr);
		    printf("q_total = %.4g\n", q_total);
		    printf("q_total_prev = %.4g\n", q_total_prev);
		    printf("dTdy = %.4g\n", dTdy);
		    printf("k = %.4g\n", k_lam_wall);
		    printf("\n");
		    // Update your wall temp with your final value
    	    IFace.fs.gas.Ttr = Twall;
		    IFace.fs.gas.p = cell.fs.gas.p;
		    gmodel.update_thermo_from_pT(IFace.fs.gas);
		    gmodel.update_trans_coeffs(IFace.fs.gas);
	    	break;
		}
	
    	// What wall temperature would reach radiative equilibrium at current q_total
    	Twall = pow( q_total / ( emissivity * SB_sigma_SI ), 0.25 );

	    // Determine new guess as a weighted average so we don't take too much of a step.
        Twall = f_relax * Twall + ( 1.0 - f_relax ) * Twall_prev;

    	Twall_prev = Twall;
    	q_total_prev = q_total;

	    // Add to counter
	    counter++;
	}
return;
    } // end solve_for_wall_temperature()
} // end class BIE_EnergyBalanceThermionic
