/++ 
 Quasi-One Dimensional Adjoint Solver (QUODAS)

 Flow Solver 
 + Quasi-1D Euler equations.
 + Finite-volume, cell-centered formulation.
 + MUSCL-type reconstruction with Limiting via a Modified Van Albada method.
 + Van Leer or AUSMDV flux calculators.

 Adjoint Solver
 + Finite-difference sensitivities.
 + Adjoint system formulated in terms of primitive, and conservative variables.

 Optimiser
 + Conjugate gradient method.
 
 author: Kyle Damm, July 2017
 location: South Korea
++/

import std.stdio;
import std.math;
import std.string;
import std.file;
import std.conv;
import std.algorithm;
import std.exception;

//-----------------------------------------------------
// Set global simulations parameters
//-----------------------------------------------------

immutable size_t np = 3;                      // number of primitive variables
immutable size_t nd = 3;                      // number of design variables
immutable double length = 10.0;               // length of nozzle in metres
immutable size_t ncells = 300;                // number of cells
immutable size_t nghost = 4;                  // number of ghost cells
immutable size_t ninterfaces = ncells+1;      // number of interfaces
immutable string data_dir = "dat-files/";      // input/output files directory
immutable double eps0 = 1.0e-10;              // adjoint finite difference perturbation
immutable double eps1 = 1.0e-07;              // dRdD finite difference perturbation 
immutable double eps2 = 1.0e-06;              // brute-force finite difference perturbation

//-----------------------------------------------------
// Main Body
//-----------------------------------------------------

void main() {
    //-----------------------------------------------------
    // Select Simulation type
     //-----------------------------------------------------
    string simulation_type = "nozzle"; // options: Sod's Shocktube, Nozzle

    //-----------------------------------------------------
    // Set simulations parameters
    //-----------------------------------------------------
    string solver = "optimisation"; // options: simulation, optimisation, verification
   
    // flow solver
    double dt = 1.0e-05;                                // flow solver time step, s
    size_t max_step = 100000000;                        // maximum number of flow solver steps
    double final_time;                                  // final simulation time, s
    if (simulation_type == "sod") final_time = 6.0e-04; 
    else if (simulation_type == "nozzle") final_time = 2.0; 
    string flux_calc = "ausmdv"; // options: van_leer, ausmdv
    size_t interpolation_order = 1;
    size_t outflow_bc_order = 1;
    double outflow_bc_back_pressure_factor = 0.5;
    
    // adjoint solver
    string adjoint_form = "primitive";           // options: conservative form, primtive form
    double tol = 1.0e-04;                        // optimisation tolerance

    // geometry parameters for nozzle

    // nozzle shock removal/addition: with bp factor 0.5
    // with shock: b = 0.5, c = 1.0, d = 3.8 
    // without shock: b = 0.347, c = 0.8, d = 4.0

    double b = 0.347;
    double c = 0.8;
    double d = 4.0;
    double yo = 0.105;
    double scale = 1.0;

    // gas properties
    double gamma = 1.4; // ideal air
    
    //-----------------------------------------------------
    // Create directory for simulation result files
    //-----------------------------------------------------
    if (exists(data_dir) && isDir(data_dir)) {} //do nothing;
    else {
	try {
	    mkdirRecurse(data_dir);
	} catch (FileException e) {
	    string msg = text("Failed to ensure directory is present: ", data_dir);
	    throw new Error(msg);
	}
    }

    //-----------------------------------------------------
    // Construct Mesh
    //-----------------------------------------------------

    fvcell[ncells+nghost] global_cells;
    double dx = length/ncells; // cell width

    // construct ghost cells ---------------------

    global_cells[0] = new fvcell(10000001, -1.5*dx, dx);
    global_cells[1] = new fvcell(10000002, -0.5*dx, dx);
    global_cells[ncells+2] = new fvcell(10000003, length+0.5*dx, dx);
    global_cells[ncells+3] = new fvcell(10000004, length+1.5*dx, dx);

    // construct inner domain cells --------------------

    double temp_pos = 0.5*dx;
    foreach(i; 0 .. ncells) {
	global_cells[i+2] = new fvcell(i, temp_pos, dx);
	temp_pos += dx ;
    }
    
    // construct interfaces ----------------------
    
    fvinterface[ninterfaces] global_interfaces;
    fvinterface[ninterfaces] perturbed_interfaces; // for implicit solver
    temp_pos = 0.0;
    foreach(i; 0 .. ninterfaces) {
	global_interfaces[i] = new fvinterface(i, temp_pos);
	perturbed_interfaces[i] = new fvinterface(i, temp_pos);
	temp_pos += dx;
    }

    // compute interface area, and set cell volume
    compute_geometry(simulation_type, global_interfaces, perturbed_interfaces, global_cells, b, c, d, yo, scale);

    //-----------------------------------------------------
    // Initial Condiions
    //-----------------------------------------------------
    set_initial_and_inflow_conditions(simulation_type, global_cells, gamma);
    
    // ---------------------------------------------------------------------------------- //
    // Now we're ready to perform some finite volume calculations on the cells
    // ---------------------------------------------------------------------------------- //
    // set some internal parameters
    double sim_time = 0.0;                              // current simulation time, s
    size_t step = 0;
    size_t opt_iter = 0;
    double dLdB_old = 0.0; double dLdC_old = 0.0; double dLdD_old = 0.0;
    double dLdB = 1e10; double dLdC = 1e10; double dLdD = 1e10;
    double Jref = 0.0;
    double[np*ncells] psi;
    
    // optimisation routine arrays
    double[nd] pk;
    double[nd] pk_old;
    double[] tol_hist;
    
    if (solver == "simulation") {
	flow_solver(global_cells, global_interfaces, perturbed_interfaces, dt, sim_time, final_time, step, max_step,
		    flux_calc, simulation_type, interpolation_order, outflow_bc_order, outflow_bc_back_pressure_factor, dx, gamma);
	// write out geometry, and flow-field
	write_out_solution(global_cells, global_interfaces, tol_hist, psi, "target_contour.dat", "target_output.dat");
    }
    else if(solver == "verification") {

	// make sure the target pressure file exists
	if (exists(data_dir~"/target_output.dat")) {} //do nothing;
	else {
		string msg = text("target pressure file (naming format: target_output.dat) does not exist in user specified data directory: ", data_dir);
		throw new Error(msg);
	}
	// read in target pressure
	double[ncells] p_target;
	read_target_pressure_file(p_target);
	
	writeln("-- computing gradient via adjoint method");
	flow_solver(global_cells, global_interfaces, perturbed_interfaces, dt, sim_time, final_time, step, max_step,
		    flux_calc, simulation_type, interpolation_order, outflow_bc_order, outflow_bc_back_pressure_factor, dx, gamma);
	adjoint_solver(global_cells, global_interfaces, perturbed_interfaces, dt, sim_time, final_time, step, max_step,
		       flux_calc, simulation_type, interpolation_order, outflow_bc_order, outflow_bc_back_pressure_factor,
		       psi, dLdB, dLdC, dLdD, dx, gamma, b, c, d, yo, scale, solver, adjoint_form, p_target);
	double adjoint_dLdB = dLdB; double adjoint_dLdC = dLdC; double adjoint_dLdD = dLdD;

	writeln("-- computing gradient via finite-difference method");
	double Jminus; double Jplus;

	// perturb b
	writeln("-- -- perturb b");
	double fd_dLdB;
	double b_orig = b;
	b = b_orig + (b_orig*eps2 + eps2);
	compute_geometry(simulation_type, global_interfaces, perturbed_interfaces, global_cells, b, c, d, yo, scale);
	flow_solver(global_cells, global_interfaces, perturbed_interfaces, dt, sim_time, final_time, step, max_step,
		    flux_calc, simulation_type, interpolation_order, outflow_bc_order, outflow_bc_back_pressure_factor, dx, gamma);
	Jplus = objective_function_evaluation(global_cells, p_target);
	b = b_orig - (b_orig*eps2 + eps2);
	compute_geometry(simulation_type, global_interfaces, perturbed_interfaces, global_cells, b, c, d, yo, scale);
	flow_solver(global_cells, global_interfaces, perturbed_interfaces, dt, sim_time, final_time, step, max_step,
		    flux_calc, simulation_type, interpolation_order, outflow_bc_order, outflow_bc_back_pressure_factor, dx, gamma);
	Jminus = objective_function_evaluation(global_cells, p_target);
	fd_dLdB = (Jplus - Jminus)/(2.0 * (b_orig*eps2 + eps2));
	b = b_orig;

	// perturb c
	writeln("-- -- perturb c");
	double fd_dLdC;
	double c_orig = c;
	c = c_orig + (c_orig*eps2 + eps2);
	compute_geometry(simulation_type, global_interfaces, perturbed_interfaces, global_cells, b, c, d, yo, scale);
	flow_solver(global_cells, global_interfaces, perturbed_interfaces, dt, sim_time, final_time, step, max_step,
		    flux_calc, simulation_type, interpolation_order, outflow_bc_order, outflow_bc_back_pressure_factor, dx, gamma);
	Jplus = objective_function_evaluation(global_cells, p_target);
        c = c_orig - (c_orig*eps2 + eps2);
	compute_geometry(simulation_type, global_interfaces, perturbed_interfaces, global_cells, b, c, d, yo, scale);
	flow_solver(global_cells, global_interfaces, perturbed_interfaces, dt, sim_time, final_time, step, max_step,
		    flux_calc, simulation_type, interpolation_order, outflow_bc_order, outflow_bc_back_pressure_factor, dx, gamma);
	Jminus = objective_function_evaluation(global_cells, p_target);
	fd_dLdC = (Jplus - Jminus)/(2.0 * (c_orig*eps2 + eps2));
	c = c_orig;

	// perturb d
	writeln("-- -- perturb d");
	double fd_dLdD;
	double d_orig = d;
	d = d_orig + (d_orig*eps2 + eps2);
	compute_geometry(simulation_type, global_interfaces, perturbed_interfaces, global_cells, b, c, d, yo, scale);
	flow_solver(global_cells, global_interfaces, perturbed_interfaces, dt, sim_time, final_time, step, max_step,
		    flux_calc, simulation_type, interpolation_order, outflow_bc_order, outflow_bc_back_pressure_factor, dx, gamma);
	Jplus = objective_function_evaluation(global_cells, p_target);
        d = d_orig - (d_orig*eps2 + eps2);
	compute_geometry(simulation_type, global_interfaces, perturbed_interfaces, global_cells, b, c, d, yo, scale);
	flow_solver(global_cells, global_interfaces, perturbed_interfaces, dt, sim_time, final_time, step, max_step,
		    flux_calc, simulation_type, interpolation_order, outflow_bc_order, outflow_bc_back_pressure_factor, dx, gamma);
	Jminus = objective_function_evaluation(global_cells, p_target);
	fd_dLdD = (Jplus - Jminus)/(2.0 * (d_orig*eps2 + eps2));
	d = d_orig;

	writeln("--------------------");
	writeln("Verification Results:");
	writeln("--------------------");
	writeln("fd_b = ", fd_dLdB);
	writeln("adjoint_b = ", adjoint_dLdB);
	writeln("% error b = ", (adjoint_dLdB-fd_dLdB)/adjoint_dLdB * 100.0);
	writeln("fd_c = ", fd_dLdC);
	writeln("adjoint_c = ", adjoint_dLdC);
	writeln("% error c = ", (adjoint_dLdC-fd_dLdC)/adjoint_dLdC * 100.0);
	writeln("fd_d = ", fd_dLdD);
	writeln("adjoint_d = ", adjoint_dLdD);
	writeln("% error d = ", (adjoint_dLdD-fd_dLdD)/adjoint_dLdD * 100.0);
    }
    else if (solver == "optimisation") {

	// make sure the target pressure file exists
	if (exists(data_dir~"/target_output.dat")) {} //do nothing;
	else {
	    string msg = text("target pressure file (naming format: target_output.dat) does not exist in user specified data directory: ", data_dir);
	    throw new Error(msg);
	}
	// read in target pressure
	double[ncells] p_target;
	read_target_pressure_file(p_target);

	writeln("BASELINE GEOMETRY (b, c, d): ", b, ", ", c, ", ", d);
	while (opt_iter < 1e6) {
	    //-----------------------------------------------------
	    // conjugate gradient update
	    //-----------------------------------------------------
	    	    	    
	    // 1. compute flow solution, gradients, and objective function for current geometry
	    flow_solver(global_cells, global_interfaces, perturbed_interfaces, dt, sim_time, final_time, step, max_step,
			flux_calc, simulation_type, interpolation_order, outflow_bc_order, outflow_bc_back_pressure_factor, dx, gamma);
	    adjoint_solver(global_cells, global_interfaces, perturbed_interfaces, dt, sim_time, final_time, step,
			   max_step, flux_calc, simulation_type, interpolation_order, outflow_bc_order, outflow_bc_back_pressure_factor,
			   psi, dLdB, dLdC, dLdD, dx, gamma, b, c, d, yo, scale, solver, adjoint_form, p_target);

	    if (opt_iter == 0) {
		// write out optimisaed geometry, and flow-field
		write_out_solution(global_cells, global_interfaces, tol_hist, psi, "initial_contour.dat", "initial_output.dat");
	    }
	    double J = objective_function_evaluation(global_cells, p_target);


	    // if on the first step, then set the reference J
	    if (opt_iter == 0) Jref = J;

	    tol_hist ~= J/Jref;
	    writeln("--------------------------------------------------------------------------");
	    writeln(opt_iter, ". J/Jref = ", J/Jref, ", b = ", b, ", c = ", c, ", d = ", d);
	    writeln("gradients: ", ", dLdB = ", dLdB,  ", dLdC = ", dLdC,  ", dLdD = ", dLdD);
	    writeln("--------------------------------------------------------------------------");
	    
	    // if the gradients and cost function are small enough then we are at a minimum
	    if (J/Jref < tol) break;
	    
	    // 2. compute search direction
	    double beta;
	    if (opt_iter == 0) { // first step uses the steepest direction method
		double norm = sqrt(dLdB*dLdB + dLdC*dLdC + dLdD*dLdD);
		pk[0] = -dLdB/norm;
		pk[1] = -dLdC/norm;
		pk[2] = -dLdD/norm;
	    }
	    else if( ( opt_iter % 3 ) == 0) { // we perform restarts for robustness
		double norm = sqrt(dLdB*dLdB + dLdC*dLdC + dLdD*dLdD);
		pk[0] = -dLdB/norm;
		pk[1] = -dLdC/norm;
		pk[2] = -dLdD/norm;
	    }
	    else { // else let's use the previous gradient knowledge to improve our step
		double norm = sqrt(dLdB*dLdB + dLdC*dLdC + dLdD*dLdD);
		double beta_numer = dLdB*dLdB + dLdC*dLdC + dLdD*dLdD;
		double beta_denom = dLdB_old*dLdB_old + dLdC_old*dLdC_old + dLdD_old*dLdD_old;
		beta = beta_numer/beta_denom;
		pk[0] = -dLdB/norm + beta*pk_old[0]; 
		pk[1] = -dLdC/norm + beta*pk_old[1]; 
		pk[2] = -dLdD/norm + beta*pk_old[2]; 
	    }
	    
	    // 3. perform line search to find step size ak (this is the expenseive stage)
	    double rho = 0.5;
	    double mu = 1.0e-3; // this has a strong effect on convergence
	    double ak = 1.0;
	    double bb = b;
	    double cc = c;
	    double dd = d;
	    double Jk = J;
	    
	    // update geometry
	    compute_geometry(simulation_type, global_interfaces, perturbed_interfaces, global_cells, b, c, d, yo, scale);
	    flow_solver(global_cells, global_interfaces, perturbed_interfaces, dt, sim_time, final_time, step, max_step,
			flux_calc, simulation_type, interpolation_order, outflow_bc_order, outflow_bc_back_pressure_factor, dx, gamma);
	    double Jk1 = objective_function_evaluation(global_cells, p_target);

	    bool are_design_vars_negative = true; // assume negative to begin with
	    while (Jk1 >= Jk + ak*mu*(dLdB*pk[0]+dLdC*pk[1]+dLdD*pk[2])) {

	      // add in constraint to ensure positive design variables
	      if (are_design_vars_negative) {
		while (are_design_vars_negative) {
		  ak = ak*rho;
		  b = bb + ak*pk[0];
		  c = cc + ak*pk[1];
		  d = dd + ak*pk[2];

		  compute_geometry(simulation_type, global_interfaces, perturbed_interfaces, global_cells, b, c, d, yo, scale);

		  size_t num_negative_areas = 0;
		  foreach (iface; global_interfaces) {
		    if (iface.area < 0.0) num_negative_areas += 1;
		  }
		  if (b > 0.0 && num_negative_areas == 0 &&  global_interfaces[0].area <  global_interfaces[$-1].area) are_design_vars_negative  = false; 
		}
	      }
	      else {
		ak = ak*rho;
		b = bb + ak*pk[0];
		c = cc + ak*pk[1];
		d = dd + ak*pk[2];
	      }
	      
	      // update geometry
	      compute_geometry(simulation_type, global_interfaces, perturbed_interfaces, global_cells, b, c, d, yo, scale);
	      flow_solver(global_cells, global_interfaces, perturbed_interfaces, dt, sim_time, final_time, step, max_step,
			  flux_calc, simulation_type, interpolation_order, outflow_bc_order, outflow_bc_back_pressure_factor, dx, gamma);
	      Jk1 = objective_function_evaluation(global_cells, p_target);
	    }

	    // update geometry
	    compute_geometry(simulation_type, global_interfaces, perturbed_interfaces, global_cells, b, c, d, yo, scale);
	    
	    dLdB_old = dLdB; 
	    dLdC_old = dLdC;
	    dLdD_old = dLdD; 
	    pk_old[0] = pk[0];
	    pk_old[1] = pk[1];
	    pk_old[2] = pk[2];
	    opt_iter += 1;
	}
	writeln("FINAL GEOMETRY (b, c, d):", b, ", ", c, ", ", d);

	// write out optimisaed geometry, and flow-field
	write_out_solution(global_cells, global_interfaces, tol_hist, psi, "optimised_contour.dat", "optimised_output.dat");
    }
}

//-----------------------------------------------------
// Classes and Functions
//-----------------------------------------------------

void write_out_solution(fvcell[] global_cells, fvinterface[] global_interfaces, double[] tol_hist, double[] psi, string contourFilename, string solutionFilename) {
    foreach(i; 0 .. ncells) {
        fvcell cell = global_cells[i+2];
        auto writer = format("%f %f %f %f %f %f %f \n", cell.xpos, cell.rho, cell.p, cell.u, psi[i*np], psi[i*np+1], psi[i*np+2]);
        append(data_dir ~ solutionFilename, writer);
    }
    foreach(i; 0 .. ninterfaces) {
        fvinterface f = global_interfaces[i];
        auto writer = format("%f %f \n", f.xpos, f.area/2.0);
        append(data_dir ~ contourFilename, writer);
    }
    int count = 0;
    foreach(i; 0 .. tol_hist.length) {
        auto writer = format("%d %f \n", count, tol_hist[i]);
        append(data_dir ~ "tol_hist.dat", writer);
	count += 1;
    }
}

struct flowstate
{
    double p;     // cell pressure
    double rho;   // cell density
    double u;     // cell velocity
}

class fvcell {
public:
    size_t id;        // id number
    double xpos;      // position
    double dx;        // width
    double vol;       // volume
    double p;         // pressure
    double rho;       // density
    double u;         // velocity
    double ru;        // momentum
    double rE;        // energy
    
    this(size_t id, double xpos, double dx) {
	this.id = id;
	this.xpos = xpos;
	this.dx = dx;
    }
}

class fvinterface {
public:
    size_t id;           // id number
    double xpos;         // position
    double area;         // area
    double mass;         // mass flux
    double momentum;     // momentum flux
    double energy;       // energy flux
    flowstate left_fs;   // left flowstate
    flowstate right_fs;  // right flowstate
    
    this(size_t id, double xpos) {
	this.id = id;
	this.xpos = xpos;
    }
}

void read_target_pressure_file(double[] p_target) {
    // target pressure distribution saved in file target_output.dat
    auto file = File(data_dir ~ "/target_output.dat", "r");
    foreach(i; 0 .. ncells) {
	auto lineContent = file.readln().strip();
	auto tokens = lineContent.split();
	p_target[i] = to!double(tokens[2]);
    }
} // end read_target_pressure_file

double objective_function_evaluation(fvcell[] global_cells, double[] p_target) {
    // J(Q) = 0.5*integral[0->l] (p-p*)^2 dx
    double J = 0.0;
    foreach (i; 0..ncells) {
	J += 0.5*(global_cells[i+2].p - p_target[i])*(global_cells[i+2].p - p_target[i]);
    }
    return J;
} // end objective_function_evaluation

void compute_geometry(string simulation_type, fvinterface[] global_interfaces, fvinterface[] perturbed_interfaces,
		      fvcell[] global_cells, double b, double c, double d, double yo, double scale) {
    // compute interface area
    double inlet_area; double exit_area;
    if (simulation_type == "sod") {
	// shock tube geometry has no area variaton
	foreach(i; 0 .. ninterfaces) {
	    global_interfaces[i].area = 1.0;
	    perturbed_interfaces[i].area = 1.0;
	}
    }
    else if (simulation_type == "nozzle") {
	// symmetric nozzle is defined by hyperbolic function
	double a = yo - b*tanh(-d/scale);
	foreach(i; 0 .. ninterfaces) {
	    double height = a + b*tanh((c*global_interfaces[i].xpos -d)/scale);
	    global_interfaces[i].area =  2.0*height;
	    perturbed_interfaces[i].area = 2.0*height;
	}
    }
    // compute cell volume ---------------------------
    foreach(i; 0 .. ncells) {
	global_cells[i+2].vol = 0.5*global_cells[i+2].dx*(global_interfaces[i].area+global_interfaces[i+1].area);
    }
} // end compute_geometry

void set_initial_and_inflow_conditions(string simulation_type, fvcell[] global_cells, double gamma) {
    if (simulation_type == "sod" ) {
	// left fill state
	double p_init = 1.0e05; // Pa
	double rho_init = 1.0; // kg/m^3
	double u_init = 0.0; // m/s
	foreach(i; 0 .. to!int(0.5*ncells)) {
	    global_cells[i+2].p = p_init;
	    global_cells[i+2].rho = rho_init;
	    global_cells[i+2].u = u_init;
	}	
	// right fill state
	p_init = 1.0e04; // Pa
	rho_init = 0.125; // kg/m^3
	u_init = 0.0; // m/s
	foreach(i; to!int(0.5*ncells) .. ncells) {
	    global_cells[i+2].p = p_init;
	    global_cells[i+2].rho = rho_init;
	    global_cells[i+2].u = u_init;
	}
    }
    else if (simulation_type == "nozzle") {
	// set inflow conditions
	double Mach = 1.5;
	double p_inflow = 101.325e3;                          // Pa
	double rho_inflow = 1.0;                              // kg/m^3
	double a_inflow = sqrt((p_inflow*gamma)/rho_inflow);  // m/s
	double u_inflow = Mach * a_inflow;                    // m/s

	foreach(i; 0 .. ncells) {
	    global_cells[i+2].p = p_inflow;
	    global_cells[i+2].rho = rho_inflow;
	    global_cells[i+2].u = 0.6 * a_inflow;
	}

	// Nozzle inflow state ----------------------------------------------
	global_cells[0].p = p_inflow; global_cells[1].p = p_inflow;
	global_cells[0].rho = rho_inflow; global_cells[1].rho = rho_inflow;
	global_cells[0].u = u_inflow; global_cells[1].u = u_inflow;

	global_cells[ncells+2].p = global_cells[ncells+1].p; global_cells[ncells+3].p = global_cells[ncells+1].p;
	global_cells[ncells+2].rho =  global_cells[ncells+1].rho; global_cells[ncells+3].rho =  global_cells[ncells+1].rho;
	global_cells[ncells+2].u =  global_cells[ncells+1].u; global_cells[ncells+3].u =  global_cells[ncells+1].u;	
    }
}

void compute_flux(string flux_calc, fvcell[] global_cells, fvinterface[] global_interfaces, size_t interpolation_order, size_t ninterfaces, double gamma) {
    // Interpolate cell centered values to interface
    foreach(i; 0 .. ninterfaces) {
	if (interpolation_order == 1)
	    {first_order_interpolation(global_cells[i+1],global_cells[i+2],global_interfaces[i]);}
	else
	    {second_order_interpolation_with_van_albada_limiting(global_cells[i+1],global_cells[i],global_cells[i+2],global_cells[i+3],global_interfaces[i]);}
    }
    // Compute Interface Flux
    foreach(i; 0 .. ninterfaces) {
	if (flux_calc == "van_leer") van_leer_flux_calculator(global_interfaces[i], gamma);
	else if (flux_calc == "ausmdv") ausmdv_flux_calculator(global_interfaces[i], gamma);
    }
} // end compute_flux

void apply_post_flux_bcs(string simulation_type, fvinterface[] global_interfaces, fvcell[] global_cells, double gamma) {
    if (simulation_type == "nozzle") {
	// for the nozzle simulation we should overwrite the inflow interface flux for the supersonic bc
	fvcell cell = global_cells[0];
	double e = cell.p / (cell.rho*(gamma - 1.0));
	double ke = 0.5*cell.u*cell.u;
	global_interfaces[0].mass = cell.rho*cell.u;
	global_interfaces[0].momentum = cell.p+cell.rho*cell.u*cell.u;
	global_interfaces[0].energy = (cell.rho*e + cell.rho*ke +cell.p)*cell.u;
    }
} // end apply_post_flux_bcs

void apply_pre_flux_bcs(string simulation_type, size_t outflow_bc_order, double outflow_bc_back_pressure_factor, fvcell[] global_cells, double dx) {
    if (simulation_type == "sod") {
	// left boundary -- wall
	global_cells[0].rho = global_cells[2].rho;
	global_cells[1].rho = global_cells[2].rho;
	global_cells[0].u = -global_cells[2].u;
	global_cells[1].u = -global_cells[2].u;
	global_cells[0].p = global_cells[2].p;
	global_cells[1].p = global_cells[2].p;
	// right boundary -- wall
	global_cells[ncells+2].rho = global_cells[ncells+1].rho;
	global_cells[ncells+3].rho = global_cells[ncells+1].rho;
	global_cells[ncells+2].u = -global_cells[ncells+1].u;
	global_cells[ncells+3].u = -global_cells[ncells+1].u;
	global_cells[ncells+2].p = global_cells[ncells+1].p;
	global_cells[ncells+3].p = global_cells[ncells+1].p;
    }
    else if (simulation_type == "nozzle") {
	// left boundary -- inflow
	global_cells[0].rho = global_cells[0].rho;
	global_cells[1].rho = global_cells[1].rho;
	global_cells[0].u = global_cells[0].u;
	global_cells[1].u = global_cells[1].u;
	global_cells[0].p = global_cells[0].p;
	global_cells[1].p = global_cells[1].p;
	// right boundary - outflow
	if (outflow_bc_order == 1) { // first order
	    global_cells[ncells+2].rho = global_cells[ncells+1].rho;
	    global_cells[ncells+3].rho = global_cells[ncells+1].rho;
	    global_cells[ncells+2].u = global_cells[ncells+1].u;
	    global_cells[ncells+3].u = global_cells[ncells+1].u;
	    global_cells[ncells+2].p = outflow_bc_back_pressure_factor*global_cells[0].p;
	    global_cells[ncells+3].p = outflow_bc_back_pressure_factor*global_cells[0].p;

	}
	else { // second order
	    double drhodx; double dudx; double dpdx;
	    drhodx = (global_cells[ncells+1].rho - global_cells[ncells].rho)/dx;
	    dudx = (global_cells[ncells+1].u - global_cells[ncells].u)/dx;
	    dpdx = (global_cells[ncells+1].p - global_cells[ncells].p)/dx;
	    
	    global_cells[ncells+2].rho = global_cells[ncells+1].rho+drhodx*dx;
	    global_cells[ncells+3].rho = global_cells[ncells+1].rho+drhodx*2.0*dx;
	    global_cells[ncells+2].u = global_cells[ncells+1].u+dudx*dx;
	    global_cells[ncells+3].u = global_cells[ncells+1].u+dudx*2.0*dx;
	    global_cells[ncells+2].p = outflow_bc_back_pressure_factor*global_cells[0].p;
	    global_cells[ncells+3].p = outflow_bc_back_pressure_factor*global_cells[0].p;
	}		    
    }
} // apply_pre_flux_bcs

void first_order_interpolation(fvcell L0, fvcell R0, fvinterface f) {
    // copies data from cell center to interface left, and right flowstates.
    // left flow state
    f.left_fs.p = L0.p; f.left_fs.rho = L0.rho; f.left_fs.u = L0.u;
    // right flow state
    f.right_fs.p = R0.p; f.right_fs.rho = R0.rho; f.right_fs.u = R0.u;
} // end first_order_interpolation

void second_order_interpolation_with_van_albada_limiting(fvcell L0, fvcell L1, fvcell R0, fvcell R1, fvinterface f) {
    // Second order interpolation with Van albada limiting (as outlined in Johnston's thesis [1999]).
    double delta_minus; double delta_plus; double S; double eps = 1.0e-12; double k = 0.0;
    // left fs state
    // pressure
    delta_minus = (L0.p - L1.p)/(0.5*(L0.dx+L1.dx));
    delta_plus = (R0.p - L0.p)/(0.5*(R0.dx+L0.dx));
    S = (2.0*delta_plus*delta_minus+eps)/(delta_plus*delta_plus+delta_minus*delta_minus+eps);
    f.left_fs.p = L0.p + (L0.dx*S)/4.0 * ((1.0-S*k)*delta_minus+(1.0+S*k)*delta_plus);
    // density
    delta_minus = (L0.rho - L1.rho)/(0.5*(L0.dx+L1.dx));
    delta_plus = (R0.rho - L0.rho)/(0.5*(R0.dx+L0.dx));
    S = (2.0*delta_plus*delta_minus+eps)/(delta_plus*delta_plus+delta_minus*delta_minus+eps);
    f.left_fs.rho = L0.rho+ (L0.dx*S)/4.0 * ((1.0-S*k)*delta_minus+(1.0+S*k)*delta_plus);
    // velocity
    delta_minus = (L0.u - L1.u)/(0.5*(L0.dx+L1.dx));
    delta_plus = (R0.u - L0.u)/(0.5*(R0.dx+L0.dx));
    S = (2.0*delta_plus*delta_minus+eps)/(delta_plus*delta_plus+delta_minus*delta_minus+eps);
    f.left_fs.u = L0.u + (L0.dx*S)/4.0 * ((1.0-S*k)*delta_minus+(1.0+S*k)*delta_plus);
    // right flow state
    // pressure
    delta_minus = (R0.p - L0.p)/(0.5*(R0.dx+L0.dx));
    delta_plus = (R1.p - R0.p)/(0.5*(R0.dx+R1.dx));
    S = (2.0*delta_plus*delta_minus+eps)/(delta_plus*delta_plus+delta_minus*delta_minus+eps);
    f.right_fs.p = R0.p + (R0.dx*S)/4.0 * ((1.0-S*k)*delta_minus+(1.0+S*k)*delta_plus);
    // density
    delta_minus = (R0.rho - L0.rho)/(0.5*(R0.dx+L0.dx));
    delta_plus = (R1.rho - R0.rho)/(0.5*(R0.dx+R1.dx));
    S = (2.0*delta_plus*delta_minus+eps)/(delta_plus*delta_plus+delta_minus*delta_minus+eps);
    f.right_fs.rho = R0.rho + (R0.dx*S)/4.0 * ((1.0-S*k)*delta_minus+(1.0+S*k)*delta_plus);
    // velocity
    delta_minus = (R0.u - L0.u)/(0.5*(R0.dx+L0.dx));
    delta_plus = (R1.u - R0.u)/(0.5*(R0.dx+R1.dx));
    S = (2.0*delta_plus*delta_minus+eps)/(delta_plus*delta_plus+delta_minus*delta_minus+eps);
    f.right_fs.u = R0.u + (R0.dx*S)/4.0 * ((1.0-S*k)*delta_minus+(1.0+S*k)*delta_plus);
} // end second_order_interpolation_with_van_albada_limiting

void van_leer_flux_calculator(fvinterface f, double gamma) {
    // note lft == plus, rght == minus
    double F_lft; double M_lft; double a_lft;
    double F_rght; double M_rght; double a_rght;
    double lft_rho = f.left_fs.rho; double rght_rho = f.right_fs.rho;
    double lft_u = f.left_fs.u; double rght_u = f.right_fs.u;
    double lft_p = f.left_fs.p; double rght_p = f.right_fs.p;
    // left state
    a_lft = sqrt( (gamma*lft_p)/lft_rho );
    M_lft = lft_u/a_lft;
    // right state
    a_rght = sqrt( (gamma*rght_p)/rght_rho );
    M_rght = rght_u/a_rght;
    double M = 0.5*(M_lft+M_rght); // average Mach number
    if (M >= 1.0) {
	// mass flux
	f.mass = lft_rho * a_lft * M_lft;
	// momentum flux
	f.momentum = lft_rho * a_lft*a_lft*(M_lft*M_lft+1.0/gamma);
	// energy flux
	f.energy = lft_rho * a_lft*a_lft*a_lft * M_lft * (0.5*M_lft*M_lft + 1.0/(gamma-1.0));
    }
    else if (M <= -1.0) {
	// mass flux
	f.mass = rght_rho * a_rght * M_rght;
	// momentum flux
	f.momentum = rght_rho * a_rght*a_rght*(M_rght*M_rght+1.0/gamma);
	// energy flux
	f.energy = rght_rho * a_rght*a_rght*a_rght * M_rght * (0.5*M_rght*M_rght + 1.0/(gamma-1.0));
    }
    else { 
	// mass flux
	F_lft = 0.25*lft_rho*a_lft*pow( (1.0+M_lft), 2.0 );
	F_rght = -0.25*rght_rho*a_rght*pow( (1.0-M_rght), 2.0 );
	f.mass = F_lft+F_rght;
	// momentum flux
	F_lft = 0.25*lft_rho*a_lft*pow( (1.0+M_lft), 2.0 ) * ( 2.0*a_lft/gamma * ( (gamma-1.0)/2.0 * M_lft + 1.0) );
	F_rght = -0.25*rght_rho*a_rght*pow( (1.0-M_rght), 2.0 ) * ( 2.0*a_rght/gamma * ( (gamma-1.0)/2.0 * M_rght - 1.0) );
	f.momentum = F_lft+F_rght;
	// energy flux
	F_lft = 0.25*lft_rho*a_lft*pow( (1.0+M_lft), 2.0 ) *
	    ((2.0*a_lft*a_lft)/(gamma*gamma - 1.0) * pow( (gamma-1.0)/2.0 * M_lft + 1.0, 2.0 ));
	F_rght = -0.25*rght_rho*a_rght*pow( (1.0-M_rght), 2.0 ) *
	    ((2.0*a_rght*a_rght)/(gamma*gamma - 1.0) * pow( (gamma-1.0)/2.0 * M_rght - 1.0, 2.0 ));
	f.energy = F_lft+F_rght;
    }
} // end van_leer_flux_calculator

void ausmdv_flux_calculator(fvinterface f, double gamma) {
    // implementation taken from Eilmer4 flow solver (http://cfcfd.mechmining.uq.edu.au/eilmer/)
    double K_SWITCH = 10.0;
    double C_EFIX = 0.125;
    double rL, rR;
    double pL, pR;
    double uL, uR;
    double aL, aR;
    double HL, HR;
    double pLrL, pRrR;
    double ML, MR;
    double eL, eR;
    double keL, keR;
    double alphaL, alphaR, am;
    double pLplus, pRminus;
    double uLplus, uRminus;
    double duL, duR;
    double p_half, ru_half, ru2_half;
    double dp, s, ru2_AUSMV, ru2_AUSMD;
    int caseA, caseB;
    double d_ua;
    /*
     * Unpack the flow-state vectors for either side of the interface.
     * Store in work vectors, those quantities that will be neede later.
     */
    rL = f.left_fs.rho;
    pL = f.left_fs.p;
    pLrL = pL / rL;
    uL = f.left_fs.u;
    eL = pL/((gamma-1.0)*rL); 
    aL = sqrt((pL*gamma)/rL);
    keL = 0.5 * (uL * uL);
    HL = eL + pLrL + keL;

    rR = f.right_fs.rho;
    pR = f.right_fs.p;
    pRrR = pR / rR;
    uR = f.right_fs.u;
    eR = pR/((gamma-1.0)*rR); 
    aR = sqrt((pR*gamma)/rR);
    keR = 0.5 * (uR * uR);
    HR = eR + pRrR + keR;
    /*
     * This is the main part of the flux calculator.
     */
    /*
     * Weighting parameters (eqn 32) for velocity splitting.
     */
    alphaL = 2.0 * pLrL / (pLrL + pRrR);
    alphaR = 2.0 * pRrR / (pLrL + pRrR);
    /*
     * Common sound speed (eqn 33) and Mach numbers.
     */
    am = fmax(aL, aR);
    ML = uL / am;
    MR = uR / am;
    /*
     * Left state: 
     * pressure splitting (eqn 34) 
     * and velocity splitting (eqn 30)
     */
    duL = 0.5 * (uL + fabs(uL));
    if (fabs(ML) <= 1.0) {
	pLplus = pL * (ML + 1.0) * (ML + 1.0) * (2.0 - ML) * 0.25;
	uLplus =
	    alphaL * ((uL + am) * (uL + am) / (4.0 * am) - duL) +
	    duL;
    } else {
	pLplus = pL * duL / uL;
	uLplus = duL;
    }
    /*
     * Right state: 
     * pressure splitting (eqn 34) 
     * and velocity splitting (eqn 31)
     */
    duR = 0.5 * (uR - fabs(uR));
    if (fabs(MR) <= 1.0) {
	pRminus = pR * (MR - 1.0) * (MR - 1.0) * (2.0 + MR) * 0.25;
	uRminus =
	    alphaR * (-(uR - am) * (uR - am) / (4.0 * am) - duR) +
	    duR;
    } else {
	pRminus = pR * duR / uR;
	uRminus = duR;
    }
    /*
     * Mass Flux (eqn 29)
     */
    // The mass flux is relative to the moving interface.
    ru_half = uLplus * rL + uRminus * rR;
    /*
     * Pressure flux (eqn 34)
     */
    p_half = pLplus + pRminus;
    /*
     * Momentum flux: normal direction
     *
     * Compute blending parameter s (eqn 37),
     * the momentum flux for AUSMV (eqn 21) and AUSMD (eqn 21)
     * and blend (eqn 36).
     */
    dp = pL - pR;
    dp = K_SWITCH * fabs(dp) / fmin(pL, pR);
    s = 0.5 * fmin(1.0, dp);

    ru2_AUSMV = uLplus * rL * uL + uRminus * rR * uR;
    ru2_AUSMD = 0.5 * (ru_half * (uL + uR) -
		       fabs(ru_half) * (uR - uL));

    ru2_half = (0.5 + s) * ru2_AUSMV + (0.5 - s) * ru2_AUSMD;

    /*
     * Assemble components of the flux vector.
     */
    f.mass = ru_half;
    f.momentum = ru2_half + p_half;
    if (ru_half >= 0.0) {
	/* Wind is blowing from the left */
	f.energy = ru_half * HL;
    } else {
	/* Wind is blowing from the right */
	f.energy = ru_half * HR;
    }
    /*
     * Apply entropy fix (section 3.5 in Wada and Liou's paper)
     */
    caseA = ((uL - aL) < 0.0) && ((uR - aR) > 0.0);
    caseB = ((uL + aL) < 0.0) && ((uR + aR) > 0.0);

    d_ua = 0.0;
    if (caseA && !caseB)
	d_ua = C_EFIX * ((uR - aR) - (uL - aL));
    if (caseB && !caseA)
	d_ua = C_EFIX * ((uR + aR) - (uL + aL));

    if (d_ua != 0.0) {
	f.mass -= d_ua * (rR - rL);
	f.momentum -= d_ua * (rR * uR - rL * uL);
	f.energy -= d_ua * (rR * HR - rL * HL);
    }   /* end of entropy fix (d_ua != 0) */
} // end ausmdv_flux_calculator

void encode_conserved_variables(fvcell cell, double gamma) {
    double e = cell.p / (cell.rho*(gamma - 1.0));
    double ke = 0.5*cell.u*cell.u;
    cell.ru = cell.rho * cell.u;
    cell.rE = cell.rho*(e + ke);
} // end encode_conserved_variables

void decode_conserved_variables(fvcell cell, double gamma) {
    cell.u = cell.ru/cell.rho;
    double ke = 0.5*cell.u*cell.u;
    double e = cell.rE/cell.rho - ke;
    cell.p = cell.rho*e*(gamma-1.0);
} // end decode_conserved_variables

void matrix_mult(ref double[np*ncells][np*ncells] A, ref double[np*ncells][np*ncells] B, ref double[np*ncells][np*ncells] C) {
    for (int i = 0; i < np*ncells; i++) {
        for (int j = 0; j < np*ncells; j++) {
            C[i][j] = 0;
            for (int k = 0; k < np*ncells; k++) {
                C[i][j] += A[i][k]*B[k][j];
	    }
	}
    }
} // end matrix_mult

void matrix_inv(ref double[np*ncells][np*ncells] matrix, ref double[np*ncells][np*ncells] inverse) {
    // gauss jordan elimination of [A|I] where the output is [I|B] where B is the inverse of A.
    int N = np*ncells;
    static double[2*np*ncells][np*ncells] c;
    foreach(i; 0 .. N) {
	foreach(j; 0.. N) {
	    c[i][j] = matrix[i][j];
	}
    }
    foreach(i; 0 .. N) {
	foreach(j; N .. 2*N) {
	    if (N+i == j) {
		c[i][j] = 1.0;
	    } else {
		c[i][j] = 0.0;
	    }
	}
    }    
    foreach(j; 0 .. N) {
	// Select pivot.
	size_t p = j;
	foreach(i; j+1 .. N) {
	    if ( abs(c[i][j]) > abs(c[p][j]) ) p = i;
	}
	//if (abs(c[p][j]) <= very_small_value) return -1; // singular
	if ( p != j ) { // Swap rows
	    foreach(col; 0 .. 2*N) {
		double tmp = c[p][col]; c[p][col] = c[j][col]; c[j][col] = tmp;
	    }
	}
	// Scale row j to get unity on the diagonal.
	double cjj = c[j][j];
	foreach(col; 0 .. 2*N) c[j][col] /= cjj;
	// Do the elimination to get zeros in all off diagonal values in column j.
	foreach(i; 0 .. N) {
	    if ( i == j ) continue;
	    double cij = c[i][j];
	    foreach(col; 0 .. 2*N) c[i][col] -= cij * c[j][col]; 
	}
    } // end foreach j
    foreach(i; 0 .. N) {
	foreach(j; N .. 2*N) {
	    inverse[i][j-N] = c[i][j];
	}
    }
} // end matrix_inv

void flow_solver(ref fvcell[ncells+nghost] global_cells, ref fvinterface[ninterfaces] global_interfaces, ref fvinterface[ninterfaces] perturbed_interfaces,
		 double dt, double sim_time, double final_time, size_t step, size_t max_step, string flux_calc, string simulation_type, size_t interpolation_order,
		 size_t outflow_bc_order, double outflow_bc_back_pressure_factor, double dx, double gamma) {
    set_initial_and_inflow_conditions(simulation_type, global_cells, gamma);
    
    // begin by computing conserved quantities at cell centers
    foreach(i; 0 .. ncells+nghost) {
	encode_conserved_variables(global_cells[i], gamma);	
    }
    
    while (sim_time <= final_time && step < max_step) {
	//-----------------------------------------------------
	// Apply Boundary Conditions
	//-----------------------------------------------------
	apply_pre_flux_bcs(simulation_type, outflow_bc_order, outflow_bc_back_pressure_factor, global_cells, dx);

	//-----------------------------------------------------
	// Compute Flux
	//-----------------------------------------------------
	compute_flux(flux_calc, global_cells, global_interfaces, interpolation_order, ninterfaces, gamma);

	//-----------------------------------------------------
	// Apply Boundary Conditions
	//-----------------------------------------------------
	apply_post_flux_bcs(simulation_type, global_interfaces, global_cells, gamma);
		
	//-----------------------------------------------------
	// Integrate flux and update cells
	//-----------------------------------------------------
	
	foreach(i; 0 .. ncells) {
	    fvcell cell = global_cells[i+2];
	    fvinterface fin = global_interfaces[i];
	    fvinterface fout = global_interfaces[i+1];
	    double delta_rho; double delta_momentum; double delta_energy;
		
	    // update mass
	    delta_rho = dt/cell.vol * (fin.mass*fin.area - fout.mass*fout.area);
	    cell.rho = delta_rho + cell.rho;
	    
	    // update momentum
	    delta_momentum = dt/cell.vol * (fin.momentum*fin.area - fout.momentum*fout.area+cell.p*(fout.area-fin.area));
	    cell.ru = cell.ru + delta_momentum;
	    
	    // update total energy
	    delta_energy = dt/cell.vol * (fin.energy*fin.area - fout.energy*fout.area);
	    cell.rE = cell.rE + delta_energy;
	    
	    decode_conserved_variables(global_cells[i+2], gamma);
	}
	
	// update time and step
	sim_time += dt;
	step += 1;
    } // end while loop
}

void adjoint_solver(ref fvcell[ncells+nghost] global_cells,
		    ref fvinterface[ninterfaces] global_interfaces,
		    ref fvinterface[ninterfaces] perturbed_interfaces,
		    double dt, double sim_time, double final_time, size_t step, size_t max_step, string flux_calc, string simulation_type,
		    size_t interpolation_order, size_t outflow_bc_order, double outflow_bc_back_pressure_factor, double[np*ncells] psi,
		    ref double dLdB, ref double dLdC, ref double dLdD, double dx, double gamma,
		    double b, double c, double d, double yo, double scale, string solver, string adjoint_form, double[] p_target) {

    // intialise some arrays  ------------------------
    static double[np*ncells][np*ncells] JV;        // flow Jacobian w.r.t. primitive variables
    static double[np*ncells][np*ncells] JU;        // flow Jacobian w.r.t. conserved variables
    static double[np*ncells][np*ncells] transform; // transform matrix (JV to JU)
    static double[np*ncells] R;                    // R.H.S. resiudal vector
    static double[np*ncells] dJdV;                 // sensitivty of cost function (J) w.r.t. primitive variables
    static double[np*ncells] dJdU;                 // sensitivty of cost function (J) w.r.t. conserved variables
    static double[np*ncells][np*ncells] JUT;       // JU transpose
    static double[np*ncells][np*ncells] invJUT;    // JU transpose inverse
    static double[np*ncells][np*ncells] JVT;       // JV transpose
    static double[np*ncells][np*ncells] invJVT;    // JV transpose inverse
    
    //-------------------------------------------------------------------------------
    // Construct flow Jacobian w.r.t. primitive variables (JV) via finite-differences
    //-------------------------------------------------------------------------------
    
    foreach(i; 0 .. ncells) {
	//-----------------------------------------------
	// Save copy of current cells primitive variables
	//-----------------------------------------------
	double orig_rho = global_cells[i+2].rho;
	double orig_u = global_cells[i+2].u;
	double orig_p = global_cells[i+2].p;
	
	//--------------------------------
	// perturb density in current cell
	//--------------------------------
	global_cells[i+2].rho += (global_cells[i+2].rho*eps0+eps0);
	
	//-----------------------------------------------------
	// perform residual vector computation
	//-----------------------------------------------------
	apply_pre_flux_bcs(simulation_type, outflow_bc_order, outflow_bc_back_pressure_factor, global_cells, dx);
	compute_flux(flux_calc, global_cells, perturbed_interfaces, interpolation_order, ninterfaces, gamma);
	apply_post_flux_bcs(simulation_type, perturbed_interfaces, global_cells, gamma);
	
	//------------------------
	// Fill column of Jacobian
	//------------------------
	foreach(j; 0 .. ncells) {
	    double resd0;
	    double resd1;
	    fvinterface fin = global_interfaces[j];
	    fvinterface fout = global_interfaces[j+1];
	    fvinterface pfin = perturbed_interfaces[j];
	    fvinterface pfout = perturbed_interfaces[j+1];
	    fvcell cell = global_cells[j+2];
	    
	    // mass flux
	    resd0 = 1.0/cell.vol * (fin.mass*fin.area - fout.mass*fout.area);
	    resd1 = 1.0/cell.vol * (pfin.mass*pfin.area - pfout.mass*pfout.area);
	    JV[j*np][i*np] = (resd1-resd0)/(orig_rho*eps0+eps0);
	    
	    // fill transform matrix
	    if (i==j) transform[j*np][i*np] = 1.0;
	    else transform[j*np][i*np] = 0.0;
	    
	    // momentum flux
	    resd0 = 1.0/cell.vol * (fin.momentum*fin.area - fout.momentum*fout.area + cell.p*(fout.area-fin.area));
	    resd1 = 1.0/cell.vol * (pfin.momentum*pfin.area-pfout.momentum*pfout.area+ cell.p*(pfout.area-pfin.area));
	    JV[j*np+1][i*np] = (resd1-resd0)/(orig_rho*eps0+eps0);
	    
	    // fill transform matrix
	    if (i==j) transform[j*np+1][i*np] = -orig_u/orig_rho;
	    else transform[j*np+1][i*np] = 0.0;
	    
	    // total energy flux
	    resd0 = 1.0/cell.vol * (fin.energy*fin.area - fout.energy*fout.area);
	    resd1 = 1.0/cell.vol * (pfin.energy*pfin.area - pfout.energy*pfout.area);
	    JV[j*np+2][i*np] = (resd1-resd0)/(orig_rho*eps0+eps0);

	    // fill transform matrix
	    if (i==j) transform[j*np+2][i*np] = 0.5*(gamma-1.0)*orig_u*orig_u;
	    else transform[j*np+2][i*np] = 0.0;
	}
	
	//--------------------------------
	// restore density in current cell
	//--------------------------------
	global_cells[i+2].rho = orig_rho;
	
	//--------------------------------
	// perturb velocity in current cell
	//--------------------------------
	global_cells[i+2].u += (global_cells[i+2].u*eps0+eps0);

	//-----------------------------------------------------
	// perform residual vector computation
	//-----------------------------------------------------
	apply_pre_flux_bcs(simulation_type, outflow_bc_order, outflow_bc_back_pressure_factor, global_cells, dx);
	compute_flux(flux_calc, global_cells, perturbed_interfaces, interpolation_order, ninterfaces, gamma);
	apply_post_flux_bcs(simulation_type, perturbed_interfaces, global_cells, gamma);
		
	//------------------------
	// Fill column of Jacobian
	//------------------------	
	foreach(j; 0 .. ncells) {
	    double resd0;
	    double resd1;
	    fvinterface fin = global_interfaces[j];
	    fvinterface fout = global_interfaces[j+1];
	    fvinterface pfin = perturbed_interfaces[j];
	    fvinterface pfout = perturbed_interfaces[j+1];
	    fvcell cell = global_cells[j+2];
	    
	    // mass flux
	    resd0 = 1.0/cell.vol * (fin.mass*fin.area - fout.mass*fout.area);
	    resd1 = 1.0/cell.vol * (pfin.mass*pfin.area - pfout.mass*pfout.area);
	    JV[j*np][i*np+1] = (resd1-resd0)/(orig_u*eps0+eps0);
	    
	    // fill transform matrix
	    if (i==j) transform[j*np][i*np+1] = 0.0;
	    else transform[j*np][i*np+1] = 0.0;
	    
	    // momentum flux
	    resd0 = 1.0/cell.vol * (fin.momentum*fin.area - fout.momentum*fout.area + cell.p*(fout.area-fin.area));
	    resd1 = 1.0/cell.vol * (pfin.momentum*pfin.area-pfout.momentum*pfout.area+cell.p*(pfout.area-pfin.area));
	    JV[j*np+1][i*np+1] = (resd1-resd0)/(orig_u*eps0+eps0);
	    
	    // fill transform matrix
	    if (i==j) transform[j*np+1][i*np+1] = 1.0/orig_rho;
	    else transform[j*np+1][i*np+1] = 0.0;
	    
	    // total energy flux
	    resd0 = 1.0/cell.vol * (fin.energy*fin.area - fout.energy*fout.area);
	    resd1 = 1.0/cell.vol * (pfin.energy*pfin.area - pfout.energy*pfout.area);
	    JV[j*np+2][i*np+1] = (resd1-resd0)/(orig_u*eps0+eps0);

	    if (i==j) transform[j*np+2][i*np+1] = -(gamma-1.0)*orig_u;
	    else transform[j*np+2][i*np+1] = 0.0;
	}
	
	//--------------------------------
	// restore velocity in current cell
	//--------------------------------
	global_cells[i+2].u = orig_u;
	
	//--------------------------------
	// perturb pressure in current cell
	//--------------------------------
	global_cells[i+2].p += (global_cells[i+2].p*eps0+eps0);

	//-------------------------------------
	//  perform residual vector computation
	//------------------------------------
	apply_pre_flux_bcs(simulation_type, outflow_bc_order, outflow_bc_back_pressure_factor, global_cells, dx);
	compute_flux(flux_calc, global_cells, perturbed_interfaces, interpolation_order, ninterfaces, gamma);
	apply_post_flux_bcs(simulation_type, perturbed_interfaces, global_cells, gamma);
		
	//------------------------
	// Fill column of Jacobian
	//------------------------	
	foreach(j; 0 .. ncells) {
	    double resd0;
	    double resd1;
	    fvinterface fin = global_interfaces[j];
	    fvinterface fout = global_interfaces[j+1];
	    fvinterface pfin = perturbed_interfaces[j];
	    fvinterface pfout = perturbed_interfaces[j+1];
	    fvcell cell = global_cells[j+2];
	    
	    // mass flux
	    resd0 = 1.0/cell.vol * (fin.mass*fin.area - fout.mass*fout.area);
	    resd1 = 1.0/cell.vol * (pfin.mass*pfin.area - pfout.mass*fout.area);
	    JV[j*np][i*np+2] = (resd1-resd0)/(orig_p*eps0+eps0);
	    
	    // fill transform matrix
	    if (i==j) transform[j*np][i*np+2] = 0.0;
	    else transform[j*np][i*np+2] = 0.0;
	    
	    // momentum flux
	    if (i == j) resd0 = 1.0/cell.vol *(fin.momentum*fin.area-fout.momentum*fout.area+orig_p*(fout.area-fin.area));
	    else resd0 = 1.0/cell.vol * (fin.momentum*fin.area - fout.momentum*fout.area+cell.p*(fout.area-fin.area));
	    resd1 = 1.0/cell.vol * (pfin.momentum*pfin.area-pfout.momentum*pfout.area+cell.p*(pfout.area-pfin.area));
	    JV[j*np+1][i*np+2] = (resd1-resd0)/(orig_p*eps0+eps0);
	    
	    // fill transform matrix
	    if (i==j) transform[j*np+1][i*np+2] = 0.0;
	    else transform[j*np+1][i*np+2] = 0.0;
	    
	    // total energy flux
	    resd0 = 1.0/cell.vol * (fin.energy*fin.area - fout.energy*fout.area);
	    resd1 = 1.0/cell.vol * (pfin.energy*pfin.area - pfout.energy*pfout.area);
	    JV[j*np+2][i*np+2] = (resd1-resd0)/(orig_p*eps0+eps0);
	    
	    // fill transform matrix
	    if (i==j) transform[j*np+2][i*np+2] = gamma-1.0;
	    else transform[j*np+2][i*np+2] = 0.0;
	}
	
	//--------------------------------
	// restore pressure in current cell
	//--------------------------------
	global_cells[i+2].p = orig_p;
	
    }

    if (adjoint_form == "conservative") matrix_mult(JV, transform, JU); // transform JV to JU
    
    //-------------------
    // Transpose Jacobian
    //-------------------
    if (adjoint_form == "conservative") {
	foreach (i; 0 .. JU.length) {
	    foreach (j; 0 .. JU.length) {
		JUT[i][j] = JU[j][i];
	    }
	}
    } else { // primitive variable
	foreach (i; 0 .. JV.length) {
	    foreach (j; 0 .. JV.length) {
		JVT[i][j] = JV[j][i];
	    }
	}
    }

    //----------------
    // Invert Jacobian
    //----------------
    if (adjoint_form == "conservative") matrix_inv(JUT, invJUT);
    else matrix_inv(JVT, invJVT);
    
    //-------------------------------------
    // Construct dJdV vector (analytically)
    //-------------------------------------
    if (adjoint_form == "conservative") {
	foreach (i; 0..ncells) {
	    double U1 = global_cells[i+2].rho; double U2 = global_cells[i+2].ru; double U3 = global_cells[i+2].rE;
	    double pstar = p_target[i];
	    dJdU[i*np] = 1.0/2.0 * (U3*(gamma-1.0)*(gamma-1.0)*U2*U2/(U1*U1)
				    -U2*U2*U2*U2*(gamma-1.0)*(gamma-1.0)/(2.0*U1*U1*U1)
				    -U2*U2/(U1*U1)*(gamma-1.0)*pstar);
	    dJdU[i*np+1] = 1.0/2.0 * (-2.0*U3*(gamma-1.0)*(gamma-1.0)*U2/U1
				      +U2*U2*U2*(gamma-1.0)*(gamma-1.0)/(U1*U1)
				      +2.0*U2/U1*(gamma-1.0)*pstar);
	    dJdU[i*np+2] = 1.0/2.0 * (2.0*U3*(gamma-1.0)*(gamma-1.0)
				      -(gamma-1.0)*(gamma-1.0)*U2*U2/(U1)
				      -2.0*(gamma-1.0)*pstar);
	}
    } else { // primitive
	foreach (i; 0..ncells) {
	    foreach (j; 0..ncells) {
		dJdV[i*np] = 0.0;
		dJdV[i*np+1] = 0.0;
		dJdV[i*np+2] = 0.5*(2.0*global_cells[i+2].p-2.0*p_target[i]);
	    }
	}
    }
    
    //--------------------------------------------------------
    // Compute adjoint variables via inv[Jacobian]^T * -dJdQ
    //--------------------------------------------------------
    foreach (i; 0 .. ncells) {
	double psi_rho = 0.0; double psi_ru = 0.0; double psi_rE = 0.0;
	foreach (j; 0 .. invJUT.length) {
	    if (adjoint_form == "conservative") {
		psi_rho += invJUT[i*np][j] * -dJdU[j];
		psi_ru += invJUT[i*np+1][j] * -dJdU[j];
		psi_rE += invJUT[i*np+2][j] * -dJdU[j];
	    } else { // primitive
		psi_rho += invJVT[i*np][j] * -dJdV[j];
		psi_ru += invJVT[i*np+1][j] * -dJdV[j];
		psi_rE += invJVT[i*np+2][j] * -dJdV[j];
	    }
	}
	psi[i*np] = psi_rho;
	psi[i*np+1] = psi_ru;
	psi[i*np+2] = psi_rE;
    }
    
    //---------------------------------------
    // Construct dR/dD via finite-differences
    //---------------------------------------
    static double[nd][np*ncells] dRdD;

    // first compute current residual vector
    foreach(i; 0 .. ncells) {
	fvinterface fin = global_interfaces[i];
	fvinterface fout = global_interfaces[i+1];
	fvcell cell = global_cells[i+2];
	R[i*np] = 1.0/cell.vol * (fin.mass*fin.area - fout.mass*fout.area);
	R[i*np+1] = 1.0/cell.vol * (fin.momentum*fin.area - fout.momentum*fout.area + cell.p*(fout.area-fin.area));
	R[i*np+2] = 1.0/cell.vol * (fin.energy*fin.area - fout.energy*fout.area);
    }

    double b_orig = b;
    b += b*eps1+eps1;
    compute_geometry(simulation_type, global_interfaces, perturbed_interfaces, global_cells, b, c, d, yo, scale);

    // now construct dR/dD using finite-differences
    foreach(j; 0 .. ncells) {
	double R_perturbed;
	fvinterface fin = perturbed_interfaces[j];
	fvinterface fout = perturbed_interfaces[j+1];
	fvcell cell = global_cells[j+2];
	
	// mass flux
	R_perturbed = 1.0/cell.vol * (fin.mass*fin.area - fout.mass*fout.area);
	dRdD[j*np][0] = (R_perturbed-R[j*np])/(b*eps1+eps1);
	
	// momentum flux
	R_perturbed = 1.0/cell.vol * (fin.momentum*fin.area - fout.momentum*fout.area + cell.p*(fout.area-fin.area));
	dRdD[j*np+1][0] = (R_perturbed-R[j*np+1])/(b*eps1+eps1);
	
	// total energy flux
	R_perturbed = 1.0/cell.vol * (fin.energy*fin.area - fout.energy*fout.area);
	dRdD[j*np+2][0] = (R_perturbed-R[j*np+2])/(b*eps1+eps1);
    }
    b = b_orig;

    double c_orig = c;
    c += c*eps1+eps1;
    compute_geometry(simulation_type, global_interfaces, perturbed_interfaces, global_cells, b, c, d, yo, scale);

    // now construct dR/dD using finite-differences 
    foreach(j; 0 .. ncells) {
	double R_perturbed;
	fvinterface fin = perturbed_interfaces[j];
	fvinterface fout = perturbed_interfaces[j+1];
	fvcell cell = global_cells[j+2];
	
	// mass flux
	R_perturbed = 1.0/cell.vol * (fin.mass*fin.area - fout.mass*fout.area);
	dRdD[j*np][1] = (R_perturbed-R[j*np])/(c*eps1+eps1);
	
	// momentum flux
	R_perturbed = 1.0/cell.vol * (fin.momentum*fin.area - fout.momentum*fout.area + cell.p*(fout.area-fin.area));
	dRdD[j*np+1][1] = (R_perturbed-R[j*np+1])/(c*eps1+eps1);
	
	// total energy flux
	R_perturbed = 1.0/cell.vol * (fin.energy*fin.area - fout.energy*fout.area);
	dRdD[j*np+2][1] = (R_perturbed-R[j*np+2])/(c*eps1+eps1);
    }
    c = c_orig;

    double d_orig = d;
    d += d*eps1+eps1;
    compute_geometry(simulation_type, global_interfaces, perturbed_interfaces, global_cells, b, c, d, yo, scale);

    // now construct dR/dD using finite-differences 
    foreach(j; 0 .. ncells) {
	double R_perturbed;
	fvinterface fin = perturbed_interfaces[j];
	fvinterface fout = perturbed_interfaces[j+1];
	fvcell cell = global_cells[j+2];
	
	// mass flux
	R_perturbed = 1.0/cell.vol * (fin.mass*fin.area - fout.mass*fout.area);
	dRdD[j*np][2] = (R_perturbed-R[j*np])/(d*eps1+eps1);
	
	// momentum flux
	R_perturbed = 1.0/cell.vol * (fin.momentum*fin.area - fout.momentum*fout.area + cell.p*(fout.area-fin.area));
	dRdD[j*np+1][2] = (R_perturbed-R[j*np+1])/(d*eps1+eps1);
	
	// total energy flux
	R_perturbed = 1.0/cell.vol * (fin.energy*fin.area - fout.energy*fout.area);
	dRdD[j*np+2][2] = (R_perturbed-R[j*np+2])/(d*eps1+eps1);
    }
    d = d_orig;
    
    //-----------------------------------------------------
    // Compute gradients via  dLdD = dRdD * psi
    //-----------------------------------------------------
    dLdB = 0.0;
    dLdC = 0.0;
    dLdD = 0.0;
    foreach (j; 0 .. psi.length) {
	dLdB += dRdD[j][0] * psi[j];
	dLdC += dRdD[j][1] * psi[j];
	dLdD += dRdD[j][2] * psi[j];
    }
}
