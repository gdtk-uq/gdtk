// block.d
module block;
import std.stdio;
import std.conv;
import std.math;
import std.file;
import std.format;
import std.string;
import config;
import finite_volume;
import number;
import complexify;
import derivative;
import flux;
import linalg;

class Block {
public:
    string label;
    string flux_calc;
    size_t ncells;
    size_t nifaces;
    int interpolation_order;
    number gamma;
    number length;
    number dx;
    number bc_factor;    
    double eps;

    Derivative derivative;
    FVInterface[] ifaces;  // interafces
    FVInterface[] pifaces; // perturbed interfaces
    FVCell[] cells;        // cells
    double[] adjoint_vars;
    
    this(Params params) {
	this.label = params.simulation_type;
        this.length = params.length;
        this.ncells = params.ncells;
        this.nifaces = ncells+1;
        this.dx = length/ncells;
        this.gamma = params.gamma;
        this.bc_factor = params.nozzle_back_pressure_factor;
        this.derivative = new Derivative(params, this.ncells, this.nifaces);
        this.eps = params.eps;
        this.flux_calc = params.flux_calc;
        this.interpolation_order = params.interpolation_order;
    }
    
    abstract void compute_geometry();    
    abstract void set_initial_flowstate_condition();
    abstract void set_inflow_and_outflow_conditions();
    abstract void apply_boundary_conditions();

    void read_flow_solution_from_file(ref number[] target) {
        auto file = File("target_flow.dat", "r");
        foreach(i; 0 .. ncells) {
            auto lineContent = file.readln().strip();
            auto tokens = lineContent.split();
            target ~= to!number(tokens[2]);
        }
    } // end read_flow_solution_from_file()
    
    void write_flow_solution_to_file() {
        string flow_filename = "flow.dat";
        remove_old_files(flow_filename);
        foreach(cell; cells) {
            if (cell.is_ghost) { continue; }
            auto writer = format("%f %f %f %f\n", cell.xpos.re, cell.fs.rho.re, cell.fs.p.re, cell.fs.vel.re);
            append(flow_filename, writer);
        }
        string contour_filename = "contour.dat";        
        remove_old_files(contour_filename);
        foreach(f; ifaces) {
            auto writer = format("%f %f \n", f.xpos.re, f.area.re/2.0);
            append(contour_filename, writer);
        }
    } // end write_flow_solution_to_file()

    void construct_geometry() {

        // cell interfaces
        number tmp_pos = 0.0; // intiial iface position
        foreach(i; 0 .. nifaces) {
            ifaces ~= new FVInterface(i, tmp_pos);
            pifaces  ~= new FVInterface(i, tmp_pos);
            tmp_pos += dx;
        }
        
        // east-end ghost cells
        cells ~= new FVCell(10000001, -1.5*dx, dx, true);
        cells ~= new FVCell(10000002, -0.5*dx, dx, true);
        
        // interior cells
        tmp_pos = 0.5*dx; // initial cell position
        foreach(i; 0 .. ncells) {
            cells ~= new FVCell(i, tmp_pos, dx);
            cells[nghost+i].left_iface = ifaces[i];
            cells[nghost+i].right_iface = ifaces[i+1];
            tmp_pos += dx ;
        }
        
        // south-end ghost cells
        cells ~= new FVCell(10000003, length+0.5*dx, dx, true);
        cells ~= new FVCell(10000004, length+1.5*dx, dx, true);   
        
        // store interface cell cloud
        foreach(i, iface; ifaces) {
            iface.cell_cloud ~= cells[i];
            iface.cell_cloud ~= cells[i+1];
            iface.cell_cloud ~= cells[i+2];
            iface.cell_cloud ~= cells[i+3];
        }
    }

    void determine_timestep(ref number dt, number cfl) {
        dt = double.max;
        foreach(i; 0 .. ncells ) {
            FVCell cell = cells[nghost+i];
            number dt_suggest = cell.calculate_stable_timestep(cfl, gamma);
            dt = fmin(dt, dt_suggest);
        }
    } // end determine_timestep()
    
    void compute_fluxes() {
        foreach (iface; ifaces) {
            iface.reconstruct_flowstate(interpolation_order);
            compute_flux(iface, gamma, flux_calc);
        }
    } // end compute_fluxes()
    
    void update_conserved_quantities(number dt) {
        foreach(cell; cells) {
            if (cell.is_ghost) { continue; }
            cell.update_conserved_quantities(dt);
        }
    } // end update_conserved_quantities()

    void encode_conserved_quantities() {
        foreach(cell; cells) {
            if (cell.is_ghost) { continue; }
            cell.encode_conserved(gamma);	
        }
    } // end encode_conserved_quantities()

    void decode_conserved_quantities() {
        foreach(cell; cells) {
            if (cell.is_ghost) { continue; }
            cell.decode_conserved(gamma);	
        }
    } // end decode_conserved_quantities()

    void evaluate_transpose_jacobian() {
        import std.stdio;
        size_t np = nprimitive;        
        // unperturbed state
        foreach (k; 0..nifaces) { derivative.ifaces_cpy[k].copy_from(ifaces[k]); }
        foreach (k; 0..nghost+ncells+nghost) { derivative.cells_cpy[k].copy_from(cells[k]); }
        
        foreach(pcell; cells) {
            if (pcell.is_ghost) { continue; }
            size_t i = pcell.id;
            mixin(evaluate_block_row_of_jacoboian("rho", "0"));
            mixin(evaluate_block_row_of_jacoboian("vel", "1"));
            mixin(evaluate_block_row_of_jacoboian("p", "2"));
        }
        jacobian_boundary_conditions();
    } // end evaluate_transpose_jacobian()

    
    void jacobian_boundary_conditions() {
        import std.stdio;
        size_t np = nprimitive;
        size_t bndary_cell_idx = ncells + 1;

        // there are two ghost cells attached to the east boundary
        // west boundary is supersonic inflow - no effect from ghost cell perturbations
        size_t[2] ghost_cell_idxs = [bndary_cell_idx+1, bndary_cell_idx+2];
        size_t[2] internal_cell_idxs = [bndary_cell_idx-1, bndary_cell_idx];
        foreach (cell_idx; internal_cell_idxs) {
            foreach (ghost_cell_idx; ghost_cell_idxs) {
                mixin(evaluate_dqdQ("rho", "0"));
                mixin(evaluate_dqdQ("vel", "1"));
                mixin(evaluate_dqdQ("p", "2"));
                
                mixin(evaluate_dRdq("rho", "0"));
                mixin(evaluate_dRdq("vel", "1"));
                mixin(evaluate_dRdq("p", "2"));
                
                matrix_multiply(derivative.dRdq, derivative.dqdQ, derivative.Lext);
                // correct boundary term
                size_t idx = bndary_cell_idx-2;
                size_t idxR = cell_idx-2;
                foreach (i; 0..3) {
                    foreach (j; 0..3) {
                        derivative.jac_VT[idx*np+i][idxR*np+j] += derivative.Lext[j][i];
                    }
                }
            }
        }
    } // end jacobian_boundary_conditions()
} // end Block class

version(complex_numbers) {
    string evaluate_dqdQ(string pvar, string idx)
    {
        return `
        cells[bndary_cell_idx].fs.`~pvar~` += complex(0.0, eps);
        apply_boundary_conditions();
        derivative.dqdQ[0][`~idx~`] = cells[ghost_cell_idx].fs.rho.im/eps;
        derivative.dqdQ[1][`~idx~`] = cells[ghost_cell_idx].fs.vel.im/eps;
        derivative.dqdQ[2][`~idx~`] = cells[ghost_cell_idx].fs.p.im/eps;
        foreach (k; 0..nifaces) { ifaces[k].copy_from(derivative.ifaces_cpy[k]); }
        foreach (k; 0..nghost+ncells+nghost) { cells[k].copy_from(derivative.cells_cpy[k]); }
        `;
    }
    
    string evaluate_dRdq(string pvar, string idx)
    {
        return `
        cells[ghost_cell_idx].fs.`~pvar~` += complex(0.0, eps);
        compute_fluxes();
        update_conserved_quantities(to!number(1.0));
        derivative.dRdq[0][`~idx~`] = cells[cell_idx].U.r.im/eps;
        derivative.dRdq[1][`~idx~`] = cells[cell_idx].U.ru.im/eps;
        derivative.dRdq[2][`~idx~`] = cells[cell_idx].U.rE.im/eps;
        foreach (k; 0..nifaces) { ifaces[k].copy_from(derivative.ifaces_cpy[k]); }
        foreach (k; 0..nghost+ncells+nghost) { cells[k].copy_from(derivative.cells_cpy[k]); }        
        `;
    }
    
    string evaluate_block_row_of_jacoboian(string pvar, string idx)
    {
        return `
        pcell.fs.`~pvar~` += complex(0.0, eps);
        compute_fluxes();
        update_conserved_quantities(to!number(1.0));
        foreach(cell; cells) {
            if (cell.is_ghost) { continue; }
            size_t j = cell.id;
            derivative.jac_VT[i*np+`~idx~`][j*np+0] = cell.U.r.im/eps;
            derivative.jac_VT[i*np+`~idx~`][j*np+1] = cell.U.ru.im/eps;
            derivative.jac_VT[i*np+`~idx~`][j*np+2] = cell.U.rE.im/eps;
        }
        foreach (k; 0..nifaces) { ifaces[k].copy_from(derivative.ifaces_cpy[k]); }
        foreach (k; 0..nghost+ncells+nghost) { cells[k].copy_from(derivative.cells_cpy[k]); }
        `;
    }
} else { // real_numbers
    string evaluate_dqdQ(string pvar, string idx)
    {
        return `
        cells[bndary_cell_idx].fs.`~pvar~` += eps;
        apply_boundary_conditions();
        derivative.dqdQ[0][`~idx~`] = cells[ghost_cell_idx].fs.rho;
        derivative.dqdQ[1][`~idx~`] = cells[ghost_cell_idx].fs.vel;
        derivative.dqdQ[2][`~idx~`] = cells[ghost_cell_idx].fs.p;
        foreach (k; 0..nifaces) { ifaces[k].copy_from(derivative.ifaces_cpy[k]); }
        foreach (k; 0..nghost+ncells+nghost) { cells[k].copy_from(derivative.cells_cpy[k]); }

        cells[bndary_cell_idx].fs.`~pvar~` -= eps;
        apply_boundary_conditions();
        derivative.dqdQ[0][`~idx~`] -= cells[ghost_cell_idx].fs.rho;
        derivative.dqdQ[1][`~idx~`] -= cells[ghost_cell_idx].fs.vel;
        derivative.dqdQ[2][`~idx~`] -= cells[ghost_cell_idx].fs.p;
        foreach (k; 0..nifaces) { ifaces[k].copy_from(derivative.ifaces_cpy[k]); }
        foreach (k; 0..nghost+ncells+nghost) { cells[k].copy_from(derivative.cells_cpy[k]); }

        derivative.dqdQ[0][`~idx~`] /= (2.0*eps);
        derivative.dqdQ[1][`~idx~`] /= (2.0*eps);
        derivative.dqdQ[2][`~idx~`] /= (2.0*eps);
        `;
    }
    
    string evaluate_dRdq(string pvar, string idx)
    {
        return `
        cells[ghost_cell_idx].fs.`~pvar~` += eps;
        compute_fluxes();
        update_conserved_quantities(to!number(1.0));
        derivative.dRdq[0][`~idx~`] = cells[cell_idx].U.r;
        derivative.dRdq[1][`~idx~`] = cells[cell_idx].U.ru;
        derivative.dRdq[2][`~idx~`] = cells[cell_idx].U.rE;
        foreach (k; 0..nifaces) { ifaces[k].copy_from(derivative.ifaces_cpy[k]); }
        foreach (k; 0..nghost+ncells+nghost) { cells[k].copy_from(derivative.cells_cpy[k]); }        

        cells[ghost_cell_idx].fs.`~pvar~` -= eps;
        compute_fluxes();
        update_conserved_quantities(to!number(1.0));
        derivative.dRdq[0][`~idx~`] -= cells[cell_idx].U.r;
        derivative.dRdq[1][`~idx~`] -= cells[cell_idx].U.ru;
        derivative.dRdq[2][`~idx~`] -= cells[cell_idx].U.rE;
        foreach (k; 0..nifaces) { ifaces[k].copy_from(derivative.ifaces_cpy[k]); }
        foreach (k; 0..nghost+ncells+nghost) { cells[k].copy_from(derivative.cells_cpy[k]); }        

        derivative.dRdq[0][`~idx~`] /= (2.0*eps);
        derivative.dRdq[1][`~idx~`] /= (2.0*eps);
        derivative.dRdq[2][`~idx~`] /= (2.0*eps);
        `;
    }
    
    string evaluate_block_row_of_jacoboian(string pvar, string idx)
    {
        return `
        pcell.fs.`~pvar~` += eps;
        compute_fluxes();
        update_conserved_quantities(to!number(1.0));
        foreach(cell; cells) {
            if (cell.is_ghost) { continue; }
            size_t j = cell.id;
            derivative.jac_VT[i*np+`~idx~`][j*np+0] = cell.U.r;
            derivative.jac_VT[i*np+`~idx~`][j*np+1] = cell.U.ru;
            derivative.jac_VT[i*np+`~idx~`][j*np+2] = cell.U.rE;
        }
        foreach (k; 0..nifaces) { ifaces[k].copy_from(derivative.ifaces_cpy[k]); }
        foreach (k; 0..nghost+ncells+nghost) { cells[k].copy_from(derivative.cells_cpy[k]); }

        pcell.fs.`~pvar~` -= eps;
        compute_fluxes();
        update_conserved_quantities(to!number(1.0));
        foreach(cell; cells) {
            if (cell.is_ghost) { continue; }
            size_t j = cell.id;
            derivative.jac_VT[i*np+`~idx~`][j*np+0] -= cell.U.r;
            derivative.jac_VT[i*np+`~idx~`][j*np+1] -= cell.U.ru;
            derivative.jac_VT[i*np+`~idx~`][j*np+2] -= cell.U.rE;
        }
        foreach (k; 0..nifaces) { ifaces[k].copy_from(derivative.ifaces_cpy[k]); }
        foreach (k; 0..nghost+ncells+nghost) { cells[k].copy_from(derivative.cells_cpy[k]); }
        
        foreach(cell; cells) {
            if (cell.is_ghost) { continue; }
            size_t j = cell.id;
            derivative.jac_VT[i*np+`~idx~`][j*np+0] /= (2.0*eps);
            derivative.jac_VT[i*np+`~idx~`][j*np+1] /= (2.0*eps);
            derivative.jac_VT[i*np+`~idx~`][j*np+2] /= (2.0*eps);
        }
        `;
    }
}

class Shocktube : Block {
public:
    this(Params params) {
        super(params);
        construct_geometry();
        compute_geometry();
    }
    
    override void compute_geometry() {
        // compute interface area
        number inlet_area;
        number exit_area;
        size_t nfaces = ncells+1;

        // shock tube geometry has no area variaton
        foreach(i; 0 .. nifaces) {
            ifaces[i].area = 1.0;
            pifaces[i].area = 1.0;
        }

        // compute cell volume
        foreach(i; 0 .. ncells) {
            cells[nghost+i].vol = 0.5*cells[nghost+i].dx*(ifaces[i].area+ifaces[i+1].area);
        }
    } // end compute_geometry()

    override void set_initial_flowstate_condition() {
        size_t ncells = ncells;
        number half_length = length/2.0;
        
        foreach(cell; cells) {
            if (cell.xpos <= half_length) {
                // left fill state (including west boundary ghost cells)
                number p_init = 1.0e05; // Pa
                number rho_init = 1.0; // kg/m^3
                number u_init = 0.0; // m/s
                cell.fs.p = p_init;
                cell.fs.rho = rho_init;
                cell.fs.vel = u_init;
            } else {
                // right fill state (including east boundary ghost cells)
                number p_init = 1.0e04; // Pa
                number rho_init = 0.125; // kg/m^3
                number u_init = 0.0; // m/s
                cell.fs.p = p_init;
                cell.fs.rho = rho_init;
                cell.fs.vel = u_init;
            }
        }
    } // end set_initial_flowstate_condition()

    override void set_inflow_and_outflow_conditions() {
        /*
          A shocktube has no inflow and outflow conditions.
        */
    } // end set_inflow_and_outflow_conditions()
    
    override void apply_boundary_conditions() {
	// west boundary is a wall
	cells[0].fs.rho = cells[2].fs.rho;
	cells[1].fs.rho = cells[2].fs.rho;
	cells[0].fs.vel = -cells[2].fs.vel;
	cells[1].fs.vel = -cells[2].fs.vel;
	cells[0].fs.p = cells[2].fs.p;
	cells[1].fs.p = cells[2].fs.p;
        
        // east boundary is a wall
	cells[ncells+2].fs.rho = cells[ncells+1].fs.rho;
	cells[ncells+3].fs.rho = cells[ncells+1].fs.rho;
	cells[ncells+2].fs.vel = -cells[ncells+1].fs.vel;
	cells[ncells+3].fs.vel = -cells[ncells+1].fs.vel;
	cells[ncells+2].fs.p = cells[ncells+1].fs.p;
	cells[ncells+3].fs.p = cells[ncells+1].fs.p;
    } // apply_boundary_conditions()   
} // end Shocktube class

class Nozzle : Block {
public:
    number b;
    number c;
    number d;
    number yo;
    number scale;
    
    this(Params params) {
        super(params);
        this.b = params.b;
	this.c = params.c;
	this.d = params.d;
        this.yo = params.yo;
        this.scale = params.scale;
        construct_geometry();
        compute_geometry();
    }
    
    override void compute_geometry() {
        // compute interface area
        number a = (yo - b*tanh(-d/scale));
        foreach(i; 0 .. nifaces) {
            number height = a + b*tanh((c*ifaces[i].xpos-d)/scale);
            ifaces[i].area = 2.0*height;
            pifaces[i].area = 2.0*height;
        }
        
        // compute cell volume
        foreach(i; 0 .. ncells) {
            cells[nghost+i].vol = 0.5*cells[nghost+i].dx*(ifaces[i].area+ifaces[i+1].area);
        }
    } // end compute_geometry()

    override void set_initial_flowstate_condition() {

        size_t ncells = ncells;
        number half_length = length/2.0;
        
        // set inflow conditions
        number Mach = 1.5;
        number p_init = 101.325e3;                          // Pa
        number rho_init = 1.0;                              // kg/m^3
        number a_init = sqrt((p_init*gamma)/rho_init);      // m/s
        number u_init = Mach * a_init;                      // m/s
        
        
        foreach(i; 0 .. ncells) {
            cells[nghost+i].fs.p = p_init;
            cells[nghost+i].fs.rho = rho_init;
            cells[nghost+i].fs.vel = 0.6 * a_init;
        }
        
    } // end set_initial_flowstate_condition()

    override void set_inflow_and_outflow_conditions() {
        
        // set inflow conditions
        number Mach = 1.5;
        number p_inflow = 101.325e3;                          // Pa
        number rho_inflow = 1.0;                              // kg/m^3
        number a_inflow = sqrt((p_inflow*gamma)/rho_inflow);  // m/s
        number u_inflow = Mach * a_inflow;                    // m/s
        
        // first two cells are ghost cells
        cells[0].fs.p = p_inflow; cells[1].fs.p = p_inflow;
        cells[0].fs.rho = rho_inflow; cells[1].fs.rho = rho_inflow;
        cells[0].fs.vel = u_inflow; cells[1].fs.vel = u_inflow;
        
        // set outflow conditions
        // last two cells are ghost cells
        cells[ncells+2].fs.p = 0.0*cells[ncells+1].fs.p; cells[ncells+3].fs.p = 0.0*cells[ncells+1].fs.p;
        cells[ncells+2].fs.rho =  0.0*cells[ncells+1].fs.rho; cells[ncells+3].fs.rho =  0.0*cells[ncells+1].fs.rho;
        cells[ncells+2].fs.vel =  0.0*cells[ncells+1].fs.vel; cells[ncells+3].fs.vel =  0.0*cells[ncells+1].fs.vel;	
        
    } // end set_inflow_and_outflow_conditions()
    
    override void apply_boundary_conditions() {
        // west boundary is an inlet
	cells[0].fs.rho = cells[0].fs.rho;
	cells[1].fs.rho = cells[1].fs.rho;
	cells[0].fs.vel = cells[0].fs.vel;
	cells[1].fs.vel = cells[1].fs.vel;
	cells[0].fs.p = cells[0].fs.p;
	cells[1].fs.p = cells[1].fs.p;
        
        // east boundary is a first-order outlet
        cells[ncells+2].fs.rho = cells[ncells+1].fs.rho;
        cells[ncells+3].fs.rho = cells[ncells+1].fs.rho;
        cells[ncells+2].fs.vel = cells[ncells+1].fs.vel;
        cells[ncells+3].fs.vel = cells[ncells+1].fs.vel;
        cells[ncells+2].fs.p = bc_factor*cells[0].fs.p;
        cells[ncells+3].fs.p = bc_factor*cells[0].fs.p;
    } // apply_boundary_conditions()
} // end Nozzle class
