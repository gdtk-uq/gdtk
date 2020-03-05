// optimizer.d
module optimizer;
import std.stdio;
import std.conv;
import std.math;
import config;
import number;
import complexify;
import flowsolve;
import finite_volume;
import block;
import derivative;
import linalg;

class Optimizer {
public:
    bool adjoint = false;
    int max_opt_iters;
    number Jref = 0.0;
    number[] obj_func_hist;
    number eps; 
    number opt_tol;
    number mu = 1.0e-02;
    number[] target;    
    number[] pk;        // current search direction
    number[] pk_old;    // previous search direction
    double[] dLdD;      // design sensitivities
    double[] dLdD_old;  // previous design sensitivities
    double[] dvar0;     // original design variables
    Nozzle block;
    Solver solver;
    
    this(Params params, Nozzle block, Solver solver) {
        this.adjoint = params.adjoint_method;
        this.max_opt_iters = params.max_opt_iters;
        this.opt_tol = params.opt_tol;  
        this.eps = params.eps;
        this.dLdD_old.length = ndesign;
        this.dLdD.length = ndesign;
        this.dvar0.length = ndesign;
        this.pk.length = 3;
        this.pk_old.length = 3;
        this.block = block;
        size_t np = nprimitive;
        size_t nc = block.ncells;
        prep_vector(block.adjoint_vars, np*nc);        
        this.solver = solver;

    }
    
    void optimize() {   
        /*
          Main optimization loop
         */
        
        writeln("BASELINE GEOMETRY (b, c, d): ", block.b, ", ", block.c, ", ", block.d);
        int iter = 0;
        number J = 0;         // objective function value
        number Jref = 0;      // reference value
        number a_star = 0.0;  // step size
        
        while (iter < max_opt_iters) {
            dvar0 = [block.b.re, block.c.re, block.d.re];

            // determine if converged
            solver.integrate_in_time(block);
            J = eval_obj_func();
            if (iter == 0) {
                // set the objective function reference value on the initial iteration
                Jref = J;
            }
            obj_func_hist ~= J/Jref;
	    // if the drop in the objective function value is sufficient then stop the optimizer
	    if ( (J/Jref) < opt_tol) break;
            
            eval_grads(adjoint);
            compute_search_direction(iter);
            a_star = line_search(J);
            update_shape_params(a_star);
            block.compute_geometry();

            // store the current steps direction for use in the next step
            foreach (i; 0..ndesign) {
                dLdD_old[i] = dLdD[i]; 
                pk_old[i] = pk[i];
            }
            
            iter += 1;
            writeln("--------------------------------------------------------------------------");
            writeln(iter, ". J/Jref = ", J/Jref, ", b = ", block.b, ", c = ", block.c, ", d = ", block.d);
            writeln("gradients: ", ", dLdB = ", dLdD[0],  ", dLdC = ", dLdD[1],  ", dLdD = ", dLdD[2]);
            writeln("--------------------------------------------------------------------------");
        } // end while()
    } // end optimize()
    
    void eval_grads(bool adjoint = false) {        
        (adjoint) ? adjoint_method() :  numerical_gradients();
    } // end eval_grads()
    
    void numerical_gradients() {
        /*
          Calculates numerical approximations of the shape gradients using:
              - finite-differences for the real-variable solver
              - complex-step differentation for the complex-variable solver.

          Ref.:
              A Coupled-Adjoint Method for High-Fidelity Aero-Structural Optimization,
              J.R.R.A. Martins, 2002, Stanford University.
        */
        number J_plus; number J_minus;
        mixin(numerical_gradient_code("b", "0"));
        mixin(numerical_gradient_code("c", "1"));
        mixin(numerical_gradient_code("d", "2"));
    } // end numerical_gradients()
    
    void adjoint_method() {
        /*
          Calculates semi-analytic gradients using the adjoint method.
          Ref.:
              Adjoint-Based Aerodynamic Design Optimisation in Hypersonic Flow,
              K. Damm, 2020, University of Queensland

              Efficient Construction of Discrete Adjoint Operators on Unstructured Grids
              Using Complex Variables,
              E.J. Nielsen & W.L. Kleb, AIAA Journal, 2006
        */
        // evaluate transpose Jacobian
        block.evaluate_transpose_jacobian();
        // evaluate objective function sensitivity w.r.t. flow state variables
        eval_obj_sens_wrt_flowstate();
        // solve for adjoint variables
        matrix_inverse(block.derivative.jac_VT);
        matrix_vector_multiply(block.derivative.jac_VT, block.derivative.obj_sens, block.adjoint_vars);
        // evaluate residual sensitivity w.r.t. to design variables
        eval_resid_sens_wrt_dvars();
        // calculate total derivative
        matrix_vector_multiply(block.derivative.resid_sens, block.adjoint_vars, dLdD);
    } // end adjoint_method()

    number eval_obj_func() {
        /* 
           Objective function evaluation: 
           inverse design, matching pressure distribution.       
        */
        number J = 0.0;
        size_t ncells = block.ncells;
        foreach (i; 0..ncells) {
            number p = block.cells[i+nghost].fs.p;
            number pstar = target[i];
            J += 0.5*(p-pstar)^^2;
        }
        return J;
    } // end eval_obj_func()

    void eval_obj_sens_wrt_flowstate() {
        /*
          Objective function sensitivity w.r.t. flow state variables/
          Hand differentiated function defined in eval_obj_func() routine.
        */
        size_t ncells = block.ncells;
        size_t np = nprimitive;
        foreach (i; 0..ncells) {
            block.derivative.obj_sens[i*np] = 0.0;
            block.derivative.obj_sens[i*np+1] = 0.0;
            block.derivative.obj_sens[i*np+2] = (0.5*(2.0*block.cells[i+2].fs.p-2.0*target[i])).re;
        }
        block.derivative.obj_sens[] = -1*block.derivative.obj_sens[];
    } // end eval_obj_sens_wrt_flowstate()

    void eval_resid_sens_wrt_dvars() {
        /*
          Residual sensitivity w.r.t. to design variables.
          Evaluated using:
              - finite-differences for the real-variable solver
              - complex-step differentation for the complex-variable solver.
         */
        size_t ncells = block.ncells;
        size_t nifaces = block.nifaces;
        size_t np = nprimitive;
        // unperturbed state
        foreach (k; 0..nifaces) { block.derivative.ifaces_cpy[k].copy_from(block.ifaces[k]); }
        foreach (k; 0..nghost+ncells+nghost) { block.derivative.cells_cpy[k].copy_from(block.cells[k]); }
       
        mixin(eval_resid_sens_row_code("b", "0"));
        mixin(eval_resid_sens_row_code("c", "1"));
        mixin(eval_resid_sens_row_code("d", "2"));
    } // end eval_resid_sens_wrt_dvars()
    
    void compute_search_direction(int iter) {
        /*
          Currently use Conjugate Gradient Method with a restart every ndesign iterations
          for robustness as suggested in
          Ref.:
              Multidisciplinary Design Optimization,
              J.R.R.A. Martins, A. Ning, & J. Hicken, 2017
        */
        number norm = sqrt(dLdD[0]*dLdD[0] + dLdD[1]*dLdD[1] + dLdD[2]*dLdD[2]);
        if (iter == 0 || (iter % ndesign) == 0 ) {
            // first step (and every n steps): use only the current gradients
            foreach (i; 0..ndesign) { pk[i] = -dLdD[i]/norm; }
        } else {
            // else let's also use the past gradient knowledge
            number beta_numer = dLdD[0]*dLdD[0] + dLdD[1]*dLdD[1] + dLdD[2]*dLdD[2];
            number beta_denom = dLdD_old[0]*dLdD_old[0] + dLdD_old[1]*dLdD_old[1] + dLdD_old[2]*dLdD_old[2];
            number beta = beta_numer/beta_denom;
            foreach (i; 0..ndesign) { pk[i] = -dLdD[i]/norm + beta*pk_old[i]; }
        }
    } // end compute_search_direction()
    
    void update_shape_params(number a) {
        /*
          Updates the block geometric parameters.
        */
        block.b = dvar0[0] + a*pk[0];
        block.c = dvar0[1] + a*pk[1];
        block.d = dvar0[2] + a*pk[2];
    } // end update_shape_params()
    
    number line_search(number j_0) {
        /*
          Simple backtracking line search algorithm
          Ref.: 
              Numerical Optimization,
              Nocedal, 2006, pg. 37.
        */
        number a = 1.0;
        number rho = 0.5;
        number j_a = double.max;
        number j_prime_0 = dLdD[0]*pk[0] + dLdD[1]*pk[1] + dLdD[2]*pk[2];
        
        while ( j_a > j_0 + mu*a*j_prime_0 ) {
            // evaluate f(x + apk)
            block.compute_geometry();
            solver.integrate_in_time(block);
            j_a = eval_obj_func();
            a *= rho;
            update_shape_params(a);
            while (block.b < 0.1 || block.c < 0.1 || block.d < 0.1) {
                a *= rho;
                update_shape_params(a);   
            }
        } // end while ()

        return a;
    } // end line_search()

} // end Optimizer class

version(complex_numbers) {
    string numerical_gradient_code(string dvar, string idx)
    {
        return `
        block.`~dvar~` += complex(0.0, eps.re);
        block.compute_geometry();
        solver.integrate_in_time(block);
        J_plus = eval_obj_func();
        block.`~dvar~` = dvar0[`~idx~`];
    
        dLdD[`~idx~`] = ((J_plus.im)/eps).re;`;
    }

    string eval_resid_sens_row_code(string dvar, string idx)
    {
        return `
        block.`~dvar~` += complex(0.0, eps.re);
        block.compute_geometry();
        block.apply_boundary_conditions();
        block.compute_fluxes();
        block.update_conserved_quantities(to!number(1.0));
        foreach(cell; block.cells) {
            if (cell.is_ghost) { continue; }
            size_t i = cell.id;
            block.derivative.resid_sens[`~idx~`][i*np+0] = (cell.U.r.im/eps).re;
            block.derivative.resid_sens[`~idx~`][i*np+1] = (cell.U.ru.im/eps).re;
            block.derivative.resid_sens[`~idx~`][i*np+2] = (cell.U.rE.im/eps).re;
        }
        block.`~dvar~` = dvar0[`~idx~`];
        foreach (k; 0..nifaces) { block.ifaces[k].copy_from(block.derivative.ifaces_cpy[k]); }
        foreach (k; 0..nghost+ncells+nghost) { block.cells[k].copy_from(block.derivative.cells_cpy[k]); }
        `;
    }
} else { //real_number
    string numerical_gradient_code(string dvar, string idx)
    {
        return `
        block.`~dvar~` += eps;
        block.compute_geometry();
        solver.integrate_in_time(block);
        J_plus = eval_obj_func();
        block.`~dvar~` = dvar0[`~idx~`];
    
        block.`~dvar~` -= eps;
        block.compute_geometry();
        solver.integrate_in_time(block);
        J_minus = eval_obj_func();
        block.`~dvar~` = dvar0[`~idx~`];
        block.compute_geometry();
   
        dLdD[`~idx~`] = (J_plus-J_minus)/(2.0*eps);`;
    }

    string eval_resid_sens_row_code(string dvar, string idx)
    {
        return `
        block.`~dvar~` += eps;
        block.compute_geometry();
        block.apply_boundary_conditions();
        block.compute_fluxes();
        block.update_conserved_quantities(to!number(1.0));
        foreach(cell; block.cells) {
            if (cell.is_ghost) { continue; }
            size_t i = cell.id;
            block.derivative.resid_sens[`~idx~`][i*np+0] = cell.U.r;
            block.derivative.resid_sens[`~idx~`][i*np+1] = cell.U.ru;
            block.derivative.resid_sens[`~idx~`][i*np+2] = cell.U.rE;
        }
        block.`~dvar~` = dvar0[`~idx~`];
        foreach (k; 0..nifaces) { block.ifaces[k].copy_from(block.derivative.ifaces_cpy[k]); }
        foreach (k; 0..nghost+ncells+nghost) { block.cells[k].copy_from(block.derivative.cells_cpy[k]); }

        block.`~dvar~` -= eps;
        block.compute_geometry();
        block.apply_boundary_conditions();
        block.compute_fluxes();
        block.update_conserved_quantities(to!number(1.0));
        foreach(cell; block.cells) {
            if (cell.is_ghost) { continue; }
            size_t i = cell.id;
            block.derivative.resid_sens[`~idx~`][i*np+0] -= cell.U.r;
            block.derivative.resid_sens[`~idx~`][i*np+1] -= cell.U.ru;
            block.derivative.resid_sens[`~idx~`][i*np+2] -= cell.U.rE;
        }
        block.`~dvar~` = dvar0[`~idx~`];
        foreach (k; 0..nifaces) { block.ifaces[k].copy_from(block.derivative.ifaces_cpy[k]); }
        foreach (k; 0..nghost+ncells+nghost) { block.cells[k].copy_from(block.derivative.cells_cpy[k]); }

        foreach(cell; block.cells) {
            if (cell.is_ghost) { continue; }
            size_t i = cell.id;
            block.derivative.resid_sens[`~idx~`][i*np+0] /= (2.0*eps);
            block.derivative.resid_sens[`~idx~`][i*np+1] /= (2.0*eps);
            block.derivative.resid_sens[`~idx~`][i*np+2] /= (2.0*eps);
        }
        `;
    }
}
