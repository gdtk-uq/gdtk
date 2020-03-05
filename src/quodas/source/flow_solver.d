// flow_solver.d
module flowsolve;
import complexify;
import number;
import config;
import block;

class Solver {
private:    
    int _step;
    int _max_step;
    number _sim_time;
    number _max_time;
    number _cfl;
    number _dt;

public:
    this(Params params) {
        this._max_time = params.max_flow_time;
        this._max_step = params.max_flow_steps;
        this._cfl = params.cfl;
    }
    
    void integrate_in_time(Block block) {
        _sim_time = 0.0;
        _step = 0;

        block.set_initial_flowstate_condition();
        block.set_inflow_and_outflow_conditions();
        block.encode_conserved_quantities();	
        // integrate in time
        while ( (_sim_time <= _max_time) && (_step < _max_step) ) {
            block.determine_timestep(_dt, _cfl);
            block.apply_boundary_conditions();
            block.compute_fluxes();
            block.update_conserved_quantities(_dt);
            block.decode_conserved_quantities();
            _sim_time += _dt;
            _step += 1;
            //writef("sim_time: %.8e    dt: %.8e\n", sim_time, dt);
        } // end while loop
    } // end integrate_in_time()

} // end Solver class
