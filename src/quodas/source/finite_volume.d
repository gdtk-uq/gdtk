// finite_volume.d
module finite_volume;
import std.stdio;
import std.conv;
import std.math;
import number;
import complexify;

struct FlowState
{
    number p;      // pressure Pa
    number rho;    // density  kg/m^3
    number vel;    // velocity m/s
} // end FlowState

struct ConservedQuantities
{
    number r;      // mass
    number ru;     // momentum
    number rE;     // total energy
} // end ConservedQuantities

class FVCell {
public:
    size_t id;        // unique id number
    number xpos;      // position
    number dx;        // width
    number vol;       // volume
    FlowState fs;
    ConservedQuantities U;
    FVInterface left_iface;
    FVInterface right_iface;
    bool is_ghost;   // denotes if a cell is a ghost cell
    
    this(size_t id, number xpos, number dx, bool is_ghost = false) {
	this.id = id;
	this.xpos = xpos;
	this.dx = dx;
        this.is_ghost = is_ghost;
    }

    void copy_from(FVCell other) {
        id = other.id;
        xpos = other.xpos;
        dx = other.dx;
        vol = other.vol;
        fs = other.fs;
        U = other.U;
        is_ghost = other.is_ghost;
    } // end copy_from()
    
    void encode_conserved(number gamma) {
        /*
          converts primitive flow state variables to conserved flow state variables.
        */
        number e = fs.p / (fs.rho*(gamma - 1.0)); // internal energy
        number ke = 0.5*fs.vel*fs.vel; // kinetic energy
        U.r = fs.rho;
        U.ru = fs.rho * fs.vel;
        U.rE = fs.rho*(e + ke);
    } // end encode_conserved_variables()

    void decode_conserved(number gamma) {
        /*
          converts conserved flow state variables to primitive flow state variables.
        */
        fs.rho = U.r;
        fs.vel = U.ru/fs.rho;
        number ke = 0.5*fs.vel*fs.vel; // kinetic energy
        number e = U.rE/fs.rho - ke; // internal energy
        fs.p = fs.rho*e*(gamma-1.0);
    } // end decode_conserved_variables()

    number calculate_stable_timestep(number cfl, number gamma) {
        /*
          calculates stable timestep based on the CFL stability condition.
        */
        number a = sqrt(gamma * fs.p/fs.rho); // speed of sound
        number lambda = fabs(fs.vel) + a;
        number dt = cfl*dx/lambda;
        return dt;
    } // end calculate_stable_timestep()
    
    void update_conserved_quantities(number dt) {
        /*
          updates conserved quantities at cell centers given fluxes across interfaces.
          Ref.:
              Transient Hypervelocity Flow in an Axisymmetric Nozzle,
              P.A. Jacobs, ICASE Report No. 91-1, 1991.
         */
        FVInterface fin = left_iface;
        FVInterface fout = right_iface;
        number d_mass;
        number d_momentum;
        number d_energy;
        // update mass
        d_mass = dt/vol * (fin.Fmass*fin.area - fout.Fmass*fout.area);
        U.r = U.r + d_mass;
        // update momentum
        d_momentum = dt/vol * (fin.Fmom*fin.area - fout.Fmom*fout.area+fs.p*(fout.area-fin.area));
        U.ru = U.ru + d_momentum;
        // update total energy
        d_energy = dt/vol * (fin.Fenergy*fin.area - fout.Fenergy*fout.area);
        U.rE = U.rE + d_energy;        
    } // end update_conserved_quantities()
}

class FVInterface {
public:
    size_t id;           // unique id number
    number xpos;         // position
    number area;         // area
    number Fmass;        // mass flux
    number Fmom;         // momentum flux
    number Fenergy;      // energy flux
    FlowState left_fs;   
    FlowState right_fs;  
    FVCell[] cell_cloud; // cells for flowstate reconstruction
                         // assumed order: [L1, L0, R0, R1]
    
    this(size_t id, number xpos) {
	this.id = id;
	this.xpos = xpos;
    }

    void copy_from(FVInterface other) {
        id = other.id;
        xpos = other.xpos;
        area = other.area;
        Fmass = other.Fmass;
        Fmom = other.Fmom;
        Fenergy = other.Fenergy;
        left_fs = other.left_fs;
        right_fs = other.right_fs;
    } // end copy_from()
    
    void reconstruct_flowstate(int interp_order) {
        if (interp_order == 1) { first_order(); }
        else if (interp_order == 2) { second_order(); }
        else {
            string msg = text("interpolation_order is not valid. Select either 1 or 2.");
            throw new Error(msg);
        }
    } // end reconstruct_flowstate()
    
    void first_order() {
        /*
          first order interpolation is a just a copy of the neighbour cell flow states
        */

        FVCell L0 = cell_cloud[1];
        FVCell R0 = cell_cloud[2];

        left_fs.p = L0.fs.p;
        left_fs.rho = L0.fs.rho;
        left_fs.vel = L0.fs.vel;

        right_fs.p = R0.fs.p;
        right_fs.rho = R0.fs.rho;
        right_fs.vel = R0.fs.vel;

    } // end first_order()
    
    void second_order() {
        /* 
           Second order interpolation with modified Van albada limiting
           Ref.:
               Simulation of flow around hypersonic blunt-nosed vehicles for the calibration of air data
               Johnston, University of Queensland, 1999
        */
        
        FVCell L0 = cell_cloud[1];
        FVCell L1 = cell_cloud[0];
        FVCell R0 = cell_cloud[2];
        FVCell R1 = cell_cloud[3];

        number delta_minus;
        number delta_plus;
        number S;
        number eps = 1.0e-12;
        number k = 0.0;
        
        mixin(interp_left_code("p"));
        mixin(interp_left_code("rho"));
        mixin(interp_left_code("vel"));

        mixin(interp_right_code("p"));
        mixin(interp_right_code("rho"));
        mixin(interp_right_code("vel"));
        
    } // end second_order()

}

string interp_left_code(string var)
{
    return `
    delta_minus = (L0.fs.`~var~` - L1.fs.`~var~`)/(0.5*(L0.dx+L1.dx));
    delta_plus = (R0.fs.`~var~` - L0.fs.`~var~`)/(0.5*(R0.dx+L0.dx));
    S = (2.0*delta_plus*delta_minus+eps)/(delta_plus*delta_plus+delta_minus*delta_minus+eps);
    left_fs.`~var~` = L0.fs.`~var~` + (L0.dx*S)/4.0 * ((1.0-S*k)*delta_minus+(1.0+S*k)*delta_plus);`;
}

string interp_right_code(string var)
{
    return `
    delta_minus = (R0.fs.`~var~` - L0.fs.`~var~`)/(0.5*(R0.dx+L0.dx));
    delta_plus = (R1.fs.`~var~` - R0.fs.`~var~`)/(0.5*(R0.dx+R1.dx));
    S = (2.0*delta_plus*delta_minus+eps)/(delta_plus*delta_plus+delta_minus*delta_minus+eps);
    right_fs.`~var~` = R0.fs.`~var~` + (R0.dx*S)/4.0 * ((1.0-S*k)*delta_minus+(1.0+S*k)*delta_plus);`;
}
