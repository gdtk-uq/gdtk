/**
 * two_temperature_argon_kinetics.d
 *
 * Two-temperature reacting argon based off of:
 * "Quasi-One-Dimensional, Nonequilibrium Gas Dynamics of Partially Ionised Two-Temperature Argon"
 * Martin I. Hoffert and Hwachii Lien
 * 
 * The Physics of Fluids 10, 1769 (1967); doi 10.1063/1.1762356
 *
 *
 * Authors: Daniel Smith and Rory Kelly
 * Version: 19-July-2017: initial cut.
 */

module kinetics.two_temperature_argon_kinetics;

import std.stdio : writeln;
import std.format;
import std.math;
import std.conv : to;
import nm.complex;
import nm.number;
import nm.bbla;

import gas;
import util.lua;
import util.lua_service;
import kinetics.thermochemical_reactor;

final class UpdateArgonFrac : ThermochemicalReactor {
    
    this(string fname, GasModel gmodel)
    {
        super(gmodel); // hang on to a reference to the gas model
        // We need to pick a number of pieces out of the gas-model file, again.
        // Although they exist in the GasModel object, they are private.
        auto L = init_lua_State();
        doLuaFile(L, fname);
        lua_getglobal(L, "TwoTemperatureReactingArgon");
        _mol_masses.length = 3;
        _mol_masses[0] = 39.948e-3; // Units are kg/mol
        _mol_masses[2] = 5.485799e-7; // Units are kg/mol
        _mol_masses[1] = _mol_masses[0] - _mol_masses[2]; // Units are kg/mol
        _m_Ar = 6.6335209e-26; //mass of argon (kg)
        _m_e = 9.10938e-31; //mass of electron (kg)
        _Kb = Boltzmann_constant;
        _Av = Avogadro_number;
        _Rgas = 208.0;
        _theta_ion = 183100.0;
        _theta_A1star = 135300.0;
        _ion_tol = getDouble(L, -1, "ion_tol");
        _chem_dt = getDouble(L, -1, "chem_dt");
        _integration_method = getString(L, -1, "integration_method");
        _Newton_Raphson_tol = getDouble(L, -1, "Newton_Raphson_tol");
        lua_close(L);
    }
    
    Matrix!number F(Matrix!number y, GasState Q)
    {
        auto y_dash = new Matrix!number([[0.0],
                                         [0.0]]);
        number Q_ea;
        number Q_ei;

        number v_ea;
        number v_ei;

        //unpack dependant variables
        number n_e = y[0,0];
        number n_Ar = _n_total-n_e;        
        Q.u = y[1,0];
        Q.u_modes[0] = _u_total-Q.u;

        //Convenient variables
        number alpha = n_e/(_n_total);
        number Te;
        number T = Q.T;

        //reconstruct the flow state from state vector
        if (n_e <= 0.0) {
            // Do not let the number of electrons go negative.
            // Force the numbers back to something physically realizable.
            Q.u = _u_total;
            Q.u_modes[0] = 0.0;
            //mass
            n_Ar = _n_total;
            n_e = 0.0;
            //energy
            //temperature
            Q.T =  2.0/3.0*Q.u/_Rgas; T = Q.T;
            Q.T_modes[0] = Q.T;
        } else if (alpha < _ion_tol) {
            Te = Q.T;
        } else {
            Te = Q.T_modes[0];
        }

        if (Q.u_modes[0] == 0.0) {Te = T;}
        
        //Must update the mass fractions before updating the gas state from rho and u...
        Q.massf[0] = n_Ar/_Av/Q.rho*_mol_masses[0];
        Q.massf[1] = n_e/_Av/Q.rho*_mol_masses[1];
        Q.massf[2] = n_e/_Av/Q.rho*_mol_masses[2];// Number density of Argon+ is the same as electron number density.
        _gmodel.update_thermo_from_rhou(Q);
        _gmodel.update_sound_speed(Q);

        //=====================================================================
        //rate constants
        number kfA = 1.68e-26*T*sqrt(T)*(_theta_A1star/T+2)*exp(-_theta_A1star/T);
        number kfe = 3.75e-22*Te*sqrt(Te)*(_theta_A1star/Te+2)*exp(-_theta_A1star/Te);
        number krA = 5.8e-49*(_theta_A1star/T+2)*exp((_theta_ion-_theta_A1star)/T);
        number kre = 1.29e-44*(_theta_A1star/Te+2)*exp((_theta_ion-_theta_A1star)/Te);

        //determine the current rate of ionisation
        number n_dot_A = kfA*pow(n_Ar,2) - krA*n_Ar*pow(n_e,2);         //production rate of electrons due to argon collisions
        number n_dot_e = kfe*n_Ar*n_e - kre*pow(n_e,3);                 //production rate of electrons due to electron collisions
        number n_dot = n_dot_A + n_dot_e;

        y_dash[0,0] = n_dot;

        //Energy moving to electronic mode due to reactions:
        number alpha_dot = n_dot/_n_total;
        number u_dot_reac = 3.0/2.0*_Rgas*alpha_dot*Te+alpha_dot*_Rgas*_theta_ion - n_dot_e*_Kb*_theta_ion/Q.rho;

        //=====================================================================
        //now need temperatures...
        if (alpha > _ion_tol) {

            //calculate Qea
            if (Te < 1.0e4) {
                    Q_ea = (0.39 - 0.551e-4*Te + 0.595e-8*pow(Te,2))*1.0e-20;
            } else if (Te<1.0e5) {
                    Q_ea = (-0.35 + 0.775e-4*Te)*1.0e-20;
            } else {
                    Q_ea = (-0.35 + 0.775e-4*50000)*1.0e-20;
            }

            //Calculate Q_ei
            Q_ei = 1.95e-10*pow(Te,-2)*log(1.53e8*pow(Te,3)/(n_e/1.0e6));
            if (Q_ei < 0.0) {Q_ei = 0.0;}

            //find the collision frequencies
            v_ea = (1.0-alpha)*Q.rho/_m_Ar*sqrt(8*_Kb*Te/to!number(PI)/_m_e)*Q_ea; //electron-Ar collisions
            v_ei = alpha*Q.rho/_m_Ar*sqrt(8*_Kb*Te/to!number(PI)/_m_e)*Q_ei; //electron-Ar+ collisions

            alpha_dot = n_dot/(n_e+n_Ar);
            number u_dot = 3*n_e*_m_e/_m_Ar*(v_ea+v_ei)*_Kb*(T-Te)/Q.rho;// energy transferred to electron mode through collisions
            y_dash[1,0] = - u_dot - u_dot_reac;
        } else {
            y_dash[1,0] = 0.0 - u_dot_reac;            
        }
        return y_dash;
    }

    Matrix!number BackEuler_F(Matrix!number y, Matrix!number y_prev, GasState Q)
    {
        auto output = y_prev - y + _chem_dt*F(y, Q);
        return output;
    }

    Matrix!number Jacobian(Matrix!number y, Matrix!number y_prev, Matrix!number h, GasState Q){

        auto J = new Matrix!number([[0.0,0.0],
                                    [0.0,0.0]]);

        auto yph0 = new Matrix!number([[y[0,0] + h[0,0]],[y[1,0]]]); //y plus increment in index 0
        auto col1 = (BackEuler_F(yph0,y_prev, Q) - BackEuler_F(y, y_prev, Q))/h[0,0];
        auto yph1 = new Matrix!number([[y[0,0]],[y[1,0] + h[1,0]]]); //y plus increment in index 1
        auto col2 = (BackEuler_F(yph1,y_prev, Q) - BackEuler_F(y, y_prev, Q))/h[1,0];

        J[0,0] = col1[0,0]; J[1,0] = col1[1,0]; J[0,1] = col2[0,0]; J[1,1] = col2[1,0];

        return J;
    }

    @nogc
    override void opCall(GasState Q, double tInterval,
                         ref double dtChemSuggest, ref double dtThermSuggest, 
                         ref number[] params)
    {
        debug {
        // Implement in debug context because of @nogc requirement.
        //
        if (Q.T > 3000.0) {
            double chem_dt_start = _chem_dt;
            int NumberSteps = to!int(tInterval/_chem_dt);
            if (NumberSteps < 1) {NumberSteps = 1;}
            _chem_dt = tInterval/NumberSteps;
            //writeln("Number of Steps ", NumberSteps);
            // Determine the current number densities.
            number n_e = Q.rho/_mol_masses[2]*Q.massf[2]*_Av; // number density of electrons
            number n_Ar = Q.rho/_mol_masses[0]*Q.massf[0]*_Av; // number density of Ar
            number alpha = n_e/(n_e + n_Ar);
            _n_total = n_e + n_Ar;
            _u_total = 3.0/2.0*_Rgas*(Q.T+alpha*Q.T_modes[0])+alpha*_Rgas*_theta_ion;
            //writeln("_u_total = ", _u_total);

            auto y = new Matrix!number([[n_e],
                                        [Q.u]]);
            auto y_prev = y;
            number norm_error = 1.0e50;
            number prev_norm_error;
            
            //Initialise variables for RK4
            auto k1 = zeros!number(2,1);
            auto k2_in = zeros!number(2,1);
            auto k2 = zeros!number(2,1);
            auto k3_in = zeros!number(2,1);
            auto k3 = zeros!number(2,1);
            auto k4_in = zeros!number(2,1);
            auto k4 = zeros!number(2,1);
            auto h = new Matrix!number([[1.0e10],[1.0e-5]]); //working: auto h = new Matrix!number([[1.0e10],[1.0e0]]); 

            foreach (n; 1 .. NumberSteps) {
                //writeln("y = ", y);
                if (_integration_method == "Forward_Euler") {
                    y = y + _chem_dt*F(y, Q);
                } else if (_integration_method == "Backward_Euler"){
                    if (y[0,0] < 0.0) {
                        y = y + _chem_dt*F(y, Q);
                    } else {
                        do {
                            //writeln(" ");
                            //prev_norm_error = norm_error;
                            auto J = Jacobian(y,y_prev,h,Q); 
                            auto J_inv = inverse!number(J);
                            y = y - dot!number(J_inv,BackEuler_F(y,y_prev,Q));
                            auto error = BackEuler_F(y,y_prev,Q);
                            //writeln("error = ", error);
                            norm_error = sqrt(pow(error[0,0]/1.0e22,2) + pow(error[1,0]/1.0e7,2));
                            //writeln("norm_error = ", norm_error);
                        } while (norm_error > _Newton_Raphson_tol);//while (norm_error < prev_norm_error); // 
                    }

                    y_prev = y;
                } else if (_integration_method == "RK4"){
                    k1 = _chem_dt*F(y, Q);

                    k2_in = y + k1/2.0;
                    k2 = _chem_dt*F(k2_in, Q);

                    k3_in = y+k2/2.0;
                    k3 = _chem_dt*F(k3_in, Q);

                    k4_in = y+k3;
                    k4 = _chem_dt*F(k4_in, Q);

                    y = y + 1.0/6.0*(k1+2.0*k2+2.0*k3+k4);
                }
            }

            //Unpack the results
            n_e = y[0,0];
            Q.u = y[1,0];

            //reconstruct the flow state
            n_Ar = _n_total - n_e;
            Q.u_modes[0] = _u_total-Q.u;    

            if (n_e <= 0.0) {
                // Do not let the number of electrons go negative.
                // Force the numbers back to something physically realizable.
                Q.u = _u_total;
                Q.u_modes[0] = 0.0;
                //mass
                n_Ar = _n_total;
                n_e = 0.0;
                alpha = n_e/(n_e + n_Ar);
                //energy
                //temperature
                Q.T =  2.0/3.0*Q.u/_Rgas;
                Q.T_modes[0] = Q.T;
            }
            
            //Must update the mass fractions before updating the gas state from rho and u...
            Q.massf[0] = n_Ar/_Av/Q.rho*_mol_masses[0];
            Q.massf[1] = n_e/_Av/Q.rho*_mol_masses[1];
            Q.massf[2] = n_e/_Av/Q.rho*_mol_masses[2];// Number density of Argon+ is the same as electron number density.       
            _gmodel.update_thermo_from_rhou(Q);
            _gmodel.update_sound_speed(Q);

            _chem_dt = chem_dt_start; // return _chem_dt back to its original value\
        }
        } else {
            // We are not in a debug build context, so...
            throw new Error("2T argon only available in debug build.");
        }
    } // end opCall()
    
private:
    double[] _mol_masses;
    double _ion_tol;
    double _chem_dt;
    double _m_Ar; //mass of argon (kg)
    double _m_e; //mass of electron (kg)
    double _Kb;
    double _Av;
    double _Rgas;
    double _theta_ion;
    double _theta_A1star;
    string _integration_method;
    double _Newton_Raphson_tol;
public:
    number _n_total;
    number _u_total;

} // end class UpdateAB


version(two_temperature_argon_kinetics_test) {
    import std.stdio;
    import util.msg_service;
    import std.math : approxEqual;
    import gas.two_temperature_reacting_argon;
    void main() {
        lua_State* L = init_lua_State();
        doLuaFile(L, "sample-input/two-temperature-reacting-argon-model.lua");
        auto gm = new TwoTemperatureReactingArgon(L);
        lua_close(L);
        auto gd = new GasState(3, 1);
        gd.p = 1.0e5;
        gd.T = 300.0;
        gd.T_modes[0] = 300;
        gd.massf[0] = 1.0; gd.massf[1] = 0.0; gd.massf[2] = 0.0;

        auto reactor = new UpdateArgonFrac("sample-input/two-temperature-reacting-argon-model.lua", gm);
        double[] params;
        double dtThermSuggest;
        double dtSuggest;

        //need molar masses to determine alpha
        double[3] _mol_masses;
        _mol_masses[0] = 39.948e-3; // Units are kg/mol
        _mol_masses[2] = 5.485799e-7; // Units are kg/mol
        _mol_masses[1] = _mol_masses[0] - _mol_masses[2]; // Units are kg/mol
        double theta_ion = 183100.0;
        double alpha;
        double new_vel; 

        //This is my own condition (Daniel Smith Condition 1 Experimental Campaign 1)
        //pre-shock conditions
        //==============================
        //double rho1 = 0.03334;
        //double T1 = 260.4;
        //double p1 = 180.6;
        //double u1 = 5.7e3;
        //double M1 = 18.96;

        ////inital gas properties
        //auto GS = new GasState(3, 1);
        //GS.p = p1*449.6;
        //GS.T = T1*113.33;
        //GS.T_modes[0] = 10000;
        //GS.massf[0] = 1.0; GS.massf[1] = 0.0; GS.massf[2] = 0.0;
        //double vel = 1436; // m/s
        //double x = 0.0;
        //==============================

        //This is the benchmark case presented in the paper
        //pre-shock conditions
        //==============================
        double rho1 = 0.0213657;
        double T1 = 300;
        double p1 = 1333.22;
        double u1 = 6.1e3;
        double M1 = 18.91;

        //inital gas properties
        auto GS = new GasState(3, 1);
        GS.p = p1*446.735124;
        GS.T = T1* 112.620756;
        GS.T_modes[0] = GS.T;//10000;
        GS.massf[0] = 1.0; GS.massf[1] = 0.0; GS.massf[2] = 0.0;
        double vel = 1548.98; // m/s
        double x = 0.0;
        //==============================

        double e_new;

        //update the gas model based on the properties specified above
        gm.update_thermo_from_pT(GS); // update gas state
        gm.update_sound_speed(GS); // upate sound speed (not necessary)

        // some time stepping information
        double maxtime = 4.0e-6; // max time for the simulation
        double dt = 1.0e-9;//1.312e-11; // time step for chemistry update
        int maxsteps = to!int(maxtime/dt + 1); // number of steps in which to iterate through
        int printfreq = maxsteps/10;
        int writefreq = maxsteps/1000;
        double Ru = R_universal;
        double m_Ar = 6.6335209e-26; //mass of argon (kg)        
        double M_Ar = Avogadro_number*m_Ar;            

        // initialise the storage arrays
        double[] t_list;
        double[] x_list;
        double[] T_list;
        double[] T_modes_list;
        double[] P_list;
        double[] u_list;
        double[] alpha_list;

        //initial values for storage arrays
        t_list ~= 0.0;
        x_list ~= x;
        T_list ~= GS.T;
        T_modes_list ~= GS.T_modes[0];
        P_list ~= GS.p;
        u_list ~= vel;
        alpha_list ~= 0.0;

        //main for loop for chemistry
        foreach (i; 1 .. maxsteps) {
            // perform the chemistry update
            reactor.opCall(GS,dt,dtSuggest,dtThermSuggest,params);
            // new alpha
            alpha = (GS.massf[2]/_mol_masses[2]) / ((GS.massf[2]/_mol_masses[2])+(GS.massf[0]/_mol_masses[0]));
            // update x position
            new_vel = u1/8.*(5 + 3.0/pow(M1,2) - sqrt(9*pow(1-1.0/pow(M1,2),2) + 
                                                      96*alpha/5./pow(M1,2)*theta_ion/T1));
            //writeln("new vel = ", new_vel);
            x += (vel + new_vel)/2*dt;
            vel = new_vel;
            //update rho, P and u.
            GS.rho = rho1*u1/vel;
            GS.p = rho1*((Ru/M_Ar)*T1 + pow(u1,2)) - GS.rho*pow(vel,2);
            //new internal energy
            e_new = 5*(Ru/M_Ar)*T1/2. + pow(u1,2)/2 - pow(vel,2)/2 - GS.p/GS.rho;
            //give all of the new energy to the heavy particles
            GS.u = GS.u + (e_new - (GS.u + GS.u_modes[0]));
            //update the temperature of the heavy paritcles based on this
            gm.update_thermo_from_rhou(GS); 
            //
            if (writefreq == 0) {writefreq = 1;}
            if ((i % writefreq) == 0) { // only save 1000 points
                t_list ~= i*dt;
                x_list ~= x;
                T_list ~= GS.T;
                T_modes_list ~= GS.T_modes[0];
                P_list ~= GS.p;
                u_list ~= vel;
                alpha_list ~= alpha;
            }
        } // end foreach i

        // writeln("writing to Data file... please wait");
        double[][] collateddata = [t_list,x_list,T_list,T_modes_list,P_list,u_list,alpha_list];
        //
        File file = File("two_temperature_argon_kinetics_test_results.data","w");
        foreach (i; 0..t_list.length) {
            file.writeln(collateddata[0][i], " ", collateddata[1][i], " ", collateddata[2][i], " ",
                         collateddata[3][i], " ", collateddata[4][i], " ", collateddata[5][i], " ",
                         collateddata[6][i]);
        }
    assert(approxEqual(GS.T, 14730, 1.0e2), failedUnitTest());
    assert(approxEqual(GS.p, 705738, 1.0e2), failedUnitTest());
    assert(approxEqual(vel, 700.095, 1.0e0), failedUnitTest());
    } // end main()
} // end two_temperature_argon_kinetics_test



