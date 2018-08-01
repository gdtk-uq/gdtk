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
import std.algorithm;
import nm.complex;
import nm.number;

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
        _T_modes_ref = getDouble(L, -1, "T_modes_ref");
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
        lua_pop(L, 1); // dispose of the table
        lua_close(L);
    }
    
    override void opCall(GasState Q, double tInterval,
                         ref double dtChemSuggest, ref double dtThermSuggest, 
                         ref number[] params)
    {
        if (Q.T > 3000.0) {
            // Determine the current number densities.
            number n_e = Q.rho/_mol_masses[2]*Q.massf[2]*_Av; // number density of electrons
            number n_Ar = Q.rho/_mol_masses[0]*Q.massf[0]*_Av; // number density of Ar
            number alpha = n_e/(n_e + n_Ar);
            number orig_n = n_e + n_Ar;
            number n_sum;
            number internal_energy_initial = 3.0/2.0*_Rgas*(Q.T+alpha*Q.T_modes[0])+alpha*_Rgas*_theta_ion;
            //
            // This is a simple euler integration with a number of steps over the full interval.
            int nSteps = min(to!int(tInterval/_chem_dt), 1);
            double dt = tInterval/nSteps;
            foreach (j; 1 .. nSteps) {
                // CHEMISTRY PART ---------------------------------------
                // Determine the current rate coefficients
                number kfA = 1.68e-26*Q.T*sqrt(Q.T)*(_theta_A1star/Q.T+2)*exp(-_theta_A1star/Q.T);
                number kfe = 3.75e-22*Q.T_modes[0]*sqrt(Q.T_modes[0])*
                    (_theta_A1star/Q.T_modes[0]+2)*exp(-_theta_A1star/Q.T_modes[0]);
                number krA = 5.8e-49*(_theta_A1star/Q.T+2)*exp((_theta_ion-_theta_A1star)/Q.T);
                number kre = 1.29e-44*(_theta_A1star/Q.T_modes[0]+2)*exp((_theta_ion-_theta_A1star)/Q.T_modes[0]);
                //
                // Determine the current rate of ionisation.
                // Production rate of electrons due to argon collisions.
                number ne_dot_A = kfA*pow(n_Ar,2) - krA*n_Ar*pow(n_e,2);
                // Production rate of electrons due to electron collisions.
                number ne_dot_e = kfe*n_Ar*n_e - kre*pow(n_e,3);
                number ne_dot = ne_dot_A + ne_dot_e;
                //
                // Determine the new number densities (using a simple Euler update).
                n_e += ne_dot*dt;
                n_Ar -= ne_dot*dt;
                //
                alpha = n_e/(n_e + n_Ar);
                //
                if (alpha <= _ion_tol) {
                    Q.T_modes[0] = _T_modes_ref;
                    Q.u_modes[0] = 3.0/2.0*_Rgas*alpha*Q.T_modes[0] + alpha*_Rgas*_theta_ion;
                } else {
                    Q.u_modes[0] = 3.0/2.0*_Rgas*alpha*Q.T_modes[0] + alpha*_Rgas*_theta_ion -
                        ne_dot_e*_Kb*_theta_ion*dt/Q.rho;
                    Q.T_modes[0] = (Q.u_modes[0]/alpha-_Rgas*_theta_ion)*2.0/3.0/_Rgas;
                }
                //
                Q.u = internal_energy_initial - Q.u_modes[0];
                Q.T = 2.0/3.0*Q.u/_Rgas;
                //
                if (n_e <= 0.0) {
                    // Do not let the number of electrons go negative.
                    // Force the numbers back to something physically realizable.
                    Q.u = internal_energy_initial;
                    Q.u_modes[0] = 0.0;
                    //mass
                    n_Ar = orig_n;
                    n_e = 0.0;
                    alpha = n_e/(n_e + n_Ar);
                    //energy
                    //temperature
                    Q.T =  2.0/3.0*Q.u/_Rgas;
                    Q.T_modes[0] = _T_modes_ref;
                } else {
                    // Positive electron numbers; proceed with relaxation process.
                    n_sum = orig_n/(n_e + n_Ar);
                    n_e = n_e*n_sum;
                    n_Ar = n_Ar*n_sum;
                    // THERMAL RELAXATION PART
                    // find the collision cross sectional area
                    if (alpha > _ion_tol) {
                        // There is no thermal relaxation if below ion tolerance.
                        number Q_ea; 
                        if (Q.T_modes[0] < 1.0e4) {
                            Q_ea = (0.39 - 0.551e-4*Q.T_modes[0] + 0.595e-8*pow(Q.T_modes[0],2))*1.0e-20;
                        } else if (Q.T_modes[0]<1.0e5) {
                            Q_ea = (-0.35 + 0.775e-4*Q.T_modes[0])*1.0e-20;
                        } else {
                            Q_ea = (-0.35 + 0.775e-4*50000)*1.0e-20;
                        }
                        number Q_ei = 1.95e-10*pow(Q.T_modes[0],-2)*log(1.53e8*pow(Q.T_modes[0],3)/(n_e/1.0e6));
                        if (Q_ei < 0.0) { Q_ei = 0.0; }
                        // Find the collision frequencies.
                        // electron-Ar collisions
                        number v_ea = (1-alpha)*Q.rho/_m_Ar*sqrt(8*_Kb*Q.T_modes[0]/to!double(PI)/_m_e)*Q_ea;
                        // electron-Ar+ collisions
                        number v_ei = alpha*Q.rho/_m_Ar*sqrt(8*_Kb*Q.T_modes[0]/to!double(PI)/_m_e)*Q_ei;
                        // update the energy of each state
                        // energy transferred to electron mode through collisions
                        number u_trans_collisions = 3*n_e*_m_e/_m_Ar*(v_ea+v_ei)*_Kb*(Q.T-Q.T_modes[0])*dt/Q.rho;
                        // writeln("u_trans_collisions = ", u_trans_collisions);
                        // Other energy states...
                        // number u_trans_ionisation;
                        // number u_trans_ionisation_heavy;
                        // number u_trans_ionisation_electron;
                        Q.u -= u_trans_collisions;
                        Q.u_modes[0] += u_trans_collisions;
                        // update thermo properties based on energy transfer
                        Q.T = 2.0/3.0*Q.u/_Rgas;
                        Q.T_modes[0] = (Q.u_modes[0]/alpha-_Rgas*_theta_ion)*2.0/3.0/_Rgas;
                    } // end if alpha > ion tol
                } // end if statement regarding number density of electrons
            } // end foreach j
            //
            // Convert back to mass fractions.
            // Note that density has not changed since we are updating an isolated system.
            Q.massf[0] = n_Ar/_Av/Q.rho*_mol_masses[0];
            // Number density of Argon+ is the same as electron number density.
            Q.massf[1] = n_e/_Av/Q.rho*_mol_masses[1];
            Q.massf[2] = n_e/_Av/Q.rho*_mol_masses[2];
            //
            // Since the internal energy and density in the (isolated) reactor is fixed,
            // we need to evaluate the new temperature, pressure, etc.
            _gmodel.update_thermo_from_rhou(Q);
            _gmodel.update_sound_speed(Q);  
        } // end if Q.T > 3000.0
    } // end opCall()
    
private:
    // Reaction rate constant
    double[] _mol_masses;
    double _ion_tol;
    double _T_modes_ref;
    double _chem_dt;
    double _m_Ar; //mass of argon (kg)
    double _m_e; //mass of electron (kg)
    double _Kb;
    double _Av;
    double _Rgas;
    double _theta_ion;
    double _theta_A1star;

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

        //pre-shock conditions
        double rho1 = 0.03334;
        double T1 = 260.4;
        double p1 = 180.6;
        double u1 = 5.7e3;
        double M1 = 18.96;

        //inital gas properties
        auto GS = new GasState(3, 1);
        GS.p = p1*449.6;
        GS.T = T1*113.33;
        GS.T_modes[0] = 10000;
        GS.massf[0] = 1.0; GS.massf[1] = 0.0; GS.massf[2] = 0.0;
        double vel = 1436; // m/s
        double x = 0.0;

        double e_new;

        //update the gas model based on the properties specified above
        gm.update_thermo_from_pT(GS); // update gas state
        gm.update_sound_speed(GS); // upate sound speed (not necessary)

        // some time stepping information
        double maxtime = 4.0e-6; // max time for the simulation
        double dt = 1.0e-10;//1.312e-11; // time step for chemistry update
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
        //
        // [TODO] Daniel, we need some assert statements here to check particular values.
        //
    } // end main()
} // end two_temperature_argon_kinetics_test



