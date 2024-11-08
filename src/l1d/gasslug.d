// gasslug.d for the Lagrangian 1D Gas Dynamics, also known as L1d4.
//
// The GasSlug is the principal dynamic component of a simulation.
//
// PA Jacobs
// 2020-04-08
//
module gasslug;

import std.conv;
import std.stdio;
import std.string;
import std.json;
import std.format;
import std.range;
import std.math;
import std.algorithm;

import util.json_helper;
import geom;
import gas;
import kinetics;
import gasdyn.gasflow;
import config;
import lcell;
import piston;
import valve;
import endcondition;
import simcore; // has the core data arrays
import misc;

class GasSlug {
public:
    size_t indx;
    string label;
    size_t gmodel_id;
    GasModel gmodel;
    size_t ncells;
    int viscous_effects;
    bool adiabatic;
    int ecL_id;
    int ecR_id;
    EndCondition ecL;
    EndCondition ecR;
    size_t nhcells;
    size_t[] hcells;

    LFace[] faces;
    LCell[] cells;

    this(size_t indx, JSONValue jsonData)
    {
        if (L1dConfig.verbosity_level >= 3) {
            writeln("construct slug[", indx, "] from json=", jsonData);
        }
        this.indx = indx;
        label = getJSONstring(jsonData, "label", "");
        gmodel_id = getJSONint(jsonData, "gmodel_id", 0);
        gmodel = gmodels[gmodel_id];
        ncells = getJSONint(jsonData, "ncells", 0);
        viscous_effects = getJSONint(jsonData, "viscous_effects", 0);
        adiabatic = getJSONbool(jsonData, "adiabatic", false);
        ecL_id = getJSONint(jsonData, "ecL_id", -1);
        ecR_id = getJSONint(jsonData, "ecR_id", -1);
        nhcells = getJSONint(jsonData, "nhcells", 1);
        hcells = to!(size_t[])(getJSONintarray(jsonData, "hcells", [0,]));
        if (L1dConfig.verbosity_level >= 1) {
            writeln("GasSlug[", indx, "]:");
            writefln("  label= \"%s\"", label);
            writeln("  gmodel_id= ", gmodel_id);
            writeln("  ncells= ", ncells);
            writeln("  viscous_effects= ", viscous_effects);
            writeln("  adiabatic= ", adiabatic);
            writeln("  ecL_id= ", ecL_id);
            writeln("  ecR_id= ", ecR_id);
            writeln("  hcells= ", hcells);
        }
        //
        foreach (i; 0 .. ncells+1) { faces ~= new LFace(); }
        foreach (i; 0 .. ncells) { cells ~= new LCell(gmodel); }
        //
        // Private workspace for quadratic reconstruction of flow data.
        gsL = GasState(gmodel);
        gsR = GasState(gmodel);
    } // end constructor

    void read_face_data(File fp, int tindx)
    {
        skip_to_data_at_tindx(fp, tindx);
        foreach (i; 0 .. ncells+1) {
            string txt = fp.readln().chomp();
            txt.formattedRead!"%e %e"(faces[i].x, faces[i].area);
        }
    } // end read_face_data

    void write_face_data(File fp, int tindx, bool write_header)
    {
        if (write_header) { fp.writeln("#   x   area"); }
        fp.writeln(format("# tindx %d", tindx));
        foreach (i; 0 .. ncells+1) {
            fp.writeln(format("%e %e", faces[i].x, faces[i].area));
        }
        fp.writeln("# end");
    } // end write_face_data()

    void read_cell_data(File fp, int tindx)
    {
        skip_to_data_at_tindx(fp, tindx);
        int nsp = gmodel.n_species;
        int nmodes = gmodel.n_modes;
        foreach (j; 0 .. ncells) {
            LCell c = cells[j];
            string txt = fp.readln().chomp();
            string[] items = txt.split();
            int k = 0;
            c.xmid = to!double(items[k]); k++;
            c.volume = to!double(items[k]); k++;
            c.vel = to!double(items[k]); k++;
            c.L_bar = to!double(items[k]); k++;
            c.gas.rho = to!double(items[k]); k++;
            c.gas.p = to!double(items[k]); k++;
            c.gas.T = to!double(items[k]); k++;
            c.gas.u = to!double(items[k]); k++;
            c.gas.a = to!double(items[k]); k++;
            c.shear_stress = to!double(items[k]); k++;
            c.heat_flux = to!double(items[k]); k++;
            foreach (i; 0 .. nsp) {
                c.gas.massf[i] = to!double(items[k]); k++;
            }
            if (nsp > 1) { c.dt_chem = to!double(items[k]); k++; }
            foreach (i; 0 .. nmodes) {
                c.gas.T_modes[i] = to!double(items[k]); k++;
                c.gas.u_modes[i] = to!double(items[k]); k++;
            }
            if (nmodes > 0) { c.dt_therm = to!double(items[k]); k++; }
            gmodel.update_thermo_from_pT(c.gas);
            gmodel.update_sound_speed(c.gas);
            gmodel.update_trans_coeffs(c.gas);
        }
    } // end read_cell_data()

    void write_cell_data(File fp, int tindx, bool write_header)
    {
        int nsp = gmodel.n_species;
        int nmodes = gmodel.n_modes;
        if (write_header) {
            fp.write("# xmid  volume  vel  L_bar  rho  p  T  u  a");
            fp.write("  shear_stress  heat_flux");
            foreach (i; 0 .. nsp) { fp.write(format("  massf[%d]", i)); }
            if (nsp > 1) { fp.write("  dt_chem"); }
            foreach (i; 0 .. nmodes) {
                fp.write(format("  T_modes[%d]  u_modes[%d]", i, i));
            }
            if (nmodes > 0) { fp.write("  dt_therm"); }
            fp.write("\n");
        }
        fp.writeln(format("# tindx %d", tindx));
        foreach (j; 0 .. ncells) {
            LCell c = cells[j];
            fp.write(format("%e %e %e %e", c.xmid, c.volume, c.vel, c.L_bar));
            fp.write(format(" %e %e %e %e", c.gas.rho, c.gas.p, c.gas.T, c.gas.u));
            fp.write(format(" %e %e %e", c.gas.a, c.shear_stress, c.heat_flux));
            foreach (i; 0 .. nsp) { fp.write(format(" %e", c.gas.massf[i])); }
            if (nsp > 1) { fp.write(format(" %e", c.dt_chem)); }
            foreach (i; 0 .. nmodes) {
                fp.write(format(" %e %e", c.gas.T_modes[i], c.gas.u_modes[i]));
            }
            if (nmodes > 0) { fp.write(format(" %e", c.dt_therm)); }
            fp.write("\n");
        }
        fp.writeln("# end");
    } // end write_cell_data()

    void write_history_loc_data(File fp, double t, double x)
    {
        // Write a subset of the data that will be common to all slugs,
        // no matter what gas model is associated with this slug.
        // It may be that nothing is written, if there is not a cell
        // that covers the x-location.
        foreach (j; 0 .. ncells) {
            if ((x >= faces[j].x) && (x <= faces[j+1].x)) {
                LCell c = cells[j];
                fp.write(format("%e %e %e", t, c.vel, c.L_bar));
                fp.write(format(" %e %e %e %e", c.gas.rho, c.gas.p, c.gas.T, c.gas.u));
                fp.write(format(" %e %e %e", c.gas.a, c.shear_stress, c.heat_flux));
                double[] massf; massf.length = overall_species_count;
                foreach (ref mf; massf) { mf = 0.0; }
                foreach (i, mf; c.gas.massf) { massf[overall_species_index[gmodel_id][i]] = mf; }
                foreach (mf; massf) { fp.write(format(" %e", mf)); }
                double[] Tmodes; Tmodes.length = overall_modes_count;
                foreach (ref tm; Tmodes) { tm = 0.0; }
                foreach (i, tm; c.gas.T_modes) { Tmodes[overall_modes_index[gmodel_id][i]] = tm; }
                foreach (tm; Tmodes) { fp.write(format(" %e", tm)); }
                fp.write("\n");
                break;
            }
        }
    } // end write_history_loc_data()

    @nogc @property
    double energy()
    {
        double e = 0.0;
        foreach (c; cells) {
            e += c.volume*c.gas.rho*(gmodel.internal_energy(c.gas) + 0.5*c.vel*c.vel);
        }
        return e;
    }

    @nogc
    void compute_areas_and_volumes()
    {
        double[6] daKT;
        foreach (f; faces) {
            daKT = tube1.eval(f.x);
            f.area = daKT[1];
        }
        foreach (i; 0 .. cells.length) {
            double xL = faces[i].x;
            double xR = faces[i+1].x;
            LCell c = cells[i];
            c.L = xR-xL;
            if (c.L <= 0.0) {
                string msg = "Oops, adjacent faces have crossed over.";
                debug { msg ~= format(" xL=%g, xR=%g", xL, xR); }
                throw new Exception(msg);
            }
            c.volume = 0.5*(faces[i].area+faces[i+1].area)*(c.L);
            c.xmid = 0.5*(xR+xL);
            daKT = tube1.eval(c.xmid);
            c.D = daKT[0];
            c.K_over_L = daKT[2];
            c.Twall = daKT[3];
            c.vf = daKT[4];
            c.htcf = daKT[5];
        }
        return;
    } // end compute_areas_and_volumes()

    @nogc
    void encode_conserved()
    {
        foreach (c; cells) { c.encode_conserved(gmodel); }
        return;
    }

    @nogc
    void decode_conserved()
    {
        foreach (c; cells) { c.decode_conserved(gmodel); }
        return;
    }

    @nogc
    void record_state()
    {
        foreach (f; faces) { f.record_state(); }
        foreach (c; cells) { c.record_state(); }
        return;
    }

    @nogc
    void restore_state()
    {
        foreach (f; faces) { f.restore_state(); }
        foreach (c; cells) { c.restore_state(gmodel); }
        compute_areas_and_volumes();
        foreach (c; cells) { c.decode_conserved(gmodel); }
        return;
    }

    @nogc
    void time_derivatives(int level, double t)
    {
        // For the interface between each of the control-mass cells,
        // compute face motion as Riemann subproblems.
        // There will be a single pressure and a single velocity
        // as a result of solving the Riemann problem.
        foreach (i, f; faces) {
            // Need to consider the end conditions.
            if (i == 0) {
                // Left-most face.
                if (ecL) {
                    if (cast(Diaphragm) ecL) {
                        auto dia = cast(Diaphragm) ecL;
                        if (dia.state == DiaphragmState.open) {
                            // Slugs interact as per GasInterface.
                            LCell cR = cells[0];
                            auto slugL =  dia.slugL;
                            LCell cL = (dia.slugL_end == End.L) ?
                                slugL.cells[0] : slugL.cells[$-1];
                            lrivp(cL.gas, cR.gas, cL.vel, cR.vel, slugL.gmodel, gmodel, f.dxdt[level], f.p);
                        } else {
                            // Diaphragm acts as a fixed wall.
                            LCell cR = cells[0];
                            piston_at_left(cR.gas, cR.vel, gmodel, 0.0, f.p);
                            f.dxdt[level] = 0.0;
                        }
                    }
                    if (cast(GasInterface) ecL) {
                        auto my_ecL = cast(GasInterface) ecL;
                        LCell cR = cells[0];
                        auto slugL =  my_ecL.slugL;
                        LCell cL = (my_ecL.slugL_end == End.L) ?
                            slugL.cells[0] : slugL.cells[$-1];
                        lrivp(cL.gas, cR.gas, cL.vel, cR.vel, slugL.gmodel, gmodel, f.dxdt[level], f.p);
                    }
                    if (cast(FreeEnd) ecL) {
                        LCell cR = cells[0];
                        f.p = cR.gas.p;
                        f.dxdt[level] = cR.vel;
                    }
                    if (cast(VelocityEnd) ecL) {
                        auto my_ecL = cast(VelocityEnd) ecL;
                        LCell cR = cells[0];
                        piston_at_left(cR.gas, cR.vel, gmodel, my_ecL.vel, f.p);
                        f.dxdt[level] = my_ecL.vel;
                    }
                    if (cast(PistonFace) ecL) {
                        auto my_ecL = cast(PistonFace) ecL;
                        Piston piston = my_ecL.pistonL;
                        LCell cR = cells[0];
                        piston_at_left(cR.gas, cR.vel, gmodel, piston.vel, f.p);
                        f.dxdt[level] = piston.vel;
                    }
                } else {
                    throw new Error("Left end of gas slug does not have an end condition.");
                }
            } else if (i+1 == faces.length) {
                // Right-most face.
                if (ecR) {
                    if (cast(Diaphragm) ecR) {
                        auto dia = cast(Diaphragm) ecR;
                        if (dia.state == DiaphragmState.open) {
                            // Slugs interact as per GasInterface.
                            LCell cL = cells[$-1];
                            auto slugR =  dia.slugR;
                            LCell cR = (dia.slugR_end == End.L) ?
                                slugR.cells[0] : slugR.cells[$-1];
                            lrivp(cL.gas, cR.gas, cL.vel, cR.vel, gmodel, slugR.gmodel, f.dxdt[level], f.p);
                        } else {
                            // Diaphragm acts as a fixed wall.
                            LCell cL = cells[$-1];
                            piston_at_right(cL.gas, cL.vel, gmodel, 0.0, f.p);
                            f.dxdt[level] = 0.0;
                        }
                    }
                    if (cast(GasInterface) ecR) {
                        auto my_ecR = cast(GasInterface) ecR;
                        LCell cL = cells[$-1];
                        auto slugR =  my_ecR.slugR;
                        LCell cR = (my_ecR.slugR_end == End.L) ?
                            slugR.cells[0] : slugR.cells[$-1];
                        lrivp(cL.gas, cR.gas, cL.vel, cR.vel, gmodel, slugR.gmodel, f.dxdt[level], f.p);
                    }
                    if (cast(FreeEnd) ecR) {
                        LCell cL = cells[$-1];
                        f.p = cL.gas.p;
                        f.dxdt[level] = cL.vel;
                    }
                    if (cast(VelocityEnd) ecR) {
                        auto my_ecR = cast(VelocityEnd) ecR;
                        LCell cL = cells[$-1];
                        piston_at_right(cL.gas, cL.vel, gmodel, my_ecR.vel, f.p);
                        f.dxdt[level] = my_ecR.vel;
                    }
                    if (cast(PistonFace) ecR) {
                        auto my_ecR = cast(PistonFace) ecR;
                        Piston piston = my_ecR.pistonR;
                        LCell cL = cells[$-1];
                        piston_at_right(cL.gas, cL.vel, gmodel, piston.vel, f.p);
                        f.dxdt[level] = piston.vel;
                    }
                } else {
                    throw new Error("Right end of gas slug does not have an end condition.");
                }
            } else {
                // Interior face i is between cell i-1 and i.
                if (L1dConfig.x_order == 2) {
                    // Do quadratic reconstruction.
                    if (i == 1) {
                        // Only one cell to the left.
                        LCell cL0 = cells[i-1];
                        LCell cR0 = cells[i];
                        LCell cR1 = cells[i+1];
                        velL = cL0.vel;
                        gsL.copy_values_from(cL0.gas);
                        interpR_prepare(cL0.L, cR0.L, cR1.L);
                        interpR(cL0, cR0, cR1, gsR, velR);
                    } else if (i == faces.length-2) {
                        // Only one cell to the right.
                        LCell cL1 = cells[i-2];
                        LCell cL0 = cells[i-1];
                        LCell cR0 = cells[i];
                        interpL_prepare(cL1.L, cL0.L, cR0.L);
                        interpL(cL1, cL0, cR0, gsL, velL);
                        velR = cR0.vel;
                        gsR.copy_values_from(cR0.gas);
                    } else {
                        // Assume at least two cells to the left and two to the right.
                        LCell cL1 = cells[i-2];
                        LCell cL0 = cells[i-1];
                        LCell cR0 = cells[i];
                        LCell cR1 = cells[i+1];
                        interpL_prepare(cL1.L, cL0.L, cR0.L);
                        interpL(cL1, cL0, cR0, gsL, velL);
                        interpR_prepare(cL0.L, cR0.L, cR1.L);
                        interpR(cL0, cR0, cR1, gsR, velR);
                    }
                    lrivp(gsL, gsR, velL, velR, gmodel, gmodel, f.dxdt[level], f.p);
                } else {
                    // Just use cell-centre values directly.
                    LCell cL0 = cells[i-1];
                    LCell cR0 = cells[i];
                    lrivp(cL0.gas, cR0.gas, cL0.vel, cR0.vel, gmodel, gmodel, f.dxdt[level], f.p);
                } // end if(x_order
            }
            // At this point pLstar and pRstar should be the same value.
            if (isNaN(f.dxdt[level]) || isNaN(f.p) || f.p < 0.0) {
                string msg = "Bad Riemann solve.";
                debug {
                    msg ~= text(" i=", i, " f.x=", f.x, " level=", level);
                    msg ~= text(" f.dxdt=", f.dxdt[level], " f.p=", f.p);
                }
                throw new Exception(msg);
            }
        } // end foreach f
        //
        // After we compute the normal Riemann solutions,
        // we apply the closing-off effect of any valves.
        foreach (f; faces) {
            // Start by assuming that the tube is completely open and
            // that there are zero pressures associated with any valve effect.
            f.fopen = 1.0; f.pLstar = 0.0; f.pRstar = 0.0;
        }
        //
        // Presently, this is limited to internal faces within the gas slug,
        // so we should be careful to set up simulations that do not have multiple
        // gas slugs running through a partially-closed valve.
        // Once a valve is fully open, this restriction no longer matters.
        //
        foreach (valve; simcore.valves) {
            // Note that we need the current value of time
            // in order to evaluate fopen for the valve.
            double fopen = valve.fopen(t);
            if (fopen > 0.999) { continue; } // Fully open so no need to do more.
            if ((faces[0].x < valve.x) && (faces[$-1].x > valve.x)) {
                // Do something to the internal face that is nearest to the valve.
                // We limit the interaction to internal faces because we do not
                // wish to have to deal with all of the special cases of an end face.
                size_t closest_f = 1;
                double closest_distance = fabs(faces[1].x - valve.x);
                foreach (i; 2 .. faces.length-1) {
                    double distance = fabs(faces[i].x - valve.x);
                    if (distance < closest_distance) {
                        closest_f = i;
                        closest_distance = distance;
                    }
                }
                // Do the two one-sided calculations as if the valve is closed.
                LFace f = faces[closest_f];
                f.fopen = fopen;
                LCell cL = cells[closest_f-1];
                LCell cR = cells[closest_f];
                piston_at_right(cL.gas, cL.vel, gmodel, 0.0, f.pLstar);
                piston_at_left(cR.gas, cR.vel, gmodel, 0.0, f.pRstar);
                // Apply the effect of the valve to the velocity.
                f.dxdt[level] *= fopen;
                // Below, the other effects of the valve will be fed
                // into the conservation equations for the gas cells.
            }
        } // end foreach valve
        //
        // Now that we have the face velocities and pressures,
        // we are ready to focus on the cell properties.
        foreach (c; cells) {
            c.source_terms(viscous_effects, adiabatic, gmodel);
        }
        // Conservation equations determine our time derivatives.
        //      faces[i]       cells[i]       faces[i+1]
        //         +----------------------------+
        //         |                            |
        // pLstar fL pRstar       c      pLstar fR pRstar
        //         |                            |
        //         +----------------------------+
        foreach (i, c; cells) {
            LFace fL = faces[i];
            LFace fR = faces[i+1];
            // Mass.
            c.dmassdt[level] = c.Q_mass;
            // Momentum -- force on cell
            // Note that the valve faces have a contribution to momentum.
            c.dmomdt[level] = fL.pRstar*fL.area*(1.0-fL.fopen) - fR.pLstar*fR.area*(1.0-fR.fopen) +
                fL.p*fL.area*fL.fopen - fR.p*fR.area*fR.fopen +
                c.gas.p*(fR.area - fL.area) + c.Q_moment;
            // Energy -- work done on cell.
            // Only the fraction-open will contribute to energy transfer across faces.
            c.dEdt[level] = fL.p*fL.area*fL.fopen*fL.dxdt[level] -
                fR.p*fR.area*fR.fopen*fR.dxdt[level] + c.Q_energy;
            // Particle distance travelled.
            c.dL_bardt[level] = fabs(c.vel);
        }
        return;
    } // end time_derivatives()

    @nogc
    void predictor_step(double dt)
    {
        foreach (f; faces) { f.predictor_step(dt); }
        foreach (c; cells) { c.predictor_step(dt, gmodel); }
        compute_areas_and_volumes();
        foreach (c; cells) { c.decode_conserved(gmodel); }
        return;
    }

    @nogc
    void corrector_step(double dt)
    {
        foreach (f; faces) { f.corrector_step(dt); }
        foreach (c; cells) { c.corrector_step(dt, gmodel); }
        compute_areas_and_volumes();
        foreach (c; cells) { c.decode_conserved(gmodel); }
        return;
    }

    @nogc
    void chemical_increment(double dt)
    {
        ThermochemicalReactor reactor = reactors[gmodel_id];
        if (!reactor) return;
        if (!(gmodel is reactor._gmodel)) {
            throw new Error("Gas model objects do not match.");
        }
        foreach (c; cells) {
            c.chemical_increment(dt, gmodel, reactor);
            if (viscous_effects > 0) {
                gmodel.update_trans_coeffs(c.gas);
            }
        }
        return;
    }

    @nogc
    int bad_cells()
    {
        int bad_cell_count = 0;
        foreach (i, c; cells) {
            bool cell_is_bad = false;
            if (c.L <= 0.0) { cell_is_bad = true; }
            if (c.mass <= 0.0) { cell_is_bad = true; }
            if (c.gas.rho <= 0.0) { cell_is_bad = true; }
            if (c.gas.T <= 0.0) { cell_is_bad = true; }
            if (c.gas.p <= 0.0) { cell_is_bad = true; }
            if (cell_is_bad) {
                bad_cell_count++;
                debug { writeln("Bad cell at xmid=", c.xmid); }
            }
        }
        return bad_cell_count;
    }

    @nogc
    double suggested_time_step(double cfl_value)
    {
        LCell c = cells[0];
        double signal_time = c.L / (c.gas.a + fabs(c.vel));
        double smallest_transit_time = signal_time;
        foreach (i; 1 .. cells.length) {
            c = cells[i];
            signal_time = c.L / (c.gas.a + fabs(c.vel));
            smallest_transit_time = min(smallest_transit_time, signal_time);
        }
        return smallest_transit_time * cfl_value;
    }

    // Quadratic reconstruction functions adapted from Eilmer, 2020-05-22.
    // For a sketch of the idea, see PJ workbook notes Jan 2001.

    immutable double epsilon_van_albada = 1.0e-12;

    @nogc
    double interpL_scalar(double qL1, double qL0, double qR0)
    {
        double qL = qL0;
        double delLminus = (qL0 - qL1) * two_over_lenL0_plus_lenL1;
        double del = (qR0 - qL0) * two_over_lenR0_plus_lenL0;
        double sL = (delLminus*del + fabs(delLminus*del) + epsilon_van_albada) /
            (delLminus*delLminus + del*del + epsilon_van_albada);
        qL = qL0 + sL * aL0 * (del * two_lenL0_plus_lenL1 + delLminus * lenR0);
        qL = clip_to_limits(qL, qL0, qR0);
        return qL;
    }

    @nogc
    void interpL_prepare(double lenL1, double lenL0, double lenR0)
    {
        this.lenL0 = lenL0;
        this.lenR0 = lenR0;
        this.aL0 = 0.5 * lenL0 / (lenL1 + 2.0*lenL0 + lenR0);
        this.two_over_lenL0_plus_lenL1 = 2.0 / (lenL0 + lenL1);
        this.two_over_lenR0_plus_lenL0 = 2.0 / (lenR0 + lenL0);
        this.two_lenL0_plus_lenL1 = (2.0*lenL0 + lenL1);
    }

    @nogc
    double interpR_scalar(double qL0, double qR0, double qR1)
    {
        double del = (qR0 - qL0) * two_over_lenR0_plus_lenL0;
        double delRplus = (qR1 - qR0) * two_over_lenR1_plus_lenR0;
        double sR = (del*delRplus + fabs(del*delRplus) + epsilon_van_albada) /
            (del*del + delRplus*delRplus + epsilon_van_albada);
        double qR = qR0 - sR * aR0 * (delRplus * lenL0 + del * two_lenR0_plus_lenR1);
        qR = clip_to_limits(qR, qL0, qR0);
        return qR;
    }

    @nogc
    void interpR_prepare(double lenL0, double lenR0, double lenR1)
    {
        this.lenL0 = lenL0;
        this.lenR0 = lenR0;
        this.aR0 = 0.5 * lenR0 / (lenL0 + 2.0*lenR0 + lenR1);
        this.two_over_lenR0_plus_lenL0 = 2.0 / (lenR0 + lenL0);
        this.two_over_lenR1_plus_lenR0 = 2.0 / (lenR1 + lenR0);
        this.two_lenR0_plus_lenR1 = (2.0*lenR0 + lenR1);
    }

    @nogc
    double clip_to_limits(double q, double A, double B)
    // Returns q if q is between the values A and B, else
    // it returns the closer limit of the range [min(A,B), max(A,B)].
    {
        double lower_limit = (A <= B) ? A : B;
        double upper_limit = (A > B) ? A : B;
        double qclipped = (q > lower_limit) ? q : lower_limit;
        return (qclipped <= upper_limit) ? qclipped : upper_limit;
    }

    @nogc
    void interpL(ref LCell cL1, ref LCell cL0, ref LCell cR0,
                 ref GasState gasL, ref double velL)
    {
        interpL_prepare(cL1.L, cL0.L, cR0.L);
        gasL.copy_values_from(cL0.gas);
        gasL.rho = interpL_scalar(cL1.gas.rho, cL0.gas.rho, cR0.gas.rho);
        gasL.T = interpL_scalar(cL1.gas.T, cL0.gas.T, cR0.gas.T);
        velL = interpL_scalar(cL1.vel, cL0.vel, cR0.vel);
        gmodel.update_thermo_from_rhoT(gasL);
    }

    @nogc
    void interpR(ref LCell cL0, ref LCell cR0, ref LCell cR1,
                 ref GasState gasR, ref double velR)
    {
        interpR_prepare(cL0.L, cR0.L, cR1.L);
        gasR.copy_values_from(cR0.gas);
        gasR.rho = interpR_scalar(cL0.gas.rho, cR0.gas.rho, cR1.gas.rho);
        gasR.T = interpR_scalar(cL0.gas.T, cR0.gas.T, cR1.gas.T);
        velR = interpR_scalar(cL0.vel, cR0.vel, cR1.vel);
        gmodel.update_thermo_from_rhoT(gasR);
    }

private:
    // Workspace for interpolation.
    double lenL0, lenR0, aL0, aR0;
    double two_over_lenL0_plus_lenL1;
    double two_over_lenR0_plus_lenL0;
    double two_over_lenR1_plus_lenR0;
    double two_lenL0_plus_lenL1;
    double two_lenR0_plus_lenR1;
    double velL, velR;
    GasState gsL, gsR;
} // end class GasSlug
