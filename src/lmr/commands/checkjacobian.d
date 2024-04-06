/**
 * Module for running a cross-check on the Jacobian.
 *
 * This command is not intended for users.
 * We'll put it in the development/diagnostics bucket.
 *
 * Authors: KAD and RJG
 * Date: 2023-07-18
 * History:
 *   2024-03-21 moved from the lmr commands to its own separate executable to allow for complex number version
 */

module checkjacobian;

import std.stdio;
import std.algorithm;
import std.random;
import std.getopt;
import std.range : empty;
import std.format : format;
import std.conv : to;
import std.math;
import std.json : JSONValue;

import nm.number : number;
import nm.smla;
import ntypes.complex;

import globalconfig;
import globaldata : localFluidBlocks;
import init : initConfiguration;
import newtonkrylovsolver;
import cmdhelper;
import command;
import json_helper;
import lmrconfig;
import blockio;
import flowsolution;

string cmdName;
Command checkJacCmd;

static this()
{
    version(complex_numbers) {
        cmdName = "lmrZ-check-jacobian";
    } else {
        cmdName = "lmr-check-jacobian";
    }
    checkJacCmd.type = LmrCmdType.dev;
    checkJacCmd.main = &main;
    checkJacCmd.description = `
Check formation of the Jacobian via two internal methods:
1. The sparse matrix representation multiplied by a test vector; and
2. The Frechet derivative approach which gives Jv directly.
`;
    checkJacCmd.shortDescription = "Check the formation of the Jacobian.";
    checkJacCmd.helpMsg = format(
`%s [options]

Check the formation of the Jacobian for a snapshot of the field variables.

If no snapshot option is passed, then process the final snapshot by default.

options ([+] can be repeated):

 -s, --snapshot
     a single snapshot to use for Jacobian check
     examples:
       --snapshot=5 : process snapshot 5
       -s 2 : process snapshot 2
     default: none

 -o, --output
     write output to a file
     example:
       --output=norms-data.txt
     default: none (just write to STDOUT)

 -rfl, --read-frozen-limiter-values
     read limiter values associated with snapshot
     default: false, do NOT read limiter values

 -v, --verbose [+]
     Increase verbosity.
`, cmdName);
}

int main(string[] args)
{
    bool helpWanted = false;
    int verbosity = 0;
    int snapshot = -1;
    string outFilename;
    bool readFrozenLimiterValues = false;
    try {
        getopt(args,
               config.bundling,
               "h|help", &helpWanted,
               "v|verbose+", &verbosity,
               "s|snapshot", &snapshot,
               "o|output", &outFilename,
               "rfl|read-frozen-limiter-values", &readFrozenLimiterValues);
    } catch (Exception e) {
        writefln("Eilmer %s program quitting.", cmdName);
        writeln("There is something wrong with the command-line arguments/options.");
        writeln(e.msg);
        return 1;
    }

    if (helpWanted || canFind(args, "help")) {
        writeln(checkJacCmd.helpMsg);
        return 0;
    }

    if (verbosity > 0) writefln("%s: Begin program.", cmdName);

    // Use stdout if no output filename is supplied,
    // or open a file ready for use if one is.
    File outfile = outFilename.empty() ? stdout : File(outFilename, "w");

    auto availSnapshots = determineAvailableSnapshots();
    int finalSnapshot = to!int(availSnapshots.length - 1);
    if (snapshot < 0) { // assume nothing set, so select final snapshot
        snapshot = finalSnapshot;
    }
    if (snapshot > finalSnapshot) {
        // Rather than handle this error, let's define it away
        // A snapshot request beyond the number available gets set to the final snapshot
        // We'll warn the user though.
        writefln("%s: WARNING", cmdName);
        writefln("W: The snapshot requested {%d} is not available. Instead, selecting final snapshot {%d}.", snapshot, finalSnapshot);
        snapshot = finalSnapshot;
    }

    if (verbosity > 1) writefln("%s: Performing check with snapshot {%d}", cmdName, snapshot);

    /*
     * 0. Initialize Newton-Krylov simulation
     */
    initNewtonKrylovSimulation(snapshot, 1, 1, "1000000");
    JSONValue jsonData = readJSONfile(lmrCfg.nkCfgFile);
    nkCfg.readValuesFromJSON(jsonData);
    // Allocate space and configure phases
    nkPhases.length = nkCfg.numberOfPhases;
    foreach (i, ref phase; nkPhases) {
        string key = "NewtonKrylovPhase_" ~ to!string(i);
        phase.readValuesFromJSON(jsonData[key]);
    }
    activePhase = nkPhases[$-1]; // set to the last phase from the simulation
    evalResidual(0); // we perform a residual evaluation here to populate the ghost-cells
    alias cfg = GlobalConfig;
    size_t nConserved = cfg.cqi.n;

    if (verbosity > 1) writefln("%s: Performing check with Jacobian perturbation parameter %.4e and Frechet perturbation parameter %.4e",
                                cmdName, nkCfg.preconditionerPerturbation, nkCfg.frechetDerivativePerturbation);

    // prepare fluidblock (note that this program only operates on a single fluidblock)
    auto blk = localFluidBlocks[0];
    if (readFrozenLimiterValues) readLimiterValues(snapshot);

    /*
     * 1. Do calculation of sparse-matrix Jacobian by test vector.
     */
    // 1a. Populate Jacobian
    blk.initialize_jacobian(cfg.interpolation_order, nkCfg.preconditionerPerturbation, 0);
    blk.evaluate_jacobian();
    // 1b. Prepare a test vector
    double[] testVec;
    testVec.length = blk.flowJacobian.local.ia.length-1;
    Mt19937 rng;
    double mag = 0.0;
    foreach (ref val; testVec) {
        val = rng.uniform01();
        mag += val*val;
    }
    testVec[] = (1./mag)*testVec[]; // we have observed that a normalized test vector behaves better

    // 1c. Compute c0 = J*testVec
    double[] c0;
    c0.length = testVec.length;
    multiply(blk.flowJacobian.local, testVec, c0);

    /*
     * 2. Do calculation of Jv with Frechet derivative.
     */
    // 2a. Populate z with test vector
    allocateGlobalGMRESWorkspace();
    blk.allocate_GMRES_workspace(1, nkCfg.useRealValuedFrechetDerivative);

    // 2b. Do Frechet step (note for the real-valued variant we need to fill R0)
    version(complex_numbers) {
        // do nothing
    } else {
        int startIdx = 0;
        foreach (cell; blk.cells) {
            foreach (ivar; 0 .. nConserved) { blk.R0[startIdx+ivar] = cell.dUdt[0][ivar].re; }
            startIdx += nConserved;
        }
    }
    blk.zed[] = testVec[];
    evalJacobianVectorProduct(nkCfg.frechetDerivativePerturbation);

    /*
     * 3. Compare vectors and write result to file.
     */
    double c0_2norm = 0.0;
    foreach (val; c0) c0_2norm += val*val;
    c0_2norm = sqrt(c0_2norm);

    double c1_2norm = 0.0;
    foreach (val; blk.zed) c1_2norm += val*val;
    c1_2norm = sqrt(c1_2norm);

    writefln("c0, 2-norm= %.16e", c0_2norm.re);
    writefln("c1, 2-norm= %.16e", c1_2norm.re);

    outfile.writeln("# array_idx cell_idx c0   c1   c0-c1   abs(c0-c1)");
    foreach (i; 0 .. c0.length) {
        auto delta = c0[i] - blk.zed[i];
        size_t cell_id = i/nConserved;
        outfile.writefln("%d %d %20.16e %20.16e %20.16e %20.16e", i, cell_id, c0[i].re, blk.zed[i].re, delta.re, fabs(delta).re);
    }
    outfile.close();

    return 0;
}

void readLimiterValues(int snapshot)
{
    /*
    Reads limiter values from file and then freezes limiter.
    localFluidBlocks must be initialised prior to calling this function.
    */
    alias cfg = GlobalConfig;
    double[][] data;
    string[] variables;
    string fileFmt;

    fileFmt = cfg.field_format;
    variables = readVariablesFromMetadata(lmrCfg.limiterMetadataFile);
    size_t[string] variableIndex;
    foreach (i, var; variables) variableIndex[var] = i;
    auto soln = new FlowSolution(to!int(snapshot), cfg.nFluidBlocks);
    foreach (i, blk; localFluidBlocks) {
        readValuesFromFile(data, limiterFilename(to!int(snapshot), to!int(i)),
                           variables, soln.flowBlocks[0].ncells, fileFmt);
        foreach (j, cell; blk.cells) {
            cell.gradients.velxPhi = data[j][variableIndex["vel.x"]];
            cell.gradients.velyPhi = data[j][variableIndex["vel.y"]];
            if (cfg.dimensions == 3) {
                cell.gradients.velzPhi = data[j][variableIndex["vel.z"]];
            } else {
                cell.gradients.velzPhi = 0.0;
            }
            final switch (cfg.thermo_interpolator) {
                case InterpolateOption.pt:
                    cell.gradients.pPhi = data[j][variableIndex["p"]];
                    cell.gradients.TPhi = data[j][variableIndex["T"]];
                    break;
                case InterpolateOption.rhou:
                    cell.gradients.rhoPhi = data[j][variableIndex["rho"]];
                    cell.gradients.uPhi = data[j][variableIndex["e"]];
                    break;
                case InterpolateOption.rhop:
                    cell.gradients.rhoPhi = data[j][variableIndex["rho"]];
                    cell.gradients.pPhi = data[j][variableIndex["p"]];
                    break;
                case InterpolateOption.rhot:
                    cell.gradients.rhoPhi = data[j][variableIndex["rho"]];
                    cell.gradients.TPhi = data[j][variableIndex["T"]];
                    break;
            }
            foreach (isp; 0 .. cfg.gmodel_master.n_species) {
                cell.gradients.rho_sPhi[isp] = data[j][variableIndex["massf-" ~ cfg.gmodel_master.species_name(isp)]];
            }
            foreach (imode; 0 .. cfg.gmodel_master.n_modes) {
                if (cfg.thermo_interpolator == InterpolateOption.rhou ||
                    cfg.thermo_interpolator == InterpolateOption.rhop ) {
                    cell.gradients.u_modesPhi[imode] = data[j][variableIndex["e-" ~ cfg.gmodel_master.energy_mode_name(imode)]];
                }
                else {
                    cell.gradients.T_modesPhi[imode] = data[j][variableIndex["T-" ~ cfg.gmodel_master.energy_mode_name(imode)]];
                }
            }
            if (cfg.turbulence_model_name != "none") {
                foreach (iturb; 0 .. cfg.turb_model.nturb) {
                    cell.gradients.turbPhi[iturb] = data[j][variableIndex["tq-" ~ cfg.turb_model.primitive_variable_name(iturb)]];
                }
            }
            if (cfg.MHD) {
                writeln("WARNING: lmr-check-jacobian currently does not support reading MHD limiter values, proceeding with re-evaluated MHD limiter values.");
            }
        }
    }
    cfg.frozen_limiter = true;
    activePhase.frozenLimiterForJacobian = true;
    activePhase.frozenLimiterForResidual = true;
}
