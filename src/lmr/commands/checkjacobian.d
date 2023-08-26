/**
 * Module for running a cross-check on the Jacobian.
 *
 * This command is not intended for users.
 * We'll put it in the development/diagnostics bucket.
 *
 * Authors: KAD and RJD
 * Date: 2023-07-18
 */

module checkjacobian;

import std.stdio;
import std.random;
import std.getopt;
import std.range : empty;
import std.format : format;
import std.conv : to;
import std.math;
import std.json : JSONValue;

import nm.number : number;
import nm.smla;
import nm.complex;

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

Command checkJacCmd;

string cmdName = "check-jacobian";

static this()
{
    checkJacCmd.type = LmrCmdType.dev;
    checkJacCmd.main = &main_;
    checkJacCmd.description = `
Check formation of the Jacobian via two internal methods:
1. The sparse matrix representation multiplied by a test vector; and
2. The Frechet derivative approach which gives Jv directly.
`;
    checkJacCmd.shortDescription = "Check the formation of the Jacobian.";
    checkJacCmd.helpMsg = format(
`lmr %s [options]

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

 -v, --verbose [+]
     Increase verbosity.
`, cmdName);
}

void main_(string[] args)
{
    int verbosity = 0;
    int snapshot = -1;
    string outFilename;
    bool readFrozenLimiterValues = false;
    getopt(args,
           config.bundling,
           "v|verbose+", &verbosity,
           "s|snapshot", &snapshot,
           "o|output", &outFilename,
           "r|read-frozen-limiter-values", &readFrozenLimiterValues);

    if (verbosity > 0) writefln("lmr %s: Begin program.", cmdName);

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
        writefln("lmr %s: WARNING", cmdName);
        writefln("W: The snapshot requested {%d} is not available. Instead, selecting final snapshot {%d}.", snapshot, finalSnapshot);
        snapshot = finalSnapshot;
    }

    if (verbosity > 1) writefln("lmt %s: Performing check with snapshot {%d}", cmdName, snapshot);

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

    /*
     * 1. Do calculation of sparse-matrix Jacobian by test vector.
     */
    // 1a. Populate Jacobian
    alias cfg = GlobalConfig;
    auto blk = localFluidBlocks[0];
    if (readFrozenLimiterValues) readLimiterValues(snapshot);
    blk.initialize_jacobian(cfg.interpolation_order, 1.0e-250, 0);
    blk.evaluate_jacobian();
    // 1b. Prepare a test vector
    number[] testVec;
    testVec.length = blk.flowJacobian.local.ia.length-1;
    Mt19937 rng;
    number mag = 0.0;
    foreach (ref val; testVec) {
        val = rng.uniform01();
        mag += val*val;
    }
    testVec[] = (1./mag)*testVec[]; // we have observed that a normalized test vector behaves better

    // 1c. Compute c0 = J*testVec
    number[] c0;
    c0.length = testVec.length;
    multiply(blk.flowJacobian.local, testVec, c0);

    /*
     * 2. Do calculation of Jv with Frechet derivative.
     */
    // 2a. Populate z with test vector
    allocateGlobalGMRESWorkspace();
    blk.allocate_GMRES_workspace(1);

    // 2b. Do Frechet step
    blk.zed[] = testVec[];
    evalJacobianVectorProduct(1.0e-250);

    /*
     * 3. Compare vectors and write result to file.
     */
    number c0_2norm = 0.0;
    foreach (val; c0) c0_2norm += val*val;
    c0_2norm = sqrt(c0_2norm);

    number c1_2norm = 0.0;
    foreach (val; blk.zed) c1_2norm += val*val;
    c1_2norm = sqrt(c1_2norm);

    writefln("c0, 2-norm= %.16e", c0_2norm.re);
    writefln("c1, 2-norm= %.16e", c1_2norm.re);

    outfile.writeln("# array_idx cell_idx c0   c1   c0-c1   abs(c0-c1)");
    size_t nConserved = cfg.cqi.n;
    foreach (i; 0 .. c0.length) {
        auto delta = c0[i] - blk.zed[i];
        size_t cell_id = i/nConserved;
        outfile.writefln("%d %d %20.16e %20.16e %20.16e %20.16e", i, cell_id, c0[i].re, blk.zed[i].re, delta.re, fabs(delta).re);
    }
    outfile.close();

    return;
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

    fileFmt = cfg.flow_format;
    variables = readVariablesFromMetadata(lmrCfg.limiterMetadataFile);
    auto soln = new FlowSolution(to!int(snapshot), cfg.nFluidBlocks);
    foreach (i, blk; localFluidBlocks) {
        readValuesFromFile(data, limiterFilename(to!int(snapshot), to!int(i)), variables, 
            soln.flowBlocks[0].ncells, fileFmt);
        foreach (j, cell; blk.cells) {
            cell.gradients.rhoPhi = data[j][0];
            cell.gradients.pPhi = data[j][1];
            cell.gradients.velxPhi = data[j][2];
            cell.gradients.velyPhi = data[j][3];
            if (blk.myConfig.dimensions == 3) {
                cell.gradients.velzPhi = data[j][4];
            } else {
                cell.gradients.velzPhi = 0.0;
            }
            foreach(it; 0 .. blk.myConfig.turb_model.nturb){
                cell.gradients.turbPhi[it] = data[j][5+it];
	    }
        }
    }
    cfg.frozen_limiter = true;
}
