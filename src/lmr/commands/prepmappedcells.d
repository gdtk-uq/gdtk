/**
 * Module for preparing a list of mapped cells.
 *
 * Author: RJG
 * Date: 2023-11-18
 */

module prepmappedcells;

import std.getopt;
import std.stdio : writeln, writefln, File;
import std.string : toStringz;
import std.conv : to;
import std.path : dirName;
import std.file : thisExePath;

import lmrconfig : lmrCfg;
import command;
import globalconfig;
import globaldata;
import fluidblock;
import init;
import geom;
import bc;

Command prepMappedCellsCmd;
string cmdName = "prep-mapped-cells";

static this()
{
    prepMappedCellsCmd.main = &main_;
    prepMappedCellsCmd.description = "Prepare mapped cells list via search for an Eilmer simulation.";
    prepMappedCellsCmd.shortDescription = "Prepare mapped cells list.";
    prepMappedCellsCmd.helpMsg = 
`lmr prep-mapped-cells [options]

Prepare a list of mapped cells via a search over unstructured blocks.

You do NOT need to use this command if you have partitioned an unstructured grid
with one of the 3rd part partitioning tools. This command is intended for a 
particular use case where Eilmer native grids have been prepared but the mapped
cells information is not available.

options ([+] can be repeated):

 -v, --verbose [+]
     Increase the verbosity during mapped cells list creation

`;

}

int main_(string[] args)
{
    int verbosity = 0;
    try {
        getopt(args,
               config.bundling,
               "v|verbose+", &verbosity,
               );
    } catch (Exception e) {
        writefln("Eilmer %s program quitting.", cmdName);
        writeln("There is something wrong with the command-line arguments/options.");
        writeln(e.msg);
        return 1;
    }

    if (verbosity > 0) {
        writeln("lmr prep-mapped-cells: Begin preparation of mapped cells list.");
    }

    /* Initialise the flow solver up to the point where we could do a mapped cells search.
     * Adapted set up from initNewtonKrylovSimulation up to the call to: initMappedCellDataExchange()
     */

    alias cfg = GlobalConfig;

    initConfiguration();
    if (cfg.nFluidBlocks == 0 && cfg.is_master_task) {
        throw new Error("No FluidBlocks; no point in continuing with mapped cell search.");
    }
    cfg.n_flow_time_levels = 2;

    initLocalBlocks();

    initThreadPool(1, 1);

    initFluidBlocksBasic();
    initFluidBlocksMemoryAllocation();
    initFluidBlocksGlobalCellIDStarts();
    initFluidBlocksZones();
    initFluidBlocksFlowField(0);

    BlockAndCellId[string][size_t] mappedCellsList;

    foreach (blk; localFluidBlocks) {
        foreach (bc; blk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce = cast(GhostCellMappedCellCopy) gce;
                if (mygce) {
                    findMappedCells(blk, mygce, mappedCellsList);
                }
            }
        }
    }


    auto of = File(lmrCfg.mappedCellsFile, "w");
    foreach (ib; mappedCellsList.keys()) {
        bool[size_t] nghbrBlks;
        foreach (faceTag; mappedCellsList[ib].keys()) {
            nghbrBlks[mappedCellsList[ib][faceTag].blkId] = true;
        }
        of.writef("MappedBlocks in BLOCK[%d]= ", ib);
        foreach (nblkId; nghbrBlks.keys()) {
            of.writef("%d ", nblkId);
        }
        of.writef("\n");
    }
    foreach (ib, blk; localFluidBlocks) {
        bool headerWritten = false;
        foreach (bc; blk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce = cast(GhostCellMappedCellCopy) gce;
                if (mygce) {
                    if (!headerWritten) {
                        of.writefln("NMappedCells in BLOCK[%d]= %d", ib, mappedCellsList[ib].length);
                        headerWritten = true;
                    }
                    foreach (face; bc.faces) {
                        size_t[] vtxs; foreach (vtx; face.vtx) vtxs ~= vtx.id;
                        string faceTag = makeFaceTag(vtxs);
                        of.writefln("%-20s %4d %6d", faceTag, mappedCellsList[ib][faceTag].blkId, mappedCellsList[ib][faceTag].cellId);
                    }
                }
            }
        }
    }
    of.close();
    return 0;
}



void findMappedCells(FluidBlock blk, GhostCellMappedCellCopy gcmc, ref BlockAndCellId[string][size_t] mappedCellsList)
{
    string[] faceTags;
    final switch(blk.grid_type) {
    case Grid_t.unstructured_grid:
        BoundaryCondition bc = blk.bc[gcmc.which_boundary];
        foreach (i, face; bc.faces) {
            gcmc.ghost_cells ~= (bc.outsigns[i] == 1) ? face.right_cell : face.left_cell;
            size_t[] vtxs; foreach (vtx; face.vtx) vtxs ~= vtx.id;
            faceTags ~= makeFaceTag(vtxs);
        }
        break;
    case Grid_t.structured_grid:
        throw new Error("This command makes no sense on a structured-grid block.");
    }

    foreach (ic, gc; gcmc.ghost_cells) {
        Vector3 ghostpos = gc.pos[0];

        // First, attempt to find via enclosing cell.
        bool found = false;
        foreach (ib, bblk; localFluidBlocks) {
            found = false;
            size_t indx = 0;
            bblk.find_enclosing_cell(ghostpos, indx, found);
            if (found) {
                mappedCellsList[blk.id][faceTags[ic]] = BlockAndCellId(ib, indx);
                break;
            }
        }
    }
}

