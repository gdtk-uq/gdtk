/*
 * Module for deforming lmr grids.
 * Many of these routines are ported from eilmer/shape_sensitivity_core.d
 *
 * Authors: RBO, KAD
 * Date: 2025-06-05
 */

module lmr.commands.deformgrid;

import core.memory : GC;
import std.algorithm.searching : canFind;
import std.array : split;
import std.conv : to;
import std.format : format;
import std.getopt;
import std.json : JSONValue;
import std.math : sqrt;
import std.range;
import std.stdio;
import std.string : strip;

import geom.elements;

import lmr.commands.cmdhelper;
import lmr.commands.command;
import lmr.lmrconfig : gridFilename;
import lmr.fluidblock : FluidBlock;
import lmr.globalconfig;
import lmr.globaldata : localFluidBlocks;
import lmr.init;
import lmr.newtonkrylovsolver;

Command deformGridCmd;
string cmdName = "deform-grid";

void loadFluidBlocks(int snapshotStart) {
    alias cfg = GlobalConfig;

    JSONValue cfgData = initConfiguration();
    if (cfg.nFluidBlocks == 0 && cfg.is_master_task) {
        throw new Error("No FluidBlocks; no point in continuing with simulation initialisation.");
    }

    initLocalBlocks();
    initThreadPool(1, 1);
    initFluidBlocksBasic(cfgData);
    initFluidBlocksMemoryAllocation();
    initFluidBlocksGlobalCellIDStarts();
    initFluidBlocksZones();
    initFluidBlocksFlowField(snapshotStart);

    initFullFaceDataExchange(cfgData);
    initMappedCellDataExchange();
    initGhostCellGeometry();

    orderBlocksBySize();
    initMasterLuaState();
    initCornerCoordinates();

    GC.collect();
    GC.minimize();
}

void inverseDistanceWeighting(FluidBlock blk, Vector3[] bndryVtxInitPos, 
    Vector3[] bndryVtxNewPos, size_t gtl, string[] nonFixedBoundaryList=[])
{
    // It is assumed that the boundary vertices have already been perturbed,
    // and have their new positions stored at the current gtl.
    // This function computes the updated positions for all internal vertices.
    double r, w, dx, dy, denomSum; 
    double[2] numerSum;
    Vector3 ds;
    size_t nBndryVtx;

    // ensure init and new boundary pos array has equal length
    nBndryVtx = bndryVtxInitPos.length;
    if (bndryVtxInitPos.length != bndryVtxNewPos.length) {
        string msg = "Inverse Distance Weighting: initial positions array ";
        msg ~= "length does not equal new positions array length.";
        throw new Error(msg);
    }
    
    foreach(vtx; blk.vertices) {
        if (!blk.boundaryVtxIndexList.canFind(vtx.id)) {
            numerSum[0] = 0.0; 
            numerSum[1] = 0.0; 
            denomSum = 0.0;

            // compute numerator and denominator of weighted average
            foreach (i; 0 .. nBndryVtx) {
                ds.x = bndryVtxNewPos[i].x - bndryVtxInitPos[i].x;
                ds.y = bndryVtxNewPos[i].y - bndryVtxInitPos[i].y;
                dx = vtx.pos[0].x - bndryVtxInitPos[i].x; 
                dy = vtx.pos[0].y - bndryVtxInitPos[i].y; 
                r = sqrt(dx*dx + dy*dy);
                w = 1.0 / (r*r);
                numerSum[0] += ds.x * w; 
                numerSum[1] += ds.y * w;
                denomSum += w;
            }

            // update vertex position
            vtx.pos[gtl].x = vtx.pos[0].x + numerSum[0]/denomSum;
            vtx.pos[gtl].y = vtx.pos[0].y + numerSum[1]/denomSum;
        }
    }
}

static this()
{
    deformGridCmd.main = &main_;
    deformGridCmd.description = "Deform a grid from updated boundaries.";
    deformGridCmd.shortDescription = deformGridCmd.description;
    deformGridCmd.helpMsg = format(
`lmr %s [options]

Deform an eilmer grid from updated boundary vertices.

options ([+] can be repeated):

 -f, --filename
     Name of file containing new vertex positions.
     
     It is assumed this file contains the following fields:
     blk.id vtx.id pos.x pos.y pos.z
     
     default:
       --filename=new-bndry-verts.txt

 -g, --group
     Name of boundary group assigned to deforming surface.
     default:
       --group=wall

 -v, --verbose [+]
     Increase verbosity during preparation and writing of VTK files.

`, cmdName);

}

int main_(string[] args)
{
    string bndryGroup = "wall";
    string newBndryFilename = "new-bndry-verts.txt";
    int verbosity = 0;
    try {
        getopt(args,
               config.bundling,
               "f|filename", &newBndryFilename,
               "g|group", &bndryGroup,
               "v|verbose+", &verbosity);
    } catch (Exception e) {
        writefln("Eilmer %s program quitting.", cmdName);
        writeln("There is something wrong with the command-line arguments/options.");
        writeln(e.msg);
        return 1;
    }

    if (verbosity > 0) writefln("lmr %s: Begin program.", cmdName);

    // load in fluid blocks
    // TODO: load from a specified sim directory
    // TODO: raise error if sim is transient mode
    loadFluidBlocks(0);

    // store boundary vertices and their ids
    // TODO: consider internal boundaries - may need to handle these separately
    Vector3[] bndryVtxInitPos;
    foreach (blk; localFluidBlocks) {
        foreach (bndry; blk.bc) {
            if (bndry.group == bndryGroup) {
                foreach(face; bndry.faces) {
                    foreach (vtx; face.vtx) {
                        if (!blk.boundaryVtxIndexList.canFind(vtx.id)) { 
                            blk.boundaryVtxIndexList ~= vtx.id;
                            bndryVtxInitPos ~= vtx.pos[0];
                        }
                    }
                }
            }
        }
    }

    // make copy of grid vertices in gtl=1
    size_t gtl = 1;
    foreach (blk; localFluidBlocks) {
        foreach(vtx; blk.vertices) {
            vtx.pos.length = 2;
            vtx.pos[gtl].x = vtx.pos[0].x;
            vtx.pos[gtl].y = vtx.pos[0].y;
            vtx.pos[gtl].z = vtx.pos[0].z;
        }
    }

    // load updated boundary vertices into gtl=1
    auto file = File(newBndryFilename);
    auto range = file.byLine();
    int blkId, vtxId;
    double x, y, z;
    Vector3 newPos;
    Vector3[] bndryVtxNewPos;
    foreach (line; range) {
        if (!line.empty && line[0] != '#') {
            auto vtxData = strip(line).split;
            blkId = to!int(vtxData[0]);
            vtxId = to!int(vtxData[1]);
            x = to!double(vtxData[2]);
            y = to!double(vtxData[3]);
            z = to!double(vtxData[4]);

            if (localFluidBlocks[blkId].boundaryVtxIndexList.canFind(vtxId)) {
                newPos = Vector3(x,y,z);
                localFluidBlocks[blkId].vertices[vtxId].pos[gtl] = newPos;
                bndryVtxNewPos ~= newPos;
            } else {
                string msg = format("Vertex %d of block %d is not in boundary group %s.", vtxId, blkId, bndryGroup);
                throw new Error(msg);
            }
        }
    }

    // compute new internal vertices with IDW
    foreach (blk; localFluidBlocks) {
        inverseDistanceWeighting(blk, bndryVtxInitPos, bndryVtxNewPos, gtl);
    }

    // write deformed grid to specified output directory
    foreach (blk; localFluidBlocks) {
        foreach(vtx; blk.vertices) {
            vtx.pos[0].x = vtx.pos[gtl].x;
            vtx.pos[0].y = vtx.pos[gtl].y;
        }
        
        // write out grid
        blk.sync_vertices_to_underlying_grid(0);
        // ensure_directory_is_present(); // make new directory for deformed grid
        auto filename = gridFilename(0, blk.id);
        blk.write_underlying_grid(filename);
    }
    
    return 0;
}
