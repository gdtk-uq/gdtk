/**
 * This command can be used to convert an entirely structured grid setup
 * to an unstructured one. It's intended to be used after a prep-grid stage
 * but before a prep-flow stage.
 *
 * Authors: RJG, PAJ, KAD
 * Date: 2023-11-29
 */

module structured2unstructured;

import std.getopt;
import std.stdio : writeln, writefln, File;
import std.json;
import std.format;
import std.conv : to;
import std.file : rename, exists, rmdirRecurse, mkdirRecurse;

import geom;
import bc;

import lmrconfig;
import command;
import json_helper : readJSONfile;


Command structured2unstructuredCmd;
string cmdName = "structured2unstructured";

static this()
{
    structured2unstructuredCmd.main = &main_;
    structured2unstructuredCmd.description = "Convert Eilmer structured grids to unstructured grids.";
    structured2unstructuredCmd.shortDescription = "Convert structured grids to unstructured grids.";
    structured2unstructuredCmd.helpMsg = format(
`lmr %s [options]

Convert Eilmer structured grids to unstructured grids.

NOTE: If command is successful, structured grids will be moved to directory: "%s".
      If command is unsuccessful, the structured grids remain in place in: "%s".
You may remove the entire directory of the original structured grids that were moved into "%s"
if no longer required. (This command won't delete these for you.)

Intended use for this command is during grid preparation before flow preparation.
After preparing structured grids with "prep-grid", you can convert these to unstructured grids
using this command. The new unstructured grids will have all the same boundary information as
their original structured grid counterparts except that (structured) full-face exchanges will
be replaced with (unstructured) mapped cell exchanges. A mapped-cells file will also be created.

options ([+] can be repeated):

 -v, --verbose [+]
     Increase the verbosity during grid conversion process.
 -t, --grid-type=gziptext|rawbinary
     Set the grid type. Default: gziptext

`, cmdName, lmrCfg.savedSgridDir, lmrCfg.gridDirectory, lmrCfg.savedSgridDir);
}

enum GridType { gziptext, rawbinary }

int main_(string[] args)
{
    int verbosity = 0;
    GridType gridType;
    try {
        getopt(args,
               config.bundling,
               "v|verbose+", &verbosity,
               config.noBundling,
               "t|grid-type", &gridType
               );
    } catch (Exception e) {
        writefln("Eilmer %s program quitting.", cmdName);
        writeln("There is something wrong with the command-line arguments/options.");
        writeln(e.msg);
        return 1;
    }

    if (verbosity >= 1) {
        writefln("lmr %s: Begin conversion of structured grids to unstructured grids.", cmdName);
    }

    /* 0. Pick up metadata and see if we can do anything with it. */
    if (verbosity >= 2) {
        writefln("   Reading grid metadata from '%s'", lmrCfg.gridMetadataFile);
    }
    JSONValue gridMetadata = readJSONfile(lmrCfg.gridMetadataFile);
    int ngrids = to!int(gridMetadata["ngrids"].integer);
    if (verbosity >= 2) {
        writefln("   Number of grids to convert: %d", ngrids);
    }
    // Pass over all grids and check that they are structured.
    foreach (ig; 0 .. ngrids) {
        JSONValue thisGridMetadata = readJSONfile(gridMetadataName(ig));
        if (thisGridMetadata["type"].str != "structured_grid") {
            string errMsg = format("Grid %d does not have type 'structured_grid'.\n", ig);
            errMsg ~= format("%s command will not proceed.\n", cmdName);
            throw new Error(errMsg);
        }
    }
    // If we get this far, we at least have structured grids (and only structured grids) to work with.

    /* 1. Let's read in the grids. */
    if (verbosity >= 2) {
        writefln("   Reading in structured grids.");
    }
    StructuredGrid[] sgrids;
    JSONValue[] sgridsMetadata;
    foreach (ig; 0 .. ngrids) {
        if (verbosity >= 3) {
            writefln("      Reading in structured grid: %04d", ig);
        }
        sgridsMetadata ~= readJSONfile(gridMetadataName(ig));
        if (gridType == GridType.gziptext)
            sgrids ~= new StructuredGrid(gridName(ig, gridType), "gziptext");
        else // presently only other choice is 'rawbinary'
            sgrids ~= new StructuredGrid(gridName(ig, gridType), "rawbinary");
    }

    /* 2. Now, form unstructured grids from the structured grids. */
    if (verbosity >= 2) {
        writeln("   Creating unstructured grids.");
    }
    UnstructuredGrid[] ugrids;
    foreach (ig, sgrid; sgrids) {
        if (verbosity >= 3) {
            writefln("      Creating unstructured grid: %04d", ig);
        }
        ugrids ~= new UnstructuredGrid(sgrid);
    }

    /* 3. Work through connections creating mapped cells. */
    BlockAndCellId[string][size_t] mappedCellsList;
    string[size_t][size_t] newTags;
    string[CellAndFaceId][size_t] globalFaceTags;
    BlockAndBoundary[BlockAndBoundary] connMap;
    convertConnections(gridMetadata, sgrids, mappedCellsList, newTags, globalFaceTags, connMap);

    /* 4. Move structured grid out of the way */
    if (exists(lmrCfg.savedSgridDir)) {
        if (verbosity >= 3) {
            writefln("      Making room for saving sgrids. Removing: %s", lmrCfg.savedSgridDir);
        }
        rmdirRecurse(lmrCfg.savedSgridDir);
    }
    if (verbosity >= 2) {
       writefln("   Moving structured grids into: %s", lmrCfg.savedSgridDir);
    }
    rename(lmrCfg.gridDirectory, lmrCfg.savedSgridDir);

    /* 5. Write out mapped cells file. */
    if (ngrids > 1) {
        if (verbosity >= 2) {
            writefln("   Writing %s file.", lmrCfg.mappedCellsFile);
        }
        writeMappedCellsFile(sgrids, mappedCellsList, newTags, globalFaceTags, connMap);
    }

    /* 6. Write out new grids. */
    if (verbosity >= 2) {
        writeln("   Writing unstructured grids.");
    }
    mkdirRecurse(lmrCfg.gridDirectory);
    foreach (ig, ugrid; ugrids) {
        if (verbosity >= 3) {
            writefln("      Writing unstructured grid: %04d", ig);
        }
        if (gridType == GridType.gziptext) {
            ugrid.write_to_gzip_file(gridName(ig, gridType));
        }
        else {
            ugrid.write_to_raw_binary_file(gridName(ig, gridType));
        }
    }

    /* 7. Write out new metadata. */
    // We'll build a local metadata writer at the moment.
    // TODO: consolidate this in one place because presently there is a
    //       Lua implementation and this implementation.
    if (verbosity >= 2) {
        writeln("   Writing metadata for unstructured grids.");
    }
    writeUnstructuredGridMetadata(ugrids, sgridsMetadata, newTags);

    if (verbosity >= 1) {
        writefln("lmr %s: Done.", cmdName);
    }
    return 0;
}

struct CellAndFacePos {
    size_t cellId;
    Vector3 facePos;
}
struct CellAndFaceId {
    size_t cellId;
    size_t faceId;
}

struct BlockAndBoundary {
    size_t blkId;
    size_t bndryId;
}

void createFaceMap(size_t bid, StructuredGrid sgrid, size_t face, ref CellAndFacePos[string] faceMap, ref string[CellAndFaceId][size_t] globalFaceTags) {

   auto bcells = sgrid.get_list_of_boundary_cells(face);

   foreach (id; bcells) {
       auto vtx_list = sgrid.get_vtx_id_list_for_cell(id);
       size_t[] vtx_on_face;
       foreach (ivtx; vtx_list) {
           auto ijk = sgrid.ijk_indices(ivtx);
           final switch(face) {
           case Face.north: if (ijk[1] == sgrid.njv-1) vtx_on_face ~= ivtx; break;
           case Face.south: if (ijk[1] == 0) vtx_on_face ~= ivtx; break;
           case Face.west: if (ijk[0] == 0) vtx_on_face ~= ivtx; break;
           case Face.east: if (ijk[0] == sgrid.niv-1) vtx_on_face ~= ivtx; break;
           case Face.top: if (ijk[2] == sgrid.nkv-1) vtx_on_face ~= ivtx; break;
           case Face.bottom: if (ijk[2] == 0) vtx_on_face ~= ivtx; break;
           }
       }
       Vector3 ctr = Vector3(*(sgrid[vtx_on_face[0]]));
       foreach (i; 1 .. vtx_on_face.length) ctr += *(sgrid[vtx_on_face[i]]);
       ctr /= vtx_on_face.length;
       auto faceTag = makeFaceTag(vtx_on_face);
       faceMap[faceTag] = CellAndFacePos(id, ctr);
       globalFaceTags[bid][CellAndFaceId(id, face)] = faceTag;
   }
}

void convertConnections(JSONValue gridMetadata, StructuredGrid[] sgrids, ref BlockAndCellId[string][size_t] mappedCellsList,
                        ref string[size_t][size_t] newTags, ref string[CellAndFaceId][size_t] globalFaceTags,
                        ref BlockAndBoundary[BlockAndBoundary] connMap)
{
    auto gridConns = gridMetadata["grid-connections"].array;
    foreach (conn; gridConns) {
        auto idA = conn["idA"].integer;
        auto faceA = conn["faceA"].str;
        auto idB = conn["idB"].integer;
        auto faceB = conn["faceB"].str;

        auto fAidx = face_index(faceA);
        auto fBidx = face_index(faceB);

        CellAndFacePos[string] facesA;
        createFaceMap(idA, sgrids[idA], fAidx, facesA, globalFaceTags);
        CellAndFacePos[string] facesB;
        createFaceMap(idB, sgrids[idB], fBidx, facesB, globalFaceTags);

        foreach (tagA, infoA; facesA) {
            string tagFound;
            bool found = false;
            foreach (tagB, infoB; facesB) {
                if (approxEqualVectors(infoA.facePos, infoB.facePos, 1.0e-6, 1.0e-10)) {
                    mappedCellsList[idA][tagA] = BlockAndCellId(idB, infoB.cellId);
                    mappedCellsList[idB][tagB] = BlockAndCellId(idA, infoA.cellId);
                    tagFound = tagB;
                    found = true;
                    break;
                }
            }
            if (!found) {
                string errMsg = format("There was a problem finding a match for face: %s at position %s\n", tagA, infoA.facePos);
                errMsg ~= format("This is in block %d, associated with cell %d\n", idA, infoA.cellId);
                throw new Error(errMsg);
            }
            // Assume we found one
            facesB.remove(tagFound);
        }

        // Assuming we were successful, let's record the new tags
        newTags[idA][fAidx] = "mapped_cells";
        newTags[idB][fBidx] = "mapped_cells";
        connMap[BlockAndBoundary(idA, fAidx)] = BlockAndBoundary(idB, fBidx);
        connMap[BlockAndBoundary(idB, fBidx)] = BlockAndBoundary(idA, fAidx);
    }
}

void writeMappedCellsFile(StructuredGrid[] sgrids, ref BlockAndCellId[string][size_t] mappedCellsList,
                          ref string[size_t][size_t] newTags, ref string[CellAndFaceId][size_t] globalFaceTags,
                          ref BlockAndBoundary[BlockAndBoundary] connMap)
{

    auto of = File(lmrCfg.mappedCellsFile, "w");
    foreach (ig; mappedCellsList.keys()) {
        bool[size_t] nghbrBlks;
        foreach (faceTag; mappedCellsList[ig].keys()) {
            nghbrBlks[mappedCellsList[ig][faceTag].blkId] = true;
        }
        of.writef("MappedBlocks in BLOCK[%d]= ", ig);
        foreach (nblkId; nghbrBlks.keys()) {
            of.writef("%d ", nblkId);
        }
        of.writef("\n");
    }
    foreach (ig, sgrid; sgrids) {
        int nboundaries = (sgrid.dimensions == 3) ? 6 : 4;
        foreach (ib; 0 .. nboundaries) {
            if (ib in newTags[ig]) {
                auto bcells = sgrid.get_list_of_boundary_cells(ib);
                // find src block and boundary
                auto src = connMap[BlockAndBoundary(ig, ib)];
                of.writefln("block-%d-boundary-%d : %d : %d : %d ", ig, ib, src.blkId, src.bndryId, bcells.length);
                foreach (cellId; bcells) {
                    auto faceTag = globalFaceTags[ig][CellAndFaceId(cellId,ib)];
                    of.writefln("%-20s %6d", faceTag, mappedCellsList[ig][faceTag].cellId);
                }
            }
        }
    }
    of.close();
}

/*
 * Keep in synch with lua-modules/grid.lua -- RegisteredGrid:tojson()
 *
 * RJG, 2024-03-06
 * There is a maintenance burden in keeping this in sync.
 * However, this technical debt will need to be paid later.
 * I'm hoping grid metadata is relatively stable over time.
 * We shouldn't need to revisit this function that often.
 */
void writeUnstructuredGridMetadata(UnstructuredGrid[] ugrids, JSONValue[] sgridsMetadata, string[size_t][size_t] newTags)
{
    auto of = File(lmrCfg.gridMetadataFile, "w");
    of.writeln("{");
    of.writefln("  \"ngrids\": %d,", ugrids.length);
    of.writeln("  \"ngridarrays\": 0,");
    of.writeln("  \"grid-connections\": [ ]");
    of.writeln("}");
    of.close();

    size_t ngrids = ugrids.length;
    foreach (ig, ugrid; ugrids) {
        of = File(gridMetadataName(ig), "w");
        of.writeln("{");
        of.writefln("  \"tag\": \"%s\",", sgridsMetadata[ig]["tag"].str);
        of.writefln("  \"fsTag\": \"%s\",", sgridsMetadata[ig]["fsTag"].str);
        of.writefln("  \"type\": \"unstructured_grid\",");
        of.writefln("  \"fieldType\": \"%s\",", sgridsMetadata[ig]["fieldType"].str);
        of.writefln("  \"dimensions\": %d,", ugrid.dimensions);
        of.writefln("  \"nvertices\": %d,", ugrid.nvertices);
        of.writefln("  \"ncells\": %d,", ugrid.ncells);
        of.writefln("  \"nfaces\": %d,", ugrid.nfaces);
        of.writefln("  \"nboundaries\": %d,", ugrid.nboundaries);
        of.writeln("  \"bcTags\": {");
        foreach (itag; 0 .. ugrid.nboundaries) {
            string tag = "";
            if (ngrids > 1) {
                if (itag in newTags[ig]) {
                    tag = newTags[ig][itag];
                }
                else {
                    tag = sgridsMetadata[ig]["bcTags"][face_name[itag]].str;
                }
            } else {
                tag = sgridsMetadata[ig]["bcTags"][face_name[itag]].str;
            }
            of.writefln("    \"%d\": \"%s\",", itag, tag);
        }
        of.writeln("    \"dummy_entry_without_trailing_comma\": \"xxxx\"");
        of.writeln("   },");
        of.writefln("  \"ssTag\": \"%s\",", sgridsMetadata[ig]["ssTag"].str);
        of.writefln("  \"solidModelTag\": \"%s\",", sgridsMetadata[ig]["solidModelTag"].str);
        // [TODO] RJG, 2024-03-06
        // We don't actually have unstructured grids for solid domains,
        // so I don't yet know how to write out the boundary information.
        // On that point, the above two entries are presently superfluous, but won't hurt.
        of.writeln("   \"gridArrayId\": -1");
        of.writeln("}");
        of.close();
    }
}

string gridName(size_t ig, GridType gt)
{
    string name = lmrCfg.gridDirectory ~ "/" ~ lmrCfg.gridPrefix ~ "-" ~ format(lmrCfg.blkIdxFmt, ig);
    if (gt == GridType.gziptext) {
        name ~= lmrCfg.gzipExt;
    }
    return name;
}

string gridMetadataName(size_t ig)
{
    return lmrCfg.gridDirectory ~ "/" ~ lmrCfg.gridPrefix ~ "-" ~ format(lmrCfg.blkIdxFmt, ig) ~ lmrJSONCfg["metadata-extension"].str;
}


