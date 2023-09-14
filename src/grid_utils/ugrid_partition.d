/++    -----------------------------------------------------------------------
       ------------------Unstructured Mesh Partitioner------------------------
       -----------------------------------------------------------------------
++/
/++
 + filename: ugrid_partitioner.d
 + Unstructured mesh partioner -- for use with Metis
 + author: KD
 + version date: 07/11/2016
 +
 + This is a stand alone executable that partitions a global mesh (su2 format)
 + into a number of blocks, and writes them to different
 + su2 formatted files, for use with the Eilmer4 flow solver.
 +
 + INSTRUCTIONS:
 +
 + To use this program, first download and install Metis 5.1.0:
 + http://glaros.dtc.umn.edu/gkhome/metis/metis/download
 +
 + To compile the partitioner:
 + $ dmd partition_core.d
 +
 + To partition a mesh:
 + $ ugrid_partition mesh_fileName mappedCells_fileName nPartitions nDim
 +
 + where,
 + mesh_fileName is the mesh file name, i.e. mesh.su2
 + mappedCells_fileName is the user defined name given to the file used to store
 +                      the ghost cell connectivity along METIS_INTERIOR boundaries
 + nPartitions is the number of partitions
 + nDim is the minimum number of nodes that constitutes a face,
 +      i.e. 2 for a line (2D), 3 for a triangle (3D)
 +
 + example:
 + $ ./ugrid_partition mesh_file.su2 mapped_cells 4 2
++/

import std.stdio;
import std.conv;
import std.file;
import std.algorithm;
import std.string;
import std.array;
import std.format;
import std.math;
import std.process;
import std.getopt;
import std.path;

void printHelp()
{
    writeln("Usage: ugrid_partition");
    writeln("");
    writeln(" > ugrid_parition grid-file.su2 mapped-cells.txt nParitions nDim");
    writeln("");
    writeln("   where:");
    writeln("   grid-file.su2      : name of grid file in SU2 format.");
    writeln("   mapped-cells.txt   : name of output file for mapped cells.");
    writeln("   nPartitions        : integer number of desired partitions.");
    writeln("   nDim               : integer (2 or 3) for dimensionality of grid.");
    writeln("");
}

// -----------------
// Class definitions
// -----------------

class Domain {

public:

    size_t dimensions;          // dimensionality of the grid (e.g. 2D/3D)
    size_t nblocks;             // number of blocks in the domain
    Cell[] cells;               // collection of cells in the domain
    Node[] nodes;               // collection of nodes in the domain
    Boundary[] boundary;        // collection of boundaries that surround the domain

    this (size_t nblocks, size_t dimensions) {
        this.nblocks = nblocks;
        this.dimensions = dimensions;
    }

} // end class Block

// -------------------------------------------------------------------------------------

class Block {

public:

    size_t id;                  // block identifier
    Cell[] cells;               // collection of cells in the block
    MappedCells[] mapped_cells; // collection of cells along any internal boundaries
    Node[] nodes;               // collection of nodes in the block
    size_t[size_t] global2local_node_transform; // a dictionary which uses the global node id
                                                // as an index to retrieve the local id
    size_t[size_t] global2local_cell_transform; // a dictionary which uses the global cell id
                                                // as an index to retrieve the local id
    Boundary[] boundary; // collection of boundary conditions from the global mesh
    bool[] bc_exists;    // Each entry represents a boundary condition
                         // from the global mesh. If the boundary condition also exists for this
                         // block the entry will be True, False otherwise

    this (size_t id) {
        this.id = id;
    }

} // end class Block

// -------------------------------------------------------------------------------------

class MappedCells {

public:

    size_t global_cell_id;            // the primary cells global identifier
    size_t global_neighbour_cell_id;  // the neighbour cells global cell identifier
    size_t block_id;                  // cell's block identifier
    size_t neighbour_block_id;        // neighbour cells's block identifier
    string faceTag;                   // the unique identifier for the shared face
                                      // constructed with the boundary face's
                                      // local (to the primary cell) node ids
                                      // in ascending order

    this (size_t global_cell_id,
          size_t global_neighbour_cell_id,
          size_t block_id,
          size_t neighbour_block_id,
          string faceTag) {

        this.global_cell_id = global_cell_id;
        this.global_neighbour_cell_id = global_neighbour_cell_id;
        this.block_id = block_id;
        this.neighbour_block_id = neighbour_block_id;
        this.faceTag = faceTag;
    }

} // end class MappedCells

// -------------------------------------------------------------------------------------

class Cell {

public:

    size_t id;                 // cell identifier
    size_t type;               // cell type (i.e. quad, tri etc. -- for numbering refer to su2 mesh format)
    size_t[] cell_cloud_ids;   // list of cell ids for the cells that share an interface with this cell
    size_t[] node_ids;         // list of node ids that make up the cell
    size_t[][] face_node_ids;  // list of the node ids that make up each face attached to the cell
    string[] face_tag_list;    // the unique identifiers for the faces that make up the cell
    size_t block_id;           // the unique identifier for what partition/block this cell is a part of
    size_t nfaces;             // number of faces attached to a cell

    this(size_t id,
         size_t type,
         size_t[] node_ids,
         size_t dimensions) {

        this.id = id;
        this.type = type;
        this.node_ids = node_ids.dup();
        this.nfaces = get_nfaces(dimensions);
    }

    void add_face_to_cell(size_t[] corners) {

        size_t[] my_face_node_ids;
        foreach(i; corners) { my_face_node_ids ~= node_ids[i]; }
        string faceTag = makeFaceTag(my_face_node_ids);
        face_node_ids ~= my_face_node_ids;
        face_tag_list ~= faceTag;
        return;

    } // end add_face_to_cell()

    void construct_faces(size_t dimensions) {
        // Attach the faces to each cell. In 2D, faces are defined as lines.
        // As we progress along the line the face normal is pointing to the right.
        // In 3D, a counter-clockwise cycles of points plus the right-hand rule
        // define the face normal. Whether a face points out of or into a cell
        // will be determined and remembered when we add the face to the cell.
        if (dimensions == 2) {
            switch(type) {
            case 5: // triangle
                add_face_to_cell([0,1]);
                add_face_to_cell([1,2]);
                add_face_to_cell([2,0]);
                break;
            case 9: // quad
                add_face_to_cell([2,3]); // north
                add_face_to_cell([1,2]); // east
                add_face_to_cell([0,1]); // south
                add_face_to_cell([3,0]); // west
                break;
            default:
                throw new Exception("invalid cell type in 2D");
            }
        } else {
            assert(dimensions == 3, "invalid dimensions");
            switch(type) {
            case 10: // tetra
                add_face_to_cell([0,1,2]);
                add_face_to_cell([0,1,3]);
                add_face_to_cell([1,2,3]);
                add_face_to_cell([2,0,3]);
                break;
            case 12: // hexa
                add_face_to_cell([2,3,7,6]); // north
                add_face_to_cell([1,2,6,5]); // east
                add_face_to_cell([1,0,4,5]); // south
                add_face_to_cell([0,3,7,4]); // west
                add_face_to_cell([4,5,6,7]); // top
                add_face_to_cell([0,1,2,3]); // bottom
                break;
            case 13: // wedge
                add_face_to_cell([0,1,2]);
                add_face_to_cell([3,4,5]);
                add_face_to_cell([0,2,5,3]);
                add_face_to_cell([0,3,4,1]);
                add_face_to_cell([1,4,5,2]);
                break;
            case 14: // pyramid
                add_face_to_cell([0,1,2,3]);
                add_face_to_cell([0,1,4]);
                add_face_to_cell([1,2,4]);
                add_face_to_cell([2,3,4]);
                add_face_to_cell([3,0,4]);
                break;
            default:
                throw new Exception("invalid cell type in 3D");
            }
        }
        return;
    } // end add_faces_to_cell()

    size_t get_nfaces(size_t dimensions) {
        size_t nfaces;
        if (dimensions == 2) {
            switch(type) {
            case 5: // triangle
                nfaces = 3;
                break;
            case 9: // quad
                nfaces = 4;
                break;
            default:
                throw new Exception("invalid cell type in 2D");
            }
        } else {
            assert(dimensions == 3, "invalid dimensions");
            switch(type) {
            case 10: // tetra
                nfaces = 4;
                break;
            case 12: // hexa
                nfaces = 6;
                break;
            case 13: // wedge
                nfaces = 5;
                break;
            case 14: // pyramid
                nfaces = 5;
                break;
            default:
                throw new Exception("invalid cell type in 3D");
            }
        }
        return nfaces;
    } // end get_nfaces()

} // end class Cell

// -------------------------------------------------------------------------------------

class Node {

public:

    size_t id;                    // node identifier
    double[3] pos;                // node position, (x,y,z-Coordinates)
    size_t[] partitions;          // tracks what blocks the node has
                                  // been added to already

    this(size_t id,
         double[3] pos,
         size_t nparts) {

        this.id = id;
        this.pos = pos.dup();
    }

} // end class Node

// -------------------------------------------------------------------------------------

class Boundary {

public:

    size_t id;               // boundary indentifier (note: METIS_INTERIOR is a special boundary
                             // type that will always be assigned the 0 id
    string tag;              // boundary condition type name
    string[] face_tag_list;  // the unique identifiers for the faces that make up the boundary
    size_t[][] face_node_ids;     // list of node ids that make up each face along the boundary
    size_t nfaces;           // number of faces along the boundary

    this(size_t id,
         string tag,
         ref size_t[][] face_node_ids,
         ref string[] face_tag_list) {

        this.id = id;
        this.tag = tag;
        this.face_node_ids = face_node_ids.dup();
        this.face_tag_list = face_tag_list.dup();
    }

} // end class Boundary

// --------------------
// Function definitions
// --------------------

bool binarySearch(string x, ref string[] container) {
    // A simple binary search algorithm : O(log n)
    int low = 0;
    int high = to!int(container.length) - 1;
    while (low <= high) {
        int mid = (low + high)/2;
        string item = container[mid];
        if (x == item) { return true; }         // x exists!
        else if (x < item) { high = mid - 1; }
        else { low = mid + 1; }
    }
    return false;                               // x is not there
} // end binarySearch

void insertion_sort(int[][] arr, int n) {
    // a helper function that applies insertion sort to a
    // list of tuples, the sorting is based on the first index
    int[2] key;
    int j;
    foreach (i; 1..n) {
        key = arr[i];
        j = i - 1;
        while (j >= 0 && arr[j][1] > key[1]) {
            arr[j + 1] = arr[j];
            j = j - 1;
        }
        arr[j + 1] = key.dup;
    }
}


string makeFaceTag(const size_t[] node_id_list) {
    // We make a tag for this face out of the vertex id numbers
    // sorted in ascending order, to be used as a key in an
    // associative array of indices.  Any cycle of the same vertices
    // should define the same face, so we will use this tag as
    // a way to check if we have made a particular face before.
    // Of course, permutations of the same values that are not
    // correct cycles will produce the same tag, however,
    // we'll not worry about that for now because we should only
    // ever be providing correct cycles.
    size_t[] my_id_list = node_id_list.dup();
    sort(my_id_list);
    string tag = "";
    size_t n = my_id_list.length;
    foreach(i; 0 .. n) {
        if (i > 0) { tag ~= "-"; }
        tag ~= format("%d", my_id_list[i]);
    }
    return tag;
} // end makeFaceTag

void cleanDir(string dualFile, string partitionFile, string metisFormatFile) {
    // Removes the files generated by Metis

    remove(dualFile);
    remove(partitionFile);
    remove(metisFormatFile);

} // end CleanDir

int metisCheck() {
    // Checks if Metis is available

    string command = "gpmetis";
    auto output = executeShell(command);
    int status = output[0];
    return status;

} // end metisCheck

string mesh2dual(string fileName, int ncommon) {
    // Calls Metis to convert the mesh into Metis' dual format

    writeln("-- Converting mesh to dual format");
    string outputFileName = "dual_"~fileName;
    string ncommon_arg = to!string(ncommon);
    string command = "m2gmetis " ~ " " ~ fileName ~ " " ~ outputFileName ~ " -gtype=dual -ncommon=" ~ ncommon_arg;
    executeShell(command);
    return outputFileName;

} // end mesh2dual

string partitionDual(string fileName, int nparts) {
    // Calls Metis to partition the mesh stored in Metis' dual format

    writeln("-- Partitioning dual format");
    string outputFile = fileName ~ ".part." ~ to!string(nparts);
    string nparts_arg = to!string(nparts);
    string command = "gpmetis " ~ " " ~ fileName ~ " " ~ nparts_arg;
    executeShell(command);
    return outputFile;

} // end partitionDual

string SU2toMetisMeshFormat(string fileName, int ncommon) {
    // Preprocesses the mesh stored in SU2 format for use with Metis
    // input: SU2 mesh stored as a .su2 file
    // output: file storing mesh in Metis format

    writeln("-- Converting SU2 mesh to Metis mesh format");
    auto f = File(fileName, "r");

    // Check that the nDim command line argument matches the file dimensionality (NNG 07/21)
    // Users setting nDim incorrectly can result in some very strange (and sometimes silent) bugs.
    int ndim=-1;
    while (!f.eof) {
        auto line = f.readln().strip();
        if (canFind(line, "NDIME")) {
            auto tokens = line.split("=");
            ndim = to!int(tokens[1].strip());
            break;
        }
    }
    if (ndim==-1) throw new Error(format("NDIME field not found in .su2 file %s",  fileName));
    if (ndim!=ncommon) {
        string msg = format("ugrid_parition dimensionality %d does not match file %s with NDIME= %d",
                             ncommon, fileName, ndim);
        throw new Error(msg);
    }

    size_t ncells;
    while (!f.eof) {
        auto line = f.readln().strip();
        if (canFind(line, "NELEM")) {
            auto tokens = line.split("=");
            ncells = to!size_t(tokens[1].strip());
            break;
        }
    }
    string[] cells;
    cells.length = ncells;
    foreach(i; 0 .. ncells) {
        auto lineContent = f.readln().strip();
        auto tokens = lineContent.split();
        size_t idx = to!size_t(tokens[tokens.length-1]);
        string[] vtxids;
        foreach(j; 1 .. tokens.length-1) { vtxids ~= format("%d \t", to!int(tokens[j])+1); }
        cells[idx] = vtxids.join();
    }

    // Strip out just the filename, it may be in another directory. (NNG 22/06/22)
    string baseFileName = baseName(fileName);
    string outputFileName = "metisFormat_"~baseFileName;
    File outFile;
    outFile = File(outputFileName, "w");
    outFile.writeln(ncells);
    foreach(cell; cells){
        outFile.writeln(cell);
    }
    return outputFileName;

} // end SU2toMetisFormat

void readSU2grid(string meshFile, ref Domain grid) {
    // Proceed through SU2 file and read in the cell, node, and boundary information
    // input: SU2 mesh stored as a .su2 file
    // output: Domain grid data filled

    writeln("-- -- Import SU2 mesh");
    auto f = File(meshFile, "r");
    string getHeaderContent(string target)
    // Helper function to proceed through file, line-by-line,
    // looking for a particular header line.
    // Returns the content from the header line and leaves the file
    // at the next line to be read, presumably with expected data.
    {
        while (!f.eof) {
            auto line = f.readln().strip();
            if (canFind(line, target)) {
                auto tokens = line.split("=");
                return tokens[1].strip();
            }
        } // end while
        return ""; // didn't find the target
    }

    // general information
    grid.dimensions = to!int(getHeaderContent("NDIME"));
    auto ncells = to!size_t(getHeaderContent("NELEM"));

    // Cells
    grid.cells.length = ncells;
    foreach(i; 0 .. ncells) {
        auto lineContent = f.readln().strip();
        auto tokens = lineContent.split();
        int cell_type = to!int(tokens[0]);
        size_t cell_id = to!size_t(tokens[$-1]);
        size_t[] my_node_id_list;
        foreach(j; 1 .. tokens.length-1) { my_node_id_list ~= to!size_t(tokens[j]); }
        grid.cells[cell_id] = new Cell(cell_id, cell_type, my_node_id_list, grid.dimensions);
    }

    // Nodes
    auto nnodes = to!size_t(getHeaderContent("NPOIN"));
    grid.nodes.length = nnodes;
    foreach(i; 0 .. nnodes) {
        auto tokens = f.readln().strip().split();
        double x=0.0; double y=0.0; double z = 0.0; size_t indx = 0;
        if (grid.dimensions == 2) {
            x = to!double(tokens[0]);
            y = to!double(tokens[1]);
            indx = to!size_t(tokens[2]);
        } else {
            assert(grid.dimensions == 3, "invalid dimensions");
            x = to!double(tokens[0]);
            y = to!double(tokens[1]);
            z = to!double(tokens[2]);
            indx = to!size_t(tokens[3]);
        }
        double[3] pos = [x, y, z];
        grid.nodes[indx] = new Node(indx, pos, grid.nblocks);
    } // end foreach i .. nnodes

    // Boundaries
    auto nboundaries = to!size_t(getHeaderContent("NMARK"));
    grid.boundary.length = nboundaries;
    foreach(i; 0 .. nboundaries) {
        string tag = getHeaderContent("MARKER_TAG");
        writeln("-- -- -- Found boundary i=", i, " tag=", tag);
        size_t[][] face_node_ids;
        string[] face_tag_list;
        size_t nelem = to!size_t(getHeaderContent("MARKER_ELEMS"));
        face_node_ids.length = nelem;
        face_tag_list.length = nelem;
        foreach(j; 0 .. nelem) {
            auto tokens = f.readln().strip().split();
            int vtk_type = to!int(tokens[0]);
            size_t[] my_face_node_id_list;
            foreach(k; 1 .. tokens.length) { my_face_node_id_list ~= to!size_t(tokens[k]); }
            string face_tag = makeFaceTag(my_face_node_id_list);
            face_tag_list[j] = face_tag;
            face_node_ids[j] ~= my_face_node_id_list;
        } // end foreach j

        // we will use a binary search on the face_tag_list, so we should order this array
        face_tag_list.sort!();
        grid.boundary[i] = new Boundary(i+1, tag, face_node_ids, face_tag_list);
    } // end foreach i
    return;
}

void constructGridBlocks(bool reorder, string meshFile, string mappedCellsFilename, string partitionFile, string dualFile, ref Domain grid, ref Block[] gridBlocks) {
    // Construct the partitioned grid blocks
    // input: Domain grid, Metis output files, empty grid blocks
    // output: grid block data filled

    writeln("-- Building ", grid.nblocks, " blocks in SU2 format");

    // Assign each cell the partition identifier it is a member of
    writeln("-- -- Reading Metis output file");
    auto f = File(partitionFile, "r");
    foreach(ref cell; grid.cells) {
        auto tokens = f.readln().strip().split();
        size_t partition_id = to!int(tokens[0]);
        cell.block_id = partition_id;
    }

    // Attach neighbour cells using cell interconnectivity
    writeln("-- -- Reading cell interconnectivity from dual format file");
    f = File(dualFile, "r");
    f.readln(); // skip first line
    foreach(ref cell; grid.cells) {
         auto lineContent = f.readln().strip();
         auto tokens = lineContent.split();
         foreach(j; 0 .. tokens.length) { cell.cell_cloud_ids ~= to!size_t(tokens[j]) - 1 ; }
    }

    // Construct blocks
    writeln("-- -- Data reading complete, proceed to building blocks");
    gridBlocks.length = grid.nblocks;
    foreach(id, ref blk; gridBlocks) {
        blk = new Block(id);
        // create an INTERIOR boundary for all blocks (we will fill this out later)
        size_t[][] face_node_id_list;
        string[] face_tag_list;
        blk.boundary ~= new Boundary(0, "METIS_INTERIOR", face_node_id_list, face_tag_list);
    }

    // Attach a copy of each global boundary condition to the blocks
    // note that only relevant boundary conditions will be populated for a block later
    foreach (ref blk; gridBlocks) {
        blk.bc_exists.length = grid.boundary.length;
        foreach( bcid, gbndry; grid.boundary) {
            blk.bc_exists[bcid] = false;
            size_t[][] face_node_id_list;
            string[] face_tag_list;
            blk.boundary ~= new Boundary(gbndry.id+1, gbndry.tag, face_node_id_list, face_tag_list);
        }
    }

    // Loop through cells and assign them to correct blocks, as well as interior boundary faces
    foreach(cid, ref cell; grid.cells) {

        // Attach cell to block
        size_t blk_id = cell.block_id;
        gridBlocks[blk_id].cells ~= cell;
        gridBlocks[blk_id].global2local_cell_transform[cid] = gridBlocks[blk_id].cells.length - 1;
        foreach(node_id;cell.node_ids) {
            if (!grid.nodes[node_id].partitions.canFind(blk_id)) {
                grid.nodes[node_id].partitions ~= blk_id;
                gridBlocks[blk_id].nodes ~= grid.nodes[node_id];
                gridBlocks[blk_id].global2local_node_transform[node_id] = gridBlocks[blk_id].nodes.length - 1;
            }
        }
    }

    // we pause the construction of the blocks here to apply an optional local reordering of the cells
    if (reorder) {
        writeln("-- -- -- applying RCM reordering");

        // we apply the reverse Cuthill-McKee (RCM) reordering from pg. 66 of
        //     A. George, J. Liu, E. Ng
        //     Computer Solution of Sparse Linear Systems
        //     Oak Ridge National Laboratory (1994)

        foreach (blk; gridBlocks) {

            // find a starting cell (we want this cell to be a minimum degree node)
            int seed_id = to!int(blk.cells[0].id);
            int min_degree = int.max;
            foreach (cell; blk.cells) {
                int degree = 0;
                // calculate the degree of the current cell (i.e. sum up it's local connections)
                foreach (ncell_id; cell.cell_cloud_ids) {
                    Cell ncell = grid.cells[ncell_id];
                    if (ncell.block_id == cell.block_id) { degree += 1; }
                }
                // check if it has a smaller degree than the current seed cell
                if (degree < min_degree) {
                    seed_id = to!int(cell.id);
                    min_degree = degree;
                }
            }

            // apply the Cuthill-McKee algorithm to determine the new ordering of the local cells
            int[] ordered_cell_ids; // note that the stored cells ids are the unique global cell ids
            int[] visited;
            ordered_cell_ids ~= to!int(seed_id); // begin with the starting node we found above
            size_t n = blk.cells.length;
            while (ordered_cell_ids.length < n) {
                foreach (cell_id; ordered_cell_ids) {
                    Cell cell = grid.cells[cell_id];
                    if (!visited.canFind(cell_id)) {
                        int[][] id_degree; // a list of cell ids and their corresponding degree
                        foreach (ncell_id; cell.cell_cloud_ids) {
                            Cell ncell = grid.cells[ncell_id];
                            if (ncell.block_id == cell.block_id && !ordered_cell_ids.canFind(ncell_id)) { // we only proceed with local cells that haven't been already added
                                int degree = 0;
                                // calculate the degree of ncell by looping through the neighbours of ncell and adding +1 for each cell in the same block
                                foreach(nncell_id; ncell.cell_cloud_ids) {
                                    Cell nncell = grid.cells[nncell_id]; // neighbour cell of ncell
                                    if (nncell.block_id == ncell.block_id) { degree += 1; }
                                }
                                id_degree ~= [to!int(ncell_id), degree];
                            }
                        }
                        // sort the list in order of ascending degree
                        insertion_sort(id_degree, to!int(id_degree.length));
                        // add the ids in the correct order
                        if (id_degree.length > 0) {
                            foreach (dat; id_degree) { ordered_cell_ids ~= dat[0]; }
                        }
                        visited ~= cell_id;
                    }
                }
            } // end while()

            // create temporary array with the new ordering
            Cell[] reordered;
            foreach (idx; ordered_cell_ids) { reordered ~= grid.cells[idx]; }

            // overwrite the cell array with the reverse order (since we apply reverse Cuthill-McKee)
            gridBlocks[blk.id].cells = [];
            int id = 0;
            foreach_reverse (cell; reordered) {
                gridBlocks[blk.id].global2local_cell_transform[cell.id] = to!int(id);
                id += 1;
                gridBlocks[blk.id].cells ~= cell;
            }

        } // end foreach (blk; gridBlocks)
    } // end if (reorder)

    // continue contructing the blocks...
    foreach(cid, ref cell; grid.cells) {
        size_t blk_id = cell.block_id;
        // Attach the cell faces to the METIS_INTERIOR boundary if applicable
        foreach(adj_cid; cell.cell_cloud_ids) {
            Cell adj_cell = grid.cells[adj_cid]; // adjacent cell
            if (adj_cell.block_id != cell.block_id) {
                // populate face_node_ids for each cell (ONLY if it hasn't already been done)
                if (adj_cell.face_node_ids.length < 1) { adj_cell.construct_faces(grid.dimensions); }
                if (cell.face_node_ids.length < 1) { cell.construct_faces(grid.dimensions); }
                foreach(i, cell_face_tag; cell.face_tag_list) {
                    if (adj_cell.face_tag_list.canFind(cell_face_tag)) {
                        gridBlocks[blk_id].boundary[0].nfaces += 1; // METIS_INTERIOR is in pos 0
                        // facetag needs to be stored with the local node ids
                        size_t[] local_node_id_list;
                        foreach(node;cell.face_node_ids[i]) { local_node_id_list ~= gridBlocks[cell.block_id].global2local_node_transform[node]; }
                        string face_tag = makeFaceTag(local_node_id_list);
                        // store mapped cell information
                        gridBlocks[blk_id].mapped_cells ~= new MappedCells(cell.id, adj_cell.id, cell.block_id, adj_cell.block_id, face_tag);
                        gridBlocks[blk_id].boundary[0].face_node_ids ~= cell.face_node_ids[i];
                    }
                }
            }
        }

        // Attach the cell faces to any other boundary conditions if applicable
        size_t nexpect = cell.nfaces - cell.cell_cloud_ids.length; // expected boundary faces
        if (nexpect > 0) {
            size_t nfound = 0;
            // populate face_node_ids (ONLY if it hasn't already been done)
            if (cell.face_node_ids.length < 1) { cell.construct_faces(grid.dimensions); }
            foreach(fi, face_tag; cell.face_tag_list) {
                if (nfound == nexpect) break;
                foreach(i; 0 .. grid.boundary.length) {
                    bool face_found = binarySearch(face_tag, grid.boundary[i].face_tag_list);
                    if (face_found) {
                        gridBlocks[blk_id].boundary[grid.boundary[i].id].face_node_ids ~= cell.face_node_ids[fi];
                        nfound += 1;
                        break;
                    }
                }
            }
        }
    }
} // end metis2su2Format

void writeGridBlockFiles(string meshFile, string mappedCellsFilename, Block[] gridBlocks, size_t dimensions) {
    // Write out the grid blocks in separate .su2 files
    // note that cells and nodes should write out  their local ids

    File outFile_mappedcells;
    outFile_mappedcells = File(mappedCellsFilename, "w");

    // If the user is partitioning more than one su2 file,
    // then we need to count the current number of partitions (su2 files) to
    // know where to start the block counting for the new partitions
    int nCurrentBlocks=0;

    string command = "ls -d block_*_* | wc -l";
    auto output = executeShell(command);
    if ( !(canFind(output[1], "ls: cannot access") | canFind(output[1], "No such file or directory")) ) {
        // Assume that we were successful in finding a number of blocks for partitioning.
        nCurrentBlocks = to!int(splitLines(output[1])[0]);
    }

    // first we add the mapped block information
    foreach(i;0..gridBlocks.length) {
        string mapped_blocks_tag = "MappedBlocks in BLOCK[" ~ to!string(i+nCurrentBlocks) ~ "]= ";
        outFile_mappedcells.writef(mapped_blocks_tag);

        // let's begin by writing out the mapped cells for the block
        size_t[] added_block_id_list;
        foreach(mapped_cell; gridBlocks[i].mapped_cells){
            size_t neighbour_block_id = mapped_cell.neighbour_block_id+nCurrentBlocks;
            if (!added_block_id_list.canFind(neighbour_block_id)) {
                outFile_mappedcells.writef("%s ", to!string(neighbour_block_id));
                added_block_id_list ~= neighbour_block_id;
            }
        }
        outFile_mappedcells.writef(" \n");
    }
    nCurrentBlocks=0; // reset block count for reparsing of data

    // Strip out the filename. The file itself may be in another directory (NNG 22/06/22)
    string meshFileName = baseName(meshFile);
    // We are ready to construct those partition files!
    foreach(i;0..gridBlocks.length) {
        writeln("-- Writing out block: #", i+nCurrentBlocks);

        // let's begin by writing out the mapped cells for the block
        string mapped_cells_tag = "NMappedCells in BLOCK[" ~ to!string(i+nCurrentBlocks) ~ "]= ";
        outFile_mappedcells.writef("%s", mapped_cells_tag);
        auto ninterior = gridBlocks[i].boundary[0].nfaces; // METIS_INTERIOR is in pos 0
        outFile_mappedcells.writef("%d \n", ninterior);
        foreach(mapped_cell; gridBlocks[i].mapped_cells){
            auto primary_cell_local_id = gridBlocks[i].global2local_cell_transform[mapped_cell.global_cell_id];
            auto secondary_cell_local_id = gridBlocks[mapped_cell.neighbour_block_id].global2local_cell_transform[mapped_cell.global_neighbour_cell_id];
            outFile_mappedcells.writef("%s \t", mapped_cell.faceTag);
            outFile_mappedcells.writef("%d \t", mapped_cell.neighbour_block_id+nCurrentBlocks);
            outFile_mappedcells.writef("%d \n", secondary_cell_local_id);
        }

        string outputFileName = "block_" ~ to!string(i+nCurrentBlocks) ~ "_" ~ meshFileName;
        File outFile;
        outFile = File(outputFileName, "w");
        outFile.writeln("%");
        outFile.writeln("% Problem dimension");
        outFile.writeln("%");
        outFile.writef("NDIME= %d \n", dimensions);
        outFile.writeln("%");
        outFile.writeln("% Inner element connectivity");
        outFile.writeln("%");
        outFile.writef("NELEM= %d \n", gridBlocks[i].cells.length);

        foreach(cell_id;0 .. gridBlocks[i].cells.length) {
            outFile.writef("%d \t", gridBlocks[i].cells[cell_id].type);
            foreach(j;0 .. gridBlocks[i].cells[cell_id].node_ids.length) {
                auto local_node_id = gridBlocks[i].global2local_node_transform[gridBlocks[i].cells[cell_id].node_ids[j]];
                outFile.writef("%d \t", local_node_id);
            }
            outFile.writef("%d \n", cell_id);
        }

        outFile.writeln("%");
        outFile.writeln("% Node coordinates");
        outFile.writeln("%");
        outFile.writef("NPOIN= %d \n", gridBlocks[i].nodes.length);
        foreach(local_node_id;0 .. gridBlocks[i].nodes.length) {
            foreach(j;0 .. dimensions) { // blocks[i].nodes[node_id].pos.length)
                outFile.writef("%.17f \t", gridBlocks[i].nodes[local_node_id].pos[j]);
            }
            outFile.writef("%d \n", local_node_id);
        }

        outFile.writeln("%");
        outFile.writeln("% Boundary elements");
        outFile.writeln("%");
        size_t nbndry = 0;
        foreach (bndry; gridBlocks[i].boundary) {
            if (bndry.face_node_ids.length > 0) { nbndry += 1; }
        }
        outFile.writef("NMARK= %d \n", nbndry); //blocks[i].boundary.length);

        foreach(bndary; gridBlocks[i].boundary) {
            if (bndary.face_node_ids.length > 0) {
                outFile.writef("MARKER_TAG= %s \n", bndary.tag);
                outFile.writef("MARKER_ELEMS= %s \n", bndary.face_node_ids.length);
                foreach(face_node_ids; bndary.face_node_ids) {
                    switch(face_node_ids.length) {
                    case 2: // line
                        outFile.writef("3 \t");
                        break;
                    case 3: // tri
                        outFile.writef("5 \t");
                        break;
                    case 4: // rectangle
                        outFile.writef("9 \t");
                        break;
                    default:
                        throw new Exception("invalid element type");
                    }
                    foreach(node_id; face_node_ids) {
                        auto local_node_id = gridBlocks[i].global2local_node_transform[node_id];
                        outFile.writef("%d \t", local_node_id);
                    }
                    outFile.writef("\n");
                }
            }
        }
    }
} // end writeSU2format

int main(string[] args){

    string inputMeshFile; string mappedCellsFilename; int nparts; int ncommon; bool reorder;
    // we will accept either 5 or 6 command line arguments
    if (args.length != 6) {
        if (args.length != 5) {
            writeln("Wrong number of command line arguments.");
            printHelp();
            return 1;
        } else {
            // we assume the user does not want the grid reordered if they didn't specify a 6th argument
            args ~= "false";
        }
    }

    // assign command line arguments to variable names
    inputMeshFile = args[1];
    mappedCellsFilename = args[2];
    nparts = to!int(args[3]);
    ncommon = to!int(args[4]);
    reorder = to!bool(args[5]);

    writeln("Begin partitioning................");
    // import entire SU2 grid file
    Domain grid = new Domain(nparts, ncommon);
    readSU2grid(inputMeshFile, grid);

    // Let's check the user doesn't have any previous partitioned su2 files in the directory
    // since this will upset the code if the user intends on partitioning more than one su2 file
    string commandCheck = "if [ -f block_0_"~inputMeshFile~" ] ; then echo 'yes' ; else echo 'no' ; fi";
    auto fileCheck = executeShell(commandCheck);
    if (fileCheck[1] != "no\n")  throw new Error("It appears the partitioned su2 files already exist. Please remove these files before running the partitioner again.");

    // check if Metis is installed correctly
    if (metisCheck() != 0)  throw new Error("Metis partition software cannot be found.");

    // convert the SU2 grid into the mesh format Metis is expecting
    string metisFormatFile = SU2toMetisMeshFormat(inputMeshFile, ncommon);

    // convert the mesh into a dual graph using Metis
    string dualFormatFile = mesh2dual(metisFormatFile, ncommon);

    // partition the dual graph using Metis
    string partitionFile;
    if (nparts > 1) {
        partitionFile = partitionDual(dualFormatFile, nparts);
    } else {
        // we need to construct our own partition file in cases when the user wants to reorder a grid without partitioning it
        partitionFile = dualFormatFile ~ ".part." ~ to!string(nparts);
        File outFile;
        outFile = File(partitionFile, "w");
        foreach(cell; 0 .. grid.cells.length) {
            outFile.writef("0 \n");
        }
    }

    // divide grid into N grid blocks
    Block[] gridBlocks;
    constructGridBlocks(reorder, inputMeshFile, mappedCellsFilename, partitionFile, dualFormatFile, grid, gridBlocks);

    // write the grid blocks out to file in SU2 format
    writeGridBlockFiles(inputMeshFile, mappedCellsFilename, gridBlocks, grid.dimensions);

    // clean up the temporary Metis files
    cleanDir(dualFormatFile, partitionFile, metisFormatFile);

    writeln("Finished partitioning................");
    return 0;
}
