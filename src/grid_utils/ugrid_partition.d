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
 + $ ./partition_core mesh_file.su2 4 2
++/

import std.stdio;
import std.conv;
import std.file;
import std.algorithm;
import std.string;
import std.array;
import std.conv;
import std.stdio;
import std.algorithm;
import std.format;
import std.math;
import std.process;
import std.getopt;

// -------------------------------------------------------------------------------------
// Classes
// -------------------------------------------------------------------------------------

class Block {
public:
    size_t id;    // block id number
    Cell[] cells; // collection of cells in block
    Node[] nodes; // collection of nodes in block
    MappedCells[] mapped_cells; // collection of mapped_cells in block
    size_t[] node_id_list; // collection of the node ids
                           // although we already have a list of the nodes
                           // having a list of the ids will make for a
                           // fast means of checking if a node is in a block
    size_t[size_t] transform_list; // a dictionary which uses the global node id
                                   // as an index to retrieve the local id
    size_t[size_t] transform_cell_list; // a dictionary which uses the global cell id
                                   // as an index to retrieve the local id
    Boundary[string] boundary; // collection of block boundaries attached to block
    string[] bc_names; // list of boundary condition types for a block
    this(size_t id) {
	this.id = id;
    }
} // end class Block
// -------------------------------------------------------------------------------------
class MappedCells {
public:
    size_t global_cell_id;            // current blocks local cell id number
    size_t global_neighbour_cell_id;  // neighbour blocks local cell id number
    size_t block_id;                  // current blocks id
    size_t neighbour_block_id;        // neighbour blocks id
    string faceTag;                   // unique face identifier constructed with the faces local node ids in ascending order
    this(size_t global_cell_id, size_t global_neighbour_cell_id, size_t block_id, size_t neighbour_block_id, string faceTag) {
	this.global_cell_id = global_cell_id;
	this.global_neighbour_cell_id = global_neighbour_cell_id;
	this.block_id = block_id;
	this.neighbour_block_id = neighbour_block_id;
	this.faceTag = faceTag;
    }
}
// -------------------------------------------------------------------------------------
class Cell {
public:
    size_t id;   // cell id number
    size_t type; // cell type (i.e. quad, tri etc. -- for numbering refer to su2 format)
    size_t[] face_id_list; // list of ids for the faces attached to a cell
    size_t[] node_id_list;  // list of ids for of the nodes attached to a cell
    size_t[] cell_neighbour_id_list; // list of cell ids of neighbouring cells
    size_t partition;
    
    this(size_t id, size_t type, size_t[] node_id_list, size_t[] face_id_list) {
	this.id = id;
	this.type = type;
	this.node_id_list = node_id_list.dup();
	this.face_id_list = face_id_list.dup();
    }
}
// -------------------------------------------------------------------------------------
class Face {
public:
    size_t id;       // Face id number
    size_t[] node_id_list;   // list of ids for the faces attached to a face
    size_t[] partition_id_list; // list of partition ids for the cells attached to face
    this(size_t[] node_id_list) {
	this.node_id_list = node_id_list.dup();
    }
}
// -------------------------------------------------------------------------------------
class Node {
public:
    size_t id;       // node id number
    double[3] pos;  // node position, (x,y,z-Coordinates)

    this(size_t id, double[3] pos) {
	this.id = id;
	this.pos = pos.dup();
    }
}
// -------------------------------------------------------------------------------------
class Boundary {
public:
    string tag;     // Boundary condition type name
    size_t[] face_id_list;  // list of ids for the faces that make up a domain boundary
    
    this(string tag, size_t[] face_id_list)
    {
	this.tag = tag;
	this.face_id_list = face_id_list.dup();
    }
}
// -------------------------------------------------------------------------------------

// -------------------------------------------------------------------------------------
// Inline Functions
// -------------------------------------------------------------------------------------

void clean_dir(string dualFile, string partitionFile, string metisFormatFile) {
    remove(dualFile);
    remove(partitionFile);
    remove(metisFormatFile);
}

int metisCheck() {
    string command = "gpmetis";
    auto output = executeShell(command);
    int status = output[0];
    return status;
}

string mesh2dual(string fileName, int ncommon) {
    writeln("-- Converting mesh to dual format");
    string outputFileName = "dual_"~fileName;
    string ncommon_arg = to!string(ncommon);
    string command = "m2gmetis " ~ " " ~ fileName ~ " " ~ outputFileName ~ " -gtype=dual -ncommon=" ~ ncommon_arg;
    executeShell(command);
    return outputFileName;
}

string partitionDual(string fileName, int nparts) {
    writeln("-- Partitioning dual format");
    string outputFile = fileName ~ ".part." ~ to!string(nparts);
    string nparts_arg = to!string(nparts);
    string command = "gpmetis " ~ " " ~ fileName ~ " " ~ nparts_arg;
    executeShell(command);
    return outputFile;
}

string SU2_to_metis_format(string fileName) {
    auto f = File(fileName, "r");
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
    writeln("-- Converting SU2 mesh to Metis input format");
    string outputFileName = "metisFormat_"~fileName;
    File outFile;
    outFile = File(outputFileName, "w");
    auto ncells = to!size_t(getHeaderContent("NELEM"));
    outFile.writeln(ncells);
    foreach(i; 0 .. ncells) {
	auto lineContent = f.readln().strip();
	auto tokens = lineContent.split();
	foreach(j; 1 .. tokens.length-1) { outFile.writef("%d \t", to!int(tokens[j])+1); }
	outFile.writef("\n");
    }
    return outputFileName;
}

void construct_blocks(string meshFile, string mappedCellsFilename, string partitionFile, string dualFile, int nparts) {
    writeln("-- Constructing ", nparts, " blocks");
    auto f = File(meshFile, "r");
    Cell[] global_cells;
    Node[] global_nodes;
    Face[] global_faces;
    Block[] blocks;
    Boundary[] global_boundaries;
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
    writeln("-- -- Reading su2 file");
    // Proceed through su2 file and collect cells
    auto dimensions = to!int(getHeaderContent("NDIME"));
    auto ncells = to!size_t(getHeaderContent("NELEM"));
    foreach(i; 0 .. ncells) {
	auto lineContent = f.readln().strip();
	auto tokens = lineContent.split();
	int cell_type = to!int(tokens[0]);
	size_t cell_id = to!size_t(tokens[$-1]);
	size_t[] node_id_list;
	foreach(j; 1 .. tokens.length-1) { node_id_list ~= to!size_t(tokens[j]); }
 	size_t[] face_id_list; // empty, so far
	global_cells ~= new Cell(cell_id, cell_type, node_id_list, face_id_list);
    }
    // Proceed through su2 file and collect nodes
    auto nnodes = to!size_t(getHeaderContent("NPOIN"));
    foreach(i; 0 .. nnodes) {
	auto tokens = f.readln().strip().split();
	double x=0.0; double y=0.0; double z = 0.0; size_t indx = 0;
	if (dimensions == 2) {
	    x = to!double(tokens[0]);
	    y = to!double(tokens[1]);
	    indx = to!size_t(tokens[2]);
	} else {
	    assert(dimensions == 3, "invalid dimensions");
	    x = to!double(tokens[0]);
	    y = to!double(tokens[1]);
	    z = to!double(tokens[2]);
	    indx = to!size_t(tokens[3]);
	}
	double[3] pos = [x, y, z];
	global_nodes ~= new Node(indx, pos);
    } // end foreach i .. nnodes
    // Collect what partition each cell belongs to
    writeln("-- -- Reading metis output file");
    f = File(partitionFile, "r");
    foreach(cell_id;0 .. ncells) {
	auto tokens = f.readln().strip().split();
	size_t partition_id = to!int(tokens[0]);
	global_cells[cell_id].partition = partition_id;
    }    
    // construct faces
    writeln("-- -- Construct mesh faces");
    //
    // Now that we have the full list of cells and vertices assigned to each cell,
    // we can construct the faces between cells and along the boundaries.
    //
    string makeFaceTag(const size_t[] node_id_list)
    {
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
    }
    size_t[string] faceIndices;
    void add_face_to_cell(ref Cell cell, size_t[] corners)
    {
	// If the face is new, we add it to the list, else we use the face
	// already stored within the list, so that it may no longer be
	// outward pointing for this cell.
	//
	size_t[] my_node_id_list;
	foreach(i; corners) { my_node_id_list ~= cell.node_id_list[i]; }
	size_t face_indx = 0;
	string faceTag = makeFaceTag(my_node_id_list);
	if (faceTag in faceIndices) {
	    // Use the face already present.
	    face_indx = faceIndices[faceTag];
	} else {
	    // Since we didn't find the face already, construct it.
	    face_indx = global_faces.length;
	    global_faces ~= new Face(my_node_id_list);
	    faceIndices[faceTag] = face_indx;
	}
	global_faces[face_indx].partition_id_list ~= cell.partition;
	cell.face_id_list ~= face_indx;
	return;
    } // end add_face_to_cell()
    foreach(cell; global_cells) {
	if (!cell) continue;
	// Attach the faces to each cell. In 2D, faces are defined as lines.
	// As we progress along the line the face normal is pointing to the right.
	// In 3D, a counter-clockwise cycles of points plus the right-hand rule
	// define the face normal. Whether a face points out of or into a cell
	// will be determined and remembered when we add the face to the cell.
	if (dimensions == 2) {
	    switch(cell.type) {
	    case 5: // triangle
		add_face_to_cell(cell, [0,1]);
		add_face_to_cell(cell, [1,2]);
		add_face_to_cell(cell, [2,0]);
		break;
	    case 9: // quad
		add_face_to_cell(cell, [2,3]); // north
		add_face_to_cell(cell, [1,2]); // east
		add_face_to_cell(cell, [0,1]); // south
		add_face_to_cell(cell, [3,0]); // west
		break;
	    default:
		throw new Exception("invalid cell type in 2D");
	    }
	} else {
	    assert(dimensions == 3, "invalid dimensions");
	    switch(cell.type) {
	    case 10: // tetra
		add_face_to_cell(cell, [0,1,2]);
		add_face_to_cell(cell, [0,1,3]);
		add_face_to_cell(cell, [1,2,3]);
		add_face_to_cell(cell, [2,0,3]);
		break;
	    case 12: // hexa
		add_face_to_cell(cell, [2,3,7,6]); // north
		add_face_to_cell(cell, [1,2,6,5]); // east
		add_face_to_cell(cell, [1,0,4,5]); // south
		add_face_to_cell(cell, [0,3,7,4]); // west
		add_face_to_cell(cell, [4,5,6,7]); // top
		add_face_to_cell(cell, [0,1,2,3]); // bottom
		break;
	    case 13: // wedge
		add_face_to_cell(cell, [0,1,2]);
		add_face_to_cell(cell, [3,4,5]);
		add_face_to_cell(cell, [0,2,5,3]);
		add_face_to_cell(cell, [0,3,4,1]);
		add_face_to_cell(cell, [1,4,5,2]);
		break;
	    case 14: // pyramid
		add_face_to_cell(cell, [0,1,2,3]);
		add_face_to_cell(cell, [0,1,4]);
		add_face_to_cell(cell, [1,2,4]);
		add_face_to_cell(cell, [2,3,4]);
		add_face_to_cell(cell, [3,0,4]);
		break;
	    default:
		throw new Exception("invalid cell type in 3D");
	    }
	}
    } // end foreach cell
    //
    // Now that we have a full set of cells and faces,
    // make lists of the boundary faces.
    //
    writeln("-- -- Fill boundary face arrays");
    f = File(meshFile, "r"); 
    auto nboundaries = to!size_t(getHeaderContent("NMARK"));
    foreach(i; 0 .. nboundaries) {
	string tag = getHeaderContent("MARKER_TAG");
	writeln("-- -- -- Found boundary i=", i, " tag=", tag);
	size_t[] face_id_list;
	size_t nelem = to!size_t(getHeaderContent("MARKER_ELEMS"));
	foreach(j; 0 .. nelem) {
	    auto tokens = f.readln().strip().split();
	    int vtk_type = to!int(tokens[0]);
	    size_t[] my_node_id_list;
	    foreach(k; 1 .. tokens.length) { my_node_id_list ~= to!size_t(tokens[k]); }
	    string faceTag = makeFaceTag(my_node_id_list);
	    if (faceTag in faceIndices) {
		size_t face_indx = faceIndices[faceTag];
		//Face my_face = faces[face_indx];
		face_id_list ~= face_indx;
	    } else {
		throw new Exception("cannot find face in collection");
	    }
	} // end foreach j
	global_boundaries ~= new Boundary(tag, face_id_list);
    } // end foreach i
    writeln("-- -- Reading cell interconnectivity from dual format file");
    // Collect cell inter-connectivity
    f = File(dualFile, "r");
    f.readln(); // skip first line
    foreach(cell_id;0 .. ncells) {
	 auto lineContent = f.readln().strip();
	 auto tokens = lineContent.split();
	 foreach(j; 0 .. tokens.length) { global_cells[cell_id].cell_neighbour_id_list ~= to!size_t(tokens[j]) - 1 ; }	 
    }
    // At this point we have all the data we need. Now construct the new blocks.
    // construct blocks
    writeln("-- -- Data reading complete, proceed to building blocks");
    foreach(i;0 .. nparts) {
	blocks ~= new Block(i);
	// create an INTERIOR boundary for all blocks (we will fill the face id list later)
	size_t[] face_id_list;
	blocks[i].boundary["METIS_INTERIOR"] = new Boundary("METIS_INTERIOR", face_id_list);
    }
    // assign cells to correct blocks, as well as interior boundary faces
    foreach(cell_id;0 .. ncells) {
	size_t partition_id = global_cells[cell_id].partition;
	blocks[partition_id].cells ~= global_cells[cell_id];
	blocks[partition_id].transform_cell_list[cell_id] = blocks[partition_id].cells.length - 1;
	foreach(node_id;global_cells[cell_id].node_id_list) {
	    if (!blocks[partition_id].node_id_list.canFind(node_id)) {
		blocks[partition_id].nodes ~= global_nodes[node_id];
		blocks[partition_id].node_id_list ~= node_id;
		blocks[partition_id].transform_list[node_id] = blocks[partition_id].nodes.length - 1; 
	    }
	}
	// assign interior boundaries
	foreach(neighbour_cell_id; global_cells[cell_id].cell_neighbour_id_list) {
	    if (global_cells[neighbour_cell_id].partition != global_cells[cell_id].partition) {
		// store shared face
		foreach(face_id; global_cells[cell_id].face_id_list) {
		    if (global_cells[neighbour_cell_id].face_id_list.canFind(face_id)) {
			blocks[partition_id].boundary["METIS_INTERIOR"].face_id_list ~= face_id;
			// facetag needs to be made with local node ids
			size_t[] local_node_id_list;
			foreach(node;global_faces[face_id].node_id_list)
			    local_node_id_list ~= blocks[global_cells[cell_id].partition].transform_list[node];
			string faceTag = makeFaceTag(local_node_id_list);
			// store mapped cell information
			blocks[partition_id].mapped_cells ~= new MappedCells(cell_id, neighbour_cell_id, global_cells[cell_id].partition, global_cells[neighbour_cell_id].partition, faceTag);
		    }
		}
	    }
	 }
    } 
    // finally assign domain boundary faces to correct blocks
    foreach(face_id;0 .. global_faces.length) {  
	if (global_faces[face_id].partition_id_list.length ==  1) {
	    foreach(i;0 .. nboundaries) {
		if (global_boundaries[i].face_id_list.canFind(face_id)) {
		    auto partition_id = global_faces[face_id].partition_id_list[0];
		    string tag = global_boundaries[i].tag;
		    if (!blocks[partition_id].bc_names.canFind(tag)) {
			size_t[] face_id_list;
			blocks[partition_id].boundary[tag] = new Boundary(tag, face_id_list);
			blocks[partition_id].bc_names ~= tag;
		    }
		    blocks[partition_id].boundary[tag].face_id_list ~= face_id;
		}                                                                 
	    }
	}
    }
    // at this point we have fully constructed the new blocks. Now write the data out to separate text files.
    // NB. cells and nodes use their local ids.
    File outFile_mappedcells;
    outFile_mappedcells = File(mappedCellsFilename, "w");

    // Let's check the user doesn't have any previous partitioned su2 files in the directory -- if there are, this will upset the code
    // if the user intends on partitioning more than one su2 file.
    int nSimilarBlocks;    
    string commandCheck = "if [ -f block_0_"~meshFile~" ] ; then echo 'yes' ; else echo 'no' ; fi";
    auto fileCheck = executeShell(commandCheck);
    if (fileCheck[1] != "no\n")  throw new Error("It appears the partitioned su2 files already exist. Please remove these files before running the partitioner again.");

    // If the user is partitioning more than one su2 file, then we need to count the current number of partitions (su2 files) to
    // know where to start the block counting for the new partitions.
    int nCurrentBlocks=0;

    string command = "ls -d block_*_* | wc -l";
    auto output = executeShell(command);
    if (output[1] != "ls: cannot access 'block_*_*': No such file or directory\n0\n") {
	nCurrentBlocks = to!int(splitLines(output[1])[0]);
    }

    // We are ready to construct those partition files!
    foreach(i;0..nparts) {
	writeln("-- Writing out block: #", i+nCurrentBlocks);
	// let's begin by writing out the mapped cells for the block
	
	string mapped_cells_tag = "NMappedCells in BLOCK[" ~ to!string(i+nCurrentBlocks) ~ "]= "; 
	outFile_mappedcells.writef("%s", mapped_cells_tag);
	auto ninterior = blocks[i].boundary["METIS_INTERIOR"].face_id_list.length;
	outFile_mappedcells.writef("%d \n", ninterior);
	
	foreach(mapped_cell; blocks[i].mapped_cells){
	    auto primary_cell_local_id = blocks[i].transform_cell_list[mapped_cell.global_cell_id];
	    auto secondary_cell_local_id = blocks[mapped_cell.neighbour_block_id].transform_cell_list[mapped_cell.global_neighbour_cell_id];
	    outFile_mappedcells.writef("%s \t", mapped_cell.faceTag);
	    outFile_mappedcells.writef("%d \t", mapped_cell.neighbour_block_id+nCurrentBlocks);
	    outFile_mappedcells.writef("%d \n", secondary_cell_local_id);
	}
	string outputFileName = "block_" ~ to!string(i+nCurrentBlocks) ~ "_" ~ meshFile;
	File outFile;
	outFile = File(outputFileName, "w");
	outFile.writeln("%");
	outFile.writeln("% Problem dimension");
	outFile.writeln("%");
	outFile.writef("NDIME= %d \n", dimensions);
	outFile.writeln("%");
	outFile.writeln("% Inner element connectivity");
	outFile.writeln("%");
	outFile.writef("NELEM= %d \n", blocks[i].cells.length);
	foreach(cell_id;0 .. blocks[i].cells.length) {
	    outFile.writef("%d \t", blocks[i].cells[cell_id].type);
	    foreach(j;0 .. blocks[i].cells[cell_id].node_id_list.length) {
		auto local_node_id = blocks[i].transform_list[blocks[i].cells[cell_id].node_id_list[j]];
		outFile.writef("%d \t", local_node_id);
	    }
	    outFile.writef("%d \n", cell_id); //blocks[i].cells[cell_id].id);
	}
	outFile.writeln("%");
	outFile.writeln("% Node coordinates");
	outFile.writeln("%");
	outFile.writef("NPOIN= %d \n", blocks[i].nodes.length);
	foreach(local_node_id;0 .. blocks[i].nodes.length) {
	    foreach(j;0 .. dimensions) { // blocks[i].nodes[node_id].pos.length)
		outFile.writef("%.17f \t", blocks[i].nodes[local_node_id].pos[j]);
	    }
	    outFile.writef("%d \n", local_node_id);
	}
	outFile.writeln("%");
	outFile.writeln("% Boundary elements");
	outFile.writeln("%");
	outFile.writef("NMARK= %d \n", blocks[i].boundary.length);
	foreach(bndary; blocks[i].boundary) {
	    outFile.writef("MARKER_TAG= %s \n", bndary.tag);
	    outFile.writef("MARKER_ELEMS= %s \n", bndary.face_id_list.length);
	    foreach(face_id; bndary.face_id_list) {
		switch(global_faces[face_id].node_id_list.length) {
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
		foreach(node_id; global_faces[face_id].node_id_list) {
		    auto local_node_id = blocks[i].transform_list[node_id];
		    outFile.writef("%d \t", local_node_id);
		}
		outFile.writef("\n");
	    }
	}
    }
}

int main(string[] args){
    // assign command line arguments to variable names
    string inputMeshFile; string mappedCellsFilename; int nparts; int ncommon;
    inputMeshFile = args[1];
    mappedCellsFilename = args[2];
    nparts = to!int(args[3]);
    ncommon = to!int(args[4]);
    writeln("Begin partitioner................");
    // check if Metis is installed correctly. 
    if (metisCheck() != 0)  throw new Error("Metis partition software cannot be found.");
    string metisFormatFile = SU2_to_metis_format(inputMeshFile);//("square_mesh.su2");
    string dualFormatFile = mesh2dual(metisFormatFile, ncommon);
    string partitionFile = partitionDual(dualFormatFile, nparts);
    construct_blocks(inputMeshFile, mappedCellsFilename, partitionFile, dualFormatFile, nparts);
    clean_dir(dualFormatFile, partitionFile, metisFormatFile);
    writeln("Finished partitioning................");
    return 0;
}
