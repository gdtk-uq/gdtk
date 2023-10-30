// Beginning of a structure model to run fluid-structure interaction simulations
module fsi.femmodel;

import std.stdio;
import std.parallelism;
import std.range;
import std.conv;
import std.algorithm;
import std.regex;
import std.format;

import util.lua_service;
import util.lua;
import geom;
import fsi.femconfig;
import nm.number;
import nm.bbla;
import globaldata;
import fluidblock;
import sfluidblock;
import fsi;

class FEMModel {
public:
    lua_State* myL;
    FEMConfig myConfig;
    int id;

    // The FEM node motion --> CFD vertex motion mapping helpers
    // What we have here is:
    // a) a 3D array, where the outermost dimensions denotes each block, second dimension for each CFD vertex
    //      in the block, last dimension for the FEM nodes to access when computing the vertex velocity.
    // b) a 3D array of same size as b), which contains the weights applied to the respective FEM node velocities.
    int[][][] nodeVelIds;
    double[][][] nodeVelWeights;

    // The CFD pressure --> FEM node force mapping
    // What we have here is (for the active plate surfaces):
    // a) A 1D array containing the block IDs that need to be referenced
    // b) a 2D array, with the first dimension corresponding to the block IDs in (a), and the second dimension
    //      the Node IDs the cell indices in (c) refer to
    // c) a 2D array, with the first dimension corresponding to the block IDs in (b), and the second dimensions
    //      the cell IDs to grab the pressure from
    int[] northSurfaceBlks, southSurfaceBlks;
    int[][] northFEMNodeIds, southFEMNodeIds;
    int[][] northCFDCellIds, southCFDCellIds;

    // At this stage, the FEM models are all classic ODEs of the form -KX + F = MA
    // Where K is stiffness matrix, X is the position vector, F is the applied force vector,
    // M is mass matrix and A is the acceleration vector.
    // Express as 2 first order ODEs: -KX + F = MVd        (d denotes time derivative)
    //                                      V = Xd
    // Assign pointers to the state vectors X, V and F, and allocate memory to them later when we
    // know how many elements in these vectors.
    number[] X, V;

    // For the moment, just use dense matrices for the M, K matrices. Likely move to sparse 
    // matrices in future once kinks in implementation are ironed out.
    Matrix!number K, M, F;

    // Rather than inverting the mass matrix to solve the first ODE, solve the linear system Ax=B
    // where A = M, x = Vd, B = (-KX + F). For this, we need the permutation matrix of M.
    size_t[2][] MpermuteList;

    // Finally, preallocate a vector to store the KX operation
    number[] KDotX;

    // Vector of node velocities in the reference frame of the plate, with the x direction being
    // normal to the plate.
    Vector3[] FEMNodeVel;

    // Orientation of the plate, to go from the plate reference frame to the global.
    Vector3 plateNormal, plateTangent1, plateTangent2;
    
    // Pressure at each node, taken from the fluid or UDF.
    number[] nodeNorthPressure, nodeSouthPressure;

    // The boundary conditions requires certain DoFs to be set to 0
    size_t[] ZeroedIndices;

    // Initialise the FEM model- follows the template used by the FluidBlocks, as we
    // may want to attach a lua interpreter for user-defined forcing.
    this(string jobName, int id) {
        myConfig = new FEMConfig(jobName);

        // Initialise the Lua interpreter which we may need for UDF pressure distributions
        if (myL) {
            writefln("We already have a nonzero pointer for the Lua interpreter for structure block %d", id);
        } else {
            myL = luaL_newstate();
        }

        if (!myL) { throw new Error("Could not allocate memory for Lua interpreter."); }
        luaL_openlibs(myL);
        lua_pushinteger(myL, id);
        lua_setglobal(myL, "blkId");
    }

    // begin plate_setup()
    void plateSetup() {
        // Set up the plate orientation
        plateNormal = Vector3(myConfig.plateNormal[0], myConfig.plateNormal[1], myConfig.plateNormal[2]);
        plateTangent1 = Vector3(myConfig.plateNormal[1], -myConfig.plateNormal[0], myConfig.plateNormal[2]);
        cross(plateTangent2, plateNormal, plateTangent1);


        // Set up the interface between the CFD and the FEM
        prepareGridMotionSetup();
        prepareNodeToVertexMap();

        int totalNodes = (myConfig.Nx + 1) * (myConfig.Nz + 1);
        nodeNorthPressure.length = totalNodes;
        nodeSouthPressure.length = totalNodes;
        FEMNodeVel.length = totalNodes;

        // The memory for the vectors and matrices is assigned in the model files, as they may have different numbers of
        // degrees of freedom therefore different numbers of elements.

        // Address boundary conditions first as they affect the formation of the matrices
        determineBoundaryConditions(myConfig.BCs);
        generateMassStiffnessMatrices();
        MpermuteList = decomp!number(M);

    } // end plate_setup()

    void finalize()
    {
        if (myL) {
            lua_close(myL);
            myL = null;
        }
    }

    // begin prepareGridMotionSetup
    void prepareGridMotionSetup() {
        // Set the size of the arrays
        nodeVelIds.length = myConfig.movingBlks.length;
        nodeVelWeights.length = myConfig.movingBlks.length;
        int nVerts;
        int[] _nodeIds;
        double[] _nodeWeights;

        // Read in the weights and indices files for each moving block
        foreach (i, blkId; myConfig.movingBlks) {
            auto byLineInd = File(format("FSI/Weights/%04d.indices", blkId), "r").byLine;
            auto byLineWeight = File(format("FSI/Weights/%04d.weights", blkId), "r").byLine();


            auto lineInd = byLineInd.front(); auto lineWeight = byLineWeight.front();
            while (!lineInd.empty) {
                nodeVelIds[i] ~= map!(to!int)(splitter(lineInd)).array;
                nodeVelWeights[i] ~= map!(to!double)(splitter(lineWeight)).array;
                byLineInd.popFront(); byLineWeight.popFront();
                lineInd = byLineInd.front(); lineWeight = byLineWeight.front();
            } 
        }
    } // end prepareGridMotionSetup

    // begin broadcastGridMotion()
    void broadcastGridMotion() {
        // Apply the grid motion to the relevant blocks
        foreach (i, blkId; parallel(myConfig.movingBlks)) {
            applyFSIMotionToBlock(to!int(i), blkId);
        }
    } // end broadcastGridMotion()

    // begin applyFSIMotionToBlock
    void applyFSIMotionToBlock(int movingBlkIndx, int blkId) {
        // Apply grid motion to a block
        FluidBlock blk = localFluidBlocks[blkId];
        Vector3 netVel;
        foreach (iv, vtx; blk.vertices) {
            netVel.clear();
            foreach (node; 0 .. 2) {
                netVel.add(FEMNodeVel[nodeVelIds[movingBlkIndx][iv][node]], nodeVelWeights[movingBlkIndx][iv][node]);
            }
            netVel.transform_to_global_frame(plateNormal, plateTangent1, plateTangent2);
            vtx.vel[0].set(netVel);
        }
    } // end applyFSIMotionToBlock

    // begin prepareNodeToVertexMap
    void prepareNodeToVertexMap() {
        // We're going to put the indices in a nice order, so we need some temporary storage space
        int[] _NodeIds, _CellIds;
        size_t[] sortIndx;
        if (myConfig.northForcing == ForcingType.Fluid) {
            auto mapFile = File("FSI/Weights/northC2N.mapping", "r").byLine();
            auto line = mapFile.front(); mapFile.popFront();
            while (!line.empty) {
                // First line is the block id
                northSurfaceBlks ~= parse!int(line);
                northFEMNodeIds.length++; northCFDCellIds.length++;
                // Next line is the NodeIds
                line = mapFile.front(); mapFile.popFront();
                _NodeIds = map!(to!int)(splitter(line)).array;
                // Then the cell Ids
                line = mapFile.front(); mapFile.popFront();
                _CellIds = map!(to!int)(splitter(line)).array;

                // Due to retrieving data from a Lua table does not have
                // a guaranteed order, we should re-order the indices to
                // make indexing more efficient
                sortIndx.length = _NodeIds.length;
                makeIndex!()(_NodeIds, sortIndx);
                foreach (i; sortIndx) {
                    northFEMNodeIds[$-1] ~= _NodeIds[i];
                    northCFDCellIds[$-1] ~= _CellIds[i];
                }

                // Grab the next line
                line = mapFile.front(); mapFile.popFront();
            }
        }

        if (myConfig.southForcing == ForcingType.Fluid) {
            auto mapFile = File("FSI/Weights/southC2N.mapping", "r").byLine();
            auto line = mapFile.front(); mapFile.popFront();
            while (!line.empty) {
                // First line is the block id
                southSurfaceBlks ~= parse!int(line);
                southFEMNodeIds.length++; southCFDCellIds.length++;
                // Next line is the NodeIds
                line = mapFile.front(); mapFile.popFront();
                _NodeIds = map!(to!int)(splitter(line)).array;
                // Then the cell Ids
                line = mapFile.front(); mapFile.popFront();
                _CellIds = map!(to!int)(splitter(line)).array;

                // Due to retrieving data from a Lua table does not have
                // a guaranteed order, we should re-order the indices to
                // make indexing more efficient
                sortIndx.length = _NodeIds.length;
                makeIndex!()(_NodeIds, sortIndx);
                foreach (i; sortIndx) {
                    southFEMNodeIds[$-1] ~= _NodeIds[i];
                    southCFDCellIds[$-1] ~= _CellIds[i];
                }

                // Grab the next line
                line = mapFile.front(); mapFile.popFront();
            }
        }
    }

    // begin retrievePressures
    void retrievePressures() {
        // Use the previously generated mapping to grab the relevant fluid pressures
        SFluidBlock blk;
        if (myConfig.northForcing == ForcingType.Fluid) {
            foreach (i, blkId; northSurfaceBlks) {
                blk = cast(SFluidBlock) globalBlocks[blkId];
                foreach (n; 0 .. northFEMNodeIds[i].length) {
                    nodeNorthPressure[northFEMNodeIds[i][n]] = blk.cells[northCFDCellIds[i][n]].fs.gas.p;
                }
            }
        }
        if (myConfig.southForcing == ForcingType.Fluid) {
            foreach (i, blkId; southSurfaceBlks) {
                blk = cast(SFluidBlock) globalBlocks[blkId];
                foreach (n; 0 .. southFEMNodeIds[i].length) {
                    nodeSouthPressure[southFEMNodeIds[i][n]] = blk.cells[southCFDCellIds[i][n]].fs.gas.p;
                }
            }
        }
    }

    // begin compute_vtx_velocities_for_FSI
    void computeVtxVelocitiesForFSI(double dt) {
        // Solve the set of ODEs
        // First, pull in the fluid pressures
        retrievePressures();

        // Then use these to fill in the force vector
        updateForceVector();

        // Now we can solve the ODE- reuse F to store the result of the linear system solution
        dt *= myConfig.couplingStep;
        dot(K, X, KDotX);
        F._data[] -= KDotX[];
        solve!number(M, F, MpermuteList);

        X[] += V[] * dt;
        V[] += F._data[] * dt;

        // Convert to the node displacement velocities used by the fluid mesh
        convertToNodeVel();
        broadcastGridMotion();

        // And then empty F
        F._data[] = 0.0;
    } // end compute_vtx_velocities_for_FSI

    // Methods that are model dependent
    abstract void generateMassStiffnessMatrices();
    abstract void updateForceVector();
    abstract void determineBoundaryConditions(string BCs);
    abstract void writeToFile(size_t tindx);
    abstract void readFromFile(size_t tindx);
    abstract void convertToNodeVel();
}
