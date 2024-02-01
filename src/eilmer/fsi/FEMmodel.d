// Master class for the solid models used in Fluid-Structure Interaction simulations
// Built by Lachlan Whyborn in late 2023/early 2024.

/+ To implement a new solid model, a few things are required, most of which are outlined
by the abstract functions defined at the end of this template.
        - model_setup(): tell the model how big the matrices/vectors should be, by defining
                1. nNodes: number of nodes in the solid model
                2. nDoF: total number of degrees of freedom (usually scalar * nNodes)
                3. nQuadPoints: how many quadrature points on the model (usually scalar * number of                     elements)
        - GenerateMassStiffnessMatrices(): build the mass and stiffness matrices for the solid
                model.
        - UpdateForceVector(): Update the force vector using the external forces (pressure or
                otherwise).
        - DetermineBoundaryConditions(string BCs): Build the vector of indices that denote which
                indices of the matrix to set to diagonal ones/elements of the vector to set to 0.
        - WriteToFile(size_t tindx): How to write the resulting position and velocity vectors
                to file.
        - ReadFromFile(size_t tindx): How to read a file containg position and velocity data to
                the relevant vectors, used for restarting simulations.
        - ConvertToNodeVel(): Conversion from the rates of change of the DoFs of the system to
                velocities applied to the mesh.

The master model currently defines methods for setting up the interfacing between the solid
and the fluid, allocating the memory required for the matrices and vectors, updating the ODEs
and broadcasting the results to the mesh.
+/

// I have tried to be consistent with formatting style, using camelCase for variables and
// UpperCamelCase for function names. The place I have not kept this consistent is with
// model_setup() and compute_vtx_velocities_for_FSI(), as they are called in the main code
// which uses primarily snake_case.

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
import gzip;
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

    // Metadata about the FEM model
    size_t nNodes, nQuadPoints, nDoF;

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
    size_t[2][] MPermuteList;

    // Preallocate memory to be used during the time stepping calculation
    number[] KDotX, XStage, VStage;
    number[][] dX, dV;
    Matrix!number rhs;

    // Vector of node velocities in the reference frame of the plate, with the x direction being
    // normal to the plate.
    Vector3[] FEMNodeVel;

    // Orientation of the plate, to go from the plate reference frame to the global.
    Vector3 plateNormal, plateTangent1, plateTangent2;
    
    // Pressure at each node, taken from the fluid or UDF, at the quadrature points.
    number[] northPressureAtQuads, southPressureAtQuads;

    // The boundary conditions requires certain DoFs to be set to 0
    size_t[] zeroedIndices;

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
    }

    // begin modelSetup()
    void model_setup() {
        // Set up the plate orientation
        plateNormal = Vector3(myConfig.plateNormal[0], myConfig.plateNormal[1], myConfig.plateNormal[2]);
        plateTangent1 = Vector3(myConfig.plateNormal[1], -myConfig.plateNormal[0], myConfig.plateNormal[2]);
        cross(plateTangent2, plateNormal, plateTangent1);

        // Set up the interface between the CFD and the FEM
        PrepareGridMotionSetup();
        PrepareNodeToVertexMap();

        // Allocate memory for the ODEs
        K = new Matrix!number(nDoF); M = new Matrix!number(nDoF); F = new Matrix!number(nDoF, 1);
        K.zeros(); M.zeros(); F.zeros();

        X.length = nDoF; V.length = nDoF; XStage.length = nDoF; VStage.length = nDoF; KDotX.length = nDoF;
        X[] = 0.0; V[] = 0.0; KDotX[] = 0.0;

        dX = new number[][](nDoF, 4); dV = new number[][](nDoF, 4);
        rhs = new Matrix!number(nDoF, 1);

        northPressureAtQuads.length = nQuadPoints; southPressureAtQuads.length = nQuadPoints;
        FEMNodeVel.length = nNodes;

        // The memory for the vectors and matrices is assigned in the model files, as they may have different numbers of
        // degrees of freedom therefore different numbers of elements.

        // Address boundary conditions first as they affect the formation of the matrices
        DetermineBoundaryConditions(myConfig.BCs);
        GenerateMassStiffnessMatrices();
        WriteMatricesToFile();
        MPermuteList = decomp!number(M);

    } // end plate_setup()

    void finalize()
    {
        if (myL) {
            lua_close(myL);
            myL = null;
        }
    }

    // begin prepareGridMotionSetup
    void PrepareGridMotionSetup() {
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
    void BroadcastGridMotion() {
        // Apply the grid motion to the relevant blocks
        foreach (i, blkId; parallel(myConfig.movingBlks)) {
            ApplyFSIMotionToBlock(to!int(i), blkId);
        }
    } // end broadcastGridMotion()

    // begin applyFSIMotionToBlock
    void ApplyFSIMotionToBlock(int movingBlkIndx, int blkId) {
        // Apply grid motion to a block
        SFluidBlock blk = cast(SFluidBlock) localFluidBlocks[blkId];
        Vector3 netVel;
        foreach (iv, vtx; blk.vertices) {
            netVel.clear();
            foreach (node; 0 .. 2) {
                if (myConfig.quasi3D) {
                    iv = iv % (blk.niv * blk.njv);
                }
                netVel.add(FEMNodeVel[nodeVelIds[movingBlkIndx][iv][node]], nodeVelWeights[movingBlkIndx][iv][node]);
            }
            netVel.transform_to_global_frame(plateNormal, plateTangent1, plateTangent2);
            vtx.vel[0].set(netVel);
        }
    } // end applyFSIMotionToBlock

    // begin prepareNodeToVertexMap
    void PrepareNodeToVertexMap() {
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
                // Then the cell Ids
                line = mapFile.front(); mapFile.popFront();
                _CellIds = map!(to!int)(splitter(line)).array;
                // Next line is the NodeIds
                line = mapFile.front(); mapFile.popFront();
                _NodeIds = map!(to!int)(splitter(line)).array;

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
                // Then the cell Ids
                line = mapFile.front(); mapFile.popFront();
                _CellIds = map!(to!int)(splitter(line)).array;
                // Next line is the NodeIds
                line = mapFile.front(); mapFile.popFront();
                _NodeIds = map!(to!int)(splitter(line)).array;

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
    void RetrievePressures() {
        // Use the previously generated mapping to grab the relevant fluid pressures
        SFluidBlock blk;
        if (myConfig.northForcing == ForcingType.Fluid) {
            foreach (i, blkId; northSurfaceBlks) {
                blk = cast(SFluidBlock) globalBlocks[blkId];
                foreach (n; 0 .. northFEMNodeIds[i].length) {
                    if (myConfig.quasi3D) {
                        northPressureAtQuads[northFEMNodeIds[i][n]] = 0.0;
                        size_t[3] twoD_indx = blk.to_ijk_indices_for_cell(northCFDCellIds[i][n]);
                        foreach (k; 0 .. blk.nkc) {
                            size_t threeD_indx = blk.cell_index(twoD_indx[0], twoD_indx[1], k);
                            northPressureAtQuads[northFEMNodeIds[i][n]] += blk.cells[threeD_indx].fs.gas.p;
                        }
                        northPressureAtQuads[northFEMNodeIds[i][n]] /= blk.nkc;
                    } else {
                        northPressureAtQuads[northFEMNodeIds[i][n]] = blk.cells[northCFDCellIds[i][n]].fs.gas.p;
                    }
                }
            }
        }
        if (myConfig.southForcing == ForcingType.Fluid) {
            foreach (i, blkId; southSurfaceBlks) {
                blk = cast(SFluidBlock) globalBlocks[blkId];
                foreach (n; 0 .. southFEMNodeIds[i].length) {
                    if (myConfig.quasi3D) {
                        southPressureAtQuads[southFEMNodeIds[i][n]] = 0.0;
                        size_t[3] twoD_indx = blk.to_ijk_indices_for_cell(southCFDCellIds[i][n]);
                        foreach (k; 0 .. blk.nkc) {
                            size_t threeD_indx = blk.cell_index(twoD_indx[0], twoD_indx[1], k);
                            southPressureAtQuads[southFEMNodeIds[i][n]] += blk.cells[threeD_indx].fs.gas.p;
                        }
                        southPressureAtQuads[southFEMNodeIds[i][n]] /= blk.nkc;
                    } else {
                        southPressureAtQuads[southFEMNodeIds[i][n]] = blk.cells[southCFDCellIds[i][n]].fs.gas.p;
                    }
                }
            }
        }
    }

    // begin compute_vtx_velocities_for_FSI (called in main code, so use snake_case).
    void compute_vtx_velocities_for_FSI(double dt) {
        // Solve the set of ODEs
        // First, pull in the fluid pressures
        RetrievePressures();

        // Then use these to fill in the force vector
        F._data[] = 0.0;
        UpdateForceVector();
        auto FFile = File("FSI/F.dat", "w+");
        FFile.writef("%1.8e", F._data[0]);
        foreach (i; 1 .. F._data.length) {
            FFile.writef(" %1.8e", F._data[i]);
        }

        // Now we can solve the ODE- reuse F to store the result of the linear system solution
        dt *= myConfig.couplingStep;

        // Attempt an RK4 step
        // First stage
        dot(K, X, KDotX);
        rhs._data[] = F._data[] - KDotX[];
        solve!number(M, rhs, MPermuteList);
        dV[][0] = rhs._data[];
        dX[][0] = V[];

        XStage[] = X[] + 0.5 * dX[0][] * dt;
        VStage[] = V[] + 0.5 * dV[0][] * dt;

        dot(K, XStage, KDotX);
        rhs._data[] = F._data[] - KDotX[];
        solve!number(M, rhs, MPermuteList);
        dV[][1] = rhs._data[];
        dX[][1] = VStage[];

        XStage[] = X[] + 0.5 * dX[1][] * dt;
        VStage[] = V[] + 0.5 * dV[1][] * dt;

        dot(K, XStage, KDotX);
        rhs._data[] = F._data[] - KDotX[];
        solve!number(M, rhs, MPermuteList);
        dV[][2] = rhs._data[];
        dX[][2] = VStage[];

        XStage[] = X[] + dX[2][] * dt;
        VStage[] = V[] + dV[2][] * dt;

        dot(K, XStage, KDotX);
        rhs._data[] = F._data[] - KDotX[];
        solve!number(M, rhs, MPermuteList);
        dV[][3] = rhs._data[];
        dX[][3] = VStage[];

        foreach (i; 0 .. X.length) {
            X[i] += (1./6.) * (dX[0][i] + 2. * dX[1][i] + 2. * dX[2][i] + dX[3][i]) * dt;
            V[i] += (1./6.) * (dV[0][i] + 2. * dV[1][i] + 2. * dV[2][i] + dV[3][i]) * dt;
        }

        // Convert to the node displacement velocities used by the fluid mesh
        ConvertToNodeVel();
        BroadcastGridMotion();

    } // end compute_vtx_velocities_for_FSI

    void WriteMatricesToFile() {
        // It will often be useful to write the mass and stiffness matrices to file,
        // which can be passed elsewhere to compute eigendecomposition, sparsity patterns etc.
        auto MFile = File("FSI/M.dat", "w+");
        auto KFile = File("FSI/K.dat", "w+");
        foreach (row; 0 .. M._nrows) {
            MFile.writef("%1.8e", M[row, 0]);
            KFile.writef("%1.8e", K[row, 0]);
            foreach (col; 1 .. M._ncols) {
                MFile.writef(" %1.8e", M[row, col]);
                KFile.writef(" %1.8e", K[row, col]);
            }
            MFile.write("\n");
            KFile.write("\n");
        }
        MFile.close(); KFile.close();
    }

    void WriteFSIToFile(size_t tindx) {
        // Loosely off the write_solution(sblk...) function in fluidblockio_old.d.
        // Open a zipfile, and pass it to the respective model writers.
        auto outfile = new GzipOut("FSI/t%04d.gzip");
        auto writer = appender!string();
        // Call the child writer method
        this.WriteToString(writer);
        // Compress the output
        outfile.compress(writer.data);
        outfile.finish();
    }

    void ReadFSIFromFile(size_t tindx) {
        // Read in the FSI file from a Gzipped file
        auto readFileByLine = new GzipByLine(format("FSI/t%04d.gzip"));
        // Call the child reader method
        this.ReadFromString(readFileByLine);
    }


    // Methods that are model dependent
    abstract void GenerateMassStiffnessMatrices();
    abstract void UpdateForceVector();
    abstract void DetermineBoundaryConditions(string BCs);
    abstract void WriteToString(Appender!string writer);
    abstract void ReadFromString(GzipByLine reader);
    abstract void ConvertToNodeVel();
}
