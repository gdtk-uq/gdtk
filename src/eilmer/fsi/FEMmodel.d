// Master class for the solid models used in Fluid-Structure Interaction simulations
// Built by Lachlan Whyborn in late 2023/early 2024.

/+ To implement a new solid model, a few things are required, most of which are outlined
by the abstract functions defined at the end of this template.
        - model_setup(): tell the model how big the matrices/vectors should be, by defining
                1. nNodes: double of nodes in the solid model
                2. nDoF: total double of degrees of freedom (usually scalar * nNodes)
                3. nQuadPoints: how many quadrature points on the model (usually scalar * double of                     elements)
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
import std.file;

import util.lua_service;
import util.lua;
import gzip;
import geom;
import fsi.femconfig;
import fsi.damping_models;
import nm.number;
import nm.bbla;
import globaldata;
import fluidblock;
import sfluidblock;
import fileutil;
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

    // At this stage, the FEM models are all classic ODEs of the form -KX - CV + F = MA
    // Where K is stiffness matrix, X is the position vector, C is the damping matrix,
    // V is the velocity vector, F is the applied force vector,
    // M is mass matrix and A is the acceleration vector.
    // Express as 2 first order ODEs: -KX - CV + F = MVdot      (dot denotes time derivative)
    //                                           V = Xdot
    // Assign pointers to the state vectors X, V and F, and allocate memory to them later when we
    // know how many elements in these vectors.
    double[] X, V, F;

    // For the moment, just use dense matrices for the M, K and C matrices. Likely move to sparse 
    // matrices in future once kinks in implementation are ironed out.
    Matrix!double K, M, C;

    // Set a boolean to track whether damping is included or not, better retrieving data from
    // an enum every step.
    bool includeDamping = false;

    // Rather than inverting the mass matrix to solve the first ODE, solve the linear system Ax=B
    // where A = M, x = Vdot, B = (-KX - CV + F). For this, we need the permutation matrix of M.
    size_t[2][] MPermuteList;

    // Preallocate memory to be used during the time stepping calculation
    double[] KDotX, CDotV, XStage, VStage;
    double[][] dX, dV;
    Matrix!double rhs;

    // Vector of node velocities in the reference frame of the plate, with the x direction being
    // normal to the plate.
    Vector3[] FEMNodeVel;

    // Orientation of the plate, to go from the plate reference frame to the global.
    Vector3 plateNormal, plateTangent1, plateTangent2;
    
    // Pressure at each node, taken from the fluid or UDF, at the quadrature points.
    double[] northPressureAtQuads, southPressureAtQuads;

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
        // The individual models have their own model_setup(), which specifies the number of nodes and DoFs
        // based on the specified number of elements. They then call this parent method, which uses those
        // computed number of nodes and DoFs to allocate memory and do other setup operations which are
        // common across all models.

        // Set up the plate orientation
        plateNormal = Vector3(myConfig.plateNormal[0], myConfig.plateNormal[1], myConfig.plateNormal[2]);
        plateTangent1 = Vector3(myConfig.plateNormal[1], -myConfig.plateNormal[0], myConfig.plateNormal[2]);
        cross(plateTangent2, plateNormal, plateTangent1);

        // Set up the interface between the CFD and the FEM
        PrepareGridMotionSetup();
        PrepareNodeToVertexMap();

        // Allocate memory for the ODEs
        // Start with the matrices. Don't include the damping matrix yet, because in most instances we want
        // to build the damping matrix based on the stiffness and mass matrices.
        K = new Matrix!double(nDoF); M = new Matrix!double(nDoF);
        K.zeros(); M.zeros();

        if (myConfig.dampingModel != DampingModel.None) {   // Only allocate for the damping if we use it
            C = new Matrix!double(nDoF); C.zeros();
            includeDamping = true;
        }

        // Allocate the vector memory for the states
        X.length = nDoF; V.length = nDoF; F.length = nDoF; XStage.length = nDoF; VStage.length = nDoF; KDotX.length = nDoF;
        X[] = 0.0; V[] = 0.0; F[] = 0.0;
        
        XStage.length = nDoF; VStage.length = nDoF; KDotX.length = nDoF; CDotV.length = nDoF;   // Don't need to zero these

        // Allocate memory for the derivate stages
        switch (myConfig.temporalScheme) {
            case FEMTemporalScheme.RK4:
                // 4 Stages
                dX = new double[][](nDoF, 4); dV = new double[][](nDoF, 4);
                break;
            case FEMTemporalScheme.euler:
                // 1 Stage
                dX = new double[][](nDoF, 1); dV = new double[][](nDoF, 1);
                break;
            default:
                throw new Error("Something went wrong in choosing the temporal scheme.");
        }

         // Allocate memory to store the results of the linear solve of the ODe
        rhs = new Matrix!double(nDoF, 1);

        // Finally, allocate memory for the interfacing vectors i.e. external pressures from fluid to solid 
        // and node velocities from solid to fluid.
        northPressureAtQuads.length = nQuadPoints; southPressureAtQuads.length = nQuadPoints;
        FEMNodeVel.length = nNodes;

        // Address boundary conditions first as they affect the formation of the matrices
        DetermineBoundaryConditions(myConfig.BCs);

        // Compute the mass and stiffness matrices by calling on the child method.
        GenerateMassStiffnessMatrices();

        // Now that we have the mass and stiffness matrices, we can call the damping matrix builder
        switch (myConfig.dampingModel) {
            case DampingModel.Rayleigh2Parameter:
                // Use the Rayleigh 2 parameter model
                RayleighDampingMatrix(C, K, M, myConfig.dampingRatios, myConfig.naturalFrequencies);
                break;
            case DampingModel.None:
                goto default;
            default:
                break;
        }

        // Pre-compute the permutation list used to solve the linear system.
        MPermuteList = decomp!double(M);

        // If the user desires, write the mass and stiffness matrices to file for later analysis.
        if (myConfig.writeMatrices) { WriteMatricesToFile(); }

        // Prepare any history files
        PrepareFSIHistoryFiles();
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
                            northPressureAtQuads[northFEMNodeIds[i][n]] += blk.cells[threeD_indx].fs.gas.p.re;
                        }
                        northPressureAtQuads[northFEMNodeIds[i][n]] /= blk.nkc;
                    } else {
                        northPressureAtQuads[northFEMNodeIds[i][n]] = blk.cells[northCFDCellIds[i][n]].fs.gas.p.re;
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
                            southPressureAtQuads[southFEMNodeIds[i][n]] += blk.cells[threeD_indx].fs.gas.p.re;
                        }
                        southPressureAtQuads[southFEMNodeIds[i][n]] /= blk.nkc;
                    } else {
                        southPressureAtQuads[southFEMNodeIds[i][n]] = blk.cells[southCFDCellIds[i][n]].fs.gas.p.re;
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
        F[] = 0.0;
        UpdateForceVector();

        // Now we can solve the ODE- reuse F to store the result of the linear system solution
        dt *= myConfig.couplingStep;

        switch (myConfig.temporalScheme) {
            case FEMTemporalScheme.RK4:
                RK4FEMUpdate(dt); break;
            case FEMTemporalScheme.euler:
                eulerFEMUpdate(dt); break;
            default:
                throw new Error("Something bad happened when choosing the FEM temporal update.");
        }

        // Convert to the node displacement velocities used by the fluid mesh
        ConvertToNodeVel();
        BroadcastGridMotion();

    } // end compute_vtx_velocities_for_FSI

    void RK4FEMUpdate(double dt) {
        // Perform a temporal update using 4th order RK.

        // First stage
        dot(K, X, KDotX);
        rhs._data[] = F[] - KDotX[];
        if (includeDamping) {
            dot(C, V, CDotV);
            rhs._data[] -= CDotV[];
        }
        solve!double(M, rhs, MPermuteList);
        dV[][0] = rhs._data[];
        dX[][0] = V[];

        XStage[] = X[] + 0.5 * dX[0][] * dt;
        VStage[] = V[] + 0.5 * dV[0][] * dt;

        // Second stage
        dot(K, XStage, KDotX);
        rhs._data[] = F[] - KDotX[];
        if (includeDamping) {
            dot(C, VStage, CDotV);
            rhs._data[] -= CDotV[];
        }
        solve!double(M, rhs, MPermuteList);
        dV[][1] = rhs._data[];
        dX[][1] = VStage[];

        XStage[] = X[] + 0.5 * dX[1][] * dt;
        VStage[] = V[] + 0.5 * dV[1][] * dt;

        // Third stage
        dot(K, XStage, KDotX);
        rhs._data[] = F[] - KDotX[];
        if (includeDamping) {
            dot(C, VStage, CDotV);
            rhs._data[] -= CDotV[];
        }
        solve!double(M, rhs, MPermuteList);
        dV[][2] = rhs._data[];
        dX[][2] = VStage[];

        XStage[] = X[] + dX[2][] * dt;
        VStage[] = V[] + dV[2][] * dt;

        // Fourth stage
        dot(K, XStage, KDotX);
        rhs._data[] = F[] - KDotX[];
        if (includeDamping) {
            dot(C, VStage, CDotV);
            rhs._data[] -= CDotV[];
        }
        solve!double(M, rhs, MPermuteList);
        dV[][3] = rhs._data[];
        dX[][3] = VStage[];

        X[] += (1./6.) * (dX[0][] + 2. * dX[1][] + 2. * dX[2][] + dX[3][]) * dt;
        V[] += (1./6.) * (dV[0][] + 2. * dV[1][] + 2. * dV[2][] + dV[3][]) * dt;
        //foreach (i; 0 .. X.length) {
            //X[i] += (1./6.) * (dX[0][i] + 2. * dX[1][i] + 2. * dX[2][i] + dX[3][i]) * dt;
            //V[i] += (1./6.) * (dV[0][i] + 2. * dV[1][i] + 2. * dV[2][i] + dV[3][i]) * dt;
        //}
    } // end RK4FEMUpdate

    void eulerFEMUpdate(double dt) {
        // Perform a simple euler temporal update
        dot(K, X, KDotX);
        rhs._data[] = F[] - KDotX[];
        if (includeDamping) {
            dot(C, V, CDotV);
            rhs._data[] -= CDotV[];
        }
        solve!double(M, rhs, MPermuteList);
        dV[][0] = rhs._data[];
        dX[][0] = V[];

        X[] += dX[0][] * dt;
        V[] += dV[0][] * dt;
    } // end eulerFEMUpdate

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
    } // end WriteMatricesToFile

    void PrepareFSIHistoryFiles() {
        // Prepare the history files

        if (myConfig.historyNodes.length > 0) {
            ensure_directory_is_present("FSI/hist");
        }

        foreach (node; myConfig.historyNodes) {
            auto histFile = File(format("FSI/hist/%04d.dat", node), "w");
            append(format("FSI/hist/%04d.dat", node), GetHistoryHeader());
        }
    } // end PrepareFSIHistoryFiles

    // Methods that are model dependent
    abstract void GenerateMassStiffnessMatrices();
    abstract void UpdateForceVector();
    abstract void DetermineBoundaryConditions(string BCs);
    abstract void WriteFSIToFile(size_t tindx);
    abstract void ReadFSIFromFile(size_t tindx);
    abstract void ConvertToNodeVel();
    abstract string GetHistoryHeader();
    abstract void WriteFSIToHistory(double t);
}
