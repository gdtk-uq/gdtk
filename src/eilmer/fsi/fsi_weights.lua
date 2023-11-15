-- Script to generate both the weighting functions that convert the FEM node
-- velocities to the CFD vertex velocities, and the mapping from cell centre
-- pressures to FEM node pressures.

-- Start with the weights. This drawing should help explain the configuration.

--[[

        -------------------------------------------------
        |               |               |               |
        |               |               |               |
Flow    |    Fluid      |    Fluid      |     Fluid     |
-->     |      0        |      1        |       2       |
        |               |               |               |
        |               |               |               |
        -------------------------------------------------
        |             plate             |     Fluid 3   |
        -------------------------------------------------
        |               |               |               |
        |               |               |               |
Flow    |    Fluid      |     Fluid     |     Fluid     |
-->     |      4        |       5       |       6       |
        |               |               |               |
        |               |               |               |
        -------------------------------------------------

-The plate can move up/down at each node. We want the top boundaries of 0, 1, 2 and
bottom of 4, 5, 6 to be unmoved. So we need distance weighting such that we get the
full velocity at the plate, and 0 at the outer boundaries. This is computed in
generateDistanceWeightings function.

Then we also need to determine which FEM nodes contribute to the vertex velocities
in a column. Using a bi-linear interpolation, so there will be contributions from
4 FEM node velocities to each vertex velocity (some contributions will be 0 for 
vertices on boundaries). This weighting can be described by 4 integers, denoting 
the FEM node indices, and 4 floats, denoting the weight of each respective node's
contribution.

The weight from each respective node can then be scaled by the previously computed 
distance weighting so get the net contribution of the given FEM node velocity to the 
CFD vertex velocity. i.e. every vertex has associated vectors

Indx = [i0, i1, i2, i3]; Weights = [w0, w1, w2, w3]

Then that vertex's velocity is

VtxVel = FEMVel[i0] * w0 + FEMVel[i1] * w1 ...

This is for the fluid blocks directly above and below the plate. At the moment,
the blocks downstream (and transverse to in 3D) i.e. 2, 6 are going to move in
lockstep with the eastmost vertices of 1, 5. All the vertices in 3 will move
according to the eastmost row of FEM velocities.

It may be useful in future to rewrite this such that the eastmost boundaries of
2 and 6 remain stationary, to make it more flexible in more complicated simulations.
--]]

local fsi_weights = {}

-- The functions that actually do the work
function fsi_weights.generateDistanceWeights(FBA, TopToBottom)
    -- Run through each column of vertices and compute the relative distance between
    -- the plate and the outer boundary
    distanceWeights = {}
    niv = FBA.niv
    njv = FBA.njv
    nkv = FBA.nkv

    -- Initialise the tables in the natural ijk order so make accessing more logical
    for k = 0, nkv-1 do
        distanceWeights[k] = {}
        for j = 0, njv-1 do
            distanceWeights[k][j] = {}
        end
    end

    -- Now run through each column
    for k = 0, nkv-1 do
        for i = 0, niv-1 do
            distances = {}
            distances[0] = 0.0
            p0 = FBA.gridArray.grid:get_vtx(i, 0, k)
            -- Doesn't actually matter which end we start from, as long as run through
            -- the correct way when we traverse back through
            for j = 1, njv-1 do
                p1 = FBA.gridArray.grid:get_vtx(i, j, k)
                distances[j] = distances[j-1] + vabs(p1 - p0)
                p0 = p1
            end
            total_distance = distances[#distances]

            -- Run back through and divide the vertex distance by the total
            for j = 0, njv-1 do
                if TopToBottom then
                    distanceWeights[k][j][i] = distances[j] / total_distance
                else
                    distanceWeights[k][j][i] = 1 - distances[j] / total_distance
                end
            end
        end
    end

    return distanceWeights
end

function fsi_weights.generateInterpolationWeights(FBA, Nx, Nz, Length, Width, Side)
    -- The linear interpolation weights for the FEM nodes --> CFD vertices
    interpWeights = {}
    interpIndices = {}

    -- The FEM nodes are always uniformly spaced, but the CFD vertices may not be so
    -- This means we can always trivially determine which FEM nodes are relevant if we
    -- can work where the CFD vertex is
    -- At the very least, we will assume the grid lines in the CFD are straight
    dxFEM = 1 / Nx;
    if Nz > 0 then
        dzFEM = 1 / Nz;
    else
        dzFEM = 1
    end

    if Side == "north" then
        j = 0
    elseif Side == "south" then
        j = FBA.njv-1
    end

    for k = 0, nkv-1 do
        interpWeights[k] = {}
        interpIndices[k] = {}
        for i = 0, niv-1 do
            interpWeights[k][i] = {}
            interpIndices[k][i] = {}
        end
    end

    niv = FBA.niv; nkv = FBA.nkv
    ref_point = FBA.gridArray.grid:get_vtx(0, j, 0)
    -- Run through i then k to simplify 2D/3D handling
    for i = 0, niv-1 do
        -- Normalised distance along plate in i direction
        xNormalised = vabs(FBA.gridArray.grid:get_vtx(i, j, 0) - ref_point) / Length

        -- Work out which FEM nodes are either side of this CFD vertex
        iIndxFEM = 1
        for _iIndxFEM = 1, Nx do
            if dxFEM * _iIndxFEM >= xNormalised then
                iIndxFEM = _iIndxFEM
                break
            end
        end

        for k = 0, nkv-1 do
            -- Normalised distance along plate in the k direction
            zNormalised = vabs(FBA.gridArray.grid:get_vtx(0, j, k) - ref_point) / Width

            kIndxFEM = 1
            -- Work out which FEM nodes
            for _kIndxFEM = 1, Nz do
                if dzFEM * _kIndxFEM >= zNormalised then
                    kIndxFEM = _kIndxFEM
                    break
                end
            end

            -- Compute the weights with bilinear interpolation
            -- [(i-1, k-1), (i-1, k+1), (i+1, k-1), (i+1, k+1)]

            -- iIndxFEM * dxFEM, kIndxFEM * dzFEM denote fractionally how far we along the plate in each direction
            -- xNormalized/zNormalized denote where the CFD vertex is along the plate
            if (config.dimensions == 3 and not FSIOptions.quasi3D) then
                interpWeights[k][i][1] = (1 - (xNormalised - (iIndxFEM - 1) * dxFEM) / dxFEM) * (1 - (zNormalised - (kIndxFEM - 1) * dzFEM) / dzFEM)
                interpIndices[k][i][1] = (kIndxFEM - 1) * (Nx + 1) + (iIndxFEM - 1)
                interpWeights[k][i][2] = (xNormalised - (iIndxFEM - 1) * dxFEM) / dxFEM * (1 - (zNormalised - (kIndxFEM - 1) * dzFEM) / dzFEM)
                interpIndices[k][i][2] = (kIndxFEM - 1) * (Nx + 1) + iIndxFEM
                interpWeights[k][i][3] = (1 - (xNormalised - (iIndxFEM - 1) * dxFEM) / dxFEM) * (zNormalised - (kIndxFEM - 1) * dzFEM) / dzFEM
                interpIndices[k][i][3] = kIndxFEM * (Nx + 1) + (iIndxFEM - 1)
                interpWeights[k][i][4] = (xNormalised - (iIndxFEM - 1) * dxFEM) / dxFEM * (zNormalised - (kIndxFEM - 1) * dzFEM) / dzFEM
                interpIndices[k][i][4] = kIndxFEM * (Nx + 1) + iIndxFEM
            else
                interpWeights[k][i][1] = 1 - (xNormalised - (iIndxFEM - 1) * dxFEM) / dxFEM
                interpIndices[k][i][1] = iIndxFEM - 1
                interpWeights[k][i][2] = (xNormalised - (iIndxFEM - 1) * dxFEM) / dxFEM
                interpIndices[k][i][2] = iIndxFEM
            end
        end
    end


    return interpWeights, interpIndices
end

function fsi_weights.writeWeightsToFile(FBA, distanceWeights, interpWeights, interpIndices, location)
    -- each orientation needs to take different elements from the weights tables
    -- keep track of the global i,j,k index
    ni_total = 0
    for ib = 1, FBA.nib do
        nj_total = 0
        for jb = 1, FBA.njb do
            nk_total = 0
            for kb = 1, FBA.nkb do
                if config.dimensions == 3 then blkId = FBA.blockArray[ib][jb][kb].id else blkId = FBA.blockArray[ib][jb].id end
                weightsFile = assert(io.open(string.format("FSI/Weights/%04d.weights", blkId), "w"))
                indicesFile = assert(io.open(string.format("FSI/Weights/%04d.indices", blkId), "w"))
                niv = FBA.gridArray.nics[ib]+1; njv = FBA.gridArray.njcs[jb]+1;
                if (config.dimensions == 3 and not FSIOptions.quasi3D) then nkv = FBA.gridArray.nkcs[kb]+1 else nkv = 1 end

                for k = 0, nkv-1 do
                    for j = 0, njv-1 do
                        for i = 0, niv-1 do
                            -- Indices in the FBA
                            ig = i + ni_total; jg = j + nj_total; kg = k + nk_total
                            if (config.dimensions == 3 and not FSIOptions.quasi3D) then nVerts = 4 else nVerts = 2 end
                            for n = 1, nVerts do
                                if location == "plate" then
                                    weightsFile:write(string.format("%.18e ", interpWeights[kg][ig][n] * distanceWeights[kg][jg][ig]))
                                    indicesFile:write(string.format("%d ", interpIndices[kg][ig][n]))
                                elseif location == "west" then
                                    weightsFile:write(string.format("%1.8e ", interpWeights[kg][0][n] * distanceWeights[kg][jg][0]))
                                    indicesFile:write(string.format("%d ", interpIndices[kg][ig][n]))
                                elseif location == "east" then
                                    weightsFile:write(string.format("%1.8e ", interpWeights[kg][#interpWeights[kg]][n] * distanceWeights[kg][jg][#distanceWeights[kg][jg]]))
                                    indicesFile:write(string.format("%d ", interpIndices[kg][#interpIndices[kg]][n]))
                                elseif location == "bottom" then
                                    weightsFile:write(string.format("%1.8e ", interpWeights[0][ig][n] * distanceWeights[0][jg][ig]))
                                    indicesFile:write(string.format("%d ", interpIndices[0][ig][n]))
                                elseif location == "top" then
                                    weightsFile:write(string.format("%1.8e ", interpWeights[#interpWeights][ig][n] * distanceWeights[#distanceWeights][jg][ig]))
                                    indicesFile:write(string.format("%d ", interpIndices[#interpIndices][ig][n]))
                                elseif location == "west-bottom" then
                                    weightsFile:write(string.format("%1.8e ", interpWeights[0][0][n] * distanceWeights[0][jg][0]))
                                    indicesFile:write(string.format("%d ", interpIndices[0][0][n]))
                                elseif location == "west-top" then
                                    weightsFile:write(string.format("%1.8e ", interpWeights[#interpWeights][0][n] * distanceWeights[#distanceWeights][jg][0]))
                                    indicesFile:write(string.format("%d ", interpIndices[#interpIndices][0][n]))
                                elseif location == "east-bottom" then
                                    weightsFile:write(string.format("%1.8e ", interpWeights[0][#interpWeights[0]][n] * distanceWeights[0][jg][#distanceWeights[0][jg]]))
                                    indicesFile:write(string.format("%d ", interpIndices[0][#interpIndices[0]][n]))
                                elseif location == "west-top" then
                                    weightsFile:write(string.format("%1.8e ", interpWeights[#interpWeights][#interpWeights[#interWeights]][n] * distanceWeights[#distanceWeights][jg][#distanceWeights[#distanceWeights][jg]]))
                                    indicesFile:write(string.format("%d ", interpIndices[#interpIndices][#interpIndices[#interpIndices]][n]))
                                elseif location == "west-adjacent" then
                                    weightsFile:write(string.format("%1.8e ", interpWeights[kg][0][n]))
                                    indicesFile:write(string.format("%d ", interpIndices[kg][0][n]))
                                elseif location == "east-adjacent" then
                                    weightsFile:write(string.format("%1.8e ", interpWeights[kg][#interpWeights[kg]][n]))
                                    indicesFile:write(string.format("%d ", interpIndices[kg][#interpIndices[kg]][n]))
                                elseif location == "bottom-adjacent" then
                                    weightsFile:write(string.format("%1.8e ", interpWeights[0][ig][n]))
                                    indicesFile:write(string.format("%d ", interpIndices[0][ig][n]))
                                elseif location == "top-adjacent" then
                                    weightsFile:write(string.format("%1.8e ", interpWeights[#interpWeights][ig][n]))
                                    indicesFile:write(string.format("%d ", interpIndices[#interpWeights][ig][n]))
                                elseif location == "west-bottom-adjacent" then
                                    weightsFile:write(string.format("%1.8e ", interpWeights[0][0][n]))
                                    indicesFile:write(string.format("%d ", interpIndices[0][0][n]))
                                elseif location == "west-top-adjacent" then
                                    weightsFile:write(string.format("%1.8e ", interpWeights[#interpWeights][0][n]))
                                    indicesFile:write(string.format("%d ", interpIndices[#interpWeights][0][n]))
                                elseif location == "east-bottom-adjacent" then
                                    weightsFile:write(string.format("%1.8e ", interpWeights[0][#interpWeights[0]][n]))
                                    indicesFile:write(string.format("%d ", interpIndices[0][#interpWeights[0]][n]))
                                elseif location == "east-top-adjacent" then
                                    weightsFile:write(string.format("%1.8e ", interpWeights[#interpWeights][#interpWeights[#interpWeights]][n]))
                                    indicesFile:write(string.format("%d ", interpIndices[#interpWeights][#interpWeights[#interpWeights]][n]))
                                end
                            end
                            weightsFile:write("\n")
                            indicesFile:write("\n")
                        end
                    end
                end

                weightsFile:close()
                indicesFile:close()

                -- Add to the global index, -1 because of shared vertices
                nk_total = nk_total + nkv - 1
            end
            nj_total = nj_total + njv - 1
        end
        ni_total = ni_total + niv - 1
    end
end

function ijkFBAIndxToFBIndx(FBA, i, j, k)
    -- Convert an (i, j, k) cell index in an FBArray to a single blkID and single index within that block
    iAccum = 0; jAccum = 0; kAccum = 0
    
    iBlk = 1; jBlk = 1; kBlk = 1
    -- Run through nics, njcs, nkcs
    for _iBlk, nic in pairs(FBA.nics) do
        if (iAccum + nic) > i then
            iInBlk = i - iAccum
            iBlk = _iBlk
            break
        end
        iAccum = iAccum + nic
    end

    for _jBlk, njc in pairs(FBA.njcs) do
        if (jAccum + njc) > j then
            jInBlk = j - jAccum
            jBlk = _jBlk
            break
        end
        jAccum = jAccum + njc
    end

    for _kBlk, nkc in pairs(FBA.nkcs) do
        if (kAccum + nkc) > k then
            kInBlk = k - kAccum
            kBlk = _kBlk
            break
        end
        kAccum = kAccum + nkc
    end

    -- We know which iBlk, kBlk and i, k in that block
    -- Convert to blkId and cellId

    if config.dimensions == 3 then blkId = FBA.blockArray[iBlk][jBlk][kBlk].id else blkId = FBA.blockArray[iBlk][jBlk].id end
    cellId = kInBlk * (FBA.nics[iBlk] * FBA.njcs[jBlk]) + jInBlk * FBA.nics[iBlk] + iInBlk

    return blkId, cellId
end

function fsi_weights.cellToQuadratureMapping(FBA, Nx, Nz, Length, Width, Side)
    -- Every element has 2 quadrature points, located at +-1/sqrt(3), in the space where an element width
    -- is (-1, 1)
    cellToQuadratureMap = {}
    if Side == "north" then
        j = 0
        j_cell = 0
    elseif Side == "south" then
        j = FBA.njv-1
        j_cell = FBA.njv-2
    end

    -- Work out where the quadrature points are along the surface
    quadLoc = (1 - 1 / math.sqrt(3)) / 2
    dxFEM = Length / Nx
    quadPoints = {}
    for i = 1, Nx do
        quadPoints[#quadPoints+1] = ((i-1) + quadLoc) * dxFEM
        quadPoints[#quadPoints+1] = ((i-1) + (1-quadLoc)) * dxFEM
    end

    ref_point = FBA.gridArray.grid:get_vtx(0, j, 0)
    cellPoints = {}
    for i = 1, FBA.niv-1 do
        xPos = vabs(FBA.gridArray.grid:get_vtx(i, j, 0) - ref_point)
        cellPoints[i] = xPos
    end

    iCell = 1
    for iQuad = 1, 2*Nx do
        while cellPoints[iCell] < quadPoints[iQuad] do
            iCell = iCell + 1
        end
        blkId, cellId = ijkFBAIndxToFBIndx(FBA, iCell-1, j, 0)
        if cellToQuadratureMap[blkId] then
            cellToQuadratureMap[blkId][iQuad] = cellId
        else
            cellToQuadratureMap[blkId] = {}
            cellToQuadratureMap[blkId][iQuad] = cellId
        end
    end

    mappingFile = assert(io.open("FSI/Weights/"..Side.."C2N.mapping", "w"))
    for blkId, blkTable in pairs(cellToQuadratureMap) do
        mappingFile:write(string.format("%d\n", blkId))
        for _, cellId in pairs(blkTable) do
            mappingFile:write(string.format("%d ", cellId))
        end
        mappingFile:write("\n")
        for quadId, _ in pairs(blkTable) do
            mappingFile:write(string.format("%s ", quadId-1))
        end
        mappingFile:write("\n")
    end
end

function fsi_weights.cellToNodeMapping(FBA, Nx, Nz, Length, Width, Side)
    -- Each node in the FEM model has a cell associated with it; determine which cell this is
    -- The mapping table contains key : value pairs of blkId : {CellIndices}
    CellToNodeMap = {}

    if Side == "north" then
        j = 0
        j_cell = 0
    elseif Side == "south" then
        j = FBA.njv-1;
        j_cell = FBA.njv-2;
    end

    -- The machinery here looks quite similar to that in GenerateInterpolationWeights
    dxFEM = 1 / Nx;
    if Nz > 0 then
        dzFEM = 1 / Nz;
    else
        dzFEM = 1
    end
    niv = FBA.niv; nkv = FBA.nkv
    --ref_point = FBA.gridArray.grid:get_vtx(0, j, 0)
    
    --iNodeIndx = 0; kNodeIndx = 0;
    --for i = 1, niv-1 do
        --xNormalised = vabs(FBA.gridArray.grid:get_vtx(i, j, 0) - ref_point) / Length

        --for k = 1, nkv-1 do
            --zNormalised = vabs(FBA.gridArray.grid:get_vtx(0, j, k) - ref_point) / Width

            --if zNormalised >= (kNodeIndx * dzFEM) then
                --nodeIndx = iNodeIndx * (Nz + 1) + kNodeIndx
                --blkId, CellIndx = ijkFBAIndxToFBIndx(FBA, i-1, j_cell, k-1)
                --if CellToNodeMap[blkId] then
                    --CellToNodeMap[blkId][nodeIndx] = CellIndx
                --else
                    --CellToNodeMap[blkId] = {[nodeIndx] = CellIndx}
                --end
                --kNodeIndx = kNodeIndx + 1
            --end
        --end

        --kNodeIndx = 0

        --if xNormalised >= (iNodeIndx * dxFEM) then
            --iNodeIndx = iNodeIndx + 1
        --end
    --end

    -- Compute normalized location of each CFD vertex
    xNormalised = {}; zNormalised = {};
    ref_point = FBA.gridArray.grid:get_vtx(0, j, 0)

    for i = 0, niv-1 do
        xNormalised[i] = vabs(FBA.gridArray.grid:get_vtx(i, j, 0) - ref_point) / Length
    end

    for k = 0, nkv-1 do
        zNormalised[k] = vabs(FBA.gridArray.grid:get_vtx(0, j, k) - ref_point) / Width
    end

    -- Now run through the FEM nodes and find out between which vertices do they reside
    mappingFile = assert(io.open("FSI/Weights/"..Side.."C2N.mapping", "w"))

    for i = 0, Nx do
        iVtx = 1
        for _iVtx = 1, niv-1 do
            if dxFEM * i <= xNormalised[_iVtx] then
                iVtx = _iVtx
                break
            end
        end

        for k = 0, Nz do
            if (config.dimensions == 3 and not FSIOptions.quasi3D) then
                kVtx = 1
                for _kVtx = 1, nkv-1 do
                    if dzFEM * k <= zNormalised[_kVtx] then
                        kVtx = _kVtx
                        break
                    end
                end
            else
                kVtx = 1
            end

            nodeIndx = k * (Nx + 1) + i
            blkId, CellIndx = ijkFBAIndxToFBIndx(FBA, iVtx-1, j_cell, kVtx-1)
            if CellToNodeMap[blkId] then
                CellToNodeMap[blkId][nodeIndx] = CellIndx
            else
                CellToNodeMap[blkId] = {[nodeIndx] = CellIndx}
            end
        end
    end

    mappingFile = assert(io.open("FSI/Weights/"..Side.."C2N.mapping", "w"))
    for blkId, CellIds in pairs(CellToNodeMap) do
        mappingFile:write(string.format("%d\n", blkId))
        for nodeId, _ in pairs(CellIds) do
            mappingFile:write(string.format("%s ", nodeId)) -- Key so it's a string
        end
        mappingFile:write("\n")
        for _, CellId in pairs(CellIds) do
            mappingFile:write(string.format("%d ", CellId))
        end
        mappingFile:write("\n")
    end
end

return fsi_weights
