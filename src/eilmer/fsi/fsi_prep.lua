-- Writing the FSI config files to JSON
local fsi_prep = {}

local fsi_check = require "fsi_check"
checkForcingSetup = fsi_check.checkForcingSetup
checkPlateDimensions = fsi_check.checkPlateDimensions
checkConfigOptions = fsi_check.checkConfigOptions

local fsi_weights = require "fsi_weights"
generateDistanceWeights = fsi_weights.generateDistanceWeights
generateInterpolationWeights = fsi_weights.generateInterpolationWeights
writeWeightsToFile = fsi_weights.writeWeightsToFile
cellToNodeMapping = fsi_weights.cellToNodeMapping
cellToQuadratureMapping = fsi_weights.cellToQuadratureMapping

local fsiconfig = require "fsioptions"

function fsi_prep.write_FSI_config(job)
    fsi_prep.check_FSI_setup()
    fsi_prep.write_FSI_config_file("config/" .. job .. ".fsi")
    fsi_prep.write_FSI_weights()
end

function fsi_prep.write_FSI_config_file(fileName)
    f = assert(io.open(fileName, "w+"))
    f:write("{\n")
    f:write(string.format('"northForcing" : "%s",\n', FSIOptions.northForcing))
    f:write(string.format('"southForcing" : "%s",\n', FSIOptions.southForcing))
    f:write(string.format('"Nx" : %d,\n', FSIOptions.Nx))
    f:write(string.format('"Nz" : %d,\n', FSIOptions.Nz))
    f:write(string.format('"length" : %1.8e,\n', FSIOptions.length))
    f:write(string.format('"width" : %1.8e,\n', FSIOptions.width))
    f:write(string.format('"thickness" : %1.8e,\n', FSIOptions.thickness))
    f:write(string.format('"density" : %1.8e,\n', FSIOptions.density))
    f:write(string.format('"plateNormal" : [%1.8e, %1.8e, %1.8e],\n', FSIOptions.plateNormal[1], FSIOptions.plateNormal[2], FSIOptions.plateNormal[3]))
    f:write(string.format('"youngsModulus" : %1.8e,\n', FSIOptions.youngsModulus))
    f:write(string.format('"poissonsRatio" : %1.8e,\n', FSIOptions.poissonsRatio))
    f:write(string.format('"quasi3D" : %s,\n', FSIOptions.quasi3D))
    f:write(string.format('"BCs" : "%s",\n', FSIOptions.BCs))
    f:write(string.format('"couplingStep" : %d,\n', FSIOptions.couplingStep))
    f:write('"movingBlks" : [ ')
    movingBlks = getMovingBlocks(FSIOptions)
    for _i, _blkId in ipairs(movingBlks) do
        if _i == #movingBlks then
            f:write(string.format("%d ]\n", _blkId))
        else
            f:write(string.format("%d, ", _blkId))
        end
    end
    f:write("}")
end

function getMovingBlocks()
    -- Pull out all the blocks that may be moving
    movingBlks = {}
    possibleEntries = {"northFBA", "northWestFBA", "northEastFBA", "northBottomFBA", "northTopFBA", "southFBA", "southWestFBA", "southEastFBA", "southBottomFBA", "southTopFBA", "westAdjacentFBA", "eastAdjacentFBA"}

    for _, entry in ipairs(possibleEntries) do
        if FSIOptions[entry] then
            for _, blk in pairs(FSIOptions[entry].blockCollection) do
                movingBlks[#movingBlks+1] = blk.id
            end
        end
    end

    return movingBlks
end

function fsi_prep.check_FSI_setup()
    -- Check forcing is valid
    fsi_check.checkForcingSetup()

    -- Check plate dimensions match
    fsi_check.checkPlateDimensions()

    -- Check the config options are compatible
    fsi_check.checkConfigOptions()
end

function fsi_prep.write_FSI_weights()
    os.execute("mkdir -p FSI")
    os.execute("mkdir -p FSI/Weights")
    -- Necessary distance weightings
    print("Start writing blocks")
    if FSIOptions.northFBA then
        NorthDistanceWeights = generateDistanceWeights(FSIOptions.northFBA, false)
        NorthInterpWeights, NorthInterpIndices = generateInterpolationWeights(FSIOptions.northFBA, FSIOptions.Nx, FSIOptions.Nz, FSIOptions.length, FSIOptions.width, "north")
        writeWeightsToFile(FSIOptions.northFBA, NorthDistanceWeights, NorthInterpWeights, NorthInterpIndices, "plate")

        if FSIOptions.northWestFBA then
            writeWeightsToFile(FSIOptions.northWestFBA, NorthDistanceWeights, NorthInterpWeights, NorthInterpIndices, "west")
        end
        if FSIOptions.northEastFBA then
            writeWeightsToFile(FSIOptions.northEastFBA, NorthDistanceWeights, NorthInterpWeights, NorthInterpIndices, "east")
        end
        if FSIOptions.northBottomFBA then
            writeWeightsToFile(FSIOptions.northBottomFBA, NorthDistanceWeights, NorthInterpWeights, NorthInterpIndices, "bottom")
        end
        if FSIOptions.northNorthFBA then
            writeWeightsToFile(FSIOptions.northTopFBA, NorthDistanceWeights, NorthInterpWeights, NorthInterpIndices, "top")
        end
        if FSIOptions.northWestBottomFBA then
            writeWeightsToFile(FSIOptions.northWestBottomFBA, NorthDistanceWeights, NorthInterpWeights, NorthInterpIndices, "west-bottom")
        end
        if FSIOptions.northWestTopFBA then
            writeWeightsToFile(FSIOptions.northWestTopFBA, NorthDistanceWeights, NorthInterpWeights, NorthInterpIndices, "west-top")
        end
        if FSIOptions.northEastBottomFBA then
            writeWeightsToFile(FSIOptions.northEastBottomFBA, NorthDistanceWeights, NorthInterpWeights, NorthInterpIndices, "east-bottom")
        end
        if FSIOptions.northEastTopFBA then
            writeWeightsToFile(FSIOptions.northEastTopFBA, NorthDistanceWeights, NorthInterpWeights, NorthInterpIndices, "east-top")
        end
        if FSIOptions.westAdjacentFBA then
            writeWeightsToFile(FSIOptions.westAdjacentFBA, NorthDistanceWeights, NorthInterpWeights, NorthInterpIndices, "west-adjacent")
        end
        if FSIOptions.eastAdjacentFBA then
            writeWeightsToFile(FSIOptions.eastAdjacentFBA, NorthDistanceWeights, NorthInterpWeights, NorthInterpIndices, "east-adjacent")
        end
        if FSIOptions.bottomAdjacentFBA then
            writeWeightsToFile(FSIOptions.bottomAdjacentFBA, NorthDistanceWeights, NorthInterpWeights, NorthInterpIndices, "bottom-adjacent")
        end
        if FSIOptions.topAdjacentFBA then
            writeWeightsToFile(FSIOptions.topAdjacentFBA, NorthDistanceWeights, NorthInterpWeights, NorthInterpIndices, "top-adjacent")
        end
        if FSIOptions.westBottomAdjacentFBA then
            writeWeightsToFile(FSIOptions.westBottomAdjacentFBA, NorthDistanceWeights, NorthInterpWeights, NorthInterpIndices, "west-bottom-adjacent")
        end
        if FSIOptions.westTopAdjacentFBA then
            writeWeightsToFile(FSIOptions.westTopAdjacentFBA, NorthDistanceWeights, NorthInterpWeights, NorthInterpIndices, "west-top-adjacent")
        end
        if FSIOptions.eastBottomAdjacentFBA then
            writeWeightsToFile(FSIOptions.eastBottomAdjacentFBA, NorthDistanceWeights, NorthInterpWeights, NorthInterpIndices, "east-bottom-adjacent")
        end
        if FSIOptions.eastTopAdjacentFBA then
            writeWeightsToFile(FSIOptions.eastTopAdjacentFBA, NorthDistanceWeights, NorthInterpWeights, NorthInterpIndices, "east-top-adjacent")
        end
    end
    if FSIOptions.southFBA then
        SouthDistanceWeights = generateDistanceWeights(FSIOptions.southFBA, true)
        SouthInterpWeights, SouthInterpIndices = generateInterpolationWeights(FSIOptions.southFBA, FSIOptions.Nx, FSIOptions.Nz, FSIOptions.length, FSIOptions.width, "south")
        writeWeightsToFile(FSIOptions.southFBA, SouthDistanceWeights, SouthInterpWeights, SouthInterpIndices, "plate")

        if FSIOptions.southWestFBA then
            writeWeightsToFile(FSIOptions.southWestFBA, SouthDistanceWeights, SouthInterpWeights, SouthInterpIndices, "west")
        end
        if FSIOptions.southEastFBA then
            writeWeightsToFile(FSIOptions.southEastFBA, SouthDistanceWeights, SouthInterpWeights, SouthInterpIndices, "east")
        end
        if FSIOptions.southBottomFBA then
            writeWeightsToFile(FSIOptions.southBottomFBA, SouthDistanceWeights, SouthInterpWeights, SouthInterpIndices, "bottom")
        end
        if FSIOptions.southTopFBA then
            writeWeightsToFile(FSIOptions.southTopFBA, SouthDistanceWeights, SouthInterpWeights, SouthInterpIndices, "top")
        end
        if FSIOptions.southWestBottomFBA then
            writeWeightsToFile(FSIOptions.southWestBottomFBA, SouthDistanceWeights, SouthInterpWeights, SouthInterpIndices, "west-bottom")
        end
        if FSIOptions.southWestTopFBA then
            writeWeightsToFile(FSIOptions.southWestTopFBA, SouthDistanceWeights, SouthInterpWeights, SouthInterpIndices, "west-top")
        end
        if FSIOptions.southEastBottomFBA then
            writeWeightsToFile(FSIOptions.southEastBottomFBA, SouthDistanceWeights, SouthInterpWeights, SouthInterpIndices, "east-bottom")
        end
        if FSIOptions.southEastTopFBA then
            writeWeightsToFile(FSIOptions.southEastTopFBA, SouthDistanceWeights, SouthInterpWeights, SouthInterpIndices, "east-top")
        end
        -- Just repeat these 2 as insurance in case we dont have a top FBA
        if FSIOptions.westAdjacentFBA then
            writeWeightsToFile(FSIOptions.westAdjacentFBA, SouthDistanceWeights, SouthInterpWeights, SouthInterpIndices, "west-adjacent")
        end
        if FSIOptions.eastAdjacentFBA then
            writeWeightsToFile(FSIOptions.eastAdjacentFBA, SouthDistanceWeights, SouthInterpWeights, SouthInterpIndices, "east-adjacent")
        end
    end
    if FSIOptions.northFBA then
        cellToQuadratureMapping(FSIOptions.northFBA, FSIOptions.Nx, FSIOptions.Nz, FSIOptions.length, FSIOptions.width, "north")
    end

    if FSIOptions.southFBA then
        cellToQuadratureMapping(FSIOptions.southFBA, FSIOptions.Nx, FSIOptions.Nz, FSIOptions.length, FSIOptions.width, "south")
    end
end

return fsi_prep
