fsiOptionsHidden = {
    -- What is the forcing function for the plate
    northForcing = "Fluid",
    southForcing = "Fluid",

    -- Discretization of the plate
    Nx = 5,
    Nz = 5,

    -- The mapping of the fluid domain relative to the moving plate
    northFBA = false,
    southFBA = false,
    
    northWestFBA = false,
    northEastFBA = false,
    northBottomFBA = false,
    northTopFBA = false,

    southWestFBA = false,
    southEastFBA = false,
    southBottomFBA = false,
    southTopFBA = false,

    westAdjacentFBA = false,
    eastAdjacentFBA = false,
    bottomAdjacentFBA = false,
    topAdjacentFBA = false,
    
    -- Characteristics of the plate
    length = 1.0,
    width = 1.0,
    thickness = 0.001,
    density = 8e3,
    youngsModulus = 190e9,
    poissonsRatio = 0.33,
    
    -- Plate orientation
    plateNormal = {0.0, 1.0, 0.0},

    -- BCs
    BCs = "CFFF",

    -- Every how many time steps do we update the plate
    couplingStep = 10,

    __index = function (t, k)
       return fsiOptionsHidden[k]
    end,
    __newindex = function (t, k, v)
       if fsiOptionsHidden[k] == nil then
           print(string.format("The field '%s' cannot be set in 'FSIOptions' table.", k))
       else
           fsiOptionsHidden[k] = v
       end
    end,
    __call = function (_, t)
       for k, v in pairs(t) do
           fsiOptionsHidden.__newindex(t, k, v)
       end
    end
}

FSIOptions = {}
setmetatable(FSIOptions, fsiOptionsHidden)
