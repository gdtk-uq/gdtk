
local fsi_check = {}

-- Make sure the forcing setup is valid
function fsi_check.checkForcingSetup()
    if FSIOptions.northForcing == "Fluid" then
        assert(FSIOptions.northFBA, "If using fluid pressure as the forcing for the north of the plate, the FBArray on the plate must be specified.")
        assert(FSIOptions.northFBA.myType == "FBArray", "The north block must be a FBArray")
    elseif FSIOptions.northForcing == "UserDefined" then
        assert(FSIOptions.udf_force_file, "Must supply a lua file containing the forcing definitions.")
    else
        assert(false, "Invalid setup for the forcing, check northForcing and northFBA.")
    end

    if FSIOptions.southForcing == "Fluid" then
        assert(FSIOptions.southFBA, "If using fluid pressure as the forcing for the south of the plate, the FBArray on the plate must be specified.")
        assert(FSIOptions.southFBA.myType == "FBArray", "The south block must be a FBArray")
    elseif FSIOptions.southForcing == "UserDefined" then
        assert(FSIOptions.udf_force_file, "Must supply a lua file containing the forcing definitions.")
    else
        assert(false, "Invalid setup for the forcing, check southForcing and southFBA.")
    end
end

-- Check that the size of the plate specified by the FSI setup is the same as the physical
-- size of the fluid block arrays used as the north and south forcing
function fsi_check.checkPlateDimensions()
    -- How closely should we match the dimensions
    tol = 1e-6
    if FSIOptions.northForcing == "Fluid" then
        FBA = FSIOptions.northFBA
        niv = FBA.niv; njv = FBA.njv; nkv = FBA.nkv
        p000 = FBA.gridArray.grid:get_vtx(0, 0, 0)
        p100 = FBA.gridArray.grid:get_vtx(niv-1, 0, 0)
        p001 = FBA.gridArray.grid:get_vtx(0, 0, nkv-1)
        length = vabs(p100 - p000); width = vabs(p001 - p000)
        assert(math.abs(length - FSIOptions.length) < tol, "Specified plate length and actual plate length are different")
        if config.dimensions == 3 then assert(math.abs(width - FSIOptions.width) < tol, "Specified plate width and actual plate width are different") end
    end
    if FSIOptions.southForcing == "Fluid" then
        FBA = FSIOptions.southFBA
        niv = FBA.niv; njv = FBA.njv; nkv = FBA.nkv
        p000 = FBA.gridArray.grid:get_vtx(0, njv-1, 0)
        p100 = FBA.gridArray.grid:get_vtx(niv-1, njv-1, 0)
        p001 = FBA.gridArray.grid:get_vtx(0, njv-1, nkv-1)
        length = vabs(p100 - p000); width = vabs(p001 - p000)
        assert(math.abs(length - FSIOptions.length) < tol, "Specified plate length and actual plate length are different")
        if config.dimensions == 3 then assert(math.abs(width - FSIOptions.width) < tol, "Specified plate width and actual plate width are different") end
    end
end

function fsi_check.checkConfigOptions()
    if config.FSIModel == "EulerBernoulli" then
        assert(config.dimensions == 3, "EulerBernoulli plate model (default) is a two dimensional model only. Use KirchhoffLove for 3 dimensions.")
    elseif config.FSIModel == "KirchhoffLove" then
        assert(config.dimensions == 2, "KirchhoffLove plate model is a three dimensional model only.")
    end

    assert(((FSIOptions.couplingStep % config.cfl_count == 0) or (config.cfl_count % FSIOptions.couplingStep == 0)) or config.fixed_time_step,
            "Configure config.cfl_count and FSIOptions.couplingStep so that one is a multiple of the other.")

end

return fsi_check
