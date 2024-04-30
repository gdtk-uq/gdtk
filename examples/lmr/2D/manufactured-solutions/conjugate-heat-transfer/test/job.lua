--
-- Authors: KAD and RJG
-- Date: 2024-03-03
-- History: Ported from eilmer4 example
--

config.solver_mode = 'transient'
config.dimensions = 2

ncells = 8

-- Pull gas constants from constants.txt
dofile("constants.txt")

setGasModel('very-viscous-air.lua')
-- Load initial condition fill functions from
-- auto-generated files
dofile('fill-fn.lua')
dofile('fill-solid-fn.lua')
p_init = 1.0e5
T_init = 350.0
initial = FlowState:new{p=p_init, T=T_init, velx=0.0, vely=0.0}
flowStates = {initial=gasFillFn}
initSolidState = {
   initial_solid_field=solidFillFn
}
registerSolidModels{
   solid_properties = ConstantPropertiesModel:new{rho=rho_s, k=k_s, Cp=Cp_s}
}

pts = {
   a = {x=0.0, y=0.0},
   b = {x=L, y=0.0},
   c = {x=0.0, y=L},
   d = {x=L, y=L},
   e = {x=0.0, y=1.5*L},
   f = {x=L, y=1.5*L}
}

fluidPatch = CoonsPatch:new{p00=pts.a, p10=pts.b, p11=pts.d, p01=pts.c}
solidPatch = CoonsPatch:new{p00=pts.c, p10=pts.d, p11=pts.f, p01=pts.e}
nicell = ncells
njcell0 = ncells
njcell1 = math.floor(njcell0/2)
--
registerFluidGrid{
   grid=StructuredGrid:new{psurface=fluidPatch, niv=nicell+1, njv=njcell0+1},
   fsTag="initial",
   bcTags={west="udf", south="udf", east="udf"} 
}
--
registerSolidGrid{
   grid=StructuredGrid:new{psurface=solidPatch, niv=nicell+1, njv=njcell1+1},
   ssTag="initial_solid_field",
   modelTag="solid_properties",
   bcTags={west="udf_solid", north="udf_solid", east="udf_solid"}
}
--
identifyGridConnections()
--
-- Now we turn our attention to simulation setup
--
-- Fluid BC
udfBC = UserDefinedBC:new{fileName='udf-bc.lua'}
fluidBCs = {udf=udfBC}
-- Solid BC
udfSolidBC = SolidUserDefinedBC:new{fileName='udf-solid-bc.lua'}
solidBCs = {udf_solid=udfSolidBC}
--
makeFluidBlocks(fluidBCs, flowStates)
makeSolidBlocks(solidBCs, initSolidState, solidProps)
-- 
config.apply_bcs_in_parallel = false
config.interpolation_order = 2
config.gasdynamic_update_scheme = 'euler'
config.flux_calculator = 'ausm_plus_up'
config.M_inf = 0.1
config.viscous = true
config.spatial_deriv_calc = 'divergence'
config.spatial_deriv_locn = 'faces'
-- Do NOT use the limiters for the verification tests
config.apply_limiter = false
config.extrema_clipping = false

config.max_time = 1.0
config.max_step = 10000000
config.dt_init = 1.0e-5
config.fixed_time_step = true
config.dt_plot = config.max_time/40

config.udf_source_terms = true
config.udf_source_terms_file = 'udf-source-terms.lua'
config.udf_solid_source_terms = true
config.udf_solid_source_terms_file = 'udf-solid-source-terms.lua'
