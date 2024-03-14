print("Sod's Classic Shock tube with moving-shock boundary")
-- Adapted from Eilmer3 and Eilmer4 examples
-- Authors: KAD, 22-Dec-2015, reworked PJ 2021-08-09, 2024-03-14
--
-- We start the calculation with a single region for the driver gas.
-- The driven gas state will appear only as part of the moving-shock
-- boundary condition.
patch_r = CoonsPatch:new{p00={x=0.5, y=0.0}, p01={x=0.5, y=0.1},
                         p10={x=1.0, y=0.0}, p11={x=1.0, y=0.1}}
grid0 = registerFluidGridArray{
   grid=StructuredGrid:new{psurface=patch_r, niv=201, njv=3},
   nib=2, njb=1,
   fsTag="initial",
   bcTags={west="inflow_sf"},
   shock_fitting=true
}
--
-- Initial gas conditions.
nsp, nmodes, gm = setGasModel('ideal-air.gas')
R = 287.1           -- J/kg.K
-- Right, driver gas
p_r = 1.0e5          -- Pa
rho_r = 1.0          -- kg/m^3
T_r = p_r/(rho_r*R)  -- K
driver = FlowState:new{p=p_r, T=T_r}
-- Left, driven gas
p_l = 1.0e4          -- Pa
rho_l = 0.125        -- kg/m^3
T_l = p_l/(rho_l*R)  -- K
driven = FlowState:new{p=p_l, T=T_l}
--
flowDict = {initial=driver}
bcDict = {inflow_sf=InFlowBC_ShockFitting:new{flowCondition=driven}}
makeFluidBlocks(bcDict, flowDict)
--
-- Run parameters
config.solver_mode = "transient"
config.gasdynamic_update_scheme = "moving-grid-2-stage"
config.flux_calculator = "ausmdv"
config.max_time = 6.0e-4  -- seconds
config.max_step = 3000
config.cfl_value = 0.25
config.dt_plot = 10.0e-6
config.grid_motion = "shock_fitting"
-- Shock fitting starts immediately, by default.
-- The default scale_factor is 0.5,
-- but we want to be time-accurate with the grid motion.
config.shock_fitting_scale_factor = 1.0
