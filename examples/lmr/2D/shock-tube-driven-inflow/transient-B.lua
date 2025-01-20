print("Sod's Classic Shock tube, just the driven tube, InFlowBC_TransientProfile")
-- Adapted from InFlowBC_StaticProfile example, PJ 2025-01-20
--
-- We start the calculation with a single region for the driven gas.
-- The effect of the driver gas enters via the inflow boundary condition.
--
-- Initial gas conditions.
nsp, nmodes, gm = setGasModel('ideal-air.gas')
R = 287.1           -- J/kg.K
-- Left, driver gas (not actually used in this lmr simulation)
p_r = 1.0e5          -- Pa
rho_r = 1.0          -- kg/m^3
T_r = p_r/(rho_r*R)  -- K
driver = FlowState:new{p=p_r, T=T_r}
-- Right, driven gas
p_l = 1.0e3          -- Pa
rho_l = 0.0125       -- kg/m^3
T_l = p_l/(rho_l*R)  -- K
driven = FlowState:new{p=p_l, T=T_l}
-- The driving boundary condition.
profile = InFlowBC_TransientProfile:new{fileName="transient-profile.zip",
                                        match="AyA-to-AyA"}
--
flowDict = {initial=driven}
bcDict = {inflow=profile}
makeFluidBlocks(bcDict, flowDict)
--
-- Run parameters
config.solver_mode = "transient"
config.gasdynamic_update_scheme = "classic-rk3"
config.flux_calculator = "ausmdv"
config.max_time = 500.0e-6  -- seconds
config.max_step = 3000
config.cfl_value = 0.5
config.dt_init = 1.0e-6
config.dt_plot = 100.0e-6
