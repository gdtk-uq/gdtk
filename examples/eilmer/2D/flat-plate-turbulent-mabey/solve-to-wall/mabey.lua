-- mabey.lua : Turbulent flow over a flat plate
-- Dimir Y.X. Pot, Samuel J. Stennett, Wilson Y.K. Chan
-- Ported from Eilmer3, 2018-03-14
--    Mabey test case (AGARDograph 223 - Test series 7402)
--    (Referenced from Fernholz & Finley (1977),
--    AGARDograph No. 223, "A critical compilation of
--    compressible turbulent boundary layer data.")
-- 
config.title = "Mabey Mach 4.5 flow over a flat plate (k-omega)"
print(config.title)
config.dimensions = 2
config.turbulence_model = "k_omega"
config.viscous = true
config.flux_calculator = "ausmdv"
config.gasdynamic_update_scheme = "classic-rk3"
--
config.max_time = 2.2e-3  -- About 4 flow lengths (1 flow length ~ 0.56 ms)
config.dt_plot =  0.2e-3
config.dt_history = 1.0e-3
config.max_step = 3000000
--
config.cfl_value = 0.4
config.cfl_count = 3
config.stringent_cfl = false -- true is more robust
config.dt_init = 1.0e-9 

-- Gas model and flow conditions to match Mabey's data set 74021802
nsp, nmodes, gm = setGasModel('ideal-air-gas-model.lua')
p_inf = 3.16e3  -- Pa
u_inf = 712.9   -- m/s
T_inf = 62.16   -- K

-- Set up gas state and update thermodynamic transfer coefficients
gas_inf = GasState:new{gm}
gas_inf.p = p_inf; gas_inf.T = T_inf
gm:updateThermoFromPT(gas_inf)
gm:updateSoundSpeed(gas_inf)
gm:updateTransCoeffs(gas_inf)

-- Use updated gas properties to estimate turbulence quantities
turb_intensity = 0.01
turb_lam_viscosity_ratio = 1.0
tke_inf = 1.5 * (turb_intensity * u_inf)^2
mu_t_inf = turb_lam_viscosity_ratio * gas_inf.mu
omega_inf = gas_inf.rho * tke_inf / mu_t_inf

-- Set up flow state
inflow = FlowState:new{p=p_inf, T=T_inf, velx=u_inf, tke=tke_inf, omega=omega_inf}
print("Inflow Check\n", inflow)

-- Geometry of the flow domain
L = 0.40 -- metres
H = 0.40 * L
--
--         wall
--        c---------b
-- flow=> |         |
--        d         |
--          -\-     |
--    flow=>    -\- |
--        0         a ----> x
-- 
a = Vector3:new{x=L, y=0.0}; b = Vector3:new{x=L, y=H}; 
c = Vector3:new{x=0.0, y=H}; d = Vector3:new{x=0.0, y=3.0*H/4.0}
patch = CoonsPatch:new{p00=d, p10=a, p11=b, p01=c}
cfx = RobertsFunction:new{end0=true,end1=false,beta=1.05}
cflist = {north=cfx, east=RobertsFunction:new{end0=false,end1=true,beta=1.0014},
	  south=cfx, west=RobertsFunction:new{end0=false,end1=true,beta=1.0074}}
grd = StructuredGrid:new{psurface=patch, niv=129, njv=87, cfList=cflist}

blks = FBArray:new{grid=grd, nib=2, njb=2, fillCondition=inflow,
		       bcList={north=WallBC_NoSlip_Adiabatic:new{wall_function=false},
			       east=OutFlowBC_Simple:new{},
			       south=InFlowBC_Supersonic:new{flowCondition=inflow},
			       west=InFlowBC_Supersonic:new{flowCondition=inflow}}}

identifyBlockConnections()
