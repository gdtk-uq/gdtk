-- mabey.lua : Turbulent flow over a flat plate (3D)
-- Dimir Y.X. Pot, Samuel J. Stennett, Wilson Y.K. Chan
-- Ported from Eilmer3, 2018-03-14
--    Mabey test case (AGARDograph 223 - Test series 7402)
--    (Referenced from Fernholz & Finley (1977),
--    AGARDograph No. 223, "A critical compilation of
--    compressible turbulent boundary layer data.")
-- 
config.title = "Mabey Mach 4.5 flow over a flat plate (k-omega) with wall function - 3D"
print(config.title)
config.dimensions = 3
config.spatial_deriv_calc = "least_squares" -- for 3D simulations
config.turbulence_model = "k_omega"
config.viscous = true
config.flux_calculator = "ausmdv"
config.gasdynamic_update_scheme = "classic-rk3"
--
config.max_time = 2.4e-3  -- About 4 flow lengths (1 flow length ~ 0.56 ms)
config.dt_plot = 0.2e-3
config.dt_history = 1.0e-3
config.max_step = 4000000
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
dz = 0.005 -- depth in z-direction (also k-index direction)

-- In 3D, look down on the top surface
--           wall
--       p01--------p11
-- flow=> |         |
--       p00        |
--          -\-     |
--    flow=>    -\- |
--        0        p10 ----> x
-- 
p000 = Vector3:new{x=0.0, y=3.0*H/4.0, z=0.0}
p100 = Vector3:new{x=L,   y=0.0,       z=0.0}
p110 = Vector3:new{x=L,   y=H,         z=0.0}
p010 = Vector3:new{x=0.0, y=H,         z=0.0}
p001 = Vector3:new{x=0.0, y=3.0*H/4.0, z=dz}
p101 = Vector3:new{x=L,   y=0.0,       z=dz}
p111 = Vector3:new{x=L,   y=H,         z=dz}
p011 = Vector3:new{x=0.0, y=H,         z=dz}

clusterx = RobertsFunction:new{end0=true, end1=false, beta=1.05}
clustery_e = RobertsFunction:new{end0=false, end1=true, beta=1.1400}
clustery_w = RobertsFunction:new{end0=false, end1=true, beta=1.7400}
clusterz = LinearFunction:new{}
cflist = {edge01=clusterx, edge12=clustery_e, edge32=clusterx, edge03=clustery_w,
	  edge45=clusterx, edge56=clustery_e, edge76=clusterx, edge47=clustery_w,
	  edge04=clusterz, edge15=clusterz,   edge26=clusterz, edge37=clusterz}
vol = TFIVolume:new{vertices={p000,p100,p110,p010,p001,p101,p111,p011}}
grd = StructuredGrid:new{pvolume=vol, cfList=cflist, niv=129, njv=87, nkv=11} 

bcList={north=WallBC_NoSlip_Adiabatic:new{wall_function=true},
	east=OutFlowBC_Simple:new{},
	south=InFlowBC_Supersonic:new{flowCondition=inflow},
        west=InFlowBC_Supersonic:new{flowCondition=inflow},
	bottom=WallBC_WithSlip:new{},
	top=WallBC_WithSlip:new{}}

blks = FBArray:new{grid=grd, nib=4, njb=2, nkb=1,
		       bcList=bcList, initialState=inflow}

identifyBlockConnections()
