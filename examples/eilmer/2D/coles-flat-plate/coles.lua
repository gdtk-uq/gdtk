-- coles.lua : Turbulent flow over a flat plate
--  Anand V  2016-05-26, Indooroo
-- Ported PJ's example from Eilmer3
-- 

-- We can set individual attributes of the global data object.
config.title = "Coles Mach 3.7 flow over a flat plate (k-omega)"
print(config.title)
config.dimensions = 2

-- Gas model and flow conditions
nsp, nmodes = setGasModel('ideal-air-gas-model.lua')
p_inf = 1.358e3  -- Pa
u_inf = 677.4    -- m/s
T_inf = 83.34    -- K

-- Turbulence quantities estimate
tke_inf = 1.5 * (u_inf * 0.05)^2
rho_inf = p_inf / (287.0 * T_inf)
mu_ref = 17.89e-6; T_ref = 273.1; S = 110.4
T_T0 = T_inf / T_ref
mu_t_inf = 100.0*( mu_ref * (T_ref + S)/(T_inf + S) * T_T0 * math.sqrt(T_T0))
omega_inf = rho_inf * tke_inf / mu_t_inf
print("Inflow turbulence: tke=", tke_inf, "omega=", omega_inf)

-- Inflow Set-Up
inflow = FlowState:new{p=p_inf, T=T_inf, velx=u_inf, vely=0.0, 
		       tke=tke_inf, omega=omega_inf}
--print("Inflow Check\n", inflow)

--Geometry of the flow domain
L = 0.60 -- metres
H = 0.4 * L
NB = 4 	 -- number of blocks

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
-- north = Line(c,b); east = Line(a,b); south = Line(d,a); west = Line(d,c)
patch = makePatch{north=Line:new{p0=c, p1=b}, east=Line:new{p0=a, p1=b},
		  south=Line:new{p0=d, p1=a}, west=Line:new{p0=d, p1=c}}

grd = StructuredGrid:new{psurface=patch,
			 cfList = {north=RobertsFunction:new{end0=true,end1=false,beta=1.05},
				   east=RobertsFunction:new{end0=false,end1=true,beta=1.01},
				  south=RobertsFunction:new{end0=true,end1=false,beta=1.05},
				  west=RobertsFunction:new{end0=false,end1=true,beta=1.05}},
			 niv=100, njv=90}

-- Define the blocks, boundary conditions and set the discretisation.
blks = SBlockArray{grid=grd, nib=NB, njb=1, 
		   fillCondition=inflow,
		   bcList={north=WallBC_NoSlip_Adiabatic:new{},
			   east=OutFlowBC_Simple:new{},
			   south=InFlowBC_Supersonic:new{flowCondition=inflow},
			   west=InFlowBC_Supersonic:new{flowCondition=inflow}
}
}
identifyBlockConnections()

-- Do a little more setting of global data.
config.turbulence_model = "k_omega"
config.viscous = true
config.flux_calculator = 'adaptive'
 	
config.max_time = 5.0e-3  -- should allow a few flow lengths   
config.dt_plot =  1.0e-3
config.dt_history = 1.0e-5
config.max_step = 3000000 

config.gasdynamic_update_scheme = "classic-rk3"
config.cfl_value = 0.5	
config.cfl_count = 3
config.stringent_cfl = false -- true is more robust
config.dt_init = 1.0e-9	-- only an initial guess, the simulation will take this over
