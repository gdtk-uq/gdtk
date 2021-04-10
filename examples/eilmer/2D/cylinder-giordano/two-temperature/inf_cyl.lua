--
-- Authors: Rowan G. and Peter J.
-- Date: 2017-07-13
--
-- This is a port of the Eilmer3 example in:
-- cfcfd3/examples/eilmer3/2D/giordano
-- 
-- -----------------------------
-- Some notes from original file
-- -----------------------------
--
-- This file can be used to simulate the test case reported by:
-- Giordano, et al. (1997)
-- Vibrationally Relaxing Flow of N2 past an Infinite Cylinder
-- Journal of Thermophysics and Heat Transfer, vol 11, no 1, pp 27 - 35
--
-- Description:
-- Pure N2 flow (chemically inert)
-- Mach = 6.5
-- T_inf = 300 K
-- p_inf = 50 - 500 Pa
--
-- ------------------
-- End original notes
-- ------------------

Rc = 1.0 -- m, cylinder radius
config.dimensions = 2
config.title = "Inifinite cylinder in nitrogen flow, M = 6.5"

nsp, nmodes, gm = setGasModel("two-temp-n2.lua")
config.reacting = true
config.reactions_file = "VT-relaxation-time-selection.lua"

-- Set flow conditions
M_inf = 6.5
T_inf = 300.0 -- K
p_inf = 50.0 -- Pa, low pressure case
inflow_gas = FlowState:new{p=p_inf, T=T_inf, T_modes=T_inf}
velx_inf = M_inf*inflow_gas.a

inflow = FlowState:new{p=p_inf, T=T_inf, T_modes=T_inf, velx=velx_inf}
initial = FlowState:new{p=p_inf/3, T=T_inf, T_modes=T_inf}

-- Geometry
a = Vector3:new{x=-Rc, y=0.0}
b = Vector3:new{x=0.0, y=Rc}
c = Vector3:new{x=0.0, y=0.0}

inflowBndry = Spline2:new{filename='inflow-spline-pts.dat'}
axis = Line:new{p0=inflowBndry(0.0), p1=a}
body = Arc:new{p0=a, p1=b, centre=c}
outflowBndry = Line:new{p0=inflowBndry(1.0), p1=b}

-- Grid
nxcells = 40
nycells = 120

quad0 = CoonsPatch:new{north=outflowBndry, east=body,
		       south=axis, west=inflowBndry}
grid0 = StructuredGrid:new{psurface=quad0, niv=nxcells+1, njv=nycells+1}

-- Blocks
FBArray:new{grid=grid0, initialState=initial, nib=1, njb=4,
		bcList={north=OutFlowBC_Simple:new{},
			west=InFlowBC_Supersonic:new{flowState=inflow}}
}

setHistoryPoint{x=a.x, y=a.y}

-- Finish off with some other configuration
config.flux_calculator = 'adaptive'
config.gasdynamic_update_scheme = 'euler'
config.max_time = 15.0e-3 -- s
config.max_step = 30000
config.cfl_value = 0.4
config.dt_init = 1.0e-8
config.dt_plot = config.max_time/50.0
config.dt_history = 1.0e-5







