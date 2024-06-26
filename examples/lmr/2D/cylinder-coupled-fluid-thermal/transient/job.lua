--
-- Conjugate heat transfer simulations of hypersonic flow over a hollow cylinder with finite thickness,
-- freestream condition and experimental data digitized from ref. [1].
--
-- References:
-- [1] Experimental Study of Shock Wave Interference Heating on a Cylindrical Leading Edge
--     A. R. Wieting
--     NASA Technical Memorandum, TM-100484
--
-- History:
-- 	2023-11-09  -- Eilmer 4 version
-- 		    Author: Kyle A. Damm
--
-- 	2024-04-29  -- Eilmer 5 version
-- 		    Authors: Saxon A. James, Rowan J. Gollan, Kyle A. Damm

config.solver_mode = "transient"
config.dimensions = 2

-- ==========================================================
-- Freestream conditions
-- ==========================================================
-- select gas model
nsp, nmodes, gmodel = setGasModel('ideal-air.gas')
gs = GasState:new{gmodel}

-- run 37 conditions taken taken from Appendix D in ref. [1]
Mach       = 6.47
T_wall     = 294.444   -- K
gs.T       = 241.5     -- K
gs.p       = 648.11    -- Pa
V_inf      = 2034.235  -- m/s

-- update some gas properties
gmodel:updateThermoFromPT(gs)
gmodel:updateTransCoeffs(gs)
gmodel:updateSoundSpeed(gs)

-- set inflow and inital condition
inflow = FlowState:new{p=gs.p, T=gs.T, velx=V_inf}
initial = inflow

-- Set up solid thermal model
registerSolidModels{
   stainless_steel_321 = ConstantPropertiesModel:new{rho=8030.0, k=16.24, Cp=502.48},
}


-- Set up the geometry for defining the grid
Ro = 3.81e-02 -- outer nose radis, metres
Ri = 2.52e-02 -- inner nose radius, metres

-- define fluid grids
n0i = 100
n0j = 50

a = Vector3:new{x=0.0, y=0.0}
b = Vector3:new{x=-Ro, y=0.0}
c = Vector3:new{x=0.0, y=Ro}
d = { Vector3:new{x=-1.5*Ro, y=0.0}, Vector3:new{x=-1.5*Ro, y=Ro},
      Vector3:new{x=-Ro, y=2*Ro}, Vector3:new{x=0.0, y=3*Ro} }

sphere_edge = Arc:new{p0=b, p1=c, centre=a}
psurf0 = makePatch{north=Line:new{p0=d[#d], p1=c}, south=Line:new{p0=d[1], p1=b},
									 east=sphere_edge, west=Bezier:new{points=d}}
cf_radial0 = RobertsFunction:new{end0=false, end1=true, beta=1.2}

registerFluidGridArray{ 
  grid=StructuredGrid:new{psurface=psurf0, niv=n0i+1, njv=n0j+1,
                          cfList={north=cf_radial0, south=cf_radial0}},
  fsTag='initial_fluid',
  nib = 2,
  njb = 2,
  bcTags={north='outflow',
          south='symmetry',	
          west='inflow'},
}

-- define solid grids
n1i = 25

e = Vector3:new{x=-Ri, y=0.0}
f = Vector3:new{x=0.0, y=Ri}

sphere_outer_edge = Arc:new{p0=b, p1=c, centre=a}
sphere_inner_edge = Arc:new{p0=e, p1=f, centre=a}
psurf1 = makePatch{north=Line:new{p0=c, p1=f}, south=Line:new{p0=b, p1=e},
									 east=sphere_inner_edge, west=sphere_outer_edge}
cf_radial1 = RobertsFunction:new{end0=true, end1=false, beta=1.1}

registerSolidGridArray{
  grid=StructuredGrid:new{psurface=psurf1, niv=n1i+1, njv=n0j+1,
                          cfList={north=cf_radial1, south=cf_radial1}},
  ssTag='init_solid_temperature',
  modelTag='stainless_steel_321',
  nib = 2,
  njb = 2,
  bcTags={north='adiabatic',
          east='fixed_temperature',
          south='solid_symmetry',
          }
}

identifyGridConnections()

flowDict = {
	initial_fluid=initial
}

flowBCDict = {
	outflow=OutFlowBC_SimpleExtrapolate:new{},
	symmetry=WallBC_WithSlip:new{},
	inflow=InFlowBC_Supersonic:new{flowState=inflow}
}

makeFluidBlocks(flowBCDict, flowDict)

solidBCDict = {
	adiabatic=SolidAdiabaticBC:new{},
	fixed_temperature=SolidFixedTBC:new{T=Twall},
	solid_symmetry=SolidConstantFluxBC:new{}
}


initSolidState = {
  init_solid_temperature=T_wall
}


makeSolidBlocks(solidBCDict, initSolidState)

-- every MPI task needs a FluidBlock, so make sure you distribute the blocks accordingly
-- e.g. here we have 4 FluidBlocks + 4 SolidBlocks divided amongst 4 MPI tasks
mpiTasks = mpiDistributeBlocks{ntasks=4, dist="load-balance", preassign={[0]=1}}

-- ==========================================================
-- Solver settings
-- ==========================================================

-- invsicid flux settings
config.flux_calculator = "ausmdv"
config.interpolation_order = 2
config.apply_limiter = true
config.extrema_clipping = true

-- viscous flux settings
config.viscous = true
config.spatial_deriv_calc = "least_squares"
config.spatial_deriv_locn = "cells"
config.viscous_signal_factor = 0.1

-- solid domain settings
config.gasdynamic_update_scheme = 'pc'
config.cfl_value = 0.4
config.dt_init = 1.0e-12
config.max_time = 0.001 -- seconds
config.max_step = 8e7
config.dt_plot = config.max_time/001.0

-- loads settings (adjacent_to_solid is a special internal boundary name for coupled fluid/solid boundaries)
config.boundary_groups_for_loads = "adjacent_to_solid"
config.write_loads = true
config.dt_loads = config.max_time/500
