-- Author: RJG
-- Date: 2016-03-29
--
-- Mach 10 flow of ideal air over a 6-deg wedge.
--
-- Note: This case can be easily adapted to suit the oblique
--       shock test case given in AIAA CFD Code Verficiation Project.
--
-- Reference:
-- Ghia et al. (2010)
-- The AIAA CFD Code Verification Project -- Test Cases for CFD Code Verification
-- Paper no. AIAA-2010-125, 48th AIAA Aerospace Sciences Meeting and Exhibit,
-- Orlando, Florida
--

-- Flow conditions, free stream
M_inf = 10.0
beta = 6.0
p_inf = 1.0e5
T_inf = 300.0

-- Grid dimensions
nx = 40
ny = 40

config.title = string.format("Mach %.1f flow over a %.1f-deg wedge.", M_inf, beta)
print(config.title)

config.dimensions = 2

nsp, nmodes, gm = setGasModel('ideal-air-gas-model.lua')
-- Compute inflow speed from Mach and free-stream conditions
Q = GasState:new{gm}
Q.p = p_inf
Q.T = T_inf
Q.massf = {air=1.0}
gm:updateThermoFromPT(Q)
gm:updateSoundSpeed(Q)
a_inf = Q.a
u_inf = M_inf*a_inf
rho = Q.rho
u = u_inf

-- Set up inflow condition
inflow = FlowState:new{p=p_inf, T=T_inf, velx=u_inf}

-- Set up geometry
-- Origin is at start of wedge.
--
--
--            D        OutFlow BC            E
--            +------------------------------+
--            |                              |
-- Inflow BC  |                              | OutFlow BC
--            |                              |
--            |                              + C
--            |                             / 
--            |                          /  
--            |                       /   Impermeable Wall BC
--            |                    /  
--            |                 /   
--            +--------------+
--            A    ^^^       B
--                  |
--                Impermeable Wall BC
--

upstreamLen = 0.1
wallLen = 0.3 -- in x dimension
domainHght = 0.1

b = Vector3:new{x=0.0, y=0.0}
a = Vector3:new{x=b.x-upstreamLen, y=0.0}
c = Vector3:new{x=wallLen, y=wallLen*math.tan(math.rad(beta))}
d = Vector3:new{x=a.x, y=domainHght}
e = Vector3:new{x=c.x, y=domainHght}

abc = Polyline:new{segments={Line:new{p0=a, p1=b}, Line:new{p0=b, p1=c}}}
ad = Line:new{p0=a, p1=d}
ce = Line:new{p0=c, p1=e}
de = Line:new{p0=d, p1=e}

-- Set up grid
grid = StructuredGrid:new{psurface=makePatch{north=de, east=ce, south=abc, west=ad},
			  niv=nx+1, njv=ny+1}

-- Set up block
blks = FBArray:new{grid=grid, initialState=inflow, nib=2, njb=1,
		       bcList={north=OutFlowBC_Simple:new{},
			       east=OutFlowBC_Simple:new{},
			       south=WallBC_WithSlip:new{},
			       west=InFlowBC_Supersonic:new{flowState=inflow}}
}

-- Set simulation parameters
config.interpolation_order = 2
config.print_count = 1
config.flux_calculator = "adaptive"
SteadyStateSolver{
   use_preconditioning = false,
   -- sigma = 1.0e-6, -- presently it's computed internally
   number_pre_steps = 3,
   number_total_steps = 30,
   -- Settings for FGMRES iterative solver
   max_outer_iterations = 60,
   -- number_inner_iterations = 5, -- not needed is preconditioning is false
   max_restarts = 9,
   -- Settings for start-up phase
   number_start_up_steps = 5,
   cfl0 = 5.0,
   eta0 = 0.5,
   tau0 = 0.1,
   sigma0 = 5.0e-6,
   -- Settings for inexact Newton phase
   cfl1 = 100.0,
   sigma1 = 5.0e-6,
   eta1 = 0.01,
   eta_strategy = "constant",
   -- Settings control write-out
   snapshots_count = 1,
   number_total_snapshots = 30,
   write_diagnostics_count = 1
}
   






