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
Q.T = {T_inf}
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
blk = SBlock:new{grid=grid, fillCondition=inflow, label="block-0"}
blk.bcList[north] = OutFlowBC_Simple:new{}
blk.bcList[east] = OutFlowBC_Simple:new{}
blk.bcList[south] = WallBC_WithSlip:new{}
blk.bcList[west] = InFlowBC_Supersonic:new{flowCondition=inflow}

-- Set simulation parameters
config.interpolation_order = 2
config.print_count = 1
SteadyStateSolver{
   cfl_init = 10.0,
   eta = 0.01,
   sigma = 1.0e-8,
   tau = 0.1,
   no_low_order_iterations = 5,
   no_outer_iterations = 40,
   max_inner_iterations = 30,
   snapshots_frequency = 1,
   snapshots_count = 40
}
   






