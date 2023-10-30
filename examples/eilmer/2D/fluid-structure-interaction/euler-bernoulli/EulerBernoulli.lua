-- Demo of the fluid structure interaction module
-- Lachlan Whyborn 30/10/23
config.dimensions = 2
config.title = "EulerBernoulli"

L = 1.0; t = 0.01; h = 0.25; w = 1.0

setGasModel("ideal-air.lua")
u = 1500; pn = 1e5; ps = 1e4; T = 300
fn = FlowState:new{velx = u, T = T, p = pn}
fs = FlowState:new{velx = u, T = T, p = ps}

-- Construct the mesh
PlateNorthLeadingEdge = Vector3:new{x = 0.0, y = t / 2}
PlateSouthLeadingEdge = Vector3:new{x = 0.0, y = -t / 2}
PlateNorthTrailingEdge = Vector3:new{x = L, y = t / 2}
PlateSouthTrailingEdge = Vector3:new{x = L, y = -t / 2}
PlateNorthOutflow = Vector3:new{x = 2 * L, y = t / 2}
PlateSouthOutflow = Vector3:new{x = 2 * L, y = -t / 2}

NorthInflowCorner = Vector3:new{x = 0.0, y = t / 2 + h}
SouthInflowCorner = Vector3:new{x = 0.0, y = -t / 2 - h}
NorthMidCorner = Vector3:new{x = L, y = t / 2 + h}
SouthMidCorner = Vector3:new{x = L, y = -t / 2 - h}
NorthOutflowCorner = Vector3:new{x = 2 * L, y = t / 2 + h}
SouthOutflowCorner = Vector3:new{x = 2 * L, y = -t / 2 - h}

PlateNorth = Line:new{p0 = PlateNorthLeadingEdge, p1 = PlateNorthTrailingEdge}
PlateSouth = Line:new{p0 = PlateSouthLeadingEdge, p1 = PlateSouthTrailingEdge}
InflowNorth = Line:new{p0 = PlateNorthLeadingEdge, p1 = NorthInflowCorner}
InflowSouth = Line:new{p0 = SouthInflowCorner, p1 = PlateSouthLeadingEdge}
DomainEdgePlateNorth = Line:new{p0 = NorthInflowCorner, p1 = NorthMidCorner}
DomainEdgePlateSouth = Line:new{p0 = SouthInflowCorner, p1 = SouthMidCorner}
NorthConnection = Line:new{p0 = PlateNorthTrailingEdge, p1 = NorthMidCorner}
PlateTrailingEdge = Line:new{p0 = PlateSouthTrailingEdge, p1 = PlateNorthTrailingEdge}
SouthConnection = Line:new{p0 = SouthMidCorner, p1 = PlateSouthTrailingEdge}
DownstreamNorth = Line:new{p0 = PlateNorthTrailingEdge, p1 = PlateNorthOutflow}
DownstreamSouth = Line:new{p0 = PlateSouthTrailingEdge, p1 = PlateSouthOutflow}
DomainEdgeDownstreamNorth = Line:new{p0 = NorthMidCorner, p1 = NorthOutflowCorner}
DomainEdgeDownstreamSouth = Line:new{p0 = SouthMidCorner, p1 = SouthOutflowCorner}
OutflowNorth = Line:new{p0 = PlateNorthOutflow, p1 = NorthOutflowCorner}
OutflowPlate = Line:new{p0 = PlateSouthOutflow, p1 = PlateNorthOutflow}
OutflowSouth = Line:new{p0 = SouthOutflowCorner, p1 = PlateSouthOutflow}

PlateNorthSurf = makePatch{south = PlateNorth, west = InflowNorth, north = DomainEdgePlateNorth, east = NorthConnection}
PlateSouthSurf = makePatch{south = DomainEdgePlateSouth, west = InflowSouth, north = PlateSouth, east = SouthConnection}
DownstreamNorthSurf = makePatch{south = DownstreamNorth, west = NorthConnection, north = DomainEdgeDownstreamNorth, east = OutflowNorth}
DownstreamPlateSurf = makePatch{south = DownstreamSouth, west = PlateTrailingEdge, north = DownstreamNorth, east = OutflowPlate}
DownstreamSouthSurf = makePatch{south = DomainEdgeDownstreamSouth, west = SouthConnection, north = DownstreamSouth, east = OutflowSouth}

N = 10
PlateNorthGrid = StructuredGrid:new{psurface = PlateNorthSurf, niv = N+1, njv = N+1}
PlateSouthGrid = StructuredGrid:new{psurface = PlateSouthSurf, niv = N+1, njv = N+1}
DownstreamNorthGrid = StructuredGrid:new{psurface = DownstreamNorthSurf, niv = N+1, njv = N+1}
DownstreamPlateGrid = StructuredGrid:new{psurface = DownstreamPlateSurf, niv = N+1, njv = 5}
DownstreamSouthGrid = StructuredGrid:new{psurface = DownstreamSouthSurf, niv = N+1, njv = N+1}

PlateNorthFBA = FBArray:new{grid = PlateNorthGrid, fillCondition = fn, bcList = {west = InFlowBC_Supersonic:new{flowState = fn}}, nib = 2, njb = 2}
PlateSouthFBA = FBArray:new{grid = PlateSouthGrid, fillCondition = fs, bcList = {west = InFlowBC_Supersonic:new{flowState = fs}}, nib = 2, njb = 2}
DownstreamNorthFBA = FBArray:new{grid = DownstreamNorthGrid, fillCondition = fn, bcList = {east = OutFlowBC_Simple:new{}}, nib = 2, njb = 2}
DownstreamPlateFBA = FBArray:new{grid = DownstreamPlateGrid, fillCondition = fs, bcList = {east = OutFlowBC_Simple:new{}}, nib = 2, njb = 2}
DownstreamSouthFBA = FBArray:new{grid = DownstreamSouthGrid, fillCondition = fs, bcList = {east = OutFlowBC_Simple:new{}}, nib = 2, njb = 2}

identifyBlockConnections()

config.fixed_time_step = true
config.dt_init = 1e-8
config.print_count = 1000
config.max_step = 100000
config.grid_motion = "FSI"          -- set config.grid_motion
config.FEMModel = "eulerBernoulli"  -- set config.FEMModel
config.gasdynamic_update_scheme = "moving_grid_1_stage"     -- set a temporal scheme with grid motion
config.dt_plot = 5e-5
config.max_time = 1e-3
config.report_invalid_cells = false
config.max_attempts_for_step = 1
-- It is important to either a) set config.fixed_time_step = true or b) set config.cfl_count to be a multiple of FSIOptions.couplingStep.
-- Otherwise, the fluid time step may change in the middle of a structure update and de-sync the systems.

-- This table controls the configuration of the FSI
FSIOptions{northForcing = "Fluid",          -- Type of external forcing, currently only fluid pressures.
           southForcing = "Fluid",          -- Ideally support for Lua defined forcing in future.
           Nx = 20,                         -- Number of elements along the beam
           northFBA = PlateNorthFBA,                -- These define the mapping of the FBArrays (note: must be FBArrays, not FluidBlocks) to allow
           southFBA = PlateSouthFBA,                -- the motion of the grid around the deforming beam. See the map below for a visual description
           northEastFBA = DownstreamNorthFBA,       -- of this map. The other possible entries not defined here are northWestFBA, southWestFBA,
           eastAdjacentFBA = DownstreamPlateFBA,    -- westAdjacentFBA.
           southEastFBA = DownstreamSouthFBA,
           plateNormal = {0.0, 1.0, 0.0},           -- Orientation of the plate
           length = L,                              -- Dimensions of the plate
           width = w,
           thickness = t,
           density = 8e3,                           -- Material properties
           youngsModulus = 20e9,
           poissonsRatio = 0.3,
           BCs = "CF",                              -- Boundary conditions, west end then east end. "C" denotes clamped (0 displacement and slope),
                                                    -- "P" denotes pinned (0 displacement and non-zero slope), "F" denotes free (non-zero displacement
                                                    -- and slope)
           couplingStep = 100                       -- Every how many fluid steps do we update the beam dynamics
       }
