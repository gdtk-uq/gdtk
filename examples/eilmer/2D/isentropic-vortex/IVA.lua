-- Advection of an isentropic vortex test case
-- "A Survey of the Isentropic Euler Vortex Problem using High-Order Methods"
-- Spiegel et al. 22nd AIAA Computational Fluid Dynamics Conference

-- Lachlan Whyborn 01/09/23

dofile("run-config.lua")

config.title = "IVA"
config.dimensions = 2
nsp, nmodes, gm = setGasModel('ideal-air.lua')

-- The initial state
beta = 1 / 5
sigma = 1
R = 0.005

Rgas = 287.15
gamma = 1.4

M_inf = 0.5; p_inf = 1e5; T_inf = 300
a_inf = math.sqrt(gamma * Rgas * T_inf)
rho_inf = p_inf / (Rgas * T_inf)

function initial(x, y)
    f = -(1 / (2 * sigma^2)) * ((x / R)^2 + (y / R)^2)
    omega = beta * math.exp(f)
    
    du = -(y / R) * omega
    dv = (x / R) * omega
    dT = -((gamma - 1) / 2) * omega^2

    T = T_inf * (1 + dT)
    p = p_inf^(-1 / (gamma - 1)) * (rho_inf * Rgas * T)^(gamma / (gamma - 1))
    velx = a_inf * (M_inf + du)
    vely = a_inf * dv

    return FlowState:new{p = p, T = T, velx = velx, vely = vely}
end

L = 0.05
botLeft = Vector3:new{x = -L, y = -L}
botRight = Vector3:new{x = L, y = -L}
topLeft = Vector3:new{x = -L, y = L}
topRight = Vector3:new{x = L, y = L}

grid = StructuredGrid:new{psurface = makePatch{south = Line:new{p0 = botLeft, p1 = botRight}, west = Line:new{p0 = botLeft, p1 = topLeft},
                                               north = Line:new{p0 = topLeft, p1 = topRight}, east = Line:new{p0 = botRight, p1 = topRight}},
                          niv = N + 1, njv = N + 1}

FBA = FBArray:new{grid = grid, bcList = {}, fillCondition = initial, nib = NB / 2, njb = NB / 2}

-- Set up the periodic boundaries
for ib = 1, NB / 2 do
    FBA.blockArray[ib][1].bcList["south"] = ExchangeBC_FullFace:new{otherBlock = FBA.blockArray[ib][NB / 2].id, otherFace = "north", orientation = 0, reorient_vector_quantities = false}
    FBA.blockArray[ib][NB / 2].bcList["north"] = ExchangeBC_FullFace:new{otherBlock = FBA.blockArray[ib][1].id, otherFace = "south", orientation = 0, reorient_vector_quantities = false}
    FBA.blockArray[1][ib].bcList["west"] = ExchangeBC_FullFace:new{otherBlock = FBA.blockArray[NB / 2][ib].id, otherFace = "east", orientation = 0, reorient_vector_quantities = false}
    FBA.blockArray[NB / 2][ib].bcList["east"] = ExchangeBC_FullFace:new{otherBlock = FBA.blockArray[1][ib].id, otherFace = "west", orientation = 0, reorient_vector_quantities = false}
end

if mpi then
    mpiDistributeBlocks{ntasks = NB, dist = "load-balance"}
end
-- Set up the timestepping- needs to be precise
config.cfl_value = 0.25
config.max_step = math.ceil(N * (M_inf + beta * math.exp(-0.5) + 1) / (M_inf * config.cfl_value))
config.max_time = 2 * L / (M_inf * a_inf)
config.dt_init = config.max_time / config.max_step
config.fixed_time_step = true
config.dt_plot = config.max_time / 5
config.gasdynamic_update_scheme = "classic-rk3"
config.flux_calculator = flux_calculator
config.interpolation_order = interpolation_order
config.apply_limiter = false
config.extrema_clipping = false

