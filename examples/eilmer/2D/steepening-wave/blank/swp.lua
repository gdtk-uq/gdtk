-- The steepening wave problem from Lifschitz and Landau in Fluid Mechanics, Vol 2 pgs 378-382
-- @author: Lachlan Whyborn and Nick Gibbons with help from Daryl Bond (Dec 2023)

config.title = "Steepening Wave Problem"

config.dimensions = 2
config.axisymmetric = false
setGasModel('ideal-air.lua')

N=128
config.apply_limiter = false
config.flux_calculator = 'asf'
config.extrema_clipping = false
config.gasdynamic_update_scheme = "classic_rk3"

-- Set the flow conditions- a sine wave.
-- There will be left and right moving waves, with velocity u=u_inf +- c
-- By setting u_inf = -c, we can fix the wave in a set spatial domain
Runi = 8.31446261815324
mMass = 0.02896000
Rgas = Runi / mMass
gamma = 1.4
p_inf = 1e5
rho_inf = 1.0
T_inf = p_inf / (Rgas * rho_inf)
c0 = math.sqrt(gamma * Rgas * T_inf)
u_inf = -1.0*c0

function initial(x, y)
    -- Define using the conditions from Bond et al (2017)
    u = u_inf * math.sin(math.pi * x)
    term = (1. + (gamma - 1) * u / (2. * c0))^(2. / (gamma - 1))
    rho = rho_inf * term
    p = p_inf * term^gamma
    T = p / (Rgas * rho)
    vely = 0
    return FlowState:new{T = T, p = p, velx = u, vely = 0}
end

-- Define the domain
nProcesses = 4
blkArray = {}
xmin = 0.0; xmax = 2.0
for i = 1, nProcesses do

    a = Vector3:new{x = xmin + (i - 1) * (xmax - xmin) / nProcesses, y = -0.01}
    b = Vector3:new{x = xmin + i * (xmax - xmin) / nProcesses, y = -0.01}
    c = Vector3:new{x = a.x, y = 0.01}
    d = Vector3:new{x = b.x, y = 0.01}

    south0 = Line:new{p0 = a, p1 = b}
    west0 = Line:new{p0 = a, p1 = c}
    north0 = Line:new{p0 = c, p1 = d}
    east0 = Line:new{p0 = b, p1 = d}

    psurf0 = makePatch{north = north0, east = east0, south = south0, west = west0}

    grid0 = StructuredGrid:new{psurface = psurf0, niv = N / nProcesses + 1, njv = 4}

    blkArray[i] = FluidBlock:new{grid = grid0, initialState = initial, bcList = {}}
    connectBlocks(blkArray[i], south, blkArray[i], north, 0)
end

connectBlocks(blkArray[1], west, blkArray[#blkArray], east, 0)
identifyBlockConnections()

-- Number of timesteps required to achieve a given CFL for a fixed time step
cfl = 0.1
tmax = 0.50 * math.abs(2. / (u_inf * math.pi * (gamma + 1)))
n = tmax * N * (math.abs(u_inf) + c0) / (cfl * 4)
n = math.ceil(n)

-- Set the configuration
mpiDistributeBlocks{ntasks = nProcesses, dist = "load-balance", preassign={[0]=1}}
config.max_time = tmax
config.max_step = n
config.dt_plot = tmax / 5
config.dt_init = tmax / n
config.fixed_time_step = true
