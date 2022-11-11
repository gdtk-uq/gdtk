-- A simple 2D Kelvin-Helmholtz instability

config.title = "KH"
config.dimensions = 2

setGasModel('ideal-air.lua')

-- Configure the domain
L = 1
a = Vector3:new{x = 0, y = 0}
b = Vector3:new{x = L, y = 0}
c = Vector3:new{x = 0, y = L}
d = Vector3:new{x = L, y = L}

south = Line:new{p0 = a, p1 = b}
west = Line:new{p0 = a, p1 = c}
north = Line:new{p0 = c, p1 = d}
east = Line:new{p0 = b, p1 = d}

psurf = makePatch{south = south, west = west, north = north, east = east}

N = 256
sgrid = StructuredGrid:new{psurface = psurf, niv = N + 1, njv = N + 1}

-- The flow conditions
U1 = 50; U2 = -50; U_m = (U1 - U2) / 2
rho1 = 1; rho2 = 2; rho_m = (rho1 - rho2) / 2
p_inf = 1e5
Rgas = 287.15; gamma = 1.4
Ls = 0.025

function KH_flow(x, y)
    if y < (L / 4) then
        rho = rho1 - rho_m * math.exp((y - L / 4) / Ls)
        U = U1 - U_m * math.exp((y - L / 4) / Ls)
    elseif (y >= (L / 4) and y < (L / 2)) then
        rho = rho2 + rho_m * math.exp((-y + L / 4) / Ls)
        U = U2 + U_m * math.exp((-y + L / 4) / Ls)
    elseif (y >= (L / 2) and y < (3 * L / 4)) then
        rho = rho2 + rho_m * math.exp(-(3 * L / 4 - y) / Ls)
        U = U2 + U_m * math.exp(-(3 * L / 4 - y) / Ls)
    else
        rho = rho1 - rho_m * math.exp(-(y - 3 * L / 4) / Ls)
        U = U1 - U_m * math.exp(-(y - 3 * L / 4) / Ls)
    end

    T = p_inf / (Rgas * rho)
    vy_perturb = 1 * math.sin(4 * math.pi * x)
    return FlowState:new{T = T, p = p_inf, velx = U, vely = vy_perturb}
end

nib = 4; njb = 4
FBA = FBArray:new{grid = sgrid, fillCondition = KH_flow, bcList = {south = OutFlowBC_SimpleExtrapolate:new{}, north = OutFlowBC_SimpleExtrapolate:new{}}, nib = nib, njb = njb}

for j = 1, njb do
    connectBlocks(FBA.blockArray[1][j], "west", FBA.blockArray[nib][j], "east")
end
for i = 1, nib do
    connectBlocks(FBA.blockArray[i][1], "south", FBA.blockArray[i][njb], "north")
end

identifyBlockConnections()
mpiDistributeBlocks{ntasks = nib * njb, dist = 'load-balance'}

ref_time = L / math.abs(U1 - U2)
config.cfl_value = 0.5
config.viscous = false
config.gasdynamic_update_scheme = 'classic-rk3'
-- config.flux_calculator = 'ldfss0'
config.flux_calculator = 'ausmdv'
-- config.apply_limiter = false
config.extrema_clipping = false
config.max_time = 2.5 * ref_time
config.dt_plot = config.max_time / 25
config.max_step = 1000000
config.new_flow_format = true
config.flow_format = "eilmer4text"

