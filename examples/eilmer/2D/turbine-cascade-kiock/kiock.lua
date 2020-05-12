-- kiock.lua
-- Plane turbine cascade for validation with experimental data.
--
-- Kiock, R et al. (1986). 'The Transonic Flow Through a Plane Turbine Cascade as Mea-
-- sured in Four European Wind Tunnels'. In: Journal of Engineering for Gas Turbines
-- and Power 108, p. 277.
--
-- From this paper, experimental results in Fig 6, BS
-- (wind tunnel in Braunschweig, West Germany)
--
-- History:
-- 2011 July: Peter Blyton, Original implementation.
-- 2011 October: Viscous simulation for better validation.
-- 2015 October: Jonathan Ho, Ported Code to Lua for Eilmer4.
-- 2020 May, PAJ, refresh and return to original flow spec.
--
---------------------------------------------------
-- Global data
---------------------------------------------------
config.title = "Simulation for Kiock turbine cascade"
config.viscous = true
config.gasdynamic_update_scheme = "euler"
config.max_time = 0.100
config.dt_plot = 0.010
config.dt_init = 1.0e-7
config.max_step = 800000

---------------------------------------------------
-- Flow conditions
---------------------------------------------------
-- We start with a given stagnation condition (atmospheric)
-- and exit Mach number, then use isentropic flow relations
-- to compute an exit pressure that we will apply at the
-- out-flow boundary condition.
nsp, nmodes = setGasModel('ideal-air-gas-model.lua')
p_tot = 100e3 --Pa
T_tot = 300 --K
gma = 1.4
Rgas = 287.086  --J/kg.K
a_tot = math.sqrt(gma*Rgas*T_tot)
M_exit = 0.78
T0_T = 1 + (gma-1.0)/2.0 * M_exit * M_exit
p0_p = math.pow(T0_T,gma/(gma-1.0))
print("p0_p=", p0_p, "T0_T=", T0_T)
p_exit = p_tot / p0_p
T_exit = T_tot / T0_T
vel_exit = M_exit * a_tot / math.sqrt(T0_T)
print("p_exit=", p_exit, "T_exit=", T_exit, "vel_exit=", vel_exit)
initialCond = FlowState:new{p=p_exit, T=T_exit, velx=vel_exit}
stagCond = FlowState:new{p=p_tot, T=T_tot}

-----------------------------------------------------
-- Geometric parameters
-----------------------------------------------------
STAGGER_ANGLE = math.rad(33.3)
PITCH = 0.71

-----------------------------------------------------
-- Mesh setup parameters
-----------------------------------------------------
mrf = 4 -- Mesh refinement factor
suct_div = 0.8 -- Divisions for suction surface blocks
suct_div2 = 0.98
pres_div1 = 0.08 -- Divisions for pressure surface blocks
pres_div2 = 0.965
clust_surf = 1.01 -- clustering normal to surface

-----------------------------------------------------
-- General path and node setup
-----------------------------------------------------
-- Full suction and pressure surface paths
suct_surface = Spline2:new{filename="suct_surface.dat"}
suct_surface = RotatedAboutZAxisPath:new{original_path=suct_surface, angle=-STAGGER_ANGLE}
pres_surface = Spline2:new{filename="press_surface.dat"}
pres_surface = RotatedAboutZAxisPath:new{original_path=pres_surface,angle=-STAGGER_ANGLE}
--Chord is line of which blade "sits" on, starting at (0, 0)
chord_vector = Vector3:new{x=1.0, y=0.0}
chord_vector:rotateAboutZAxis(-STAGGER_ANGLE)
chord_normal = Vector3:new{x=1.0, y=0.0}
chord_normal:rotateAboutZAxis(-STAGGER_ANGLE+math.pi/2.0)
pitch_vector = Vector3:new{x=0.0, y=PITCH/2.0}
---------------------------------------------------
-- Block[0]: Front 1 block on suction surface
----------------------------------------------------
--Block path is section of suction surface path
suct_front_surface = SubRangedPath:new{underlying_path=suct_surface, t0 =0.0, t1=suct_div}
--Spline for north path of block
LE_up = suct_front_surface(0)- 0.04*chord_vector
mid_1 = suct_front_surface(0.1) + 0.14*pitch_vector
mid_2 = suct_front_surface(0.34) + 0.1*pitch_vector
mid_3 = suct_front_surface(0.6) + 0.04*chord_normal
TE_up = suct_front_surface(1) + 0.04*chord_normal + 0.01*chord_vector
suct_front_north = Spline:new{points={LE_up, mid_1, mid_2, mid_3, TE_up}}
suct_front1_north = SubRangedPath:new{underlying_path=suct_front_north,t0=0, t1=pres_div2-0.1}
suct_front1_surface = SubRangedPath:new{underlying_path=suct_front_surface,t0=0, t1=pres_div2-0.1}
suct_front1_east = Line:new{p0=suct_front1_surface(1), p1=suct_front1_north(1)}
suct_front_west = Line:new{p0=suct_front1_surface(0), p1=suct_front1_north(0)}
patch = CoonsPatch:new{north=suct_front1_north, east=suct_front1_east,
                       south=suct_front1_surface, west=suct_front_west}
cflist = {east=RobertsFunction:new{end0=true, end1=false,beta=clust_surf},
          west=RobertsFunction:new{end0=true,end1=false, beta=clust_surf}}
suct_front1_grid = StructuredGrid:new{psurface=patch, cfList=cflist, niv=8*mrf, njv=2*mrf}
suct_front1_block = FluidBlock:new{grid=suct_front1_grid, label="suct_front1", initialState=initialCond}
suct_front1_block.bcList[south] = WallBC_NoSlip_Adiabatic:new{}
---------------------------------------------------
-- Block[1]: Front 2 block on suction surface
----------------------------------------------------
suct_front2_north = SubRangedPath:new{underlying_path=suct_front_north, t0=pres_div2-0.1, t1=1.0}
suct_front2_surface = SubRangedPath:new{underlying_path=suct_front_surface, t0=pres_div2-0.1, t1=1.0}
suct_front2_east = Line:new{p0=suct_front2_surface(1), p1=suct_front2_north(1)}
patch = CoonsPatch:new{north=suct_front2_north, east=suct_front2_east,
                       south=suct_front2_surface, west=suct_front1_east}
cflist = {east=RobertsFunction:new{end0=true, end1=false,beta=clust_surf},
          west=RobertsFunction:new{end0=true,end1=false, beta=clust_surf}}
suct_front2_grid = StructuredGrid:new{psurface=patch, cfList=cflist, niv=mrf, njv=2*mrf}
suct_front2_block = FluidBlock:new{grid=suct_front2_grid, label="suct_front2", initialState=initialCond}
suct_front2_block.bcList[south] = WallBC_NoSlip_Adiabatic:new{}
---------------------------------------------------
-- Block[2]: Middle block on suction surface
---------------------------------------------------
suct_mid_surface = SubRangedPath:new{underlying_path=suct_surface, t0=suct_div, t1=suct_div2}
TE_out_mid = suct_mid_surface(1) + 0.04*chord_normal

suct_mid_north = Line:new{p0=TE_up, p1=TE_out_mid}
suct_mid_east = Line:new{p0=suct_mid_surface(1), p1=TE_out_mid}

patch = CoonsPatch:new{north=suct_mid_north, east=suct_mid_east,
                       south=suct_mid_surface, west=suct_front2_east}
cflist = {east=RobertsFunction:new{end0=true, end1=false, beta=clust_surf},
          west=RobertsFunction:new{end0=true, end1=false, beta=clust_surf}}
suct_mid_grid = StructuredGrid:new{psurface=patch, cfList=cflist, niv=3*mrf, njv=2*mrf}--suct_front_grid.njv}
suct_mid_block = FluidBlock:new{grid=suct_mid_grid,label="suct_mid", initialState=initialCond}
suct_mid_block.bcList[south] = WallBC_NoSlip_Adiabatic:new{}
---------------------------------------------------
-- Block[3]: Rear block on suction surface
---------------------------------------------------
suct_rear_surface = SubRangedPath:new{underlying_path=suct_surface, t0=suct_div2, t1=1.0}

--Spline for north path of block
TE_out = suct_rear_surface(1) + 0.036*chord_vector - 0.02*chord_normal
mid_1 = suct_rear_surface(0.5) + Vector3:new{x=0.033, y=0}
suct_rear_north = Spline:new{points={TE_out_mid, mid_1, TE_out}}
suct_rear_east = Line:new{p0=suct_rear_surface(1), p1=suct_rear_north(1)}

patch = CoonsPatch:new{north=suct_rear_north, east=suct_rear_east,
                       south=suct_rear_surface, west=suct_mid_east}
cflist = {east=RobertsFunction:new{end0=true, end1=false, beta=clust_surf},
          west=RobertsFunction:new{end0=true, end1=false, beta=clust_surf}}
suct_rear_grid = StructuredGrid:new{psurface=patch, cfList=cflist, niv=mrf, njv=2*mrf}--suct_front_grid.njv}
suct_rear_block = FluidBlock:new{grid=suct_rear_grid,label="suct_rear", initialState=initialCond}
suct_rear_block.bcList[south] = WallBC_NoSlip_Adiabatic:new{}
---------------------------------------------------
-- Block[4]: Front block on pressure surface
---------------------------------------------------
pres_front_surface = SubRangedPath:new{underlying_path=pres_surface, t0=pres_div1, t1=0.0}
--Spline for west path of block
LE_down = pres_front_surface(0)- 0.04*chord_normal + 0.01*chord_vector
mid_1 = pres_front_surface(0.5) + Vector3:new{x=-0.04, y=0.0}
mid_2 = pres_front_surface(0.5) + Vector3:new{x=-0.02, y=-0.05}
pres_front_west = Spline:new{points={LE_down, mid_2, mid_1, LE_up}}
pres_front_south = Line:new{p0=pres_front_west(0), p1=pres_front_surface(0)}
patch = CoonsPatch:new{north=ReversedPath:new{underlying_path=suct_front_west},
                       east=pres_front_surface,
                       south=pres_front_south,
                       west=pres_front_west}
cflist = {north=RobertsFunction:new{end0=false, end1=true, beta=clust_surf},
          south=RobertsFunction:new{end0=false, end1=true, beta=clust_surf}}
pres_front_grid = StructuredGrid:new{psurface=patch, cfList=cflist, niv=2*mrf, njv=2*mrf} --changed niv from referention front block
pres_front_block = FluidBlock:new{grid=pres_front_grid, label="pres_front", initialState=initialCond}
pres_front_block.bcList[east] = WallBC_NoSlip_Adiabatic:new{}
---------------------------------------------------
-- Block[5]: Middle block on pressure surface
---------------------------------------------------
pres_mid_surface = SubRangedPath:new{underlying_path=pres_surface,t0=pres_div1,t1=pres_div2}
--Spline for south path of block
TE_down = pres_mid_surface(1) - 0.04*chord_normal
mid_1 = pres_mid_surface(0.43) - 0.040*chord_normal
mid_2 = pres_mid_surface(0.6) - 0.030*chord_normal
pres_mid_south = Spline:new{points={LE_down, mid_1, TE_down}}
pres_mid_east = Line:new{p0=pres_mid_south(1), p1=pres_mid_surface(1)}
patch = CoonsPatch:new{north=pres_mid_surface, east=pres_mid_east,
                       south=pres_mid_south, west=pres_front_south}
cflist = {east=RobertsFunction:new{end0=false, end1=true,beta=clust_surf},
          west=RobertsFunction:new{end0=false,end1=true, beta=clust_surf}}
pres_mid_grid = StructuredGrid:new{psurface=patch, cfList=cflist, niv=8*mrf, njv=2*mrf} --njv previously linked to pres_frontg
pres_mid_block = FluidBlock:new{grid=pres_mid_grid, label="pres_mid", initialState=initialCond}
pres_mid_block.bcList[north] = WallBC_NoSlip_Adiabatic:new{}
---------------------------------------------------
-- Block[6]: Rear block on pressure surface
---------------------------------------------------
pres_rear_surface = SubRangedPath:new{underlying_path=pres_surface, t0=pres_div2, t1=1.0}
--Spline for south path of block
mid_1 = pres_rear_surface(0.56) - 0.13*pitch_vector
pres_rear_south = Spline:new{points={TE_down, mid_1, TE_out}}

pres_rear_east = ReversedPath:new{underlying_path=suct_rear_east}

patch = CoonsPatch:new{north=pres_rear_surface, east=pres_rear_east,
                       south=pres_rear_south, west=pres_mid_east}
cflist = {east=RobertsFunction:new{end0=false, end1=true,beta=clust_surf},
          west=RobertsFunction:new{end0=false,end1=true, beta=clust_surf}}
pres_rear_grid = StructuredGrid:new{psurface=patch, cfList=cflist, niv=mrf, njv=2*mrf} --njv=pres_mid_block.njv
pres_rear_block = FluidBlock:new{grid=pres_rear_grid, label="pres_rear", initialState=initialCond}
pres_rear_block.bcList[north] = WallBC_NoSlip_Adiabatic:new{}
---------------------------------------------------
-- Inflow 1 block
---------------------------------------------------
in1_top_left = Vector3:new{x=-0.8, y=PITCH/2.0}
in1_top_right = in1_top_left + Vector3:new{x=0.78, y=0}
in1_bottom_left = in1_top_left + Vector3:new{x=0, y=-0.2}
in1_north = Line:new{p0=in1_top_left, p1=in1_top_right}
in1_east = Line:new{p0=LE_up, p1=in1_top_right}
in1_south = Line:new{p0=in1_bottom_left, p1=LE_up}
in1_west = Line:new{p0=in1_bottom_left, p1=in1_top_left}
patch = AOPatch:new{north=in1_north, east=in1_east, south=in1_south, west=in1_west}
in1_grid = StructuredGrid:new{psurface=patch, niv=5*mrf, njv=2*mrf}
in1_block = FluidBlock:new{grid=in1_grid, label="in1", initialState=initialCond}
in1_block.bcList[west] = InFlowBC_FromStagnation:new{stagnationState=stagCond,
						     direction_type = "uniform",
						     direction_x=math.sqrt(3)/2,
						     direction_y = 0.5,
						     label="inflow-boundary"}
---------------------------------------------------
--Inflow 2 block
---------------------------------------------------
in2_bottom_left = in1_bottom_left + Vector3:new{x=0, y=-0.4}
in2_south = Line:new{p0=in2_bottom_left, p1=LE_down}
in2_west = Line:new{p0=in2_bottom_left, p1=in1_bottom_left}
patch = AOPatch:new{north=in1_south, east=pres_front_west, south=in2_south, west=in2_west}
in2_grid = StructuredGrid:new{psurface=patch, niv=5*mrf, njv=2*mrf} --niv = 5*mrf; njv=pres_front_block.njb
in2_block = FluidBlock:new{grid=in2_grid, label="in2", initialState=initialCond}
in2_block.bcList[west] = InFlowBC_FromStagnation:new{stagnationState=stagCond,
						     direction_type = "uniform",
						     direction_x=math.sqrt(3)/2,
						     direction_y = 0.5,
						     label="inflow-boundary"}
---------------------------------------------------
--Inflow 3 block
---------------------------------------------------
in3_bottom_left = Vector3:new{x=-0.8, y=-PITCH/2.0}
in3_bottom_right = in3_bottom_left + Vector3:new{x=0.78, y=0}
in3_east = Line:new{p0=in3_bottom_right, p1=LE_down}
in3_south = Line:new{p0=in3_bottom_left, p1=in3_bottom_right}
in3_west = Line:new{p0=in3_bottom_left, p1=in2_bottom_left}
patch = AOPatch:new{north=in2_south, east=in3_east, south=in3_south, west=in3_west}
in3_grid = StructuredGrid:new{psurface=patch, niv=5*mrf, njv=2*mrf} --niv=in2_grid.niv,njv=in1_grid.njv
in3_block = FluidBlock:new{grid=in3_grid, label="in3", initialState=initialCond}
in3_block.bcList[west] = InFlowBC_FromStagnation:new{stagnationState=stagCond,
						     direction_type = "uniform",
						     direction_x=math.sqrt(3)/2,
						     direction_y = 0.5,
						     label="inflow-boundary"}
---------------------------------------------------
--Outflow 1 block
---------------------------------------------------
out1_top_right = suct_surface(1) + Vector3:new{x=1, y=0.6*PITCH}
out1_top_left = out1_top_right + Vector3:new{x=-0.98, y=0}
out1_bottom_right = out1_top_right + Vector3:new{x=0, y=-0.2}
out1_north = Line:new{p0=out1_top_left, p1=out1_top_right}
out1_east = Line:new{p0=out1_bottom_right, p1=out1_top_right}
out1_south = Line:new{p0=TE_up, p1=out1_bottom_right}
out1_west = Line:new{p0=TE_up, p1=out1_top_left}
patch = AOPatch:new{north=out1_north, east=out1_east, south=out1_south, west=out1_west}
out1_grid = StructuredGrid:new{psurface=patch, niv=5*mrf, njv=2*mrf} --njv=in1_grid.njv
out1_block = FluidBlock:new{grid=out1_grid, label="out1", initialState=initialCond}
out1_block.bcList[east] = OutFlowBC_FixedP:new{p_outside=p_exit, label="OUTLET"}
---------------------------------------------------
--Outflow 2 block
---------------------------------------------------
out2_bottom_right = out1_bottom_right + Vector3:new{x=0, y=-0.2}
out2_east = Line:new{p0=out2_bottom_right, p1=out1_bottom_right}
out2_south = Line:new{p0=TE_out_mid, p1=out2_bottom_right}
out2_west = ReversedPath:new{underlying_path=suct_mid_north}
patch = AOPatch:new{north=out1_south, east=out2_east, south=out2_south, west=out2_west}
out2_grid = StructuredGrid:new{psurface=patch, niv=5*mrf, njv=3*mrf} --,niv=out1_grid.niv,njv=suct_mid_grid.niv}
out2_block = FluidBlock:new{grid=out2_grid, label="out2", initialState=initialCond}
out2_block.bcList[east] = OutFlowBC_FixedP:new{p_outside=p_exit, label="OUTLET"}
---------------------------------------------------
--Outflow 2low block
---------------------------------------------------
out2low_bottom_right = out2_bottom_right + Vector3:new{x=0, y=-0.1}
out2low_east = Line:new{p0=out2low_bottom_right, p1=out2_bottom_right}
out2low_south = Line:new{p0=TE_out, p1=out2low_bottom_right}
out2low_west = ReversedPath:new{underlying_path=suct_rear_north}
patch = AOPatch:new{north=out2_south, east=out2low_east, south=out2low_south, west=out2low_west}
out2low_grid = StructuredGrid:new{psurface=patch, niv=5*mrf, njv =mrf} --njv=suct_rear_grid.niv}
out2low_block = FluidBlock:new{grid=out2low_grid, label="out2low", initialState=initialCond}
out2low_block.bcList[east] = OutFlowBC_FixedP:new{p_outside=p_exit, label="OUTLET"}

---------------------------------------------------
--Outflow 3 block
---------------------------------------------------
out3_bottom_right = suct_surface(1) + Vector3:new{x=1, y=-0.4*PITCH}
out3_bottom_left = out3_bottom_right + Vector3:new{x=-0.98, y=0}
out3_east = Line:new{p0=out3_bottom_right, p1=out2low_bottom_right}
out3_south = Line:new{p0=out3_bottom_left, p1=out3_bottom_right}
out3_west = Line:new{p0=out3_bottom_left, p1=TE_out}
patch = AOPatch:new{north=out2low_south, east=out3_east, south=out3_south, west=out3_west}
out3_grid = StructuredGrid:new{psurface=patch, niv=5*mrf, njv=2*mrf} --out2_grid.niv,njv=in3_grid.njv}
out3_block = FluidBlock:new{grid=out3_grid, label="out3", initialState=initialCond}
out3_block.bcList[east] = OutFlowBC_FixedP:new{p_outside=p_exit, label="OUTLET"}
-----------------------------------------------------
-- Top 1 block
-----------------------------------------------------
mid_1 = in1_top_right + Vector3:new{x=0.44, y=0.01}
top_north = Spline:new{points={in1_top_right, mid_1, out1_top_left}}
top1_north = SubRangedPath:new{underlying_path=top_north,t0=0.0,t1=pres_div2-0.1}
top1_east = Line:new{p0=suct_front1_north(1), p1=top1_north(1)}
patch = AOPatch:new{north=top1_north, east=top1_east, south=suct_front1_north, west=in1_east}
top1_grid = StructuredGrid:new{psurface=patch, niv=8*mrf, njv=2*mrf} --niv=suct_front_grid.niv,njv=in1_grid.njv}
top1_block = FluidBlock:new{grid=top1_grid, label="top1", initialState=initialCond}
-----------------------------------------------------
-- Top 2 block
-----------------------------------------------------
mid_1 = in1_top_right + Vector3:new{x=0.44, y=0.01}
top2_north = SubRangedPath:new{underlying_path=top_north, t0=pres_div2-0.1, t1=1.0}
patch = AOPatch:new{north=top2_north, east=out1_west, south=suct_front2_north, west=top1_east}
top2_grid = StructuredGrid:new{psurface=patch, niv=mrf, njv=2*mrf} --niv=suct_front_grid.niv,njv=in1_grid.njv}
top2_block = FluidBlock:new{grid=top2_grid, label="top2", initialState=initialCond}
-----------------------------------------------------
-- Lower 1 block
-----------------------------------------------------
low1_south= TranslatedPath:new{original_path=top_north, shift=-2.0*pitch_vector}
low1_south = SubRangedPath:new{underlying_path=low1_south, t0=0.0, t1=pres_div2-0.1}
low1_east = Line:new{p0=low1_south(1), p1=TE_down}
patch = AOPatch:new{north=pres_mid_south, east=low1_east, south=low1_south, west=in3_east}
low1_grid = StructuredGrid:new{psurface=patch, niv=8*mrf, njv=2*mrf} --,niv=pres_mid_grid.niv,njv=in3_grid.njv}
low1_block = FluidBlock:new{grid=low1_grid, label="low1", initialState=initialCond}
-----------------------------------------------------
-- Lower 2 block
-----------------------------------------------------
low2_southtemp = TranslatedPath:new{original_path=top_north, shift=-2.0*pitch_vector}
low2_south = SubRangedPath:new{underlying_path=low2_southtemp, t0=pres_div2-0.1, t1=1.0}
patch = AOPatch:new{north=pres_rear_south, east=out3_west, south=low2_south, west=low1_east}
low2_grid = StructuredGrid:new{psurface=patch, niv=mrf, njv=2*mrf} --niv=pres_rear_grid.niv,njv=low1_grid.njv}
low2_block = FluidBlock:new{grid=low2_grid, label="low2", initialState=initialCond}

identifyBlockConnections()
-- This flow domain holds one blade that is part of a larger cascade.
-- Apply periodic conditions, north and south.
connectBlocks(in1_block,north, in3_block,south, 0)
connectBlocks(out1_block,north, out3_block,south, 0)
connectBlocks(top1_block,north, low1_block,south, 0)
connectBlocks(top2_block,north, low2_block,south, 0)
--
mpiDistributeBlocks{ntasks=3, dist="load-balance"}
