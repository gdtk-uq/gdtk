----------------------------------------------------------------------------

-- Plane turbine cascade for validation with experimental data.

-- Kiock, R et al. (1986). 'The Transonic Flow Through a Plane Turbine Cascade as Mea-
-- sured in Four European Wind Tunnels'. In: Journal of Engineering for Gas Turbines
-- and Power 108, p. 277.

-- From this paper, experimental results in Fig 6, BS (wind tunnel in Braunschweig, West Germany)

-- Peter Blyton
--     July 2011: Original implementation.
--     October 2011: Viscous simulation for better validation.

-- Jonathan Ho
-- 	October 2015: Ported Code to lua for eilmer4
----------------------------------------------------------------------------

---------------------------------------------------
--Global data
---------------------------------------------------
config.title = "CO2 Simulation for Kiock turbine cascade"
config.viscous = false
config.max_time = 0.05
config.max_step = 800000
config.dt_plot = 0.005
config.dt_init = 1.0e-12
config.include_quality =true
config.cfl_value = 0.5
config.thermo_interpolator = "rhoe"

--------------------------------------------------
--Flow conditions
---------------------------------------------------
nsp, nmodes = setGasModel('kiock-gas-model.lua')
print("GasModel set to CO2. nsp= ", nsp, " nmodes= ", nmodes)
p_tot =  10e6--Pa
T_tot = 450 --K
-- gma = 1.4
-- Rgas = 287.086  --J/kg.K
gma = 1.3
Rgas = 188.9241
a_tot = math.sqrt(gma*Rgas*T_tot)
M_exit = 0.78
T0_T = 1 + (gma-1.0)/2.0 * M_exit * M_exit
p0_p = math.pow(T0_T,gma/(gma-1.0))
print("p0_p=", p0_p, "T0_T=", T0_T)
p_exit = p_tot / p0_p
T_exit = T_tot / T0_T
u_exit = M_exit * a_tot / math.sqrt(T0_T)
print("p_exit=", p_exit, "T_exit=", T_exit, "u_exit=", u_exit)
--solved externally for p_tot 101.3e4, T_tot = 500K
-- p_exit = 704362
-- T_exit = 466.518
-- u_exit = 256.624
--solved externally for p_tot 20MPa, T0 = 500K
-- p_exit = 1.35398e7
-- T_exit = 459.762
-- u_exit = 250.421
--solved externally for p_tot 20MPa, T0 = 600K
-- p_exit = 1.36765e7
-- T_exit = 555.615
-- u_exit = 283.306
--solved externally for p_tot 20MPa, T0 = 475K
-- p_exit = 1.34768e7
-- T_exit = 435.598
-- u_exit = 240.566
--solved externally for p_tot 10MPa, T0 = 450K
-- p_exit = 6.88931e6
-- T_exit = 414.68
-- u_exit = 233.004
--solved externally for p_tot 10MPa, T0 = 350K
-- p_exit = 6.82001e6
-- T_exit = 319.067
-- u_exit = 178.197
--solved externally for p_tot 10MPa, T0 = 340K
-- p_exit = 6.79549e6
-- T_exit = 310.226
-- u_exit = 168.036
--solved externally for p_tot 10MPa, T0 = 335K
-- p_exit = 6.77837e6
-- T_exit = 306.133
-- u_exit = 161.708


initialCond = FlowState:new{p=p_exit, velx=u_exit/1.414, vely=-u_exit/1.414, T=T_exit}
stagCond = FlowState:new{p=p_tot, T = T_tot, velx = 0.0, vely = 0.0}
--Inlet calculated from previous simulations
initialInletCond = FlowState:new{p=985711.955, T = 497.423071481,
                velx = 0.21087*338.9806*math.sqrt(3)/2, vely = 0.21087*338.9806*0.5}

-----------------------------------------------------
-- Geometric parameters
-----------------------------------------------------
STAGGER_ANGLE = (33.3/180.0)*math.pi
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
suct_surface = RotatedAboutZAxisPath:new{suct_surface, angle=-STAGGER_ANGLE}
pres_surface = Spline2:new{filename="press_surface.dat"}
pres_surface = RotatedAboutZAxisPath:new{pres_surface,angle=-STAGGER_ANGLE}

--Chord is line of which blade "sits" on, starting at (0, 0)
chord_vector = Vector3:new{1.0, 0.0}
chord_vector:rotateAboutZAxis(-STAGGER_ANGLE)
chord_normal = Vector3:new{1.0, 0.0}
chord_normal:rotateAboutZAxis(-STAGGER_ANGLE+math.pi/2.0)
pitch_vector = Vector3:new{0.0, PITCH/2.0}

---------------------------------------------------
--Front 1 block on suction surface
----------------------------------------------------
suct_front_surface = SubRangedPath:new{suct_surface, t0 =0.0, t1=suct_div} --Block path is section of suction surface path

--Spline for north path of block
LE_up = suct_front_surface(0)- 0.04*chord_vector
mid_1 = suct_front_surface(0.1) + 0.14*pitch_vector
mid_2 = suct_front_surface(0.34) + 0.1*pitch_vector
mid_3 = suct_front_surface(0.6) + 0.04*chord_normal
TE_up = suct_front_surface(1) + 0.04*chord_normal + 0.01*chord_vector
suct_front_north = Spline:new{points={LE_up, mid_1, mid_2, mid_3, TE_up}}
suct_front1_north = SubRangedPath:new{suct_front_north,t0=0, t1=pres_div2-0.1}
suct_front1_surface = SubRangedPath:new{suct_front_surface,t0=0, t1=pres_div2-0.1}

suct_front1_east = Line:new{suct_front1_surface(1), suct_front1_north(1)}
suct_front_west = Line:new{suct_front1_surface(0), suct_front1_north(0)}

patch = makePatch{suct_front1_north, suct_front1_east, suct_front1_surface, suct_front_west}
cflist = {None, RobertsFunction:new{end0=true, end1=false,beta=clust_surf}, 
          None, RobertsFunction:new{end0=true,end1=false, beta=clust_surf}}

suct_front1_grid = StructuredGrid:new{psurface=patch, cfList=cflist,niv=8*mrf,njv=2*mrf} 
suct_front1_block = SBlock:new{grid=suct_front1_grid,label="suct_front1", fillCondition=initialCond}

suct_front1_block.bcList[south] = WallBC_NoSlip_Adiabatic:new{}
---------------------------------------------------
--Front 2 block on suction surface
----------------------------------------------------

suct_front2_north = SubRangedPath:new{suct_front_north,t0=pres_div2-0.1,t1=1.0}
suct_front2_surface = SubRangedPath:new{suct_front_surface,t0=pres_div2-0.1,t1=1.0}

suct_front2_east = Line:new{suct_front2_surface(1), suct_front2_north(1)}

patch = makePatch{suct_front2_north, suct_front2_east, suct_front2_surface, suct_front1_east}
cflist = {None, RobertsFunction:new{end0=true, end1=false,beta=clust_surf}, 
          None, RobertsFunction:new{end0=true,end1=false, beta=clust_surf}}

suct_front2_grid = StructuredGrid:new{psurface=patch, cfList=cflist,niv=mrf,njv=2*mrf} 
suct_front2_block = SBlock:new{grid=suct_front2_grid,label="suct_front2", fillCondition=initialCond}

suct_front2_block.bcList[south] = WallBC_NoSlip_Adiabatic:new{}
---------------------------------------------------
--Middle block on suction surface
---------------------------------------------------
suct_mid_surface = SubRangedPath:new{suct_surface,t0=suct_div, t1=suct_div2}
TE_out_mid = suct_mid_surface(1) + 0.04*chord_normal

suct_mid_north = Line:new{TE_up, TE_out_mid}
suct_mid_east = Line:new{suct_mid_surface(1), TE_out_mid}

patch = makePatch{suct_mid_north, suct_mid_east, suct_mid_surface, suct_front2_east}
cflist = {None, RobertsFunction:new{end0=true, end1=false,beta=clust_surf},
          None, RobertsFunction:new{end0=true,end1=false, beta=clust_surf}}


suct_mid_grid = StructuredGrid:new{psurface=patch, cfList=cflist,niv=3*mrf,njv=2*mrf}--suct_front_grid.njv} 
suct_mid_block = SBlock:new{grid=suct_mid_grid,label="suct_mid", fillCondition=initialCond}
suct_mid_block.bcList[south] = WallBC_NoSlip_Adiabatic:new{}
---------------------------------------------------
--Rear block on suction surface
---------------------------------------------------
suct_rear_surface = SubRangedPath:new{suct_surface,t0=suct_div2,t1=1.0}

--Spline for north path of block
TE_out = suct_rear_surface(1) + 0.036*chord_vector - 0.02*chord_normal
mid_1 = suct_rear_surface(0.5) + Vector3:new{0.033, 0}
suct_rear_north = Spline:new{points={TE_out_mid, mid_1, TE_out}}

suct_rear_east = Line:new{suct_rear_surface(1), suct_rear_north(1)}

patch = makePatch{suct_rear_north, suct_rear_east, suct_rear_surface, suct_mid_east}
cflist = {None, RobertsFunction:new{end0=true, end1=false,beta=clust_surf}, 
          None, RobertsFunction:new{end0=true,end1=false, beta=clust_surf}}


suct_rear_grid = StructuredGrid:new{psurface=patch, cfList=cflist,niv=mrf,njv=2*mrf}--suct_front_grid.njv} 
suct_rear_block = SBlock:new{grid=suct_rear_grid,label="suct_rear", fillCondition=initialCond}
suct_rear_block.bcList[south] = WallBC_NoSlip_Adiabatic:new{}
---------------------------------------------------
--Front block on pressure surface
---------------------------------------------------
pres_front_surface = SubRangedPath:new{pres_surface,t0=pres_div1,t1=0.0} 
--Spline for west path of block
LE_down = pres_front_surface(0)- 0.04*chord_normal + 0.01*chord_vector
mid_1 = pres_front_surface(0.5) + Vector3:new{-0.04, 0.0}
mid_2 = pres_front_surface(0.5) + Vector3:new{-0.02, -0.05}
pres_front_west = Spline:new{points={LE_down, mid_2, mid_1, LE_up}}

pres_front_south = Line:new{pres_front_west(0), pres_front_surface(0)}

patch = makePatch{ReversedPath:new{suct_front_west}, pres_front_surface, pres_front_south, pres_front_west}
cflist = {RobertsFunction:new{end0=false, end1=true,beta=clust_surf}, None,
          RobertsFunction:new{end0=false,end1=true, beta=clust_surf},None}

pres_front_grid = StructuredGrid:new{psurface=patch, cfList=cflist,niv=2*mrf,njv=2*mrf} --changed niv from referention front block
pres_front_block = SBlock:new{grid=pres_front_grid,label="pres_front", fillCondition=initialCond}
pres_front_block.bcList[east] = WallBC_NoSlip_Adiabatic:new{}
---------------------------------------------------
-- Middle block on pressure surface
---------------------------------------------------
pres_mid_surface = SubRangedPath:new{pres_surface,t0=pres_div1,t1=pres_div2}

--Spline for south path of block
TE_down = pres_mid_surface(1) - 0.04*chord_normal
mid_1 = pres_mid_surface(0.43) - 0.040*chord_normal
mid_2 = pres_mid_surface(0.6) - 0.030*chord_normal
pres_mid_south = Spline:new{points={LE_down, mid_1, TE_down}}

pres_mid_east = Line:new{pres_mid_south(1), pres_mid_surface(1)}

patch = makePatch{pres_mid_surface, pres_mid_east, pres_mid_south, pres_front_south}
cflist = {None, RobertsFunction:new{end0=false, end1=true,beta=clust_surf}, 
          None, RobertsFunction:new{end0=false,end1=true, beta=clust_surf}}

pres_mid_grid = StructuredGrid:new{psurface=patch, cfList=cflist,niv=8*mrf,njv=2*mrf} --njv previously linked to pres_frontg
pres_mid_block = SBlock:new{grid=pres_mid_grid,label="pres_mid", fillCondition=initialCond}
pres_mid_block.bcList[north] = WallBC_NoSlip_Adiabatic:new{}
---------------------------------------------------
--Rear block on pressure surface
---------------------------------------------------
pres_rear_surface = SubRangedPath:new{pres_surface,t0=pres_div2,t1=1.0}
--Spline for south path of block
mid_1 = pres_rear_surface(0.56) - 0.13*pitch_vector
pres_rear_south = Spline:new{points={TE_down, mid_1, TE_out}}

pres_rear_east = ReversedPath:new{suct_rear_east}

patch = makePatch{pres_rear_surface, pres_rear_east, pres_rear_south, pres_mid_east}
cflist = {None, RobertsFunction:new{end0=false, end1=true,beta=clust_surf}, 
          None, RobertsFunction:new{end0=false,end1=true, beta=clust_surf}}


pres_rear_grid = StructuredGrid:new{psurface=patch, cfList=cflist,niv=mrf,njv=2*mrf} --njv=pres_mid_block.njv
pres_rear_block = SBlock:new{grid=pres_rear_grid,label="pres_rear", fillCondition=initialCond}
pres_rear_block.bcList[north] = WallBC_NoSlip_Adiabatic:new{}
---------------------------------------------------
-- Inflow 1 block
---------------------------------------------------
in1_top_left = Vector3:new{-0.8, PITCH/2.0}
in1_top_right = in1_top_left + Vector3:new{0.78, 0}
in1_bottom_left = in1_top_left + Vector3:new{0, -0.2}
in1_north = Line:new{in1_top_left, in1_top_right}
in1_east = Line:new{LE_up, in1_top_right}
in1_south = Line:new{in1_bottom_left, LE_up}
in1_west = Line:new{in1_bottom_left, in1_top_left}

patch = makePatch{in1_north, in1_east, in1_south, in1_west}

in1_grid = StructuredGrid:new{psurface=patch, gridType = "ao",niv=5*mrf,njv=2*mrf} 
in1_block = SBlock:new{grid=in1_grid,label="in1", fillCondition=initialCond}
-- in1_block.bcList[west] = UserDefinedBC:new{fileName="udf-subsonic-kiock.lua"}
in1_block.bcList[west] = InFlowBC_FromStagnation:new{stagCondition=stagCond, 
               direction_type = "uniform", direction_x=math.sqrt(3)/2, direction_y = 0.5, label="inflow-boundary"}
---------------------------------------------------
--Inflow 2 block
---------------------------------------------------
in2_bottom_left = in1_bottom_left + Vector3:new{0, -0.4}
in2_south = Line:new{in2_bottom_left, LE_down}
in2_west = Line:new{in2_bottom_left, in1_bottom_left}

patch = makePatch{in1_south, pres_front_west, in2_south, in2_west}


in2_grid = StructuredGrid:new{psurface=patch, gridType = "ao",niv=5*mrf,njv=2*mrf} --niv = 5*mrf; njv=pres_front_block.njb 
in2_block = SBlock:new{grid=in2_grid,label="in2", fillCondition=initialCond}
-- in2_block.bcList[west] = UserDefinedBC:new{fileName="udf-subsonic-kiock.lua"}
in2_block.bcList[west] = InFlowBC_FromStagnation:new{stagCondition=stagCond, 
               direction_type = "uniform", direction_x=math.sqrt(3)/2, direction_y = 0.5, label="inflow-boundary"}


---------------------------------------------------
--Inflow 3 block
---------------------------------------------------
in3_bottom_left = Vector3:new{-0.8, -PITCH/2.0}
in3_bottom_right = in3_bottom_left + Vector3:new{0.78, 0}
in3_east = Line:new{in3_bottom_right, LE_down}
in3_south = Line:new{in3_bottom_left, in3_bottom_right}
in3_west = Line:new{in3_bottom_left, in2_bottom_left}

patch = makePatch{in2_south, in3_east, in3_south, in3_west}


in3_grid = StructuredGrid:new{psurface=patch, gridType = "ao",niv=5*mrf,njv=2*mrf} --niv=in2_grid.niv,njv=in1_grid.njv
in3_block = SBlock:new{grid=in3_grid,label="in3", fillCondition=initialCond}
-- in3_block.bcList[west] = UserDefinedBC:new{fileName="udf-subsonic-kiock.lua"}
in3_block.bcList[west] = InFlowBC_FromStagnation:new{stagCondition=stagCond, 
               direction_type = "uniform", direction_x=math.sqrt(3)/2, direction_y = 0.5, label="inflow-boundary"}

---------------------------------------------------
--Outflow 1 block
---------------------------------------------------
out1_top_right = suct_surface(1) + Vector3:new{1, 0.6*PITCH}
out1_top_left = out1_top_right + Vector3:new{-0.98, 0}
out1_bottom_right = out1_top_right + Vector3:new{0, -0.2}
out1_north = Line:new{out1_top_left, out1_top_right}
out1_east = Line:new{out1_bottom_right, out1_top_right}
out1_south = Line:new{TE_up, out1_bottom_right}
out1_west = Line:new{TE_up, out1_top_left}

patch = makePatch{out1_north, out1_east, out1_south, out1_west}


out1_grid = StructuredGrid:new{psurface=patch,gridType = "ao",niv=5*mrf,njv=2*mrf}--njv=in1_grid.njv
out1_block = SBlock:new{grid=out1_grid,label="out1", fillCondition=initialCond}
out1_block.bcList[east] = OutFlowBC_FixedP:new{p_outside=p_exit, label = "OUTLET"}


---------------------------------------------------
--Outflow 2 block
---------------------------------------------------
out2_bottom_right = out1_bottom_right + Vector3:new{0, -0.2}
out2_east = Line:new{out2_bottom_right, out1_bottom_right}
out2_south = Line:new{TE_out_mid, out2_bottom_right}
out2_west = ReversedPath:new{suct_mid_north}

patch = makePatch{out1_south, out2_east, out2_south, out2_west}


out2_grid = StructuredGrid:new{psurface=patch,gridType = "ao",niv=5*mrf,njv=3*mrf}--,niv=out1_grid.niv,njv=suct_mid_grid.niv} 
out2_block = SBlock:new{grid=out2_grid,label="out2", fillCondition=initialCond}
out2_block.bcList[east] = OutFlowBC_FixedP:new{p_outside=p_exit, label = "OUTLET"}

---------------------------------------------------
--Outflow 2low block
---------------------------------------------------
out2low_bottom_right = out2_bottom_right + Vector3:new{0, -0.1}
out2low_east = Line:new{out2low_bottom_right, out2_bottom_right}
out2low_south = Line:new{TE_out, out2low_bottom_right}
out2low_west = ReversedPath:new{suct_rear_north}

patch = makePatch{out2_south, out2low_east, out2low_south, out2low_west}


out2low_grid = StructuredGrid:new{psurface=patch,gridType = "ao",niv=5*mrf,njv =mrf}--njv=suct_rear_grid.niv} 
out2low_block = SBlock:new{grid=out2low_grid,label="out2low", fillCondition=initialCond}
out2low_block.bcList[east] = OutFlowBC_FixedP:new{p_outside=p_exit, label = "OUTLET"}

---------------------------------------------------
--Outflow 3 block
---------------------------------------------------
out3_bottom_right = suct_surface(1) + Vector3:new{1, -0.4*PITCH}
out3_bottom_left = out3_bottom_right + Vector3:new{-0.98, 0}
out3_east = Line:new{out3_bottom_right, out2low_bottom_right}
out3_south = Line:new{out3_bottom_left, out3_bottom_right}
out3_west = Line:new{out3_bottom_left, TE_out}

patch = makePatch{out2low_south, out3_east, out3_south, out3_west}


out3_grid = StructuredGrid:new{psurface=patch,gridType = "ao",niv=5*mrf, njv=2*mrf}--out2_grid.niv,njv=in3_grid.njv} 
out3_block = SBlock:new{grid=out3_grid,label="out3", fillCondition=initialCond}
out3_block.bcList[east] = OutFlowBC_FixedP:new{p_outside=p_exit, label = "OUTLET"}


-----------------------------------------------------
-- Top 1 block
-----------------------------------------------------
mid_1 = in1_top_right + Vector3:new{0.44, 0.01}
top_north = Spline:new{points={in1_top_right, mid_1, out1_top_left}}
top1_north = SubRangedPath:new{top_north,t0=0.0,t1=pres_div2-0.1}
top1_east = Line:new{suct_front1_north(1),top1_north(1)}
patch = makePatch{top1_north, top1_east, suct_front1_north, in1_east}

top1_grid = StructuredGrid:new{psurface=patch,gridType = "ao",niv=8*mrf,njv=2*mrf}--niv=suct_front_grid.niv,njv=in1_grid.njv} 
top1_block = SBlock:new{grid=top1_grid,label="top1", fillCondition=initialCond}
-----------------------------------------------------
-- Top 2 block
-----------------------------------------------------
mid_1 = in1_top_right + Vector3:new{0.44, 0.01}
top2_north = SubRangedPath:new{top_north,t0=pres_div2-0.1,t1=1.0}

patch = makePatch{top2_north, out1_west, suct_front2_north, top1_east}


top2_grid = StructuredGrid:new{psurface=patch,gridType = "ao",niv=mrf,njv=2*mrf}--niv=suct_front_grid.niv,njv=in1_grid.njv} 
top2_block = SBlock:new{grid=top2_grid,label="top2", fillCondition=initialCond}
-----------------------------------------------------
-- Lower 1 block
-----------------------------------------------------
low1_south= TranslatedPath:new{top_north,shift=-2.0*pitch_vector}
low1_south = SubRangedPath:new{low1_south,t0=0.0,t1=pres_div2-0.1}
low1_east = Line:new{low1_south(1), TE_down}

patch = makePatch{pres_mid_south, low1_east, low1_south, in3_east}



low1_grid = StructuredGrid:new{psurface=patch,gridType = "ao",niv=8*mrf,njv=2*mrf}--,niv=pres_mid_grid.niv,njv=in3_grid.njv} 
low1_block = SBlock:new{grid=low1_grid,label="low1", fillCondition=initialCond}
-----------------------------------------------------
-- Lower 2 block
-----------------------------------------------------
low2_southtemp = TranslatedPath:new{top_north,shift=-2.0*pitch_vector}
low2_south = SubRangedPath:new{low2_southtemp,t0=pres_div2-0.1,t1=1.0}


patch = makePatch{pres_rear_south, out3_west, low2_south, low1_east}

low2_grid = StructuredGrid:new{psurface=patch,gridType = "ao",niv=mrf,njv=2*mrf}--niv=pres_rear_grid.niv,njv=low1_grid.njv} 
low2_block = SBlock:new{grid=low2_grid,label="low2", fillCondition=initialCond}

identifyBlockConnections()
connectBlocks(in1_block,north,in3_block,south,0)
connectBlocks(out1_block,north,out3_block,south,0)
connectBlocks(top1_block,north,low1_block,south,0)
connectBlocks(top2_block,north,low2_block,south,0)
config.apply_bcs_in_parallel = false