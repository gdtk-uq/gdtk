-- Sliding_Mesh_test.lua
-- Script to simulate and test sliding meshes

-- Author: Ingo Jahn
-- Last modified: 25/04/2017 



--####################################
--### Setting up Basic Information ###
--####################################
config.title = "SLIDING MESH TEST CASE"
print(config.title)
config.dimensions = 3
config.axisymmetric = false

-- The gas model is defined via a gas-model file.
nsp, nmodes, gm = setGasModel('ideal-air-gas-model.lua')
print("GasModel set to ideal air. nsp= ", nsp, " nmodes= ", nmodes)


-- Do a little more setting of global data.
config.udf_supervisor_file="udf-process_e4.lua"
config.max_time = 40.e-4  -- seconds
config.max_step = 100000
config.dt_init = 1.0e-8
config.cfl_value = 0.5
-- config.dt_max = 10.0e-6
config.dt_plot = 5.e-5
config.dt_history = 10.0e-5

-- Call configuration script to do pre-processing for rotating simulation
dofile('udf-config_e4.lua')

-- Use following if loops to set initial and boundary conditions
if false then
    initial1 = FlowState:new{p=1e3, T=800.0, velx=0.0, vely=0.0}
    initial2 = FlowState:new{p=1e3, T=800.0, velx=0.0, vely=0.0}
    initial3 = FlowState:new{p=1e3, T=800.0, velx=0.0, vely=0.0}
    stagnation1 = FlowState:new{p=1e3, T=800.0, velx=0.0, vely=0.0}
    stagnation2 = FlowState:new{p=1e3, T=800.0, velx=0.0, vely=0.0}
    Pout = 1e3
end


if false then
    initial1 = FlowState:new{p=1e3, T=100.0, velx=0.0, vely=0.0}
    initial2 = FlowState:new{p=1e3, T=100.0, velx=0.0, vely=0.0}
    initial3 = FlowState:new{p=1e3, T=100.0, velx=0.0, vely=0.0}
    stagnation1 = FlowState:new{p=1e5, T=800.0, velx=0.0, vely=0.0}
    stagnation2 = FlowState:new{p=1e5, T=800.0, velx=0.0, vely=0.0}
    Pout = 1e3
end

if false then
    initial1 = FlowState:new{p=0.99e5, T=100.0, velx=0.0, vely=0.0}
    initial2 = FlowState:new{p=0.99e5, T=100.0, velx=0.0, vely=0.0}
    initial3 = FlowState:new{p=0.99e5, T=100.0, velx=0.0, vely=0.0}
    stagnation1 = FlowState:new{p=1e5, T=833.0, velx=0.0, vely=0.0}
    stagnation2 = FlowState:new{p=1e5, T=400.0, velx=0.0, vely=0.0}
    Pout = 0.99e5
end

if false then
    initial1 = FlowState:new{p=1e5, T=100.0, velx=0.0, vely=0.0}
    initial2 = FlowState:new{p=1e5, T=100.0, velx=0.0, vely=0.0}
    initial3 = FlowState:new{p=5e5, T=100.0, velx=0.0, vely=0.0}
    stagnation1 = FlowState:new{p=1.1e5, T=500.0, velx=0.0, vely=0.0}
    stagnation2 = FlowState:new{p=1.0e5, T=500.0, velx=0.0, vely=0.0}
    Pout = 0.99e5
end

if true then
    initial1 = FlowState:new{p=1e2, T=100.0, velx=0.0, vely=0.0}
    initial2 = FlowState:new{p=1e2, T=100.0, velx=0.0, vely=0.0}
    initial3 = FlowState:new{p=1e2, T=100.0, velx=0.0, vely=0.0}
    stagnation1 = FlowState:new{p=5e5, T=500.0, velx=-1000.0, vely=0.0}
    stagnation2 = FlowState:new{p=5e5, T=1000.0, velx=-1000.0, vely=0.0}
    Pout = 0.99e5
end


-- Create Stator and Rotor Mesh 
R1 = 0.05
R2 = 0.045 
R3 = 0.04
theta = 2.* math.pi / 16  -- 16 blades
theta_a = -0.5*theta 
theta_b = 0.0*theta
theta_c = 0.5*theta


C = Vector3:new{x=0.,y=0.,z=0.,label="C"}
R1a = Vector3:new{x=R1*math.cos(theta_a),y=R1*math.sin(theta_a),z=0.,label="R1a"}
R1b = Vector3:new{x=R1*math.cos(theta_b),y=R1*math.sin(theta_b),z=0.,label="R1b"}
R1c = Vector3:new{x=R1*math.cos(theta_c),y=R1*math.sin(theta_c),z=0.,label="R1c"}
R2a = Vector3:new{x=R2*math.cos(theta_a),y=R2*math.sin(theta_a),z=0.,label="R2a"}
R2b = Vector3:new{x=R2*math.cos(theta_b),y=R2*math.sin(theta_b),z=0.,label="R2b"}
R2c = Vector3:new{x=R2*math.cos(theta_c),y=R2*math.sin(theta_c),z=0.,label="R2c"}
R3a = Vector3:new{x=R3*math.cos(theta_a),y=R3*math.sin(theta_a),z=0.,label="R3a"}
R3b = Vector3:new{x=R3*math.cos(theta_b),y=R3*math.sin(theta_b),z=0.,label="R3b"}
R3c = Vector3:new{x=R3*math.cos(theta_c),y=R3*math.sin(theta_c),z=0.,label="R3c"}

R1_1 = Arc:new{p0=R1a,p1=R1b,centre=C}; R1_2 = Arc:new{p0=R1b,p1=R1c,centre=C}
R2_1 = Arc:new{p0=R2a,p1=R2b,centre=C}; R2_2 = Arc:new{p0=R2b,p1=R2c,centre=C}
R3_1 = Arc:new{p0=R3a,p1=R3b,centre=C}; R3_2 = Arc:new{p0=R3b,p1=R3c,centre=C}
R2R1a = Line:new{p0=R2a,p1=R1a}; R2R1b = Line:new{p0=R2b,p1=R1b}; R2R1c = Line:new{p0=R2c,p1=R1c}
R2R3a = Line:new{p0=R2a,p1=R3a}; R2R3b = Line:new{p0=R2b,p1=R3b}; R2R3c = Line:new{p0=R2c,p1=R3c}

H_1 = 0.01
H_2 = 0.02
H1 = Line:new{p0=R2a, p1=Vector3:new{x=R2a.x,y=R2a.y,z=R2a.z+H_2} }
H2 = Line:new{p0=R2b, p1=Vector3:new{x=R2b.x,y=R2b.y,z=R2b.z+H_2} }
H3 = Line:new{p0=R2a, p1=Vector3:new{x=R2a.x,y=R2a.y,z=R2a.z+H_1} }
H4 = Line:new{p0=Vector3:new{x=R2a.x,y=R2a.y,z=R2a.z+H_1}, p1=Vector3:new{x=R2a.x,y=R2a.y,z=R2a.z+H_2} }

vol0 = SweptSurfaceVolume:new{face0123=makePatch{north=R2R1b, east=R1_1, south=R2R1a, west=R2_1}, edge04=H1}
vol1 = SweptSurfaceVolume:new{face0123=makePatch{north=R2R1c, east=R1_2, south=R2R1b, west=R2_2}, edge04=H2}

vol2 = SweptSurfaceVolume:new{face0123=makePatch{north=Polyline:new{segments={R3_1, R3_2}}, east=R2R3c, south=Polyline:new{segments={R2_1, R2_2}}, west=R2R3a}, edge04=H3}

vol3 = SweptSurfaceVolume:new{face0123=makePatch{north=Polyline:new{segments={R3_1, R3_2}}, east=R2R3c, south=Polyline:new{segments={R2_1, R2_2}}, west=R2R3a}, edge04=H4}

-- Define number of cells
n1 = 20
n2 = 12
n3 = 6
n4 = 5
n5 = 4

blk0 = FluidBlock:new{grid=StructuredGrid:new{pvolume=vol0, niv=n1, njv=n2, nkv=n3}, initialState=initial1, label="block-0"}
blk1 = FluidBlock:new{grid=StructuredGrid:new{pvolume=vol1, niv=n1, njv=n2, nkv=n3}, initialState=initial2, label="block-1"}

blk2 = FluidBlock:new{grid=StructuredGrid:new{pvolume=vol2, niv=n1, njv=n1, nkv=n4}, initialState=initial3, label="block-2"}
blk3 = FluidBlock:new{grid=StructuredGrid:new{pvolume=vol3, niv=n1, njv=n1, nkv=n5}, initialState=initial3, label="block-3"}

-- Looking radially outwards at upstream patch faces
--    +------+------+
--    | blk1 | blk0 | ^
--    |      |      | | Z,k
--    +------+------+ |
--               <----+
--                  T,j

-- Looking radially outwards at downstream patch faces
--    +-------------+
--    |    blk3     |  
--    |             |  
--    +-------------+     
--    |    blk2     |  ^
--    |             |  | Z,k
--    +-------------+  |   
--               <-----+ 
--                  T,i
identifyBlockConnections()

if true then

    print("Periodic connection to blocks , blk0 and blk1")
    blk0.bcList[south] = ExchangeBC_FullFace:new{otherBlock=1, otherFace=north,
    orientation=0, reorient_vector_quantities=true,
    Rmatrix={math.cos(-theta), -math.sin(-theta), 0.0, 
             math.sin(-theta),  math.cos(-theta), 0.0, 
             0.0,              0.0,             1.0}}
    
    blk1.bcList[north] = ExchangeBC_FullFace:new{otherBlock=0, otherFace=south,
    orientation=0, reorient_vector_quantities=true,
    Rmatrix={math.cos(theta), -math.sin(theta), 0.0, 
             math.sin(theta),  math.cos(theta), 0.0, 
             0.0,              0.0,             1.0}}

    print("Periodic connection to blocks , blk2 and blk2")
    blk2.bcList[east] = ExchangeBC_FullFace:new{otherBlock=2, otherFace=west,
    orientation=0, reorient_vector_quantities=true,
    Rmatrix={math.cos(theta), -math.sin(theta), 0.0, 
             math.sin(theta),  math.cos(theta), 0.0, 
             0.0,              0.0,             1.0}}

    blk2.bcList[west] = ExchangeBC_FullFace:new{otherBlock=2, otherFace=east,
    orientation=0, reorient_vector_quantities=true,
    Rmatrix={math.cos(-theta), -math.sin(-theta), 0.0, 
             math.sin(-theta),  math.cos(-theta), 0.0, 
             0.0,              0.0,             1.0}}

    print("Periodic connection to blocks , blk3 and blk3")
    blk3.bcList[east] = ExchangeBC_FullFace:new{otherBlock=3, otherFace=west,
    orientation=0, reorient_vector_quantities=true,
    Rmatrix={math.cos(theta), -math.sin(theta), 0.0, 
             math.sin(theta),  math.cos(theta), 0.0, 
             0.0,              0.0,               1.0}}

    blk3.bcList[west] = ExchangeBC_FullFace:new{otherBlock=3, otherFace=east,
    orientation=0, reorient_vector_quantities=true,
    Rmatrix={math.cos(-theta), -math.sin(-theta), 0.0, 
             math.sin(-theta),  math.cos(-theta), 0.0, 
             0.0,              0.0,               1.0}}
end



-- set external BC

if true then
    -- supersonic flow in the radial inwards direction
    blk0.bcList[east] = InFlowBC_Supersonic:new{flowState=stagnation1}
    blk1.bcList[east] = InFlowBC_Supersonic:new{flowState=stagnation2}

    blk2.bcList[north] = OutFlowBC_Simple:new{}
    blk3.bcList[north] = OutFlowBC_Simple:new{}
end

if false then
    -- supersonic flow in the radial outwards direction
    blk0.bcList[east] = OutFlowBC_Simple:new{}
    blk1.bcList[east] = OutFlowBC_Simple:new{}

    blk2.bcList[north] = InFlowBC_Supersonic:new{flowState=stagnation1}
    blk3.bcList[north] = InFlowBC_Supersonic:new{flowState=stagnation2}
end

if false then
    -- subsonic flow in the radial inwards direction
    blk0.bcList[east] = InFlowBC_FromStagnation:new{stagnationState=stagnation1} 
    blk1.bcList[east] = InFlowBC_FromStagnation:new{stagnationState=stagnation2} 

    blk2.bcList[north] = OutFlowBC_FixedP:new{p_outside=Pout}
    blk3.bcList[north] = OutFlowBC_FixedP:new{p_outside=Pout}
end


if false then
    -- subsonic flow in the radial outwards direction
    blk0.bcList[east] = OutFlowBC_FixedP:new{p_outside=Pout}
    blk1.bcList[east] = OutFlowBC_FixedP:new{p_outside=Pout}

    blk2.bcList[north] = InFlowBC_FromStagnation:new{stagnationState=stagnation1} 
    blk3.bcList[north] = InFlowBC_FromStagnation:new{stagnationState=stagnation2} 
end



-- #######################################
-- ### Couple Meshes and Overwrite BCs ###
-- #######################################
if false then -- MappedCell_BC 
    -- Stator Outlet faces: E0.bc_list[WEST]; E1.bc_list[WEST]; E2.bc_list[WEST]; E3.bc_list[WEST]; E4.bc_list[WEST]
    blk0.bcList[west] = ExchangeBC_MappedCell:new{}  -- MappedCellBC(ghost_cell_trans_fn=lambda x, y, z: (x, y, z))
    blk1.bcList[west] = ExchangeBC_MappedCell:new{}  -- MappedCellBC(ghost_cell_trans_fn=lambda x, y, z: (x, y, z))

   
    -- Rotor Inlet faces: BLI_IT0.bc_list[SOUTH]; BLI_IT1.bc_list[SOUTH]; BLI_IB0.bc_list[SOUTH]; BLI_IB1.bc_list[SOUTH]; BLI_IC0.bc_list[SOUTH]; BLI_IC1.bc_list[SOUTH]; 
    blk2.bcList[south] = ExchangeBC_MappedCell:new{}  -- MappedCellBC(ghost_cell_trans_fn=lambda x, y, z: (x, y, z))
    blk3.bcList[south] = ExchangeBC_MappedCell:new{}  -- MappedCellBC(ghost_cell_trans_fn=lambda x, y, z: (x, y, z))
end


if true then  -- attempts at sliding interface
    print("#######################################################")
    print("Creating simulation environment for Mapped (sliding mesh)")
    print("WARNING: Work in progress")
    print("#######################################################")

    config.apply_bcs_in_parallel = false

    -- Do preprocessing of B/C data
    dofile('udf-process_e4.lua')

    blk0.bcList[west] = BoundaryCondition:new{type="user_defined",
					    ghost_cell_data_available=false,
					    convective_flux_computed_in_bc=true,
					    postConvFluxAction = {UserDefinedFlux:new{fileName='udf-stator_out_e4.lua',
                        funcName='convective_flux'}}}
    blk1.bcList[west] = BoundaryCondition:new{type="user_defined",
					    ghost_cell_data_available=false,
					    convective_flux_computed_in_bc=true,
					    postConvFluxAction = {UserDefinedFlux:new{fileName='udf-stator_out_e4.lua',
                      funcName='convective_flux'}}}

    --blk0.bcList[west] = OutFlowBC_Simple:new{}
    --blk1.bcList[west] = OutFlowBC_Simple:new{}


    blk2.bcList[south] = BoundaryCondition:new{type="user_defined",
					    ghost_cell_data_available=false,
					    convective_flux_computed_in_bc=true,
					    postConvFluxAction = {UserDefinedFlux:new{fileName='udf-rotor_in_e4.lua',
                        funcName='convective_flux'}}}
    blk3.bcList[south] = BoundaryCondition:new{type="user_defined",
					    ghost_cell_data_available=false,
					    convective_flux_computed_in_bc=true,
					    postConvFluxAction = {UserDefinedFlux:new{fileName='udf-rotor_in_e4.lua',
                        funcName='convective_flux'}}}

    -- by set_conv_flux=1 --> need convective_flux() in lua
    -- by set_visc_flux=0 --> need interface() in lua  if using viscous simulation  
end




