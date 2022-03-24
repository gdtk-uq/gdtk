-- Structured mesh around a sphere-cone.
-- Simulatiuon is based on spher-cone example used in MECH4480.
-- T = 384.698 K
-- p = 21343.397 Pa
-- v = 1509.033 m/s
--
-- Author: Ingo Jahn
-- Created: 22/3/2022
-- Modified:
--
config.dimensions = 3
nsp, nmodes, gm = setGasModel('ideal-air-gas-model.lua')
inflow = FlowState:new{p=21343.397, T=384.698, velz=-1509.033}
initial = FlowState:new{p=21343.397, T=384.698, velz=-1509.033}
--
-- Refer to README.txt for details. 
-- Domain is constructed by 5 blocks attached to sphere section (nose, side_gap)
-- plus 4 blocks wrapped around the cone section (skirt). 
--
-- Define geometry
R_sphere = 15e-3  -- radius of sphere
R_mesh = 25e-3  -- radius of mesh around sphere
alpha = 70./180*math.pi  -- cone half angle
height = 25e-3  -- height of cone section
delta_radius_rear = 10e-3  -- growth in mesh height along cone
-- 
-- Define Meshing Parameters
nose_fraction = 0.35  -- variable to define size of rectangle on nose
N_REFINE = 1.0
N_nose = 20  -- vertices for square on front
N_ring = 16  -- vertices for rectangles skirting around sides
N_cone = 25  -- vertices in axial direction along cone
N_normal = 15  -- vertices in wall normal direction
N_nose = math.ceil(N_REFINE*N_nose)
N_ring = math.ceil(N_REFINE*N_ring)
N_cone = math.ceil(N_REFINE*N_cone)
N_normal = math.ceil(N_REFINE*N_normal)
--
-- define square sitting on nose
function front_square(s, r, radius, nose_fraction)
    dr = (r-0.5)*2 * nose_fraction
    ds = (s-0.5)*2 * nose_fraction
    x_dash = ds * math.sqrt(1-dr*dr/2)
    y_dash = dr * math.sqrt(1-ds*ds/2)
    -- map to sphere
    rad = math.sqrt(x_dash*x_dash + y_dash*y_dash)
    phi = math.atan2(y_dash, x_dash)
    theta = rad * alpha
    -- conver to cartesian coordinates
    x = radius * math.sin(theta) * math.cos(phi)
    y = radius * math.sin(theta) * math.sin(phi)
    z = radius * math.cos(theta)
    return {x=x , y=y, z=z}
end
-- define rectangles wrapped around front square named by edge they attach to.
function side_gap(r, s, radius, nose_fraction, edge)
    if edge == 'east' then
        ds = -1* (s-0.5)*2 * nose_fraction
        dr = 1 * nose_fraction        
        x_dash_up = math.sin(-1*(-0.5+s)*math.pi/2)
        y_dash_up = math.cos(-1*(-0.5+s)*math.pi/2)
    elseif edge == 'west' then
        ds = 1 * (s-0.5)*2 * nose_fraction
        dr = -1 * nose_fraction        
        x_dash_up = math.sin(-1*(-0.5+s+2)*math.pi/2)
        y_dash_up = math.cos(-1*(-0.5+s+2)*math.pi/2)
    elseif edge == 'north' then
        ds = 1 * nose_fraction  
        dr = 1 * (s-0.5)*2 * nose_fraction   
        x_dash_up = math.sin(1*(0.5-s+1)*math.pi/2)
        y_dash_up = math.cos(1*(0.5-s+1)*math.pi/2)
    elseif edge == 'south' then
        ds = -1 * nose_fraction  
        dr = -1 * (s-0.5)*2 * nose_fraction   
        x_dash_up = math.sin(1*(0.5-s+3)*math.pi/2)
        y_dash_up = math.cos(1*(0.5-s+3)*math.pi/2)
    else
        print('edge option no supported')
    end
    x_dash_low = ds * math.sqrt(1-dr*dr/2)
    y_dash_low = dr * math.sqrt(1-ds*ds/2)
    -- variable t maps between 
    x_dash = x_dash_low * (1-r) + x_dash_up * r
    y_dash = y_dash_low * (1-r) + y_dash_up * r
    -- map to sphere
    rad = math.sqrt(x_dash*x_dash + y_dash*y_dash)
    phi = math.atan2(y_dash, x_dash)
    theta = rad * alpha
    -- conver to cartesian coordinates
    x = radius * math.sin(theta) * math.cos(phi)
    y = radius * math.sin(theta) * math.sin(phi)
    z = radius * math.cos(theta)
    return {x=x , y=y, z=z}
end
-- define rectangles that are wrapped around the cone
function skirt(r, s, radius_front, height, delta_radius_rear, alpha, edge)
    if edge == 'c_east' then
        x_dash = math.sin(-1*(-0.5+s)*math.pi/2)
        y_dash = math.cos(-1*(-0.5+s)*math.pi/2)
    elseif edge == 'c_west' then
        x_dash = math.sin(-1*(-0.5+s+2)*math.pi/2)
        y_dash = math.cos(-1*(-0.5+s+2)*math.pi/2)
    elseif edge == 'c_north' then
        x_dash = math.sin(1*(0.5-s+1)*math.pi/2)
        y_dash = math.cos(1*(0.5-s+1)*math.pi/2)
    elseif edge == 'c_south' then
        x_dash = math.sin(1*(0.5-s+3)*math.pi/2)
        y_dash = math.cos(1*(0.5-s+3)*math.pi/2)
    else
        print('edge option no supported')
    end
    -- map to sphere
    rad = math.sqrt(x_dash*x_dash + y_dash*y_dash)
    phi = math.atan2(y_dash, x_dash)
    theta = rad * alpha
    -- conver to cartesian coordinates
    x_front = radius_front * math.sin(theta) * math.cos(phi)
    y_front = radius_front * math.sin(theta) * math.sin(phi)
    z_front = radius_front * math.cos(theta)
    --
    x_rear = x_front + (height/math.tan(alpha) + delta_radius_rear) * math.sin(theta) * math.cos(phi)
    y_rear = y_front + (height/math.tan(alpha) + delta_radius_rear)  * math.sin(theta) * math.sin(phi)
    z_rear = z_front - height
    --
    x = x_front * (1-r) + x_rear * r
    y = y_front * (1-r) + y_rear * r
    z = z_front * (1-r) + z_rear * r
    return {x=x , y=y, z=z}
end
--
-- Create 3-D blocks
blocks = {"nose", "east", "west", "north", "south", "c_east", "c_west", "c_north", "c_south", }
vols = {}
for _, b in ipairs(blocks) do
   -- define LuaFnVolume functions locally.
   local v_name = string.format("vf_%s", b)
   if b == "nose" then
      -- do nose
      _G[v_name] = function (r, s, t)
          pos0 = front_square(r, s, R_sphere, nose_fraction)
          pos1 = front_square(r, s, R_mesh, nose_fraction)
          x = (1-t)*pos0.x + t*pos1.x
          y = (1-t)*pos0.y + t*pos1.y
          z = (1-t)*pos0.z + t*pos1.z
          return {x=x , y=y, z=z}
      end
   elseif b == "north" or b == "south" or b == "east" or b == "west" then
      -- do blocks around nose
      _G[v_name] = function (r, s, t)
          pos0 = side_gap(r, s, R_sphere, nose_fraction, b)
          pos1 = side_gap(r, s, R_mesh, nose_fraction, b)
          x = (1-t)*pos0.x + t*pos1.x
          y = (1-t)*pos0.y + t*pos1.y
          z = (1-t)*pos0.z + t*pos1.z
          return {x=x , y=y, z=z}
      end 
   else
      -- do blocks around cone
      _G[v_name] = function (r, s, t)
          pos0 = skirt(r, s, R_sphere, height, 0, alpha, b)
          pos1 = skirt(r, s, R_mesh, height, delta_radius_rear, alpha, b)
          x = (1-t)*pos0.x + t*pos1.x
          y = (1-t)*pos0.y + t*pos1.y
          z = (1-t)*pos0.z + t*pos1.z
          return {x=x , y=y, z=z}
      end 
   end
   -- create parametric volume
    vols[b] = LuaFnVolume:new{luaFnName=v_name}
end
--
-- apply clustering and create grids
cf = RobertsFunction:new{end0=true, end1=false, beta=1.1}
grids = {}
for _, b in ipairs(blocks) do
   if b == "nose" then
       -- do nose
       niv = N_nose
       njv = N_nose
       cfList = {edge04=cf, edge15=cf, edge26=cf, edge37=cf}
   elseif b == "north" or b == "south" or b == "east" or b == "west" then
       -- do blocks around nose
       niv = N_ring
       njv = N_nose
       cfList = {edge04=cf, edge15=cf, edge26=cf, edge37=cf}
   else
       -- do blocks around cone
       niv = N_cone
       njv = N_nose
       cfList = {edge04=cf, edge15=cf, edge26=cf, edge37=cf}
   end
   nkv = N_normal
   grids[b] = StructuredGrid:new{pvolume=vols[b], cfList=cfList,
                                 niv=niv, njv=njv, nkv=nkv}
end
-- and, finally, the FluidBlocks
blks = {}
nib = 1; njb = 1; nkb = 1
wall_bc = WallBC_WithSlip:new{}
inflow_bc = InOutFlowBC_Ambient:new{flowState=inflow}
for _, b in ipairs(blocks) do
   if b == "nose" or b == "east" or b == "south" or b == "top" then
      -- do nose and blocks around nose
      bcList = {top=inflow_bc, bottom=wall_bc}
   else
      -- do blockes around cone
      bcList = {top=inflow_bc, bottom=wall_bc, east=inflow_bc}
   end
   blks[b] = FBArray:new{grid=grids[b], initialState=initial,
                         bcList=bcList, nib=nib, njb=njb, nkb=nkb}
end
identifyBlockConnections()
--
config.max_time = 2.0e-3  -- seconds
config.max_step = 200000
config.dt_plot = 0.25e-3  -- seconds