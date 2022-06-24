-- Structured mesh around a sphere-cone in 3D, for parallel testing.
--
-- Author: Nick Gibbons and Ingo Jahn
-- Created: 01/06/2022


-- Define geometry
R_sphere = 5e-3  -- radius of sphere
R_mesh = 11e-3  -- radius of mesh around sphere
sphere_angle = 75./180*math.pi  -- cone half angle
mesh_angle = 55./180*math.pi  -- cone half angle
height = 35e-3  -- height of cone section
offset = 4e-3
 
-- Define Meshing Parameters
nose_fraction = 0.35  -- variable to define size of rectangle on nose
N_REFINE = 0.5
N_nose = 36  -- vertices for square on front
N_ring = 26  -- vertices for rectangles skirting around sides
N_cone = 80  -- vertices in axial direction along cone
N_normal = 60  -- vertices in wall normal direction
N_nose = math.ceil(N_REFINE*N_nose)
N_ring = math.ceil(N_REFINE*N_ring)
N_cone = math.ceil(N_REFINE*N_cone)
N_normal = math.ceil(N_REFINE*N_normal)

-- define square sitting on nose
function front_square(s, r, radius, nose_fraction, alpha, offset)
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
    z = radius * math.cos(theta) - offset
    return {x=x , y=y, z=z}
end
-- define rectangles wrapped around front square named by edge they attach to.
function side_gap(r, s, radius, nose_fraction, alpha, offset, edge)
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
    z = radius * math.cos(theta) - offset
    return {x=x , y=y, z=z}
end
-- define rectangles that are wrapped around the cone
function skirt(r, s, radius_front, height, alpha, offset, edge)
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
    -- Cone flank is specified by a length and an angle
    cone_z = radius_front*math.cos(alpha) + height - radius_front
    --cone_z = height
    cone_h = cone_z/math.tan(alpha)-- + delta_radius_rear
    x_rear = x_front + cone_h*math.sin(theta)*math.cos(phi)
    y_rear = y_front + cone_h*math.sin(theta)*math.sin(phi)
    z_rear = z_front - cone_z
    --
    x = x_front * (1-r) + x_rear * r
    y = y_front * (1-r) + y_rear * r
    z = z_front * (1-r) + z_rear * r - offset
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
          pos0 = front_square(r, s, R_sphere, nose_fraction, sphere_angle, 0.0)
          pos1 = front_square(r, s, R_mesh, nose_fraction, mesh_angle, offset)
          x = (1-t)*pos0.x + t*pos1.x
          y = (1-t)*pos0.y + t*pos1.y
          z = (1-t)*pos0.z + t*pos1.z
          return {x=x , y=y, z=z}
      end
   elseif b == "north" or b == "south" or b == "east" or b == "west" then
      -- do blocks around nose
      _G[v_name] = function (r, s, t)
          pos0 = side_gap(r, s, R_sphere, nose_fraction, sphere_angle, 0.0, b)
          pos1 = side_gap(r, s, R_mesh, nose_fraction, mesh_angle, offset, b)
          x = (1-t)*pos0.x + t*pos1.x
          y = (1-t)*pos0.y + t*pos1.y
          z = (1-t)*pos0.z + t*pos1.z
          return {x=x , y=y, z=z}
      end 
   else
      -- do blocks around cone
      _G[v_name] = function (r, s, t)
          pos0 = skirt(r, s, R_sphere, height, sphere_angle, 0.0, b)
          pos1 = skirt(r, s, R_mesh, height+(R_mesh-R_sphere)-offset, mesh_angle, offset, b)
          x = (1-t)*pos0.x + t*pos1.x
          y = (1-t)*pos0.y + t*pos1.y
          z = (1-t)*pos0.z + t*pos1.z
          return {x=x , y=y, z=z}
      end 
   end
   -- create parametric volume
    vols[b] = LuaFnVolume:new{luaFnName=v_name}
end

-- apply clustering and create grids
a_factor = 0.005
cf = GeometricFunction:new{a=a_factor, r=1.25, N=N_normal}
cfskirt = GaussianFunction:new{m=0.0001, s=0.7, ratio=0.2}

grids = {}
for _, b in ipairs(blocks) do
   if b == "nose" then
       -- do nose
       niv = N_nose
       njv = N_nose
       cfList = {edge04=cf, edge15=cf, edge26=cf, edge37=cf,
                 edge01=cfnose, edge12=cfnose, edge23=cfnose, edge34=cfnose,
                 edge45=cfnose, edge56=cfnose, edge67=cfnose, edge47=cfnose}
   elseif b == "north" or b == "south" or b == "east" or b == "west" then
       -- do blocks around nose
       niv = N_ring
       njv = N_nose
       cfList = {edge04=cf, edge15=cf, edge26=cf, edge37=cf}
   else
       -- do blocks around cone
       niv = N_cone
       njv = N_nose
       cfList = {edge04=cf, edge15=cf, edge26=cf, edge37=cf,
                 edge01=cfskirt, edge32=cfskirt, edge45=cfskirt, edge76=cfskirt}
   end
   nkv = N_normal
   grids[b] = StructuredGrid:new{pvolume=vols[b], cfList=cfList,
                                 niv=niv, njv=njv, nkv=nkv}
end

ugrids = {}
for name, sgrid in pairs(grids) do
   ugrids[name] = UnstructuredGrid:new{sgrid=sgrid}
end

grid = ugrids["nose"]
grid:joinGrid(ugrids["north"])
grid:joinGrid(ugrids["south"])
grid:joinGrid(ugrids["east"])
grid:joinGrid(ugrids["west"])
grid:joinGrid(ugrids["c_north"])
grid:joinGrid(ugrids["c_south"])
grid:joinGrid(ugrids["c_east"])
grid:joinGrid(ugrids["c_west"])
--grid:write_to_vtk_file('grid.vtk')
grid:write_to_su2_file('grid.su2')
