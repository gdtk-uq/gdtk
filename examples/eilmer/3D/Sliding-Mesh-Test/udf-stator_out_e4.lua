-- udf-supersonic-in.lua
-- Lua script for the user-defined functions 
-- called by the UserDefinedGhostCell BC.
--
-- This particular example is defining the constant supersonic inflow
-- for the cone20 test case.

--function ghost_cell(args)
   -- Function that returns the flow state for a ghost cell
   -- for use in the inviscid flux calculations.
   --
   -- args contains t, x, y, z, csX, csY, csZ, i, j, k, which_boundary
   --
   -- Sample the flow field at the current cell 
   -- which is beside the boundary.
--   cell = sample_flow(block_id, args.i, args.j, args.k)
--   return cell, cell
--end

function interface(args)
   --print("Shouldnt be here")
   -- Function that returns the conditions at the boundary 
   -- when viscous terms are active.
   --print(block_id)
   --print(sample_flow(block_id, args.i, args.j, args.k))
   return sample_flow(block_id, args.i, args.j, args.k)
end

function convective_flux(args)
   -- Function that returns the fluxes of conserved quantities.
   -- For use in the inviscid flux calculations.
   --
   -- args contains t, x, y, z, csX, csY, csZ, i, j, k, which_boundary
   --
   -- Set constant conditions across the whole boundary.
   --print("--------------------------------------")
   --print("Hello from function converctive_flux.")
   --print("You are in udf-stator_out.lua")


   --print("Available functions and variables")
   --for k,v in pairs(_G) do
   --    print("Global key:", k,"value:",v)
   --end
   --print("content of args")
   --for k,v in pairs(args) do
   --     print("args key:", k,"value:",v)
   --end


   if (args.i==2 and args.j==2 and args.k==2) then -- loading routine only needs to be run at first call (imin, jmin, kmin) 
      -- might be able to pre-allocate this
      UP,err = table.load( "UP_tbl.lua" )
      UP_Nfaces = UP[0]
      UP_face_ind = UP[1]
      UP_row_ind = UP[2]
      UP_weight = UP[3]
      UP_area = UP[4]
      T_down_pos_low = UP[5]
      T_down_pos_high = UP[6]
      T_down_col = UP[7]
      T_down_i = UP[8]
      T_down_length = UP[9]

      BLK_MAP,err = table.load("BLK_MAP_tbl.lua")

      DATA,err = table.load("CONFIG_tbl.lua")
      N_blade = DATA[0]
      OMEGA = DATA[1]
      DOWN_row_list = DATA[3]
      --downstream_list_row1 = DATA[6]
      --downstream_list_row2 = DATA[7]

      Theta_min = DATA[11]
      relax = DATA[13]   

      -- calculate angle shift
      Arc = 2. * math.pi / N_blade
      Theta_max = Theta_min + Arc

      if OMEGA > 0 then
         Rotation = ((OMEGA*args.t) / Arc - math.floor((OMEGA*args.t) / Arc)) * Arc
      elseif OMEGA < 0 then
         Rotation = ((-OMEGA*args.t) / Arc - math.floor((-OMEGA*args.t) / Arc)) * Arc * (-1)
      else
         Rotation = 0. 
         --Rotation = -Arc/2
      end
      --print("\n")
      --print("Time:", args.t)
      --print("Rotation:",Rotation)
      --print("Theta_min:",Theta_min)
      --print("Theta_max:",Theta_max)
   end

    --print("Time:", args.t)
    --print("Rotation:",Rotation)

    --print("/n")
    --print("to do LIST")
    --print("turbulence models")


   -- access vtx locations for current patch. 
   --print("\n")
   --print("Stator outlet, working on: Block,I,J,K",blkId,args.i,args.j,args.k)
   vtx3 = getVtxPosition(blkId,args.i,args.j+1,args.k)
   vtx0 = getVtxPosition(blkId,args.i,args.j,args.k)
   -- vtx4 = sample_vtx(blk_id,args.i,args.j,args.k+1)
   -- vtx7 = sample_vtx(blk_id,args.i,args.j+1,args.k+1)
   T_low = math.atan(vtx0.y/vtx0.x)
   T_high = math.atan(vtx3.y/vtx3.x)
   --print("Angles before Rotation: T_low, T_high:", T_low,T_high)

   -- adjust because of rotation
   T_low = T_low - Rotation
   T_high = T_high - Rotation
   --print("Angles after Rotation: T_low, T_high:", T_low,T_high)
   Periodic_flag = 0
   if (T_high <= Theta_max and T_low >= Theta_min) then 
      --print("nominal case")
      Angle = Rotation
   elseif (T_high > Theta_max and T_low >= Theta_max) then
      --print("both above Theta_max")
      T_low = T_low - Arc 
      T_high = T_high - Arc 
      --Angle = Rotation - Arc  --dont touch 
      Angle = Arc + Rotation
      -- Angle = Theta_min + (Rotation - Theta_max)
   elseif (T_high <= Theta_min and T_low < Theta_min) then
      --print("both below Theta_min")
      T_low = T_low + Arc 
      T_high = T_high + Arc 
      Angle = Arc - Rotation
      --Angle = Rotation + Arc  --dont touch
      -- Angle = Theta_max + (Rotation - Theta_max)
   else -- cell goes over periodic B/C
      Periodic_flag = 1
      if T_high > Theta_max then
         T2 = T_high - Arc
         T1 = Theta_min  
         T_low = T_low
         T_high = Theta_max
         Angle = Rotation
         Angle2 = Arc + Rotation 
      end
      if T_low < Theta_min then
         T2 = Theta_max  
         T1 = T_low + Arc
         T_high = T_high
         T_low = Theta_min
         Angle = Rotation
         Angle2 = Arc - Rotation 
      end
      --print("Periodic: Angles, T_low, T_high, T1, T2:", T_low,T_high, T1, T2)
      --print("Angle:", Angle, "Angle2:", Angle2)
   end 
   --if block_id == 12 then
      --print("Angles, T_low, T_high:", T_low,T_high)
      --print("Angle:", Angle)
   --end


   -- get list of overalpping cells in theta_direction
   Target_low_vtx = T_low
   Target_high_vtx = T_high
   Source_length = T_down_length
   Source_pos_low = T_down_pos_low
   Source_pos_high = T_down_pos_high
   T_count = 0
   T_faces = {}
   T_cols = {}
   T_weights = {}
   for j = 1,Source_length do 
      Source_low_vtx = Source_pos_low[j]    
      Source_high_vtx = Source_pos_high[j]
      --print("Source_angles:",Source_low_vtx,Source_high_vtx)
      --print("Target_angles:",Target_low_vtx,Target_high_vtx)
      if ( (Source_high_vtx <= Target_high_vtx) and (Source_low_vtx >= Target_low_vtx) ) then -- face fully inside
         --print(Source_high_vtx <= Target_high_vtx)
         --print(Source_low_vtx >= Target_low_vtx)
         T_count = T_count+1
         T_faces[T_count] = T_down_i[j]
         T_cols[T_count] = T_down_col[j]
         T_weights[T_count] = 1.0
         --print("fully inside")
      elseif ( (Source_high_vtx > Target_high_vtx) and (Source_low_vtx < Target_high_vtx) ) then -- only overlap at the top
         T_count = T_count+1
         T_faces[T_count] = T_down_i[j]
         T_cols[T_count] = T_down_col[j]
         T_weights[T_count] = (Target_high_vtx-Source_low_vtx) / (Source_high_vtx-Source_low_vtx)
         --print("top overlap")
      elseif ( (Source_high_vtx > Target_low_vtx) and (Source_low_vtx < Target_low_vtx) ) then -- only overlap at the bottom
         T_count = T_count+1
         T_faces[T_count] = T_down_i[j]
         T_cols[T_count] = T_down_col[j]
         T_weights[T_count] = (Source_high_vtx - Target_low_vtx) / (Source_high_vtx-Source_low_vtx)
         --print("bottom overlap")
      elseif ( (Source_high_vtx > Target_high_vtx) and (Source_low_vtx < Target_low_vtx) ) then -- overlap at both side
         T_count = T_count+1
         T_faces[T_count] = T_down_i[j]
         T_cols[T_count] = T_down_col[j]
         T_weights[T_count] = (Target_high_vtx - Target_low_vtx) / (Source_high_vtx-Source_low_vtx)
         --print("both sides")
      else -- other cases
         --print("No overlap so no action required.") 
      end
   end

   --print("Number of upstream faces overlaping with current downstream face is :",T_count)
   --for j=1,T_count do
   --   print("Face_ind, Col_ind, weight:", T_faces[j],T_cols[j],T_weights[j])
   --end

   -- use row and column to find block index
   -- Looking Radially outwards
   --   Col1  Col0 
   -- +-----+-----+
   -- |     |     |
   -- | IT0 | IT1 |  Row2
   -- |  20 |  21 |  
   -- +-----+-----+
   -- |     |     |
   -- | IC0 | IC1 |  Row1
   -- |  23 |  24 |  
   -- +-----+-----+
   -- |     |     |
   -- | IB0 | IB1 |  Row0
   -- |  26 |  27 |  ^
   -- +-----+-----+  | Z,k
   --         T <----+
   --              i  


   -- area average properties to ensure conservative exchange
   Nz = args.k - 1 -- minus one to account for ghost cells 
   Z_Nfaces = UP_Nfaces[Nz]
   Z_face_ind = UP_face_ind[Nz]
   Z_row_ind = UP_row_ind[Nz]
   Z_weight = UP_weight[Nz]
   Tot_area = 0
   rho = 0.; T = 0.; p = 0.
   u_temp = 0.; v_temp = 0.; w = 0.
   tke = 0.; omega =0.
   --print("Number of faces being sample in Z direction:",Z_Nfaces," in T direction:", T_count)
   for i=1,Z_Nfaces do
      for j=1,T_count do
         weight = Z_weight[i] * T_weights[j]

         --if Z_row_ind[i] == 0 then
         --   b_id = downstream_list_row0[T_cols[j]+1]  -- +1 as column number start from 0
         --elseif Z_row_ind[i] == 1 then
         --   b_id = downstream_list_row1[T_cols[j]+1]  -- +1 as column number start from 0
         --elseif Z_row_ind[i] == 2 then 
         --   b_id = downstream_list_row2[T_cols[j]+1]  -- +1 as column number start from 0
         --else
         --   print("error")
         --   print(Z_row_ind[i])
         --end  

         row = DOWN_row_list[ Z_row_ind[i] +1]  
         b_id = row[T_cols[j]+1]  -- +1 as column number start from 0

         -- print("Stator",b_id,T_cols[j]) 

         --print("Stator_out Sampling from, blk, i, j, k:",  b_id,T_faces[j],2,Z_face_ind[i]) -- set i = 2 as going along WEST face
         --print("weight =", weight)  
         if args.timeStep == 0 then   
            flow = sampleFlow(b_id,T_faces[j],2,Z_face_ind[i])
            flow2 = sampleFace("j",b_id,T_faces[j],2,Z_face_ind[i])

            area = flow2.area
         else 
            flow = sampleFace("j",b_id,T_faces[j],2,Z_face_ind[i])
            area = flow.area
         end

         Tot_area = Tot_area + area * weight
         rho    =  rho    + flow.rho  * weight * area
         p      =  p      + flow.p    * weight * area 
         T      =  T      + flow.T    * weight * area 
         u_temp =  u_temp + flow.velx * weight * area 
         v_temp =  v_temp + flow.vely * weight * area 
         w      =  w      + flow.velz * weight * area 
         tke    =  tke    + flow.tke  * weight * area 
         omega  =  omega  + flow.omega* weight * area 
      end 
   end

   -- normalise the properties 
   --print("Tot_weight:",Tot_weight)
   if Tot_area == 0 then
      rho = 0; T = 0; p = 0
      u_temp = 0; v_temp = 0; w = 0
   else
      rho = rho / Tot_area; T = T / Tot_area; p = p / Tot_area  
      u_temp = u_temp / Tot_area; v_temp = v_temp / Tot_area; w = w / Tot_area
   end 

   -- rotate u and v to account for relative rotation
   u = u_temp * math.cos(Angle) - v_temp * math.sin(Angle)
   v = u_temp * math.sin(Angle) + v_temp * math.cos(Angle)

   --print("Velocities:", u_temp, v_temp, u, v)
   --print("Angle:", Angle)


   if Periodic_flag == 1 then
      --print("Yuppy")
      -- repeat above for second segment at other end of range
      T_low = T1
      T_high = T2

      Target_low_vtx = T_low
      Target_high_vtx = T_high
      Source_length = T_down_length
      Source_pos_low = T_down_pos_low
      Source_pos_high = T_down_pos_high
      T_count = 0
      T_faces = {}
      T_cols = {}
      T_weights = {}
      for j = 1,Source_length do 
         Source_low_vtx = Source_pos_low[j]    
         Source_high_vtx = Source_pos_high[j]
         --print("Source_angles:",Source_low_vtx,Source_high_vtx)
         if ( (Source_high_vtx <= Target_high_vtx) and (Source_low_vtx >= Target_low_vtx) ) then -- face fully inside
            --print(Source_high_vtx <= Target_high_vtx)
            --print(Source_low_vtx >= Target_low_vtx)
            T_count = T_count+1
            T_faces[T_count] = T_down_i[j]
            T_cols[T_count] = T_down_col[j]
            T_weights[T_count] = 1.0
            --print("fully inside")
         elseif ( (Source_high_vtx > Target_high_vtx) and (Source_low_vtx < Target_high_vtx) ) then -- only overlap at the top
            T_count = T_count+1
            T_faces[T_count] = T_down_i[j]
            T_cols[T_count] = T_down_col[j]
            T_weights[T_count] = (Target_high_vtx-Source_low_vtx) / (Source_high_vtx-Source_low_vtx)
            --print("top overlap")
         elseif ( (Source_high_vtx > Target_low_vtx) and (Source_low_vtx < Target_low_vtx) ) then -- only overlap at the bottom
            T_count = T_count+1
            T_faces[T_count] = T_down_i[j]
            T_cols[T_count] = T_down_col[j]
            T_weights[T_count] = (Source_high_vtx - Target_low_vtx) / (Source_high_vtx-Source_low_vtx)
            --print("bottom overlap")
         elseif ( (Source_high_vtx > Target_high_vtx) and (Source_low_vtx < Target_low_vtx) ) then -- overlap at both side
            T_count = T_count+1
            T_faces[T_count] = T_down_i[j]
            T_cols[T_count] = T_down_col[j]
            T_weights[T_count] = (Target_high_vtx - Target_low_vtx) / (Source_high_vtx-Source_low_vtx)
            --print("both sides")
         else -- other cases
            --print("No overlap so no action required.") 
         end
      end

      --print("Number of upstream faces overlaping with current downstream face is :",T_count)
      --for j=1,T_count do
      --   print("Face_ind, Col_ind, weight:", T_faces[j],T_cols[j],T_weights[j])
      --end

      -- use row and column to find block index
      -- Looking Radially outwards
      --   Col1  Col0 
      -- +-----+-----+
      -- |     |     |
      -- | IT0 | IT1 |  Row2
      -- |  20 |  21 |  
      -- +-----+-----+
      -- |     |     |
      -- | IC0 | IC1 |  Row1
      -- |  23 |  24 |  
      -- +-----+-----+
      -- |     |     |
      -- | IB0 | IB1 |  Row0
      -- |  26 |  27 |  ^
      -- +-----+-----+  | Z,k
      --         T <----+
      --              i  

      -- area average properties to ensure conservative exchange
      --Nz = args.k - 1 -- minus one to account for ghost cells 
      --Z_Nfaces = DOWN_Nfaces[Nz]
      --Z_face_ind = DOWN_face_ind[Nz]
      --Z_row_ind = DOWN_row_ind[Nz]
      -- Z_weight = DOWN_weight[Nz]
      Tot_area2 = 0
      rho2 = 0.; T2 = 0.; p2 = 0.
      u_temp2 = 0.; v_temp2 = 0.; w2 = 0.
      tke2 = 0.; omega2 =0.
      --print("Number of faces being sample in Z direction:",Z_Nfaces," in T direction:", T_count)
      for i=1,Z_Nfaces do
         for j=1,T_count do
            weight = Z_weight[i] * T_weights[j]

            --if Z_row_ind[i] == 0 then
            --   b_id = downstream_list_row0[T_cols[j]+1]  -- +1 as column number start from 0
            --elseif Z_row_ind[i] == 1 then
            --   b_id = downstream_list_row1[T_cols[j]+1]  -- +1 as column number start from 0
            --elseif Z_row_ind[i] == 2 then 
            --   b_id = downstream_list_row2[T_cols[j]+1]  -- +1 as column number start from 0
            --else
            --   print("error")
            --   print(Z_row_ind[i])
            --end   


            row = DOWN_row_list[ Z_row_ind[i] +1]  
            b_id = row[T_cols[j]+1]  -- +1 as column number start from 0


            --print("Stator_out Periodic Sampling from, blk, i, j, k:",  b_id,T_faces[j],2,Z_face_ind[i]) -- set i = 2 as going along WEST face
            --print("weight =", weight)  
            if args.timeStep == 0 then   
               flow = sampleFlow(b_id,T_faces[j],2,Z_face_ind[i])
               flow2 = sampleFace("j",b_id,T_faces[j],2,Z_face_ind[i])
               area = flow2.area
            else 
               flow = sampleFace("j",b_id,T_faces[j],2,Z_face_ind[i])
               area = flow.area
            end
            Tot_area2 = Tot_area2 + area * weight
            rho2    =  rho2    + flow.rho  * weight * area
            p2      =  p2      + flow.p    * weight * area
            T2      =  T2      + flow.T    * weight * area 
            u_temp2 =  u_temp2 + flow.velx * weight * area 
            v_temp2 =  v_temp2 + flow.vely * weight * area 
            w2      =  w2      + flow.velz * weight * area 
            tke2    =  tke2    + flow.tke  * weight * area 
            omega2  =  omega2  + flow.omega* weight * area 
         end 
      end
      -- normalise properties
      --rho2 = rho2 / Tot_weight2; T2 = T2 / Tot_weight2; p2 = p2 / Tot_weight2  
      --u_temp2 = u_temp2 / Tot_weight2; v_temp2 = v_temp2 / Tot_weight2; w2 = w2 / Tot_weight2

      u2 = u_temp2 * math.cos(Angle2) - v_temp2 * math.sin(Angle2)
      v2 = u_temp2 * math.sin(Angle2) + v_temp2 * math.cos(Angle2)

      --print("Periodic Velocities:", u_temp2, v_temp2, u2, v2)

      rho = (rho*Tot_area + rho2) / (Tot_area+Tot_area2)
      p = (p*Tot_area + p2) / (Tot_area+Tot_area2)
      T = (T*Tot_area + T2) / (Tot_area+Tot_area2)
      u = (u*Tot_area + u2) / (Tot_area+Tot_area2)
      v = (v*Tot_area + v2) / (Tot_area+Tot_area2)
      w = (w*Tot_area + w2) / (Tot_area+Tot_area2)
      tke = (tke*Tot_area + tke2) / (Tot_area+Tot_area2)
      omega = (omega*Tot_area + omega2) / (Tot_area+Tot_area2)
   end



   -- remove very small velocities to eliminate uncontrolled run-away growth. There is a feedback loop in mixing plane between cells that don't fully cover each other.
   if math.abs(u) < 1.e-3 then
      u = 0.
   end
   if math.abs(v) < 1.e-3 then
      v = 0.
   end
   if math.abs(w) < 1.e-3 then
      w = 0.
   end

   --print("Stator_out")
   --print("Output: u,v,w,p,rho,T[0]:",u,v,w,p,rho,T)

   if math.abs(args.csZ) < 1e-10 then
      CSZ = 0
   else
      CSZ = args.csZ
   end 
   
   FLOW = sampleFace("i",blkId,args.i,args.j,args.k)
   --relax = 0.1 -- fraction to which properties are updated.


   u_old = FLOW.velx
   u_new = u
   u = u_new * relax + (1-relax) * u_old    
   v_old = FLOW.vely
   v_new = v
   v = v_new * relax + (1-relax) * v_old  
   w_old = FLOW.velz
   w_new = w
   w = w_new * relax + (1-relax) * w_old  
   rho_old = FLOW.rho
   rho_new = rho
   rho = rho_new * relax + (1-relax) * rho_old  
   p_old = FLOW.p
   p_new = p
   p = p_new * relax + (1-relax) * p_old 
   T_old = FLOW.T
   T_new = T
   T = T_new * relax + (1-relax) * T_old 
   tke_old = FLOW.tke
   tke_new = tke
   tke = tke_new * relax + (1-relax) * tke_old 
   omega_old = FLOW.omega
   omega_new = omega
   omega = omega_new * relax + (1-relax) * omega_old 

    -- get gas model information
    Q = GasState:new{gmodel}
    Cv = gmodel:Cv{Q}

   massf = Q.massf   -- mass fractions are indexed from 0 to nsp-1
   -- Assemble flux vector
   flux = {}
   flux.mass = rho * (u*args.csX + v*args.csY + w*CSZ)  -- kg/s/m**2
   --print("Fmass",F.mass)
   flux.momentum_x = p * args.csX + u * flux.mass
   flux.momentum_y = p * args.csY + v * flux.mass
   flux.momentum_z = p * CSZ + w * flux.mass
   flux.total_energy = flux.mass * (Cv*T + 0.5*(u*u+v*v+w*w) + p/rho)
   flux.tke = tke * flux.mass
   flux.omega = omega * flux.mass
   --flux.species = {}
   --flux.species[0] = flux.mass * massf[0]
   --flux.renergies = {}
   --flux.renergies[0] = flux.mass * (Cv*T)

   --if block_id == 12 then
      --print("Block,I,J,K",blkId,args.i,args.j,args.k,"Zmomentum:",flux.momentum_z, p, args.csZ, CSZ, w)
      --print(flux.mass, flux.momentum_x,flux.momentum_y,flux.momentum_z,flux.total_energy)
      --print(u,v,w)
   --end

   return flux
end


--// The Load Function
function table.load( sfile )
   local ftables,err = loadfile( sfile )
   if err then return _,err end
   local tables = ftables()
   for idx = 1,#tables do
      local tolinki = {}
      for i,v in pairs( tables[idx] ) do
         if type( v ) == "table" then
            tables[idx][i] = tables[v[1]]
         end
         if type( i ) == "table" and tables[i[1]] then
            table.insert( tolinki,{ i,tables[i[1]] } )
         end
      end
      -- link indices
      for _,v in ipairs( tolinki ) do
         tables[idx][v[2]],tables[idx][v[1]] =  tables[idx][v[1]],nil
      end
   end
   return tables[1]
end



