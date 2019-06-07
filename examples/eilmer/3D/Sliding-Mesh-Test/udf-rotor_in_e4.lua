-- udf-supersonic-in.lua
-- Lua script for the user-defined functions 
-- called by the UserDefinedGhostCell BC.
--
-- This particular example is defining the constant supersonic inflow
-- for the cone20 test case.


--function ghost_cell(args)
   -- Function that returns the flow states for a ghost cells.
   -- For use in the inviscid flux calculations.
   --
   -- args contains t, x, y, z, csX, csY, csZ, i, j, k, which_boundary
   -- but we don't happen to us any of them.
   --
   -- Set constant conditions across the whole boundary.
   -- print("Hello from function ghost_cell.")
--   ghost = {}
--   ghost.p = 95.84e3 -- pressure, Pa
--   ghost.T = {}  -- temperatures, K (as a table)
--   ghost.T[0] = 1103.0
--   ghost.u = 1000.0  -- x-velocity, m/s
--   ghost.v = 0.0     -- y-velocity, m/s
--   ghost.w = 0.0
--   ghost.massf = {} -- mass fractions to be provided as a table
--   ghost.massf[0] = 1.0 -- mass fractions are indexed from 0 to nsp-1
--   return ghost, ghost
--end



-- use dofile ([filename]) to load table saved as strings in file. I.e. autmactially recrete the table.

require 'lua_helper'  -- loads lua helper functions

function convectiveFlux(args)
   -- Function that returns the fluxes of conserved quantities.
   -- For use in the inviscid flux calculations.
   --
   -- args contains t, x, y, z, csX, csY, csZ, i, j, k, which_boundary
   --
   -- Set constant conditions across the whole boundary.
   --print("--------------------------------------")
   --print("Hello from function converctive_flux.")
   --print("You are in udf-rotor_in.lua")


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
      DOWN,err = table.load("DOWN_tbl.lua")
      DOWN_Nfaces = DOWN[0]
      DOWN_face_ind = DOWN[1]
      DOWN_row_ind = DOWN[2]
      DOWN_weight = DOWN[3]
      DOWN_area = DOWN[4]
      T_up_pos_low = DOWN[5]
      T_up_pos_high = DOWN[6]
      T_up_col = DOWN[7]
      T_up_j = DOWN[8]
      T_up_length = DOWN[9]

      BLK_MAP,err = table.load("BLK_MAP_tbl.lua")

      DATA,err = table.load("CONFIG_tbl.lua")
      N_blade = DATA[0]
      OMEGA = DATA[1]
      UP_row_list = DATA[2]
      Theta_min = DATA[11]
      relax = DATA[14]  

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
      --print("Rotation",Rotation)
   end

   --print("Time:", args.t)
   --print("Rotation",Rotation)

   -- access vtx locations for current patch. 
   --print("\n")
   --print("Rotor inlet, working on: Block,I,J,K",blkId,args.i,args.j,args.k)
   vtx0 = getVtxPosition(blkId,args.i,args.j,args.k)
   vtx1 = getVtxPosition(blkId,args.i+1,args.j,args.k)
   -- vtx5 = getVtxPosition(blkId,args.i+1,args.j,args.k+1)
   -- vtx4 = getVtxPosition(blkId,args.i,args.j,args.k+1)
   T_low = math.atan(vtx0.y/vtx0.x)
   T_high = math.atan(vtx1.y/vtx1.x)
   --print("Angles before Rotation: T_low, T_high:", T_low,T_high)

   -- adjust because of rotation
   T_low = T_low + Rotation
   T_high = T_high + Rotation
   --print("Angles after Rotation: T_low, T_high:", T_low,T_high)
   Periodic_flag = 0
   if T_high <= Theta_max and T_low >= Theta_min then 
      --print("nominal case")
      Angle = - Rotation
   elseif T_high > Theta_max and T_low >= Theta_max then
      --print("both above Theta_max")
      T_low = T_low - Arc --dont touch
      T_high = T_high - Arc --dont touch
      -- Angle = -Rotation + Arc  --dont touch
      Angle = -Arc + Rotation 
   elseif T_high <= Theta_min and T_low < Theta_min then
      --print("both below Theta_min")
      T_low = T_low + Arc --dont touch
      T_high = T_high + Arc --dont touch
      -- Angle = -Rotation - Arc  --dont touch
      Angle = -Arc - Rotation 
   else -- cell goes over periodic B/C
      Periodic_flag = 1
      if T_high > Theta_max then
         T2 = T_high - Arc
         T1 = Theta_min
         T_low = T_low
         T_high = Theta_max
         Angle = - Rotation
         Angle2 = -Arc + Rotation
      end
      if T_low < Theta_min then
         T2 = Theta_max
         T1 = T_low + Arc
         T_high = T_high
         T_low = Theta_min
         Angle = - Rotation
         Angle2 = -Arc -Rotation
      end
      --print("Periodic: Angles, T_low, T_high, T1, T2:", T_low,T_high, T1, T2)
      --print("Angle:", Angle, "Angle2:", Angle2)
   end
   --print("Angles, T_low, T_high:", T_low,T_high)
   --print("Angle:", Angle) 

   -- get list of overalpping cells in theta_direction
   D_low_vtx = T_low
   D_high_vtx = T_high
   T_count = 0
   T_faces = {}
   T_cols = {}
   T_weights = {}
   for j = 1,T_up_length do 
      U_low_vtx = T_up_pos_low[j]    
      U_high_vtx = T_up_pos_high[j]
      if ( (U_high_vtx <= D_high_vtx) and (U_low_vtx >= D_low_vtx) ) then -- face fully inside
         T_count = T_count+1
         T_faces[T_count] = T_up_j[j]
         T_cols[T_count] = T_up_col[j]
         T_weights[T_count] = 1.0
         --print("fully inside")
      elseif ( (U_high_vtx > D_high_vtx) and (U_low_vtx < D_high_vtx) ) then -- only overlap at the top
         T_count = T_count+1
         T_faces[T_count] = T_up_j[j]
         T_cols[T_count] = T_up_col[j]
         T_weights[T_count] = (D_high_vtx-U_low_vtx) / (U_high_vtx-U_low_vtx)
         --print("top overlap")
      elseif ( (U_high_vtx > D_low_vtx) and (U_low_vtx < D_low_vtx) ) then -- only overlap at the bottom
         T_count = T_count+1
         T_faces[T_count] = T_up_j[j]
         T_cols[T_count] = T_up_col[j]
         T_weights[T_count] = (U_high_vtx - D_low_vtx) / (U_high_vtx-U_low_vtx)
         --print("bottom overlap")
      elseif ( (U_high_vtx > D_high_vtx) and (U_low_vtx < D_low_vtx) ) then -- overlap at both side
         T_count = T_count+1
         T_faces[T_count] = T_up_j[j]
         T_cols[T_count] = T_up_col[j]
         T_weights[T_count] = (D_high_vtx - D_low_vtx) / (U_high_vtx-U_low_vtx)
         --print("both sides")
      else -- other cases
         --print("No overlap so no action required.") 
      end
   end

   --print("Number of upstream faces overlaping with current downstream face is :",T_count)
   --for j=1,T_count do
   --   print("Face_ind, Col_ind, weight:", T_faces[j],T_cols[j],T_weights[j])
   --end
   
   -- Example 1:   
   --    Col1  Col0 
   -- +------+------+
   -- | blk1 | blk0 | ^
   -- |      |      | | Z,k
   -- +------+------+ |
   --            <----+
   --              T,j

   -- Example 2: 
   -- use row and column to find block index
   --   Col4  Col3  Col2  Col1  Col0 
   -- +-----+-----+-----+-----+-----+
   -- |     |     |     |     |     |
   -- |  12 |  13 |  14 |  15 |  16 |  Row0
   -- |     |     |     |     |     |  ^
   -- +-----+-----+-----+-----+-----+  | Z
   --                           T <----+

   -- area average properties to ensure conservative exchange
   Nz = args.k - 1 -- minus one to account for ghost cells 
   Z_Nfaces = DOWN_Nfaces[Nz]
   Z_face_ind = DOWN_face_ind[Nz]
   Z_row_ind = DOWN_row_ind[Nz]
   Z_weight = DOWN_weight[Nz]
   Tot_area = 0
   rho = 0.; T = 0.; p = 0.
   u_temp = 0.; v_temp = 0.; w = 0.
   tke = 0.; omega =0.
   --print("Number of faces being sample in Z direction:",Z_Nfaces," in T direction:", T_count)
    
   for i=1,Z_Nfaces do
      for j=1,T_count do
         weight = Z_weight[i] * T_weights[j]

         row = UP_row_list[ Z_row_ind[i]+1]  
         b_id = row[T_cols[j]+1]  -- +1 as column number start from 0
    
         -- b_id = upstream_list_row0[T_cols[j]+1]  -- +1 as column number start from 0
         --print("Rotor",b_id,T_cols[j]) 

         --print("Rotor_in Sampling from, blk, i, j, k:",  b_id,2,T_faces[j],Z_face_ind[i])  -- set j = 2 as going along SOUTH face
         --print("weight =", weight)     
         if args.timeStep == 0 then     
            flow = sampleFluidCell(b_id,2,T_faces[j],Z_face_ind[i])
            flow2 = sampleFluidFace("i",b_id,2,T_faces[j],Z_face_ind[i])
            area = flow2.area
         else            
            flow = sampleFluidFace("i",b_id,2,T_faces[j],Z_face_ind[i])
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
      -- repeat above for second segment at other end of range
      T_low = T1
      T_high = T2

      -- get list of overalpping cells in theta_direction
      D_low_vtx = T_low
      D_high_vtx = T_high
      T_count = 0
      T_faces = {}
      T_cols = {}
      T_weights = {}
      for j = 1,T_up_length do 
         U_low_vtx = T_up_pos_low[j]    
         U_high_vtx = T_up_pos_high[j]
         if ( (U_high_vtx <= D_high_vtx) and (U_low_vtx >= D_low_vtx) ) then -- face fully inside
            T_count = T_count+1
            T_faces[T_count] = T_up_j[j]
            T_cols[T_count] = T_up_col[j]
            T_weights[T_count] = 1.0
            --print("fully inside")
         elseif ( (U_high_vtx > D_high_vtx) and (U_low_vtx < D_high_vtx) ) then -- only overlap at the top
            T_count = T_count+1
            T_faces[T_count] = T_up_j[j]
            T_cols[T_count] = T_up_col[j]
            T_weights[T_count] = (D_high_vtx-U_low_vtx) / (U_high_vtx-U_low_vtx)
            --print("top overlap")
         elseif ( (U_high_vtx > D_low_vtx) and (U_low_vtx < D_low_vtx) ) then -- only overlap at the bottom
            T_count = T_count+1
            T_faces[T_count] = T_up_j[j]
            T_cols[T_count] = T_up_col[j]
            T_weights[T_count] = (U_high_vtx - D_low_vtx) / (U_high_vtx-U_low_vtx)
            --print("bottom overlap")
         elseif ( (U_high_vtx > D_high_vtx) and (U_low_vtx < D_low_vtx) ) then -- overlap at both side
            T_count = T_count+1
            T_faces[T_count] = T_up_j[j]
            T_cols[T_count] = T_up_col[j]
            T_weights[T_count] = (D_high_vtx - D_low_vtx) / (U_high_vtx-U_low_vtx)
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
      --   Col4  Col3  Col2  Col1  Col0 
      -- +-----+-----+-----+-----+-----+
      -- |     |     |     |     |     |
      -- |  12 |  13 |  14 |  15 |  16 |  Row0
      -- |     |     |     |     |     |  ^
      -- +-----+-----+-----+-----+-----+  | Z
      --                           T <----+

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

            row = UP_row_list[ Z_row_ind[i] +1]  
            b_id = row[T_cols[j]+1]  -- +1 as column number start from 0

            --b_id = upstream_list_row0[T_cols[j]+1]  -- +1 as column number start from 0

            --print("\n") 
            --print("Rotor_in Periodic sampling from, blk, i, j, k:",  b_id,2,T_faces[j],Z_face_ind[i])
            --print("weight =", weight)
            if args.timeStep == 0 then     
               flow = sampleFluidCell(b_id,2,T_faces[j],Z_face_ind[i])
               flow2 = sampleFluidFace("i",b_id,2,T_faces[j],Z_face_ind[i])
               area = flow2.area
            else            
               flow = sampleFluidFace("i",b_id,2,T_faces[j],Z_face_ind[i])
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


      --print("Velocities2:", u_temp2, v_temp2, u2, v2)
      --print("Angle2:", Angle2)

      rho = (rho*Tot_area + rho2) / (Tot_area+Tot_area2)
      p = (p*Tot_area + p2) / (Tot_area+Tot_area2)
      T = (T*Tot_area + T2) / (Tot_area+Tot_area2)
      u = (u*Tot_area + u2) / (Tot_area+Tot_area2)
      v = (v*Tot_area + v2) / (Tot_area+Tot_area2)
      w = (w*Tot_area + w2) / (Tot_area+Tot_area2)
      tke = (tke*Tot_area + tke2) / (Tot_area+Tot_area2)
      omega = (omega*Tot_area + omega2) / (Tot_area+Tot_area2)
      --print(rho, Tot_weight, rho2, Tot_weight2)

   end


 


   -- remove very small velocities to eliminate uncontrolled run-away growth. There is a feedback loop in mixing plane between cells that don't full cover.
   if math.abs(u) < 1.e-3 then
      u = 0.
   end
   if math.abs(v) < 1.e-3 then
      v = 0.
   end
   if math.abs(w) < 1.e-3 then
      w = 0.
   end

   --print("Rotor_in")
   --print("Output: u,v,w,p,rho,T[0]:",u,v,w,p,rho,T)

   if math.abs(args.csZ) < 1e-10 then
      CSZ = 0
   else
      CSZ = args.csZ
   end 


   FLOW = sampleFluidFace("i",blkId,args.i,args.j,args.k)
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
    Q.T = T
    Q.p = p 
    gmodel:updateThermoFromPT(Q)

   massf = Q.massf       -- mass fractions to be provided as a table
   -- Assemble flux vector
   flux = {}
   flux.mass = rho * (u*args.csX + v*args.csY + w*CSZ) -- kg/s/m**2
   --print("Fmass",F.mass)
   flux.momentum_x = p * args.csX + u * flux.mass
   flux.momentum_y = p * args.csY + v * flux.mass
   flux.momentum_z = p * CSZ + w * flux.mass
   --print("Block,I,J,K",block_id,args.i,args.j,args.k,"Zmomentum:",F.momentum_z, p, args.csZ, CSZ, w)
   flux.total_energy = flux.mass * (Q.u + 0.5*(u*u+v*v+w*w) + p/rho)
   flux.tke = tke * flux.mass
   flux.omega = omega * flux.mass
   --flux.species = {}
   --flux.species[0] = F.mass * massf[0]
   --flux.renergies = {}
   --flux.renergies[0] = F.mass * (Cv*T)

   --print("Block,I,J,K",blkId,args.i,args.j,args.k,"Zmomentum:",flux.momentum_z, p, args.csZ, CSZ, w)
   --print(flux.mass, flux.momentum_x,flux.momentum_y,flux.momentum_z,flux.total_energy)
   --print(u,v,w)


   return flux
end


