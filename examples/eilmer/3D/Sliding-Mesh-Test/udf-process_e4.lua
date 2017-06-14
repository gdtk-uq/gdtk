-- udf-process.lua
-- This file sets up functions that will be called
-- from the main time-stepping loop.

-- Version 3 is based on vtx locations and implements a conservative scheme. I.e. overalpping areas are considered to set correct fluxes.


function atTimestepStart(sim_time,steps,gasBlocks)
   --print("in supervisor function --> atTimeStepStart")

   -- Rport the rotational angle at regular intervals. 
   if (steps % 50) == 0 then
      -- load DATA  
      DATA,err = table.load("CONFIG_tbl.lua")
      N_BLADE = DATA[0]
      OMEGA = DATA[1]
      UP_row_list = DATA[2]
      DOWN_row_list = DATA[3]
      LIST = DATA[4]
      theta_min = DATA[11]
      theta_max = DATA[12]

      print("Rotatation Angle = ", sim_time*OMEGA, " rad or ", sim_time*OMEGA/math.pi*180, " deg")
   end


   --print("Available functions and variables")
   --for k,v in pairs(_G) do
   --     print("Global key:", k,"value:",v)
   --end

   --A = sampleFlow(0,1,1,1)
   --print("Contents of sampleFlow")
   --for k,v in pairs(A) do
   --     print("Global key:", k,"value:",v)
   -- end

   --A = sampleFace("i",0,1,1,1)
   --print("Contents of sampleFace")
   --for k,v in pairs(A) do
   --     print("Global key:", k,"value:",v)
   -- end


   -- The followign code does some pre-allocation work for the sliding mesh interfce.       
   if steps == 0 then

      -- load DATA  
      DATA,err = table.load("CONFIG_tbl.lua")
      N_BLADE = DATA[0]
      OMEGA = DATA[1]
      UP_row_list = DATA[2]
      DOWN_row_list = DATA[3]
      LIST = DATA[4]
      theta_min = DATA[11]
      theta_max = DATA[12]

      -- check that blks are in current block list. 
      -- This ensures code is only executed by correct MPI process
      flag = 0
      for i,block in pairs(blockData) do
        if i == 0 then
            flag = 1
            break
         end
      end

    -- Must only be run by MPI process that has interface.
      if flag == 1 then

          -- check that all interface blocks are present on current Node
          A = {}
          for i,block in pairs(blockData) do
              table.insert(A,i)
          end
          for k,b in pairs(LIST) do
            for i,block in pairs(blockData) do
                if i == b then
                   break
                end
                --if i == nblks-1 then
                --    print('ERROR: Not all blocks that are part of the sliding interface appear on same MPI Node.')
                --    print('Blocks in interface: ', LIST)
                --    print('Blocks on current Node: ', A)
                --end
            end
          end

          BLK_MAP = {}
          -- Create Mapping table that allows mapping from global index to local indices.
          for i,block in pairs(blockData) do

             for k,ind in pairs(LIST) do
                if i == ind then
                   BLK_MAP[ind] = i
                   print("Mapped Global Block ", ind, " to local:",i) 
                else
                   print("Global Block:",i, "Not mapped")       
                end
             end
          end

          -- upstream face
          -- Example 1:
          -- Looking Radially outwards
          --    Col1  Col0 
          -- +------+------+
          -- | blk1 | blk0 | ^
          -- |      |      | | Z,k
          -- +------+------+ |
          --            <----+
          --              T,j
          -- Example 2:
          -- Looking Radially outwards
          --   Col4  Col3  Col2  Col1  Col0 
          -- +-----+-----+-----+-----+-----+
          -- |     |     |     |     |     |
          -- |  12 |  13 |  14 |  15 |  16 |  Row0
          -- |     |     |     |     |     |  ^
          -- +-----+-----+-----+-----+-----+  | Z,k
          --                           T <----+
          --                               j


          print("\n")
          print("Getting Z-positions of vtx on upstream face")
          Z_up_pos_low = {}
          Z_up_pos_high = {}
          Z_up_k = {}
          Z_up_row = {}
          Z_up_area = {}
          indx = 1
          row = -1
          for _, ROW in pairs(UP_row_list) do -- get Z position
             global_id = ROW[1]  
             blk_id = BLK_MAP[global_id]
             print("Block currently being explored. blk_id =",blk_id)
             row = row + 1 
             imin = blockData[blk_id].vtxImin; imax = blockData[blk_id].vtxImax
             jmin = blockData[blk_id].vtxJmin; jmax = blockData[blk_id].vtxJmax
             kmin = blockData[blk_id].vtxKmin; kmax = blockData[blk_id].vtxKmax
             -- print(imin,imax,jmin,jamx,kmin,kmax) 
             i = imin -- Summing up across WEST face
             j = jmin -- only need to follow one line
             print(kmin,kmax)
             for k=kmin,kmax-1 do  
                vtx3 = getVtxPosition(blk_id, i, j+1, k)
                vtx0 = getVtxPosition(blk_id, i, j, k)
                vtx4 = getVtxPosition(blk_id, i, j, k+1)
                vtx7 = getVtxPosition(blk_id, i, j+1, k+1)
                area = face_area(vtx0,vtx3,vtx4,vtx7)
                Z_up_pos_low[indx] = vtx0.z
                Z_up_pos_high[indx] = vtx4.z
                Z_up_k[indx] = k
                Z_up_row[indx] = row
                Z_up_area[indx] = area
                print("check creation of Z_up list: pos_low, pos_high, k, row:",Z_up_pos_low[indx],Z_up_pos_high[indx],Z_up_k[indx],Z_up_row[indx])
                indx = indx + 1
             end

          end 
          Z_up_length = indx-1 -- subtract +1 from last loop
          --store number of items in list for later access  

          print("\n")
          print("Getting Theta-positions of vtx on upstream face")
          T_up_pos_low={}
          T_up_pos_high={}
          T_up_j={}
          T_up_col={}
          T_up_area={} 
          indx = 1
          col = -1
          List = UP_row_list[1]
          for _,global_id in pairs(List) do -- get Theta position
             blk_id = BLK_MAP[global_id]
             print("Block currently being explored. blk_id =",blk_id)
             col = col+1 
             imin = blockData[blk_id].vtxImin; imax = blockData[blk_id].vtxImax
             jmin = blockData[blk_id].vtxJmin; jmax = blockData[blk_id].vtxJmax
             kmin = blockData[blk_id].vtxKmin; kmax = blockData[blk_id].vtxKmax
             i = imin -- Summing up across WEST face
             k = kmin
             for j=jmin,jmax-1 do  
                vtx3 = getVtxPosition(blk_id, i, j+1, k)
                vtx0 = getVtxPosition(blk_id, i, j, k)
                vtx4 = getVtxPosition(blk_id, i, j, k+1)
                vtx7 = getVtxPosition(blk_id, i, j+1, k+1)

                area = face_area(vtx0,vtx3,vtx4,vtx7)

                T_up_pos_low[indx] = math.atan(vtx0.y/vtx0.x)
                T_up_pos_high[indx] = math.atan(vtx3.y/vtx3.x)
                T_up_j[indx] = j
                T_up_col[indx] = col
                T_up_area[indx] = area
                print("check creation of T_up list: angle_low, angle_high, j, col:",T_up_pos_low[indx],T_up_pos_high[indx],T_up_j[indx],T_up_col[indx])
                indx = indx + 1  
             end 
          end
          T_up_length = indx - 1 -- subtract +1 from last loop


          -- downstream face
          -- Example 1:
          -- Looking Radially outwards
          --      Col0 
          -- +-------------+
          -- |    blk3     |  
          -- |             |  
          -- +-------------+     
          -- |    blk2     |  ^
          -- |             |  | Z,k
          -- +-------------+  |   
          --            <-----+ 
          --              T,i
          --
          -- Example 2:
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

          print("\n")
          print("Getting Z-positions of vtx on downstream face")
          Z_down_pos_low = {}
          Z_down_pos_high = {}
          Z_down_k = {}
          Z_down_row = {}
          Z_down_area = {}
          indx = 1
          row = - 1
          for _, ROW in pairs(DOWN_row_list) do -- get Z position
             global_id = ROW[1] 
             blk_id = BLK_MAP[global_id]
             print("Block currently being explored. blk_id =",blk_id)
             row = row + 1 
             imin = blockData[blk_id].vtxImin; imax = blockData[blk_id].vtxImax
             jmin = blockData[blk_id].vtxJmin; jmax = blockData[blk_id].vtxJmax
             kmin = blockData[blk_id].vtxKmin; kmax = blockData[blk_id].vtxKmax
             i = imin -- Summing up across WEST face
             j = jmin -- only need to follow one line
             for k=kmin,kmax-1 do  
                vtx0 = getVtxPosition(blk_id, i, j, k)
                vtx1 = getVtxPosition(blk_id, i+1, j, k)
                vtx5 = getVtxPosition(blk_id, i+1, j, k+1)
                vtx4 = getVtxPosition(blk_id, i, j, k+1)

                area = face_area(vtx0,vtx3,vtx4,vtx7)

                Z_down_pos_low[indx] = vtx0.z
                Z_down_pos_high[indx] = vtx4.z
                Z_down_k[indx] = k
                Z_down_row[indx] = row
                Z_down_area[indx] = area
                print("check creation of Zdown list: pos_low, pos_high, k, row:",Z_down_pos_low[indx],Z_down_pos_high[indx],Z_down_k[indx],Z_down_row[indx])
                indx = indx + 1
             end
          end 
          Z_down_length = indx-1 -- subtract +1 from last loop
          --store nuber of items in list for later access  

          print("\n")
          print("Getting Theta-positions of vtx on downstream face")
          T_down_pos_low = {}
          T_down_pos_high = {}
          T_down_i = {}
          T_down_col = {}
          T_down_area = {}
          indx = 1
          col = -1
          List = DOWN_row_list[1]
          for _,global_id in pairs(List) do -- get Theta position
             blk_id = BLK_MAP[global_id]
             print("Block currently being explored. blk_id =",blk_id)
             col = col+1 
             imin = blockData[blk_id].vtxImin; imax = blockData[blk_id].vtxImax
             jmin = blockData[blk_id].vtxJmin; jmax = blockData[blk_id].vtxJmax
             kmin = blockData[blk_id].vtxKmin; kmax = blockData[blk_id].vtxKmax
             i = imin -- Summing up across WEST face
             k = kmin
             for i=imin,imax-1 do  
                vtx0 = getVtxPosition(blk_id, i, j, k)
                vtx1 = getVtxPosition(blk_id, i+1, j, k)
                vtx5 = getVtxPosition(blk_id, i+1, j, k+1)
                vtx4 = getVtxPosition(blk_id, i, j, k+1)

                area = face_area(vtx0,vtx3,vtx4,vtx7)

                T_down_pos_low[indx] = math.atan(vtx0.y/vtx0.x)
                T_down_pos_high[indx] = math.atan(vtx1.y/vtx1.x)
                T_down_i[indx] = i
                T_down_col[indx] = col
                T_down_area[indx] = area
                print("check creation of Tdown list: ang_low, ang_high, j, col:",T_down_pos_low[indx],T_down_pos_high[indx],T_down_i[indx],T_down_col[indx])
                indx = indx + 1    
             end 
          end
          T_down_length = indx-1 -- subtract +1 from last loop
          --store number of items in list for later access 


          print("\n")
          print("UP UP UP")
          -- create map for z direction (for mapping DOWN --> UP)
          UP_Nfaces = {}
          UP_face_ind = {}
          UP_row_ind = {}
          UP_weight = {}
          UP_area = {}
          for i=1,Z_up_length do
             U_low_vtx = Z_up_pos_low[i]
             U_high_vtx = Z_up_pos_high[i]
             count = 0
             faces = {}
             rows = {}
             weights = {}
             areas = {}
             for j = 1,Z_down_length do 
                D_low_vtx = Z_down_pos_low[j]    
                D_high_vtx = Z_down_pos_high[j]
                if ( (D_high_vtx <= U_high_vtx) and (D_low_vtx >= U_low_vtx) ) then -- face fully inside
                    count = count+1
                    faces[count] = Z_down_k[j]
                    rows[count] = Z_down_row[j]
                    weights[count] = 1.0
                    areas[count] = Z_down_area[j]
                    --print("fully inside")
                elseif ( (D_high_vtx > U_high_vtx) and (D_low_vtx < U_high_vtx) ) then -- only overlap at the top
                    count = count+1
                    faces[count] = Z_down_k[j]
                    rows[count] = Z_down_row[j]
                    weights[count] = (U_high_vtx-D_low_vtx) / (D_high_vtx-D_low_vtx)
                    areas[count] = Z_down_area[j]
                    --print("top overlap")
                elseif ( (D_high_vtx > U_low_vtx) and (D_low_vtx < U_low_vtx) ) then -- only overlap at the bottom
                    count = count+1
                    faces[count] = Z_down_k[j]
                    rows[count] = Z_down_row[j]
                    weights[count] = (D_high_vtx - U_low_vtx) / (D_high_vtx-D_low_vtx)
                    areas[count] = Z_down_area[j]
                    --print("bottom overlap")
                elseif ( (D_high_vtx > U_high_vtx) and (D_low_vtx < U_low_vtx) ) then -- overlap at both side
                    count = count+1
                    faces[count] = Z_down_k[j]
                    rows[count] = Z_down_row[j]
                    weights[count] = (U_high_vtx - U_low_vtx) / (D_high_vtx-D_low_vtx)
                    areas[count] = Z_down_area[j]
                    --print("both sides")
                else -- other cases
                    -- No overlap so no action required. 
                end
             end
             UP_Nfaces[i] = count
             UP_face_ind[i] = faces
             UP_row_ind[i] = rows
             UP_weight[i] = weights 
             UP_area[i] = areas        
          end

          for i=1,Z_up_length do
             print("Number of downstream faces overlaping with upstream face ",  i," is :",UP_Nfaces[i])
             for j=1,UP_Nfaces[i] do
                print("Face_ind, Row_ind, weight:", UP_face_ind[i][j],UP_row_ind[i][j],UP_weight[i][j])
             end
          end

          UP = {}
          UP[0] = UP_Nfaces
          UP[1] = UP_face_ind
          UP[2] = UP_row_ind
          UP[3] = UP_weight
          UP[4] = UP_area
          UP[5] = T_down_pos_low  
          UP[6] = T_down_pos_high
          UP[7] = T_down_col
          UP[8] = T_down_i
          UP[9] = T_down_length

          --print(T_down_i[0],T_down_i[1],T_down_i[2],T_down_i[3],T_down_i[4],T_down_i[5])
          
          print("\n")
          print("DOWN DOWN DOWN")
          -- create map for z direction (for mapping UP --> DOWN)
          DOWN_Nfaces = {}
          DOWN_face_ind = {}
          DOWN_row_ind = {}
          DOWN_weight = {}
          DOWN_area = {}
          for i=1,Z_down_length do
             D_low_vtx = Z_down_pos_low[i]
             D_high_vtx = Z_down_pos_high[i]
             count = 0
             faces = {}
             rows = {}
             weights = {}
             areas = {} 
             for j = 1,Z_up_length do 
                U_low_vtx = Z_up_pos_low[j]    
                U_high_vtx = Z_up_pos_high[j]
                if ( (U_high_vtx <= D_high_vtx) and (U_low_vtx >= D_low_vtx) ) then -- face fully inside
                    count = count+1
                    faces[count] = Z_up_k[j]
                    rows[count] = Z_up_row[j]
                    weights[count] = 1.0
                    areas[count] = Z_up_area[j]
                    --print("fully inside")
                elseif ( (U_high_vtx > D_high_vtx) and (U_low_vtx < D_high_vtx) ) then -- only overlap at the top
                    count = count+1
                    faces[count] = Z_up_k[j]
                    rows[count] = Z_up_row[j]
                    weights[count] = (D_high_vtx-U_low_vtx) / (U_high_vtx-U_low_vtx)
                    areas[count] = Z_up_area[j]
                    --print("top overlap")
                elseif ( (U_high_vtx > D_low_vtx) and (U_low_vtx < D_low_vtx) ) then -- only overlap at the bottom
                    count = count+1
                    faces[count] = Z_up_k[j]
                    rows[count] = Z_up_row[j]
                    weights[count] = (U_high_vtx - D_low_vtx) / (U_high_vtx-U_low_vtx)
                    areas[count] = Z_up_area[j]
                    --print("bottom overlap")
                elseif ( (U_high_vtx > D_high_vtx) and (U_low_vtx < D_low_vtx) ) then -- overlap at both side
                    count = count+1
                    faces[count] = Z_up_k[j]
                    rows[count] = Z_up_row[j]
                    weights[count] = (D_high_vtx - D_low_vtx) / (U_high_vtx-U_low_vtx)
                    areas[count] = Z_up_area[j]
                    --print("both sides")
                else -- other cases
                    -- No overlap so no action required. 
                end
             end
             DOWN_Nfaces[i] = count
             DOWN_face_ind[i] = faces
             DOWN_row_ind[i] = rows
             DOWN_weight[i] = weights 
             DOWN_area[i] = areas        
          end

          for i=1,Z_down_length do
             print("Number of upstream faces overlaping with downstream face ",  i," is :",DOWN_Nfaces[i])
             for j=1,DOWN_Nfaces[i] do
                print("Face_ind, Row_ind, weight:", DOWN_face_ind[i][j],DOWN_row_ind[i][j],DOWN_weight[i][j])
             end
          end


          DOWN = {}
          DOWN[0] = DOWN_Nfaces
          DOWN[1] = DOWN_face_ind
          DOWN[2] = DOWN_row_ind
          DOWN[3] = DOWN_weight
          DOWN[4] = DOWN_area
          DOWN[5] = T_up_pos_low 
          DOWN[6] = T_up_pos_high
          DOWN[7] = T_up_col
          DOWN[8] = T_up_j
          DOWN[9] = T_up_length
          --print(T_up_col[1],T_up_col[2],T_up_col[3])
     
          -- write mapping tables to files  
          assert( table.save( UP, "UP_tbl.lua" ) == nil )
          -- write mapping tables to files  
          assert( table.save( DOWN, "DOWN_tbl.lua" ) == nil )

          -- PRINT NOTICE THAT FILE HAS BEEN RUN
          print("#################################################################")
          print("###       udf-process_e4.lua  Successfully Completed          ###")
          print("#################################################################")

     
       end
   end

   return
end

function face_area(vtx0,vtx1,vtx2,vtx3,Area)
-- Function to evaluate the are of a face defined by 4 vertices

    -- Center of rectangle
    Cx = 0.25* (vtx0.x + vtx1.x + vtx2.x + vtx3.x) 
    Cy = 0.25* (vtx0.y + vtx1.y + vtx2.y + vtx3.y) 
    Cz = 0.25* (vtx0.z + vtx1.z + vtx2.z + vtx3.z)

    -- a01
    a01x = (vtx0.y-Cy)*(vtx1.z-Cz) - (vtx0.z-Cz)*(vtx1.y-Cy) 
    a01y = (vtx0.x-Cx)*(vtx1.z-Cz) - (vtx0.z-Cz)*(vtx1.x-Cx)   
    a01z = (vtx0.x-Cx)*(vtx1.y-Cy) - (vtx0.y-Cy)*(vtx1.x-Cx)  
    a01 = math.sqrt( a01x*a01x + a01y*a01y + a01z*a01z)

    -- a12
    a12x = (vtx1.y-Cy)*(vtx2.z-Cz) - (vtx1.z-Cz)*(vtx2.y-Cy) 
    a12y = (vtx1.x-Cx)*(vtx2.z-Cz) - (vtx1.z-Cz)*(vtx2.x-Cx)   
    a12z = (vtx1.x-Cx)*(vtx2.y-Cy) - (vtx1.y-Cy)*(vtx2.x-Cx)  
    a12 = math.sqrt( a12x*a12x + a12y*a12y + a12z*a12z)

    -- a23
    a23x = (vtx2.y-Cy)*(vtx3.z-Cz) - (vtx2.z-Cz)*(vtx3.y-Cy) 
    a23y = (vtx2.x-Cx)*(vtx3.z-Cz) - (vtx2.z-Cz)*(vtx3.x-Cx)   
    a23z = (vtx2.x-Cx)*(vtx3.y-Cy) - (vtx2.y-Cy)*(vtx3.x-Cx)  
    a23 = math.sqrt( a23x*a23x + a23y*a23y + a23z*a23z)

    -- a30
    a30x = (vtx3.y-Cy)*(vtx0.z-Cz) - (vtx3.z-Cz)*(vtx0.y-Cy) 
    a30y = (vtx3.x-Cx)*(vtx0.z-Cz) - (vtx3.z-Cz)*(vtx0.x-Cx)   
    a30z = (vtx3.x-Cx)*(vtx0.y-Cy) - (vtx3.y-Cy)*(vtx0.x-Cx)  
    a30 = math.sqrt( a30x*a30x + a30y*a30y + a30z*a30z)

    Area = a01 + a12 + a23 + a30
    --print("A", Area)
    return Area
end

--function at_timestep_end(args)
--   --if (args.step % 100) == 0 then
--   --   -- Run every 100 timesteps
--   --   print("At end of timestep ", args.step, " t=", args.t)
--   --   sum_total_mass(args)
--   --end
--   --print("At end of timestep")
--
--
--   return
--end

function sum_total_mass(args)
   --print("running mass summing fucntion")
   mass = 0.0
   for ib=0,(nblks-1) do
      imin = blks[ib].imin; imax = blks[ib].imax
      jmin = blks[ib].jmin; jmax = blks[ib].jmax
      kmin = blks[ib].kmin; kmax = blks[ib].kmax
      --print(imin, imax,jmin,jmax,kmin,kmax)
      blk_id = blks[ib].id
      for k=kmin,kmax do
         for j=jmin,jmax do
            for i=imin,imax do
               cell = sample_flow(blk_id, i, j, k)
               --print( blk_id, i, j, k)
               -- We are only given p and T
               -- so need to compute density
               -- using gas model
               Q = create_empty_gas_table()
               Q.p = cell.p
               Q.T = cell.T
               for isp=0,(nsp-1) do Q.massf[isp] = cell.massf[isp] end
               print("Here?",Q.p,Q.T[0])
               eval_thermo_state_pT(Q)
               rho = Q.rho
               -- Now we can compute mass in cell using volume of cell
               mass = mass + rho*cell.vol
            end
         end
      end
   end
   print("At timestep, ",args.step,", Mass (kg) of gas in domain: ", mass)
   return
end

function source_vector(args, cell_data)
   -- args contains t
   -- cell_data table contains most else
   src = {}
   src.mass = 0.0
   src.momemtum_x = 0.0
   src.momentum_y = 0.0
   src.momentum_z = 0.0
   if cell_data.x > 0.2 and cell_data.x < 0.4 and 
      cell_data.y > 0.2 and cell_data.y < 0.3 then
      src.total_energy = 100.0e+6  -- J/m**3
      src.energies = {[0]=100.0e6}
   else
      src.total_energy = 0.0
      src.energies = {[0]=0.0}
   end
   src.species = {[0]=0.0,}
   return src
end

