-- Author: Flynn Hack
-- Date: 2022-02-14

function assignVtxVelocities(t, dt)

   -- ==============================================================
   -- Pull our userPad values to be used for grid motion
   -- ==============================================================    
   
   for _, blk_id in ipairs(localFluidBlockIds) do
      omegab_g = Vector3:new{x = 0, y = 0, z = userPad[2]}
      steps = userPad[3]      
       
      -- body position and translational velocity not changing.
      xb_g = Vector3:new{x=0,y=0,z=0}
      xdotb_g = Vector3:new{x=0,y=0,z=0}
        
      -- Key function here is:

      setVtxVelocitiesForRigidBlock(blk_id, omegab_g, xb_g, xdotb_g)
        
      -- Further functions for different cases (stretching etc) can be found
      -- in gdtk/src/eilmer/grid_motion_udf.d

   end    
end





