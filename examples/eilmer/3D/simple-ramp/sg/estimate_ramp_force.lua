-- estimate_ramp_force.py
-- Example postprocessing script to look at the data along the ramp
-- and compute some potentially useful information.
-- Invoke with the command line:
-- $ e4shared --custom-post --script-file=estimate_ramp_force.lua
-- PJ, 2015-10-26 adapted from Eilmer3 version
--
print("Begin estimate force on ramp surface")
nb = 2
fsol = FlowSolution:new{jobName="ramp", dir=".", tindx=5, nBlocks=nb}
print("fsol=", fsol)

-- Integrate the pressure force over the BOTTOM surface of the block.
local force = Vector3:new(0.0, 0.0, 0.0)
local ib = 1
for i = 0, fsol:get_nic(ib)-1 do
   for j = 0, fsol:get_njc(ib)-1 do
      -- The bottom cell face has p0, p1, p2, p3 as corners.
      face = quadProperties{p0=fsol:get_vtx{ib=ib, i=i, j=j, k=0},
			    p1=fsol:get_vtx{ib=ib, i=i+1, j=j, k=0},
			    p2=fsol:get_vtx{ib=ib, i=i+1, j=j+1, k=0},
			    p3=fsol:get_vtx{ib=ib, i=i, j=j+1, k=0}}
      cellData = fsol:get_cell_data{ib=ib, i=i, j=j, k=0}
      df = face.area * cellData.p * face.n
      force = force - df 
      -- negative because the unit normal of this cell face is into the volume
   end
end
print("force=", force, "Newtons")

print("done.")
