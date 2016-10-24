-- locate-bow-shock.lua
-- Invoke with the command line:
-- $ e4shared --custom-post --script-file=locate-bow-shock.lua
--
-- PJ, 2016-10-24, updated for Eilmer4
--
print("Locate a bow shock by its pressure jump.")
print("Start by reading full flow solution.")
fsol = FlowSolution:new{jobName="cyl", dir=".", tindx=10, nBlocks=4}
print("fsol=", fsol)

function locate_shock_along_strip()
   local p_max = ps[1]
   for i = 2, #ps do
      p_max = math.max(ps[i], p_max)
   end
   local p_trigger = ps[1] + 0.3 * (p_max - ps[1])
   local x_old = xs[1]; local p_old = ps[1]
   local x_new = x_old; local p_new = p_old
   for i = 2, #ps do
      x_new = xs[i]; p_new = ps[i]
      if p_new > p_trigger then break end
      x_old = x_new; p_old = p_new
   end
   local frac = (p_trigger - p_old) / (p_new - p_old)
   x_loc = x_old * (1.0 - frac) + x_new * frac
   return
end

-- Since this is a 3D simulation, the shock is not expected
-- to be flat in the k-direction (along the cylinder axis).
-- Sample the shock layer in a few places near the stagnation line.
-- Block 0 contains the stagnation point and the bottom surface is
-- the plane of symmetry that cuts the cylinder half-way along its axis.
-- The south boundary is the plane of symmetry that cuts the cylinder
-- along its axis. Supersonic flow comes in from the west boundary
-- and exits from the north boundary. The east boundary is the
-- cylinder surface.
local xshock = {}; local yshock = {}
local ib = 0
local nk = fsol:get_nkc(0)
for k = 0, nk-1 do
   xs = {}; ys = {}; ps = {}
   local j = 0
   local ni = fsol:get_nic(ib)
   for i = 0, ni-1 do
      cellData = fsol:get_cell_data{ib=ib, i=i, j=j, k=k}
      xs[#xs+1] = cellData["pos.x"]
      ps[#ps+1] = cellData["p"]
   end
   locate_shock_along_strip()
   xshock[#xshock+1] = x_loc
   yshock[#yshock+1] = y_loc
   if #xshock >= 6 then break end
end

x_sum = 0.0
for _,x in ipairs(xshock) do x_sum = x_sum + x end
x_average = x_sum / #xshock
print("Average x-location=", x_average)
D = 15.0e-3 -- cylinder diameter
delta = -x_average - D/2
print("shock displacement=", delta*1000.0, "mm")
print("delta/R=", delta/(D/2))


