-- estimate_shock_angle.lua
-- Invoke with the command line:
-- $ e4shared --custom-post --script-file=estimate_shock_angle.lua
-- PJ, 2015-10-20
--
print("Begin estimate_shock_angle")
nb = 2
fsol = FlowSolution:new{jobName="cone20", dir=".", tindx=4, nBlocks=nb}
print("fsol=", fsol)

function locate_shock_along_strip()
   local p_max = ps[1]
   for i = 2, #ps do
      p_max = math.max(ps[i], p_max)
   end
   local p_trigger = ps[1] + 0.3 * (p_max - ps[1])
   local x_old = xs[1]; local y_old = ys[1]; local p_old = ps[1]
   local x_new = x_old; local y_new = y_old; local p_new = p_old
   for i = 2, #ps do
      x_new = xs[i]; y_new = ys[i]; p_new = ps[i]
      if p_new > p_trigger then break end
      x_old = x_new; y_old = y_new; p_old = p_new
   end
   local frac = (p_trigger - p_old) / (p_new - p_old)
   x_loc = x_old * (1.0 - frac) + x_new * frac
   y_loc = y_old * (1.0 - frac) + y_new * frac
   return
end

xshock = {}; yshock = {}
local nj = fsol:get_njc(0)
for j = 0, nj-1 do
   xs = {}; ys = {}; ps = {}
   for ib = 0, nb-1 do
      local ni = fsol:get_nic(ib)
      for i = 0, ni-1 do
	 cellData = fsol:get_cell_data{ib=ib, i=i, j=j}
	 xs[#xs+1] = cellData["pos.x"]
	 ys[#ys+1] = cellData["pos.y"]
	 ps[#ps+1] = cellData["p"]
      end
   end
   locate_shock_along_strip()
   if x_loc < 0.9 then
      -- Keep only the good part of the shock.
      xshock[#xshock+1] = x_loc
      yshock[#yshock+1] = y_loc
   end
end

-- Least-squares fit of a straight line for the shock
-- Model is y = alpha0 + alpha1 * x
sum_x = 0.0; sum_y = 0.0; sum_x2 = 0.0; sum_xy = 0.0 
for j = 1, #xshock do
   sum_x = sum_x + xshock[j]
   sum_x2 = sum_x2 + xshock[j]*xshock[j]
   sum_y = sum_y + yshock[j]
   sum_xy = sum_xy + xshock[j]*yshock[j]
end
N = #xshock
alpha1 = (sum_xy/N - sum_x/N * sum_y/N) / (sum_x2/N - sum_x/N * sum_x/N)
alpha0 = sum_y/N - alpha1 * sum_x/N
shock_angle = math.atan(alpha1)
sum_y_error = 0.0
for j = 1, N do
   sum_y_error = sum_y_error + math.abs((alpha0 + alpha1 * xshock[j]) - yshock[j])
end
print("shock_angle_deg=", shock_angle*180.0/math.pi)
print("average_deviation_metres=", sum_y_error/N)
