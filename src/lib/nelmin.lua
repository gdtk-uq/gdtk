-- nelmin.lua
--
--- @module nelmin
local nelmin = {}

-----
-- Optimise using the Nelder-Mead simplex minimisation of a nonlinear (multivariate) function.
--
-- This code has been adapted from nelmin.py by Peter Jacobs.
-- Consult that file to in turn see where that code was adapted from.
--
-- References:
-- -----------
-- Nelder, J.A. and Mead, R. (1965)
-- A simplex method for function minimization.
-- Computer Journal, Volume 7, pp 308--313.
--
-- O'Neill, R. (1971)
-- Algorithm AS47. Function minimization using a simplex algorithm.
-- Applied Statistics, Volume 20, pp 338--345.
--
--
-- @author Rowan J. Gollan
-- @version 05-Mar-2009
-- @version 30-Jun-2018 made into a module by PJ
-----
function nelmin.minimise(f, x, dx, tol,
                         maxfe, n_check, delta,
                         Kreflect, Kextend, Kcontract)
   -- Some defaults
   tol = tol or 1.0e-6
   maxfe = maxfe or 300
   n_check = n_check or 20
   delta = delta or 0.001
   Kreflect = Kreflect or 1.0
   Kextend = Kextend or 2.0
   Kcontract = Kcontract or 0.5
   --
   local converged = false
   if not dx then
      dx = {}
      for i=1,#x do dx[i] = 0.1 end
   end
   --
   smplx = nelmin.NMSimplex:new{x=x, dx=dx, f=f}
   --
   while ( not converged ) and ( smplx.nfe < maxfe ) do
      -- Take some steps and then check for convergence.
      for i=1,n_check do
	 nelmin.take_a_step(smplx, Kreflect, Kextend, Kcontract)
      end
      -- Pick out the currect best vertex.
      i_best = smplx:lowest()
      x_best = smplx:get_vertex(i_best)
      f_best = smplx.f_list[i_best]
      -- Check the scatter of vertex values to see if
      -- we are close enough to call it quits.
      local mean, stddev = smplx:f_statistics()
      if stddev < tol then
	 -- All of the points are close together, but we need to
	 -- test more carefully.
	 converged = smplx:test_for_minimum(i_best, delta)
	 if not converged then
	    -- The function evaulations are all very close together,
	    -- but we are not at a true minimum; rescale the simplex.
	    smplx:rescale(delta)
	 end
      end
   end
   --
   return x_best, f_best, converged, smplx.nfe, smplx.nrestarts
end -- nelmin.minimise()

local NMSimplex = {}
nelmin.NMSimplex = NMSimplex

local function copy(src)
   dest = {}
   for i=1,#src do
      dest[i] = src[i]
   end
   return dest
end

function NMSimplex:new(o)
   o.N = #o.x
   o.vertex_list = {}
   o.f_list = {}
   o.nfe = 0
   o.nrestarts = 0
   for i=1,o.N+1 do
      local p = copy(o.x)
      if i >= 2 then
	 p[i-1] = p[i-1] + o.dx[i-1]
      end
      o.vertex_list[#o.vertex_list+1] = copy(p)
      o.f_list[#o.f_list+1] = o.f(p)
      o.nfe = o.nfe + 1
   end
   --
   setmetatable(o, self)
   self.__index = self
   return o
end

function NMSimplex:rescale(ratio)
   -- Pick out the current minimum and rebuild the simplex about that point.
   local i_min = self:lowest()
   for i=1,self.N do
      self.dx[i] = self.dx[i]*ratio
   end
   local x = self:get_vertex(i_min)
   self.vertex_list = {}
   self.f_list = {}
   for i=1,self.N+1 do
      local p = copy(x)
      if i >= 2 then
	 p[i-1] = p[i-1] + self.dx[i-1]
      end
      self.vertex_list[#self.vertex_list+1] = copy(p)
      self.f_list[#self.f_list+1] = self.f(p)
      self.nfe = self.nfe + 1
   end
   self.nrestarts = self.nrestarts + 1
end

function NMSimplex:get_vertex(i)
   return copy(self.vertex_list[i])
end

function NMSimplex:replace_vertex(i, x, fvalue)
   self.vertex_list[i] = copy(x)
   self.f_list[i] = fvalue
end

function NMSimplex:lowest(exclude)
   -- Returns the index of the lowest vertex, excluding the one specified.
   exclude = exclude or -1
   local indx
   if exclude == 1 then
      indx = 2
   else
      indx = 1
   end
   local lowest_f_value = self.f_list[indx]
   for i=1,self.N+1 do
      if i ~= exclude then
	 if self.f_list[i] < lowest_f_value then
	    lowest_f_value = self.f_list[i]
	    indx = i
	 end
      end
   end
   --
   return indx
end

function NMSimplex:highest(exclude)
   --  Returns the index of the highest vertex, excluding the one specified.
   exclude = exclude or -1
   if exclude == 1 then
      indx = 2
   else
      indx = 1
   end
   local highest_f_value = self.f_list[indx]
   for i=1,self.N+1 do
      if i ~= exclude then
	 if self.f_list[i] > highest_f_value then
	    highest_f_value = self.f_list[i]
	    indx = i
	 end
      end
   end
   --
   return indx
end

function NMSimplex:f_statistics()
   -- Returns mean and standard deviation of the vertex fn values.
   local sum = 0.0
   for i=1,self.N+1 do
      sum = sum + self.f_list[i]
   end
   local mean = sum / (self.N + 1)
   local sum = 0.0
   for i=1,self.N+1 do
      diff = self.f_list[i] - mean
      sum = sum + diff^2
   end
   local stddev = math.sqrt(sum / self.N)
   --
   return mean, stddev
end

function NMSimplex:centroid(exclude)
   -- Returns the centroid of all vertices excluding the one specified.
   exclude = exclude or -1
   local xmid = {}
   for i=1,self.N do xmid[i] = 0.0 end
   for i=1,self.N+1 do
      if i ~= exclude then
	 for j=1,self.N do
	    xmid[j] = xmid[j] + self.vertex_list[i][j]
	 end
      end
   end
   for j=1,self.N do
      xmid[j] = xmid[j] / self.N
   end
   return xmid
end

function NMSimplex:contract_about_one_point(i_con)
   -- Contract the simplex about the vertex i_con.
   local p_con = self.vertex_list[i_con]
   for i=1,self.N+1 do
      if i ~= i_con then
	 local p = self.vertex_list[i]
	 for j=1,self.N do
	    p[j] = 0.5 * (p[j] + p_con[j])
	 end
	 self.f_list[i] = self.f(p)
	 self.nfe = self.nfe + 1
      end
   end
   --
end

function NMSimplex:test_for_minimum(i_min, delta)
   -- Perturb the minimum vertex and check that it is a local minimum.
   local is_minimum = true -- Assume minimum and test for failure
   local f_min = self.f_list[i_min]
   for j=1,self.N do
      -- Check either side of the minimum, perturbing one
      -- ordinate at a time.
      local p = self:get_vertex(i_min)
      p[j] = p[j] + self.dx[j] * delta
      local f_p = self.f(p)
      self.nfe = self.nfe + 1
      if f_p < f_min then
	 is_minimum = false
	 break
      end
      p[j] = p[j] - self.dx[j] * delta * 2
      f_p = self.f(p)
      self.nfe = self.nfe + 1
      if f_p < f_min then
	 is_minimum = false
	 break
      end
   end
   --
   return is_minimum
end

function nelmin.take_a_step(smplx, Kreflect, Kextend, Kcontract)
   -- Try to move away from the worst point in the simplex.
   -- The new point will be inserted into the simplex (in plane).
   --
   local i_low = smplx:lowest()
   local i_high = smplx:highest()
   local x_high = smplx.vertex_list[i_high] -- this only points
   local f_high = smplx.f_list[i_high]
   -- Centroid of simplex, excluding worst point
   local x_mid = smplx:centroid(i_high)
   local f_mid = smplx.f(x_mid)
   smplx.nfe = smplx.nfe + 1
   --
   -- First, try moving away from worst point by
   -- reflection through centroid
   local x_refl = nelmin.create_new_point(1+Kreflect, x_mid, -Kreflect, x_high)
   local f_refl = smplx.f(x_refl)
   smplx.nfe = smplx.nfe + 1
   if f_refl < f_mid then
      -- The reflection through the centroid is good,
      -- try to extend in the same direction.
      local x_ext = nelmin.create_new_point(Kextend, x_refl, 1-Kextend, x_mid)
      local f_ext = smplx.f(x_ext)
      smplx.nfe = smplx.nfe + 1
      if f_ext < f_refl then
	 -- Keep the extension because it's best.
	 smplx:replace_vertex(i_high, x_ext, f_ext)
      else
	 -- Settle for the original reflection.
	 smplx:replace_vertex(i_high, x_refl, f_refl)
      end

   else
      -- The reflection is not going in the right direction, it seems.
      -- See how many vertices are better than the reflected point.
      local count = 0
      for i=1,smplx.N+1 do
	 if smplx.f_list[i] > f_refl then count = count + 1 end
      end
      if count <= 1 then
	 -- Not too many points are higher than the original reflection.
	 -- Try a contraction on the reflection-side of the centroid.
	 local x_con = nelmin.create_new_point(1-Kcontract, x_mid, Kcontract, x_high)
	 local f_con = smplx.f(x_con)
	 smplx.nfe = smplx.nfe + 1
	 if f_con < f_high then
	    -- At least we haven't gone uphill; accept.
	    smplx:replace_vertex(i_high, x_con, f_con)
	 else
	    -- We have not been successful in taking a single step.
	    -- Contract the simplex about the current lowest point.
	    smplx:contract_about_one_point(i_low)
	 end
      else
	 -- Retain the original reflection because there are many
	 -- vertices with higher values of the objective function.
	 smplx:replace_vertex(i_high, x_refl, f_refl)
      end
   end
end

function nelmin.create_new_point(c1, p1, c2, p2)
   -- Create a new N-dimensional point as a weighting of points p1 and p2.
   local p_new = {}
   for j=1,#p1 do
      p_new[#p_new+1] = c1*p1[j] + c2*p2[j]
   end
   return p_new
end

return nelmin
