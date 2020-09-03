-- post.lua
-- A place to put helper functions for custom post-processing activities.
--
print("Loading post.lua...")

FlowSolution.get_cell_data_along_i_index = function(self, args)
   -- Return a table of CellData tables, for cells along an i-index slice
   -- for specified j,k indices of a structured-grid.
   if not (type(self)=='userdata') then
      error("Make sure that you are using FlowSolution:get_xxxx{args} not FlowSolution.get_xxxx{args}", 2)
   end
   if (not (args.ib and args.j)) then
      error("FlowSolution:get_cell_data_along_i_index{ib=n1, j=n2, k=n3}", 2)
   end
   local ib = args.ib
   local j = args.j
   local k = args.k or 0
   local imin = 0
   local imax = self:get_nic(ib) - 1
   local flowData = {}
   for i = imin, imax do
      flowData[#flowData+1] = self:get_cell_data{ib=ib, i=i, j=j, k=k}
   end
   return flowData
end

FlowSolution.get_cell_data_along_j_index = function(self, args)
   -- Return a table of CellData tables, for cells along an j-index slice
   -- for specified i,k indices of a structured-grid.
   if not (type(self)=='userdata') then
      error("Make sure that you are using FlowSolution:get_xxxx{args} not FlowSolution.get_xxxx{args}", 2)
   end
   if (not (args.ib and args.i)) then
      error("FlowSolution:get_cell_data_along_j_index{ib=n1, i=n2, k=n3}", 2)
   end
   local ib = args.ib
   local i = args.i
   local k = args.k or 0
   local jmin = 0
   local jmax = self:get_njc(ib) - 1
   local flowData = {}
   for j = jmin, jmax do
      flowData[#flowData+1] = self:get_cell_data{ib=ib, i=i, j=j, k=k}
   end
   return flowData
end

FlowSolution.get_cell_data_along_k_index = function(self, args)
   -- Return a table of CellData tables, for cells along an k-index slice
   -- for specified i,j indices of a structured-grid.
   if not (type(self)=='userdata') then
      error("Make sure that you are using FlowSolution:get_xxxx{args} not FlowSolution.get_xxxx{args}", 2)
   end
   if (not (args.ib and args.i and args.j)) then
      error("FlowSolution:get_cell_data_along_k_index{ib=n1, i=n2, j=n3}", 2)
   end
   local ib = args.ib
   local i = args.i
   local j = args.j
   local kmin = 0
   local kmax = self:get_nkc(ib) - 1
   if kmax < 0 then kmax = 0 end
   local flowData = {}
   for k = kmin, kmax do
      flowData[#flowData+1] = self:get_cell_data{ib=ib, i=i, j=j, k=k}
   end
   return flowData
end
