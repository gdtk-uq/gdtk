-- gridarray.lua
-- A module for constructing arrays of structured grids.
--
-- PJ, 2021-10-04
--

module(..., package.seeall)


-- Class for GridArray objects.
GridArray = {
   myType = "GridArray"
}

function GridArray:new(o)
   local flag = type(self)=='table' and self.myType=='GridArray'
   if not flag then
      error("Make sure that you are using GridArray:new{} and not GridArray.new{}", 2)
   end
   o = o or {}
   local flag = checkAllowedNames(o, {"grid", "gridArray", "nib", "njb", "nkb"})
   if not flag then
      error("Invalid name for item supplied to GridArray constructor.", 2)
   end
   setmetatable(o, self)
   self.__index = self
   -- We will embed the GridArray identity in the individual grids
   -- and we would like that identity to start from 0 for the D code.
   o.id = #(gridArraysList)
   --
   return o
end -- GridArray:new

function GridArray:tojson(o)
   str = "[TODO] implement GridArray:tojson"
   return str
end
