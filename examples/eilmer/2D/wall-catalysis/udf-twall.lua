-- udf-twall.lua
-- Lua script for defining spatially variable wall temperatures
--
-- @author: Nick Gibbons (22/02/21)


function interface(args)
   -- Ensure that this geometry matches the input file!
   -- Grid Geometry Specs  --            bo   ___---* co
   ri = 0.010              --        _-- *---/      |    -- Notes: - ai is always at the origin
   thetai = math.rad(5.0)  --       /               |              - ci and co are at x=L
   ro = 0.0215             --      /    _-*---------* ci           - oo is downstream of oi by diffo
   diffo = 0.0075          --     /    /  bi           
   thetao = math.rad(45.0) --     |    |
   L = 0.018               --  ao *----*  *  *

   oix = ri
   bix = oix + ri*-math.sin(thetai) -- bix is the x position of the join between the ETC material and the insulator

   -- Hyperbolic Tangent Function Profile
   w = 0.008      --  Width of the thermal boundary 
   Thot = 2200.0
   Tcold = 300.0
   y = 2*6/w*(args.x - bix)
   Twall = (Thot-Tcold)/2*math.tanh(-y) + (Thot-Tcold)/2 + Tcold

   -- The UDF boundary interface effect is smart enough to only set things that are in this table.
   face = {}
   face.T = Twall
   face.T_modes = {Twall}

   return face
end
