-- Author: Rowan J. Gollan
-- Date: 15-Feb-2013
-- Place: The University of Queensland, Brisbane, Australia
-- 
-- History:
--   15-Feb-2013 : extracted common elements from reac.lua and exch.lua
module(..., package.seeall)

Space = lpeg.S(" \n\t")^0
Number = lpeg.R("09")^1
Underscore = lpeg.S("_")
Element = ((lpeg.R("AZ") * lpeg.R("az")^0) + lpeg.P("e"))
Solid = lpeg.P("S")
ElecLevel = (lpeg.R("az", "AZ", "09"))^-3 -- ^-3 says to match at most 3 occurrences
PM = lpeg.S("+-")
Species = lpeg.C(((Element * Number^0)^1 * PM^0)^1 * (Underscore * (Solid + ElecLevel))^0)
Tilde = lpeg.P("~")
Dash = lpeg.P("-") * Space
Comma = Space * lpeg.P(",") * Space
Colon = lpeg.P(":") * Space
Open = "(" * Space
Close = Space * ")" * Space
FArrow = lpeg.C(lpeg.P("=>")) * Space
RArrow = lpeg.C(lpeg.P("<=>")) * Space
Plus = lpeg.P("+") * Space

