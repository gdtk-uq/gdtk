-- Author: Rowan J. Gollan
-- Date: 15-Feb-2013
-- Place: The University of Queensland, Brisbane, Australia
-- 
-- History:
--   15-Feb-2013 : extracted common elements from reac.lua and exch.lua
module(..., package.seeall)

local S = lpeg.S
local R = lpeg.R
local C = lpeg.C
local P = lpeg.P

Digit = R("09")
Integer = (S("+-") ^ -1) * (Digit^1)
Fractional = (P(".")   ) * (Digit^1)
Decimal = 
     (Integer *              -- Integer
     (Fractional ^ -1)) +    -- Fractional
     (S("+-") * Fractional)  -- Completely fractional number
Scientific = 
     Decimal * -- Decimal number
     S("Ee") * -- E or e
     Integer   -- Exponent
Number = Scientific + Decimal

Space = S(" \n\t")^0
Underscore = S("_")
Element = ((R("AZ") * R("az")^0) + P("e"))
S_letter = P("S")
Star = P("*")
ElecLevel = (R("az", "AZ", "09"))^-3 -- ^-3 says to match at most 3 occurrences
PM = S("+-")
Species = C( (R("ad","fz")^0 * (Element * Digit^0)^1)^1 * (PM + (Underscore * (S_letter + ElecLevel)))^0)
Tilde = P("~")
Dash = P("-") * Space
Comma = Space * P(",") * Space
Slash = Space * P("/") * Space
Colon = P(":") * Space
Open = "(" * Space
Close = Space * ")" * Space
FArrow = C(P("=>")) * Space
FRArrow = C(P("<=>") + P("=")) * Space
Equals = C(P("=")) * Space
Plus = P("+") * Space

