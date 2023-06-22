-- Author: Rowan J. Gollan
-- Date: 15-Feb-2013
-- Place: The University of Queensland, Brisbane, Australia
-- 
-- History:
--   15-Feb-2013 : extracted common elements from reac.lua and exch.lua

local lpeg = require 'lpeg'

local S = lpeg.S
local R = lpeg.R
local C = lpeg.C
local P = lpeg.P

local Digit = R("09")
local Integer = (S("+-") ^ -1) * (Digit^1)
local Fractional = (P(".")   ) * (Digit^1)
local Decimal = 
     (Integer *              -- Integer
     (Fractional ^ -1)) +    -- Fractional
     (S("+-")^0 * Fractional)  -- Completely fractional number
local Scientific = 
     Decimal * -- Decimal number
     S("Ee") * -- E or e
     Integer   -- Exponent
local Number = Scientific + Decimal

local Space = S(" \n\t")^0
local Underscore = S("_")
local Element = ((R("AZ") * R("az")^0) + P("e"))
local S_letter = P("S")
local Star = P("*")
local ElecLevel = (R("az", "AZ", "09"))^-3 -- ^-3 says to match at most 3 occurrences
local PM = S("+-")
local Species = C( ((R("ad","fz") * P("-")^0)^0 * (Element * Digit^0)^1)^1 * (PM + (Underscore * (S_letter + ElecLevel)))^0)
local Tilde = P("~")
local Dash = P("-") * Space
local Comma = Space * P(",") * Space
local Slash = Space * P("/") * Space
local Colon = P(":") * Space
local Open = "(" * Space
local Close = Space * ")" * Space
local FArrow = C(P("=>"))
local FRArrow = C(P("<=>"))
local Equals = C(P("="))
local EqnSeparator = (FArrow + FRArrow + Equals) * Space
local Plus = P("+") * Space
local AllColliders = C(P("*all"))
local MolcColliders = C(P("*molcs"))
local IonColliders = C(P("*ions"))
local HeavyColliders = C(P("*heavy"))
local DoubleTilde = Space * P("~~") * Space

-- return whatever is needed externally
return {
   Space = Space,
   Open = Open,
   Close = Close,
   EqnSeparator = EqnSeparator,
   FArrow = FArrow,
   FRArrow = FRArrow,
   Equals = Equals,
   Plus = Plus,
   Number = Number,
   Species = Species,
   DoubleTilde = DoubleTilde,
   Comma = Comma,
   MolcColliders = MolcColliders,
   AllColliders = AllColliders,
   IonColliders = IonColliders,
   HeavyColliders = HeavyColliders
}
