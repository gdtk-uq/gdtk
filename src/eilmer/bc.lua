-- Module for collecting constructed boundary conditions.
-- 
-- Authors: PJ and RJG
-- Date: 2015-10-01
--         Extracted from prep.lua
--

module(..., package.seeall)

-- -----------------------------------------------------------------------
-- Classes for constructing boundary conditions.
-- Each "complete" boundary condition is composed of lists of actions to do
-- at specific points in the superloop of the main simulation code.

-- For the classes below, we just follow the prototype pattern
-- given in Ierusalimchy's book "Programming in Lua"

-- Base class and subclasses for GhostCellEffect
GhostCellEffect = {
   type = ""
}
function GhostCellEffect:new(o)
   o = o or {}
   setmetatable(o, self)
   self.__index = self
   return o
end

InternalCopyThenReflect = GhostCellEffect:new()
InternalCopyThenReflect.type = "internal_copy_then_reflect"
function InternalCopyThenReflect:tojson()
   local str = string.format('          {"type" : "%s"}', self.type)
   return str
end

FlowStateCopy = GhostCellEffect:new{flowCondition=nil}
FlowStateCopy.type = "flowstate_copy"
function FlowStateCopy:tojson()
   local str = string.format('          {"type": "%s",', self.type)
   str = str .. string.format(' "flowstate": %s', self.flowCondition:toJSONString())
   str = str .. '}'
   return str
end

ExtrapolateCopy = GhostCellEffect:new{xOrder=0}
ExtrapolateCopy.type = "extrapolate_copy"
function ExtrapolateCopy:tojson()
   local str = string.format('          {"type": "%s", "x_order": %d}', self.type, self.xOrder)
   return str
end

FixedPT = GhostCellEffect:new{p_outside=1.0e5, T_outside=300.0}
FixedPT.type = "fixed_pressure_temperature"
function FixedPT:tojson()
   local str = string.format('          {"type": "%s", "p_outside": %f, "T_outside": %f}',
			     self.type, self.p_outside, self.T_outside)
   return str
end

FromStagnation = GhostCellEffect:new{stagCondition=nil,
				     direction_type="normal"}
-- other options for direction type: uniform, radial, axial
FromStagnation.type = "from_stagnation_condition"
function FromStagnation:tojson()
   local str = string.format('          {"type": "%s",', self.type)
   str = str .. string.format(' "stagnation_condition": %s,',
			      self.stagCondition:toJSONString())
   str = str .. string.format(' "direction_type": "%s"', self.direction_type)
   str = str .. '}'
   return str
end

FullFaceExchangeCopy = GhostCellEffect:new{otherBlock=nil, otherFace=nil, orientation=-1}
FullFaceExchangeCopy.type = "full_face_exchange_copy"
function FullFaceExchangeCopy:tojson()
   local str = string.format('          {"type": "%s", ', self.type)
   str = str .. string.format('"other_block": %d, ', self.otherBlock)
   str = str .. string.format('"other_face": "%s", ', self.otherFace)
   str = str .. string.format('"orientation": %d', self.orientation)
   str = str .. '}'
   return str
end

UserDefinedGhostCell = GhostCellEffect:new{fileName='user-defined-bc.lua'}
UserDefinedGhostCell.type = "user_defined"
function UserDefinedGhostCell:tojson()
   local str = string.format('         {"type": "%s", ', self.type)
   str = str .. string.format('"filename": "%s"', self.fileName)
   str = str .. '}'
   return str
end

-- Base class and subclasses for BoundaryInterfaceEffect
BoundaryInterfaceEffect = {
   type = ""
}
function BoundaryInterfaceEffect:new(o)
   o = o or {}
   setmetatable(o, self)
   self.__index = self
   return o
end

CopyCellData = BoundaryInterfaceEffect:new()
CopyCellData.type = "copy_cell_data"
function CopyCellData:tojson()
   local str = string.format('          {"type" : "%s"}', self.type)
   return str
end

ZeroVelocity = BoundaryInterfaceEffect:new()
ZeroVelocity.type = "zero_velocity"
function ZeroVelocity:tojson()
   local str = string.format('          {"type" : "%s"}', self.type)
   return str
end

FixedT = BoundaryInterfaceEffect:new{Twall=nil}
FixedT.type = "fixed_temperature"
function FixedT:tojson()
   local str = string.format('          {"type": "%s",', self.type)
   str = str .. string.format(' "Twall": %f', self.Twall)
   str = str .. '}'
   return str
end

UpdateThermoTransCoeffs = BoundaryInterfaceEffect:new()
UpdateThermoTransCoeffs.type = "update_thermo_trans_coeffs"
function UpdateThermoTransCoeffs:tojson()
   local str = string.format('          {"type" : "%s"}', self.type)
   return str
end

WallKOmega = BoundaryInterfaceEffect:new()
WallKOmega.type = "wall_k_omega"
function WallKOmega:tojson()
   local str = string.format('          {"type" : "%s"}', self.type)
   return str
end

TemperatureFromGasSolidInterface = BoundaryInterfaceEffect:new{otherBlock=nil, otherFace=nil, orientation=-1}
TemperatureFromGasSolidInterface.type = "temperature_from_gas_solid_interface"
function TemperatureFromGasSolidInterface:tojson()
   local str = string.format('          {"type": "%s", ', self.type)
   str = str .. string.format('"other_block": %d, ', self.otherBlock)
   str = str .. string.format('"other_face": "%s", ', self.otherFace)
   str = str .. string.format('"orientation": %d', self.orientation)
   str = str .. '}'
   return str
end

UserDefinedInterface = BoundaryInterfaceEffect:new{fileName='user-defined-bc.lua'}
UserDefinedInterface.type = "user_defined"
function UserDefinedInterface:tojson()
   local str = string.format('         {"type": "%s", ', self.type)
   str = str .. string.format('"filename": "%s"', self.fileName)
   str = str .. '}'
   return str
end

-- Base class and subclasses for BoundaryFluxEffect
BoundaryFluxEffect = {
   type = ""
}
function BoundaryFluxEffect:new(o)
   o = o or {}
   setmetatable(o, self)
   self.__index = self
   return o
end

EnergyFluxFromAdjacentSolid = BoundaryFluxEffect:new{otherBlock=nil, otherFace=nil, orientation=-1}
EnergyFluxFromAdjacentSolid.type = "energy_flux_from_adjacent_solid"
function EnergyFluxFromAdjacentSolid:tojson()
   local str = string.format('          {"type": "%s", ', self.type)
   str = str .. string.format('"other_block": %d, ', self.otherBlock)
   str = str .. string.format('"other_face": "%s", ', self.otherFace)
   str = str .. string.format('"orientation": %d', self.orientation)
   str = str .. '}'
   return str
end

-- Class for (complete) BoundaryCondition

BoundaryCondition = {
   label = "",
   type = "",
   group = "",
   preReconAction = {},
   preSpatialDerivAction = {},
   postDiffFluxAction = {}
}
function BoundaryCondition:new(o)
   o = o or {}
   setmetatable(o, self)
   self.__index = self
   return o
end
function BoundaryCondition:tojson()
   local str = '{'
   str = str .. string.format('"label": "%s", \n', self.label)
   str = str .. string.format('        "type": "%s", \n', self.type)
   str = str .. string.format('        "group": "%s", \n', self.group)
   str = str .. '        "pre_recon_action": [\n'
   for i,effect in ipairs(self.preReconAction) do
      str = str .. effect:tojson()
      -- Extra code to deal with annoying JSON trailing comma deficiency
      if i ~= #self.preReconAction then str = str .. "," end
   end
   str = str .. '\n        ],\n'
   str = str .. '        "pre_spatial_deriv_action": [\n'
   for i,effect in ipairs(self.preSpatialDerivAction) do
      str = str .. effect:tojson()
      if i ~= #self.preSpatialDerivAction then str = str .. "," end
   end
   str = str .. '\n        ],\n'
   str = str .. '        "post_diff_flux_action": [\n'
   for i,effect in ipairs(self.postDiffFluxAction) do
      str = str .. effect:tojson()
      if i ~= #self.postDiffFluxAction then str = str .. "," end
   end
   str = str .. '\n        ]\n'
   str = str .. '    }'
   return str
end

WallBC_WithSlip = BoundaryCondition:new()
WallBC_WithSlip.type = "wall_with_slip"
function WallBC_WithSlip:new(o)
   o = BoundaryCondition.new(self, o)
   o.preReconAction = { InternalCopyThenReflect:new() }
   o.preSpatialDerivAction = { CopyCellData:new() }
   return o
end

WallBC_NoSlip_FixedT = BoundaryCondition:new()
WallBC_NoSlip_FixedT.type = "wall_no_slip_fixed_t"
function WallBC_NoSlip_FixedT:new(o)
   o = BoundaryCondition.new(self, o)
   o.preReconAction = { InternalCopyThenReflect:new() }
   o.preSpatialDerivAction = { CopyCellData:new(), ZeroVelocity:new(),
			       FixedT:new{Twall=o.Twall},
			       UpdateThermoTransCoeffs:new(),
			       WallKOmega:new() }
   return o
end

WallBC_NoSlip_Adiabatic = BoundaryCondition:new()
WallBC_NoSlip_Adiabatic.type = "wall_no_slip_adiabatic"
function WallBC_NoSlip_Adiabatic:new(o)
   o = BoundaryCondition.new(self, o)
   o.preReconAction = { InternalCopyThenReflect:new() }
   o.preSpatialDerivAction = { CopyCellData:new(), ZeroVelocity:new(),
			       WallKOmega:new() }
   return o
end

InFlowBC_Supersonic = BoundaryCondition:new()
InFlowBC_Supersonic.type = "inflow_supersonic"
function InFlowBC_Supersonic:new(o)
   o = BoundaryCondition.new(self, o)
   o.preReconAction = { FlowStateCopy:new{flowCondition=o.flowCondition} }
   o.preSpatialDerivAction = { CopyCellData:new() }
   return o
end

InFlowBC_FromStagnation = BoundaryCondition:new()
InFlowBC_FromStagnation.type = "inflow_from_stagnation_condition"
function InFlowBC_FromStagnation:new(o)
   o = BoundaryCondition.new(self, o)
   o.preReconAction = { FromStagnation:new{stagCondition=o.stagCondition,
					   direction_type="normal"} }
   o.preSpatialDerivAction = { CopyCellData:new() }
   return o
end

OutFlowBC_Simple = BoundaryCondition:new()
OutFlowBC_Simple.type = "outflow_simple_extrapolate"
function OutFlowBC_Simple:new(o)
   o = BoundaryCondition.new(self, o)
   o.preReconAction = { ExtrapolateCopy:new{xOrder = o.xOrder} }
   o.preSpatialDerivAction = { CopyCellData:new() }
   return o
end

OutFlowBC_FixedPT = BoundaryCondition:new()
OutFlowBC_FixedPT.type = "outflow_fixed_p_and_t"
function OutFlowBC_FixedPT:new(o)
   o = BoundaryCondition.new(self, o)
   o.preReconAction = { ExtrapolateCopy:new{xOrder = o.xOrder},
			FixedPT:new{p_outside=o.p_outside, T_outside=o.T_outside} }
   o.preSpatialDerivAction = { CopyCellData:new() }
   return o
end

ExchangeBC_FullFace = BoundaryCondition:new()
ExchangeBC_FullFace.type = "exchange_over_full_face"
function ExchangeBC_FullFace:new(o)
   o = BoundaryCondition.new(self, o)
   o.preReconAction = { FullFaceExchangeCopy:new{otherBlock=o.otherBlock,
						 otherFace=o.otherFace,
						 orientation=o.orientation} }
   o.preSpatialDerivAction = { UpdateThermoTransCoeffs:new() }
   return o
end

UserDefinedBC = BoundaryCondition:new()
UserDefinedBC.type = "user_defined"
function UserDefinedBC:new(o)
   o = BoundaryCondition.new(self, o)
   o.preReconAction = { UserDefinedGhostCell:new{fileName=o.fileName} }
   o.preSpatialDerivAction = { UserDefinedInterface:new{fileName=o.fileName} } 
   return o
end

WallBC_AdjacentToSolid = BoundaryCondition:new()
WallBC_AdjacentToSolid.type = "wall_adjacent_to_solid"
function WallBC_AdjacentToSolid:new(o)
   o = BoundaryCondition.new(self, o)
   o.preReconAction = { InternalCopyThenReflect:new() }
   o.preSpatialDerivAction = { CopyCellData:new(), ZeroVelocity:new(),
			       TemperatureFromGasSolidInterface:new{otherBlock=o.otherBlock,
							    otherFace=o.otherFace,
							    orientation=o.orientation},
			       WallKOmega:new() }
   o.postDiffFluxAction = { EnergyFluxFromAdjacentSolid:new{otherBlock=o.otherBlock,
							    otherFace=o.otherFace,
							    orientation=o.orientation }
   }
   return o
end

-- Retain the old BC names as aliases, for now.
-- They are deprecated.
allowOldBCNames = true
if allowOldBCNames then
   print("Old boundary condition names are available.")
   SlipWallBC = WallBC_WithSlip
   FixedTWallBC = WallBC_NoSlip_FixedT
   AdiabaticWallBC = WallBC_NoSlip_Adiabatic
   SupInBC = InFlowBC_Supersonic
   ExtrapolateOutBC = OutFlowBC_Simple
   FixedPTOutBC = OutFlowBC_FixedPT
   FullFaceExchangeBC = ExchangeBC_FullFace
   AdjacentToSolidBC = WallBC_AdjacentToSolid
end

-- ---------------------------------------------------------------------------
-- Classes related to Solid blocks and boundary conditions

SolidBoundaryInterfaceEffect = {
   type = ""
}
function SolidBoundaryInterfaceEffect:new(o)
   o = o or {}
   setmetatable(o, self)
   self.__index = self
   return o
end

SolidBIE_FixedT = SolidBoundaryInterfaceEffect:new{Twall=300.0}
SolidBIE_FixedT.type = "fixed_temperature"
function SolidBIE_FixedT:tojson()
   local str = string.format('          {"type": "%s", ', self.type)
   str = str .. string.format('"Twall": %12.6e }', self.Twall)
   return str
end

SolidBIE_UserDefined = SolidBoundaryInterfaceEffect:new{fileName='user-defined-solid-bc.lua'}
SolidBIE_UserDefined.type = "user_defined"
function SolidBIE_UserDefined:tojson()
   local str = string.format('          {"type": "%s", ', self.type)
   str = str .. string.format('"filename": "%s" }', self.fileName)
   return str
end

-- Class for SolidBoundaryCondition
-- This class is a convenience class: it translates a high-level
-- user name for the boundary condition into a sequence of
-- lower-level operators.
SolidBoundaryCondition = {
   label = "",
   type = "",
   group = "",
   setsFluxDirectly = false,
   preSpatialDerivAction = {}
}
function SolidBoundaryCondition:new(o)
   o = o or {}
   setmetatable(o, self)
   self.__index = self
   return o
end
function SolidBoundaryCondition:tojson()
   local str = '{'
   str = str .. string.format('"label": "%s", \n', self.label)
   str = str .. string.format('        "type": "%s", \n', self.type)
   str = str .. string.format('        "group": "%s", \n', self.group)
   str = str .. string.format('        "sets_flux_directly": %s,\n', tostring(self.setsFluxDirectly))
   str = str .. '        "pre_spatial_deriv_action": [\n'
   for i,effect in ipairs(self.preSpatialDerivAction) do
      str = str .. effect:tojson()
      -- Extra code to deal with annoying JSON trailing comma deficiency
      if i ~= #self.preSpatialDerivAction then str = str .. "," end
   end
   str = str .. '\n        ]\n    }'
   return str
end

SolidFixedTBC = SolidBoundaryCondition:new()
SolidFixedTBC.type = "SolidFixedT"
function SolidFixedTBC:new(o)
   o = SolidBoundaryCondition.new(self, o)
   o.preSpatialDerivAction = { SolidBIE_FixedT:new{Twall=o.Twall} }
   return o
end

SolidUserDefinedBC = SolidBoundaryCondition:new()
SolidUserDefinedBC.type = "SolidUserDefined"
function SolidUserDefinedBC:new(o)
   o = SolidBoundaryCondition.new(self, o)
   o.preSpatialDerivAction = { SolidBIE_UserDefined:new{fileName=o.fileName} }
   return o
end

SolidAdjacentToGasBC = SolidBoundaryCondition:new()
SolidAdjacentToGasBC.type = "SolidAdjacentToGas"
function SolidAdjacentToGasBC:new(o)
   o = SolidBoundaryCondition.new(self, o)
   o.setsFluxDirectly = true
   return o
end
