-- Module for collecting constructed boundary conditions.
-- 
-- Authors: PJ and RJG
-- Date: 2015-10-01
--         Extracted from prep.lua
--

module(..., package.seeall)

-- Helper function to check arguments for BC construction
function checkForInvalidArgs(tab, allowedArgs, funcName)
   -- allowedArgs is in array form build a reverse mapping.
   -- We'll put the args in as keys and set those values to true.
   argsTruthTable = {}
   for _,k in ipairs(allowedArgs) do
      argsTruthTable[k] = true
   end
   for k,_ in pairs(tab) do
      if not argsTruthTable[k] then
	 print(string.format("WARNING: %s is an INVALID argument in function: %s\n", k, funcName))
      end
   end 
end

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

FixedP = GhostCellEffect:new{p_outside=1.0e5}
FixedP.type = "fixed_pressure"
function FixedP:tojson()
   local str = string.format('          {"type": "%s", "p_outside": %f}',
			     self.type, self.p_outside)
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
				     direction_type="normal",
                                     direction_x=1.0, direction_y=0.0, direction_z=0.0,
                                     alpha=0.0, beta=0.0,
                                     mass_flux=0.0, relax_factor=0.10}
-- other options for direction type: uniform, radial, axial
FromStagnation.type = "from_stagnation_condition"
function FromStagnation:tojson()
   local str = string.format('          {"type": "%s",', self.type)
   str = str .. string.format(' "stagnation_condition": %s,',
			      self.stagCondition:toJSONString())
   str = str .. string.format(' "direction_type": "%s",', self.direction_type)
   str = str .. string.format(' "direction_x": %f, "direction_y": %f, "direction_z": %f,',
			      self.direction_x, self.direction_y, self.direction_z)
   str = str .. string.format(' "alpha": %f, "beta": %f,', self.alpha, self.beta)
   str = str .. string.format(' "mass_flux": %f, "relax_factor": %f',
			      self.mass_flux, self.relax_factor)
   str = str .. '}'
   return str
end

FullFaceExchangeCopy = GhostCellEffect:new{otherBlock=nil, otherFace=nil, orientation=-1,
					   reorient_vector_quantities=false,
					   Rmatrix={1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0}}
FullFaceExchangeCopy.type = "full_face_exchange_copy"
function FullFaceExchangeCopy:tojson()
   local str = string.format('          {"type": "%s", ', self.type)
   str = str .. string.format('"other_block": %d, ', self.otherBlock)
   str = str .. string.format('"other_face": "%s", ', self.otherFace)
   str = str .. string.format('"orientation": %d, ', self.orientation)
   str = str .. string.format('"reorient_vector_quantities": %s, ',
			      tostring(self.reorient_vector_quantities))
   str = str .. string.format('"Rmatrix": [')
   for i,v in ipairs(self.Rmatrix) do
      str = str .. string.format('%f', v)
      if i < #self.Rmatrix then str = str .. ', ' end
   end
   str = str .. ']' -- end of Rmatrix
   str = str .. '}' -- end of JSON value
   return str
end

MappedCellExchangeCopy = GhostCellEffect:new{transform_position=false,
					     c0=Vector3:new{0.0,0.0,0.0},
					     n=Vector3:new{0.0,0.0,1.0},
					     alpha=0.0,
					     delta=Vector3:new{0.0,0.0,0.0},
					     list_mapped_cells=false,
					     reorient_vector_quantities=false,
					     Rmatrix={1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0}}
MappedCellExchangeCopy.type = "mapped_cell_exchange_copy"
function MappedCellExchangeCopy:tojson()
   local str = string.format('          {"type": "%s", ', self.type)
   str = str .. string.format('"transform_position": %s, ',
			      tostring(self.transform_position))
   str = str .. string.format('"c0": [%f, %f, %f], ', self.c0:x(), self.c0:y(), self.c0:z())
   str = str .. string.format('"n": [%f, %f, %f], ', self.n:x(), self.n:y(), self.n:z())
   str = str .. string.format('"alpha": %f, ', self.alpha)
   str = str .. string.format('"delta": [%f, %f, %f], ', self.delta:x(), self.delta:y(), self.delta:z())
   str = str .. string.format('"list_mapped_cells": %s, ',
			      tostring(self.list_mapped_cells))
   str = str .. string.format('"reorient_vector_quantities": %s, ',
			      tostring(self.reorient_vector_quantities))
   str = str .. string.format('"Rmatrix": [')
   for i,v in ipairs(self.Rmatrix) do
      str = str .. string.format('%f', v)
      if i < #self.Rmatrix then str = str .. ', ' end
   end
   str = str .. ']' -- end of Rmatrix
   str = str .. '}' -- end of JSON value
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

TranslatingSurface = BoundaryInterfaceEffect:new{v_trans=nil}
-- Note that we are expecting v_trans as a table of 3 floats.
TranslatingSurface.type = "translating_surface"
function TranslatingSurface:tojson()
   local str = string.format('          {"type": "%s",', self.type)
   str = str .. string.format(' "v_trans": [%f, %f, %f]', self.v_trans:x(), 
			      self.v_trans:y(), self.v_trans:z())
   str = str .. '}'
   return str
end

RotatingSurface = BoundaryInterfaceEffect:new{centre=nil, r_omega=nil}
-- Note that we are expecting tables, each with three floats.
RotatingSurface.type = "rotating_surface"
function RotatingSurface:tojson()
   local str = string.format('          {"type": "%s",', self.type)
   str = str .. string.format(' "centre": [%f, %f, %f],', self.centre:x(), 
			      self.centre:y(), self.centre:z())
   str = str .. string.format(' "r_omega": [%f, %f, %f]', self.r_omega:x(), 
			      self.r_omega:y(), self.r_omega:z())
   str = str .. '}'
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
ConstFlux = BoundaryFluxEffect:new{flowCondition=nil}
ConstFlux.type = "const_flux"
function ConstFlux:tojson()
   local str = string.format('          {"type": "%s", ', self.type)
   str = str .. string.format(' "flowstate": %s', self.flowCondition:toJSONString())
   str = str .. '}'
   return str
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
-- BoundaryConditions consist of lists of actions to be done
-- at particular stages of the gas-dynamic update.

BoundaryCondition = {
   label = "",
   type = "",
   group = "",
   is_wall = true,
   ghost_cell_data_available = true,
   convective_flux_computed_in_bc = false,
   preReconAction = {},
   postConvFluxAction = {},
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
   str = str .. string.format('"is_wall": %s, ', tostring(self.is_wall))
   str = str .. string.format('"ghost_cell_data_available": %s, ',
			      tostring(self.ghost_cell_data_available))
   str = str .. string.format('"convective_flux_computed_in_bc": %s, ',
			      tostring(self.convective_flux_computed_in_bc))
   str = str .. '        "pre_recon_action": [\n'
   for i,effect in ipairs(self.preReconAction) do
      str = str .. effect:tojson()
      -- Extra code to deal with annoying JSON trailing comma deficiency
      if i ~= #self.preReconAction then str = str .. "," end
   end
   str = str .. '\n        ],\n'
   str = str .. '        "post_conv_flux_action": [\n'
   for i,effect in ipairs(self.postConvFluxAction) do
      str = str .. effect:tojson()
      if i ~= #self.postConvFluxAction then str = str .. "," end
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

WallBC_TranslatingSurface_FixedT = BoundaryCondition:new()
WallBC_TranslatingSurface_FixedT.type = "wall_translating_surface_fixed_t"
function WallBC_TranslatingSurface_FixedT:new(o)
   o = BoundaryCondition.new(self, o)
   o.preReconAction = { InternalCopyThenReflect:new() }
   o.preSpatialDerivAction = { CopyCellData:new(),
			       TranslatingSurface:new{v_trans=o.v_trans},
			       FixedT:new{Twall=o.Twall},
			       UpdateThermoTransCoeffs:new(),
			       WallKOmega:new() }
   return o
end

WallBC_TranslatingSurface_Adiabatic = BoundaryCondition:new()
WallBC_TranslatingSurface_Adiabatic.type = "wall_translating_surface_adiabatic"
function WallBC_TranslatingSurface_Adiabatic:new(o)
   o = BoundaryCondition.new(self, o)
   o.preReconAction = { InternalCopyThenReflect:new() }
   o.preSpatialDerivAction = { CopyCellData:new(),
			       TranslatingSurface:new{v_trans=o.v_trans},
			       WallKOmega:new() }
   return o
end

WallBC_RotatingSurface_FixedT = BoundaryCondition:new()
WallBC_RotatingSurface_FixedT.type = "wall_rotating_surface_fixed_t"
function WallBC_RotatingSurface_FixedT:new(o)
   o = BoundaryCondition.new(self, o)
   o.preReconAction = { InternalCopyThenReflect:new() }
   o.preSpatialDerivAction = { CopyCellData:new(),
			       RotatingSurface:new{r_omega=o.r_omega, centre=o.centre},
			       FixedT:new{Twall=o.Twall},
			       UpdateThermoTransCoeffs:new(),
			       WallKOmega:new() }
   return o
end

WallBC_RotatingSurface_Adiabatic = BoundaryCondition:new()
WallBC_RotatingSurface_Adiabatic.type = "wall_rotating_surface_adiabatic"
function WallBC_RotatingSurface_Adiabatic:new(o)
   o = BoundaryCondition.new(self, o)
   o.preReconAction = { InternalCopyThenReflect:new() }
   o.preSpatialDerivAction = { CopyCellData:new(),
			       TranslatingSurface:new{r_omega=o.r_omega, centre=o.centre},
			       WallKOmega:new() }
   return o
end

InFlowBC_Supersonic = BoundaryCondition:new()
InFlowBC_Supersonic.type = "inflow_supersonic"
function InFlowBC_Supersonic:new(o)
   o = BoundaryCondition.new(self, o)
   o.is_wall = false
   o.preReconAction = { FlowStateCopy:new{flowCondition=o.flowCondition} }
   o.preSpatialDerivAction = { CopyCellData:new() }
   return o
end

InFlowBC_ConstFlux = BoundaryCondition:new()
InFlowBC_ConstFlux.type = "inflow_const_flux"
function InFlowBC_ConstFlux:new(o)
   o = BoundaryCondition.new(self, o)
   o.is_wall = false
   o.convective_flux_computed_in_bc = true
   o.ghost_cell_data_available = false
   o.postConvFluxAction = { ConstFlux:new{flowCondition=o.flowCondition} }
   o.preSpatialDerivAction = { CopyCellData:new() }
   return o
end

InFlowBC_ShockFitting = BoundaryCondition:new()
InFlowBC_ShockFitting.type = "inflow_shock_fitting"
function InFlowBC_ShockFitting:new(o)
   o = BoundaryCondition.new(self, o)
   o.is_wall = false
   o.convective_flux_computed_in_bc = true
   o.ghost_cell_data_available = false
   o.postConvFluxAction = { ConstFlux:new{flowCondition=o.flowCondition} }
   o.preSpatialDerivAction = { CopyCellData:new() }
   return o
end

InFlowBC_FromStagnation = BoundaryCondition:new()
InFlowBC_FromStagnation.type = "inflow_from_stagnation_condition"
InFlowBC_FromStagnation.allowedArgs = {'type', 'label', 'group',
				       'is_wall', 'ghost_cell_data_available',
				       'convective_flux_computed_in_bc', 
				       'stagCondition',
				       'direction_type',
				       'direction_x',
				       'direction_y',
				       'direction_z',
				       'alpha',
				       'mass_flux',
				       'relax_factor'}
function InFlowBC_FromStagnation:new(o)
   o = BoundaryCondition.new(self, o)
   checkForInvalidArgs(o, self.allowedArgs, "InFlowBC_FromStagnation:new{}")
   o.is_wall = false
   o.preReconAction = { FromStagnation:new{stagCondition=o.stagCondition,
					   direction_type=o.direction_type,
					   direction_x=o.direction_x,
					   direction_y=o.direction_y,
					   direction_z=o.direction_z,
					   alpha=o.alpha, beta=o.beta,
					   mass_flux=o.mass_flux,
					   relax_factor=o.relax_factor} }
   o.preSpatialDerivAction = { CopyCellData:new() }
   return o
end

OutFlowBC_Simple = BoundaryCondition:new()
OutFlowBC_Simple.type = "outflow_simple_extrapolate"
function OutFlowBC_Simple:new(o)
   o = BoundaryCondition.new(self, o)
   o.is_wall = false
   o.preReconAction = { ExtrapolateCopy:new{xOrder = o.xOrder} }
   o.preSpatialDerivAction = { CopyCellData:new() }
   return o
end

OutFlowBC_FixedP = BoundaryCondition:new()
OutFlowBC_FixedP.type = "outflow_fixed_p"
function OutFlowBC_FixedP:new(o)
   o = BoundaryCondition.new(self, o)
   o.is_wall = false
   o.preReconAction = { ExtrapolateCopy:new{xOrder = o.xOrder},
			FixedP:new{p_outside=o.p_outside} }
   o.preSpatialDerivAction = { CopyCellData:new() }
   return o
end

OutFlowBC_FixedPT = BoundaryCondition:new()
OutFlowBC_FixedPT.type = "outflow_fixed_p_and_t"
function OutFlowBC_FixedPT:new(o)
   o = BoundaryCondition.new(self, o)
   o.is_wall = false
   o.preReconAction = { ExtrapolateCopy:new{xOrder = o.xOrder},
			FixedPT:new{p_outside=o.p_outside, T_outside=o.T_outside} }
   o.preSpatialDerivAction = { CopyCellData:new() }
   return o
end

ExchangeBC_FullFace = BoundaryCondition:new()
ExchangeBC_FullFace.type = "exchange_over_full_face"
function ExchangeBC_FullFace:new(o)
   o = BoundaryCondition.new(self, o)
   o.is_wall = false
   o.preReconAction = { FullFaceExchangeCopy:new{otherBlock=o.otherBlock,
						 otherFace=o.otherFace,
						 orientation=o.orientation,
						 reorient_vector_quantities=o.reorient_vector_quantities,
						 Rmatrix=o.Rmatrix} }
   o.preSpatialDerivAction = { UpdateThermoTransCoeffs:new() }
   return o
end

ExchangeBC_MappedCell = BoundaryCondition:new()
ExchangeBC_MappedCell.type = "exchange_using_mapped_cells"
function ExchangeBC_MappedCell:new(o)
   o = BoundaryCondition.new(self, o)
   o.is_wall = false
   o.preReconAction = { MappedCellExchangeCopy:new{transform_position=o.transform_position,
						   c0=o.c0, n=o.n, alpha=o.alpha, delta=o.delta,
						   list_mapped_cells=o.list_mapped_cells,
						   reorient_vector_quantities=o.reorient_vector_quantities,
						   Rmatrix=o.Rmatrix} }
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
   o.is_wall = true
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
