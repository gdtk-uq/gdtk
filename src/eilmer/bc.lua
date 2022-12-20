-- Module for collecting constructed boundary conditions.
--
-- Authors: PJ and RJG
-- Date: 2015-10-01
--       Extracted from prep.lua (to be imported into prep.lua via require)
--

local MASSF_ERROR_TOL = 1.0e-6

function convertSpeciesTableToArray(massfTable)
   -- This function is essentially the Lua version of:
   -- luagas_model.d:getSpeciesValsFromTable
   local gm = getGasModel()
   local nsp = gm:nSpecies()
   -- 1. Check all keys are valid species names
   for sp,_ in pairs(massfTable) do
      isp = gm:speciesIndex(sp)
      if isp == -1 then
	 errMsg = "Species name used in table does not exist: " .. sp
	 error(errMsg)
      end
   end
   -- 2. Set all values to 0.0
   massfArray = {}
   for isp=0,nsp-1 do
      massfArray[isp] = 0.0
   end
   -- 3. Then populate those that are specified explicitly
   for sp,massf in pairs(massfTable) do
      isp = gm:speciesIndex(sp)
      massfArray[isp] = massf
   end
   -- 4. Do a check and possible scale of mass fractions
   massfSum = 0.0
   for isp=0,nsp-1 do massfSum = massfSum + massfArray[isp] end
   if math.abs(massfSum - 1) > MASSF_ERROR_TOL then
      errMsg = "The given mass faction values to do sum to 1.0\n"
      errMsg = errMsg .. string.format("The sum value is: %e\n", massfSum)
      errMsg = errMsg .. string.format("The error is larger than the tolerance: %e\n",
                                       MASSF_ERR_TOL)
      error(errMsg)
   end
   -- otherwise, perform a scaling
   for isp=0,nsp-1 do massfArray[isp] = massfArray[isp]/massfSum end
   return massfArray
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

FlowStateCopy = GhostCellEffect:new{flowState=nil, x0=0.0, y0=0.0, z0=0.0, r=0.0}
FlowStateCopy.type = "flowstate_copy"
function FlowStateCopy:tojson()
   local str = string.format('          {"type": "%s",', self.type)
   str = str .. string.format(' "flowstate": %s,', self.flowState:toJSONString())
   str = str .. string.format(' "x0": %.18e,', self.x0)
   str = str .. string.format(' "y0": %.18e,', self.y0)
   str = str .. string.format(' "z0": %.18e,', self.z0)
   str = str .. string.format(' "r": %.18e', self.r)
   str = str .. '}'
   return str
end

FlowStateCopyFromProfile = GhostCellEffect:new{filename=nil, match=nil}
FlowStateCopyFromProfile.type = "flowstate_copy_from_profile"
function FlowStateCopyFromProfile:tojson()
   local str = string.format('          {"type": "%s",', self.type)
   str = str .. string.format(' "filename": "%s",', self.filename)
   str = str .. string.format(' "match": "%s"', self.match)
   str = str .. '}'
   return str
end

FlowStateCopyFromHistory = GhostCellEffect:new{filename=nil}
FlowStateCopyFromHistory.type = "flowstate_copy_from_history"
function FlowStateCopyFromHistory:tojson()
   local str = string.format('          {"type": "%s",', self.type)
   str = str .. string.format(' "filename": "%s"', self.filename)
   str = str .. '}'
   return str
end

SynthesiseFlowState = GhostCellEffect:new{filename=nil}
SynthesiseFlowState.type = "synthesise_flowstate"
function SynthesiseFlowState:tojson()
   local str = string.format('          {"type": "%s",', self.type)
   str = str .. string.format(' "filename": "%s"', self.filename)
   str = str .. '}'
   return str
end

FromUpwindCopy = GhostCellEffect:new{flowState=nil}
FromUpwindCopy.type = "from_upwind_copy"
function FromUpwindCopy:tojson()
   local str = string.format('          {"type": "%s",', self.type)
   str = str .. string.format(' "flowstate": %s', self.flowState:toJSONString())
   str = str .. '}'
   return str
end

FromUpwindCopyDualState = GhostCellEffect:new{flowState1=nil, flowState2=nil, p=nil, n=nil}
FromUpwindCopyDualState.type = "from_upwind_copy_dual_state"
function FromUpwindCopyDualState:tojson()
   local str = string.format('          {"type": "%s",', self.type)
   str = str .. string.format(' "flowstate1": %s,', self.flowState1:toJSONString())
   str = str .. string.format(' "flowstate2": %s,', self.flowState2:toJSONString())
   str = str .. string.format('"p": [%.18e, %.18e, %.18e], ', self.p.x, self.p.y, self.p.z)
   str = str .. string.format('"n": [%.18e, %.18e, %.18e] ', self.n.x, self.n.y, self.n.z)
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
   local str = string.format('          {"type": "%s", "p_outside": %.18e}',
			     self.type, self.p_outside)
   return str
end

FixedPT = GhostCellEffect:new{p_outside=1.0e5, T_outside=300.0}
FixedPT.type = "fixed_pressure_temperature"
function FixedPT:tojson()
   local str = string.format('          {"type": "%s", "p_outside": %.18e, "T_outside": %.18e}',
			     self.type, self.p_outside, self.T_outside)
   return str
end

FromStagnation = GhostCellEffect:new{stagnationState=nil, fileName="",
				     direction_type="normal",
                                     direction_x=1.0, direction_y=0.0, direction_z=0.0,
                                     alpha=0.0, beta=0.0,
                                     mass_flux=0.0, relax_factor=0.10}
-- other options for direction type: uniform, radial, axial
FromStagnation.type = "from_stagnation_condition"
function FromStagnation:tojson()
   local str = string.format('          {"type": "%s",', self.type)
   str = str .. string.format(' "stagnation_condition": %s,',
			      self.stagnationState:toJSONString())
   str = str .. string.format(' "filename": "%s", ', self.fileName)
   str = str .. string.format(' "direction_type": "%s",', self.direction_type)
   str = str .. string.format(' "direction_x": %.18e, "direction_y": %.18e, "direction_z": %.18e,',
			      self.direction_x, self.direction_y, self.direction_z)
   str = str .. string.format(' "alpha": %.18e, "beta": %.18e,', self.alpha, self.beta)
   str = str .. string.format(' "mass_flux": %.18e, "relax_factor": %.18e',
			      self.mass_flux, self.relax_factor)
   str = str .. '}'
   return str
end

FullFaceCopy = GhostCellEffect:new{otherBlock=nil, otherFace=nil, orientation=-1,
                                   reorient_vector_quantities=false,
                                   Rmatrix={1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0}}
FullFaceCopy.type = "full_face_copy"
function FullFaceCopy:tojson()
   local str = string.format('          {"type": "%s", ', self.type)
   str = str .. string.format('"other_block": %d, ', self.otherBlock)
   str = str .. string.format('"other_face": "%s", ', self.otherFace)
   str = str .. string.format('"orientation": %d, ', self.orientation)
   str = str .. string.format('"reorient_vector_quantities": %s, ',
			      tostring(self.reorient_vector_quantities))
   str = str .. string.format('"Rmatrix": [')
   for i,v in ipairs(self.Rmatrix) do
      str = str .. string.format('%.18e', v)
      if i < #self.Rmatrix then str = str .. ', ' end
   end
   str = str .. ']' -- end of Rmatrix
   str = str .. '}' -- end of JSON value
   return str
end
MappedCellCopy = GhostCellEffect:new{cell_mapping_from_file=false, fileName='mapped_cells',
                                     transform_position=false,
                                     c0=Vector3:new{x=0.0,y=0.0,z=0.0},
                                     n=Vector3:new{x=0.0,y=0.0,z=1.0},
                                     alpha=0.0,
                                     delta=Vector3:new{x=0.0,y=0.0,z=0.0},
                                     list_mapped_cells=false,
                                     reorient_vector_quantities=false,
                                     Rmatrix={1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0}}
MappedCellCopy.type = "mapped_cell_copy"
function MappedCellCopy:tojson()
   local str = string.format('          {"type": "%s", ', self.type)
   str = str .. string.format('"cell_mapping_from_file": %s, ',
			      tostring(self.cell_mapping_from_file))
   str = str .. string.format('"filename": "%s", ', self.fileName)
   str = str .. string.format('"transform_position": %s, ',
			      tostring(self.transform_position))
   str = str .. string.format('"c0": [%.18e, %.18e, %.18e], ', self.c0.x, self.c0.y, self.c0.z)
   str = str .. string.format('"n": [%.18e, %.18e, %.18e], ', self.n.x, self.n.y, self.n.z)
   str = str .. string.format('"alpha": %.18e, ', self.alpha)
   str = str .. string.format('"delta": [%.18e, %.18e, %.18e], ',
                              self.delta.x, self.delta.y, self.delta.z)
   str = str .. string.format('"list_mapped_cells": %s, ',
			      tostring(self.list_mapped_cells))
   str = str .. string.format('"reorient_vector_quantities": %s, ',
			      tostring(self.reorient_vector_quantities))
   str = str .. string.format('"Rmatrix": [')
   for i,v in ipairs(self.Rmatrix) do
      str = str .. string.format('%.18e', v)
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

FlowStateCopyToInterface = BoundaryInterfaceEffect:new{flowState=nil, x0=0.0, y0=0.0, z0=0.0, r=0.0}
FlowStateCopyToInterface.type = "flow_state_copy_to_interface"
function FlowStateCopyToInterface:tojson()
   local str = string.format('          {"type": "%s",', self.type)
   str = str .. string.format(' "flowstate": %s,', self.flowState:toJSONString())
   str = str .. string.format(' "x0": %.18e,', self.x0)
   str = str .. string.format(' "y0": %.18e,', self.y0)
   str = str .. string.format(' "z0": %.18e,', self.z0)
   str = str .. string.format(' "r": %.18e', self.r)
   str = str .. '}'
   return str
end

FlowStateCopyFromProfileToInterface = BoundaryInterfaceEffect:new{filename=nil, match=nil}
FlowStateCopyFromProfileToInterface.type = "flow_state_copy_from_profile_to_interface"
function FlowStateCopyFromProfileToInterface:tojson()
   local str = string.format('          {"type": "%s",', self.type)
   str = str .. string.format(' "filename": "%s",', self.filename)
   str = str .. string.format(' "match": "%s"', self.match)
   str = str .. '}'
   return str
end

FlowStateCopyFromHistoryToInterface = BoundaryInterfaceEffect:new{filename=nil}
FlowStateCopyFromHistoryToInterface.type = "flow_state_copy_from_history_to_interface"
function FlowStateCopyFromHistoryToInterface:tojson()
   local str = string.format('          {"type": "%s",', self.type)
   str = str .. string.format(' "filename": "%s"', self.filename)
   str = str .. '}'
   return str
end

SynthesiseFlowStateToInterface = BoundaryInterfaceEffect:new{filename=nil}
SynthesiseFlowStateToInterface.type = "synthesise_flow_state_to_interface"
function SynthesiseFlowStateToInterface:tojson()
   local str = string.format('          {"type": "%s",', self.type)
   str = str .. string.format(' "filename": "%s"', self.filename)
   str = str .. '}'
   return str
end

ZeroVelocity = BoundaryInterfaceEffect:new()
ZeroVelocity.type = "zero_velocity"
function ZeroVelocity:tojson()
   local str = string.format('          {"type" : "%s"}', self.type)
   return str
end

ZeroSlipWallVelocity = BoundaryInterfaceEffect:new()
ZeroSlipWallVelocity.type = "zero_slip_wall_velocity"
function ZeroSlipWallVelocity:tojson()
   local str = string.format('          {"type" : "%s"}', self.type)
   return str
end

TranslatingSurface = BoundaryInterfaceEffect:new{v_trans=nil}
-- Note that we are expecting v_trans as a table of 3 floats, or a Vector3 object.
TranslatingSurface.type = "translating_surface"
function TranslatingSurface:tojson()
   local str = string.format('          {"type": "%s",', self.type)
   str = str .. string.format(' "v_trans": [%.18e, %.18e, %.18e]', self.v_trans.x,
			      self.v_trans.y, self.v_trans.z)
   str = str .. '}'
   return str
end

RotatingSurface = BoundaryInterfaceEffect:new{centre=nil, r_omega=nil}
-- Note that we are expecting tables, each with three floats, or Vector3 objects.
RotatingSurface.type = "rotating_surface"
function RotatingSurface:tojson()
   local str = string.format('          {"type": "%s",', self.type)
   str = str .. string.format(' "centre": [%.18e, %.18e, %.18e],', self.centre.x,
			      self.centre.y, self.centre.z)
   str = str .. string.format(' "r_omega": [%.18e, %.18e, %.18e]', self.r_omega.x,
			      self.r_omega.y, self.r_omega.z)
   str = str .. '}'
   return str
end

FixedT = BoundaryInterfaceEffect:new{Twall=nil}
FixedT.type = "fixed_temperature"
function FixedT:tojson()
   local str = string.format('          {"type": "%s",', self.type)
   str = str .. string.format(' "Twall": %.18e', self.Twall)
   str = str .. '}'
   return str
end

FixedComposition = BoundaryInterfaceEffect:new{wall_massf_composition={}}
FixedComposition.type = "fixed_composition"
function FixedComposition:tojson()
   local gm = getGasModel()
   local nsp = gm:nSpecies()
   local str = string.format('          {"type": "%s",', self.type)
   str = str .. ' "wall_massf_composition": [ '
   for isp=0,nsp-2 do
      str = str .. string.format("%.18e, ", self.wall_massf_composition[isp])
   end
   str = str .. string.format("%.18e ]\n", self.wall_massf_composition[nsp-1])
   str = str .. '}'
   return str
end


UpdateThermoTransCoeffs = BoundaryInterfaceEffect:new()
UpdateThermoTransCoeffs.type = "update_thermo_trans_coeffs"
function UpdateThermoTransCoeffs:tojson()
   local str = string.format('          {"type" : "%s"}', self.type)
   return str
end

WallTurbulent = BoundaryInterfaceEffect:new()
WallTurbulent.type = "wall_turbulent"
function WallTurbulent:tojson()
   local str = string.format('          {"type" : "%s"}', self.type)
   return str
end

WallFunctionInterfaceEffect = BoundaryInterfaceEffect:new()
WallFunctionInterfaceEffect.type = "wall_function_interface_effect"
function WallFunctionInterfaceEffect:tojson()
   local str = string.format('          {"type" : "%s"', self.type)
   str = str .. '}'
   return str
end

AdiabaticWallFunctionInterfaceEffect = BoundaryInterfaceEffect:new()
AdiabaticWallFunctionInterfaceEffect.type = "adiabatic_wall_function_interface_effect"
function AdiabaticWallFunctionInterfaceEffect:tojson()
   local str = string.format('          {"type" : "%s"', self.type)
   str = str .. '}'
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

TemperatureFromGasSolidInterface2 = BoundaryInterfaceEffect:new{otherBlock=nil, otherFace=nil, orientation=-1}
TemperatureFromGasSolidInterface2.type = "temperature_from_gas_solid_interface2"
function TemperatureFromGasSolidInterface2:tojson()
   local str = string.format('          {"type": "%s", ', self.type)
   str = str .. string.format('"other_block": %d, ', self.otherBlock)
   str = str .. string.format('"other_face": "%s", ', self.otherFace)
   str = str .. string.format('"orientation": %d', self.orientation)
   str = str .. '}'
   return str
end

ThermionicRadiativeEquilibrium = BoundaryInterfaceEffect:new{emissivity=1.0, Ar=0.0, phi=0.0,
                          ThermionicEmissionActive=0,catalytic_type="none", wall_massf_composition={}}
ThermionicRadiativeEquilibrium.type = "thermionic_radiative_equilibrium"
function ThermionicRadiativeEquilibrium:tojson()
   local str = string.format('          {"type": "%s",', self.type)
   str = str .. string.format(' "emissivity": %.18e,', self.emissivity)
   str = str .. string.format(' "Ar": %.18e,', self.Ar)
   str = str .. string.format(' "phi": %.18e,', self.phi)
   str = str .. string.format(' "ThermionicEmissionActive": %d,', self.ThermionicEmissionActive)
   str = str .. string.format(' "catalytic_type": "%s"', self.catalytic_type)
   if self.catalytic_type == "fixed_composition" then
       if next(self.wall_massf_composition) == nil then
           error("catalytic_type='fixed_composition' requires a wall_massf_composition table")
       end

       local gm = getGasModel()
       local nsp = gm:nSpecies()
       local mymf = convertSpeciesTableToArray(self.wall_massf_composition)
       str = str .. ', "wall_massf_composition": [ '
       for isp=0,nsp-2 do
          str = str .. string.format("%.18e, ", mymf[isp])
       end
       str = str .. string.format("%.18e ]", mymf[nsp-1])
   end
   str = str .. '}'
   return str
end

EquilibriumComposition = BoundaryInterfaceEffect:new{}
EquilibriumComposition.type = "equilibrium_composition"
function EquilibriumComposition:tojson()
   local str = string.format('          {"type": "%s"}', self.type)
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

-- Base class and subclasses for BoundaryCellEffect
BoundaryCellEffect = {
   type = ""
}
function BoundaryCellEffect:new(o)
   o = o or {}
   setmetatable(o, self)
   self.__index = self
   return o
end

WallFunctionCellEffect = BoundaryCellEffect:new()
WallFunctionCellEffect.type = "wall_function_cell_effect"
function WallFunctionCellEffect:tojson()
   local str = string.format('          {"type" : "%s"}', self.type)
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

ConstFlux = BoundaryFluxEffect:new{flowState=nil, x0=0.0, y0=0.0, z0=0.0, r=0.0}
ConstFlux.type = "const_flux"
function ConstFlux:tojson()
   local str = string.format('          {"type": "%s", ', self.type)
   str = str .. string.format(' "flowstate": %s,', self.flowState:toJSONString())
   str = str .. string.format(' "x0": %.18e,', self.x0)
   str = str .. string.format(' "y0": %.18e,', self.y0)
   str = str .. string.format(' "z0": %.18e,', self.z0)
   str = str .. string.format(' "r": %.18e', self.r)
   str = str .. '}'
   return str
end

SimpleOutflowFlux = BoundaryFluxEffect:new{}
SimpleOutflowFlux.type = "simple_outflow_flux"
function SimpleOutflowFlux:tojson()
   local str = string.format('          {"type": "%s" ', self.type)
   str = str .. '}'
   return str
end

UserDefinedFlux = BoundaryFluxEffect:new{fileName='user-defined-bc.lua', funcName='convective_flux'}
UserDefinedFlux.type = "user_defined"
function UserDefinedFlux:tojson()
   local str = string.format('          {"type": "%s", ', self.type)
   str = str .. string.format('"filename": "%s", ', self.fileName)
   str = str .. string.format('"function_name": "%s" ', self.funcName)
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

UpdateEnergyWallNormalVelocity = BoundaryFluxEffect:new{}
UpdateEnergyWallNormalVelocity.type = "update_energy_wall_normal_velocity"
function UpdateEnergyWallNormalVelocity:tojson()
   local str = string.format('          {"type" : "%s"}', self.type)
   return str
end

ThermionicElectronFlux = BoundaryFluxEffect:new{Ar=nil, phi=nil}
ThermionicElectronFlux.type = "thermionic_electron_flux"
function ThermionicElectronFlux:tojson()
   local str = string.format('          {"type": "%s",', self.type)
   str = str .. string.format(' "Ar": %.18e,', self.Ar)
   str = str .. string.format(' "phi": %.18e', self.phi)
   str = str .. '}'
   return str
end

-- Data for field boundaries (NNG)
FieldBoundary = {
   name = "unspecified",
}
function FieldBoundary:new(o)
   o = o or {}
   setmetatable(o, self)
   self.__index = self
   return o
end
function FieldBoundary:tojson()
   local str = '{'
   str = str .. string.format('"name": "%s"', self.name)
   str = str .. '}'
   return str
end

ZeroNormalGradient = FieldBoundary:new()
ZeroNormalGradient.name = "ZeroNormalGradient"
function ZeroNormalGradient:new(o)
   o = FieldBoundary.new(self, o)
   return o
end

FixedField = FieldBoundary:new{value=0.0}
FixedField.name = "FixedField"
function FixedField:new(o)
   o = FieldBoundary.new(self, o)
   return o
end
function FixedField:tojson()
   local str = string.format(' {"name": "%s", ', self.name)
   str = str .. string.format('"value": %.18e', self.value)
   str = str .. '}'
   return str
end

MixedField = FieldBoundary:new{differential=1.0, xinsulator=0.0, xcollector=0.0}
MixedField.name = "MixedField"
function MixedField:new(o)
   o = FieldBoundary.new(self, o)
   return o
end
function MixedField:tojson()
   local str = string.format(' {"name": "%s", ', self.name)
   str = str .. string.format('"differential": %.18e, ', self.differential)
   str = str .. string.format('"xinsulator": %.18e, ', self.xinsulator)
   str = str .. string.format('"xcollector": %.18e', self.xcollector)
   str = str .. '}'
   return str
end

FixedGradient_Test = FieldBoundary:new()
FixedGradient_Test.name = "FixedGradient_Test"
function FixedGradient_Test:new(o)
   o = FieldBoundary.new(self, o)
   return o
end

FixedField_Test = FieldBoundary:new()
FixedField_Test.name = "FixedField_Test"
function FixedField_Test:new(o)
   o = FieldBoundary.new(self, o)
   return o
end

-- Class for (complete) BoundaryCondition
--
-- BoundaryConditions consist of lists of actions to be done
-- at particular stages of the gas-dynamic update.
--
-- We expect that the user will construct instances of this class
-- in their input script, so it is worth checking arguments in
-- the following constructors.

BoundaryCondition = {
   label = "",
   type = "",
   group = "",
   is_gas_domain_bc = true,
   is_wall_with_viscous_effects = true,
   is_configured = false,
   ghost_cell_data_available = true,
   convective_flux_computed_in_bc = false,
   is_design_surface = false,
   num_cntrl_pts = 0,
   field_bc = FieldBoundary:new{},
   preReconAction = {},
   postConvFluxAction = {},
   preSpatialDerivActionAtBndryFaces = {},
   preSpatialDerivActionAtBndryCells = {},
   postDiffFluxAction = {}
}
function BoundaryCondition:new(o)
   o = o or {}
   setmetatable(o, self)
   self.__index = self
   -- If any of the actions are set, then we assume the boundary condition is configured.
   if ( #o.preReconAction > 0 or
	#o.postConvFluxAction > 0 or
	#o.preSpatialDerivActionAtBndryFaces > 0 or
	#o.preSpatialDerivActionAtBndryCells > 0 or
	#o.postDiffFluxAction > 0 ) then
      o.is_configured = true
   end
   return o
end
function BoundaryCondition:tojson()
   local str = '{'
   str = str .. string.format('"label": "%s", \n', self.label)
   str = str .. string.format('        "type": "%s", \n', self.type)
   str = str .. string.format('        "group": "%s", \n', self.group)
   str = str .. string.format('        "is_wall_with_viscous_effects": %s, \n',
                              tostring(self.is_wall_with_viscous_effects))
   str = str .. string.format('        "ghost_cell_data_available": %s, \n',
			      tostring(self.ghost_cell_data_available))
   str = str .. string.format('        "convective_flux_computed_in_bc": %s, \n',
			      tostring(self.convective_flux_computed_in_bc))
   str = str .. string.format('        "is_design_surface": %s, \n',
                              tostring(self.is_design_surface))
   str = str .. string.format('        "num_cntrl_pts": %s, \n',
                              tostring(self.num_cntrl_pts))
   str = str .. string.format('        "field_bc": %s,\n', self.field_bc:tojson())
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
   str = str .. '        "pre_spatial_deriv_action_at_bndry_faces": [\n'
   for i,effect in ipairs(self.preSpatialDerivActionAtBndryFaces) do
      str = str .. effect:tojson()
      if i ~= #self.preSpatialDerivActionAtBndryFaces then str = str .. "," end
   end
   str = str .. '\n        ],\n'
   str = str .. '        "pre_spatial_deriv_action_at_bndry_cells": [\n'
   for i,effect in ipairs(self.preSpatialDerivActionAtBndryCells) do
      str = str .. effect:tojson()
      if i ~= #self.preSpatialDerivActionAtBndryCells then str = str .. "," end
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

WallBC_WithSlip0 = BoundaryCondition:new()
WallBC_WithSlip0.type = "wall_with_slip"
function WallBC_WithSlip0:new(o)
   local flag = type(self)=='table' and self.type=='wall_with_slip'
   if not flag then
      error("Make sure that you are using WallBC_WithSlip0:new{}"..
               " and not WallBC_WithSlip0.new{}", 2)
   end
   o = o or {}
   flag = checkAllowedNames(o, {"label", "group", "is_design_surface", "num_cntrl_pts", "field_bc"})
   if not flag then
      error("Invalid name for item supplied to WallBC_WithSlip0 constructor.", 2)
   end
   o = BoundaryCondition.new(self, o)
   -- In a turbulence model sense, a slip wall is NOT a wall
   -- We mean a wall where the speed of the gas matches the speed of the wall
   o.ghost_cell_data_available = true
   o.is_wall_with_viscous_effects = false
   o.preReconAction = { InternalCopyThenReflect:new() }
   o.preSpatialDerivActionAtBndryFaces = { CopyCellData:new() }
   o.is_configured = true
   return o
end

-- Copy for trying out the new no-ghost-cell type of boundary condition
-- that makes use of the one-sided flux calculators at walls.
WallBC_WithSlip1 = BoundaryCondition:new()
WallBC_WithSlip1.type = "wall_with_slip"
function WallBC_WithSlip1:new(o)
   local flag = type(self)=='table' and self.type=='wall_with_slip'
   if not flag then
      error("Make sure that you are using WallBC_WithSlip1:new{}"..
               " and not WallBC_WithSlip1.new{}", 2)
   end
   o = o or {}
   flag = checkAllowedNames(o, {"label", "group", "is_design_surface", "num_cntrl_pts", "field_bc"})
   if not flag then
      error("Invalid name for item supplied to WallBC_WithSlip1 constructor.", 2)
   end
   o = BoundaryCondition.new(self, o)
   -- In a turbulence model sense, a slip wall is NOT a wall
   -- We mean a wall where the speed of the gas matches the speed of the wall
   o.ghost_cell_data_available = false
   o.is_wall_with_viscous_effects = false
   o.preReconAction = {}
   o.preSpatialDerivActionAtBndryFaces = { CopyCellData:new() }
   o.is_configured = true
   return o
end

-- Copy of WallBC_WithSlip for development/trialing of wall BCs
-- for moving mesh simulations for walls with normal velocity
WallBC_WithSlip2 = BoundaryCondition:new()
WallBC_WithSlip2.type = "wall_with_slip2"
function WallBC_WithSlip2:new(o)
   local flag = type(self)=='table' and self.type=='wall_with_slip2'
   if not flag then
      error("Make sure that you are using WallBC_WithSlip2:new{}"..
               " and not WallBC_WithSlip2.new{}", 2)
   end
   o = o or {}
   flag = checkAllowedNames(o, {"label", "group", "is_design_surface", "num_cntrl_pts", "field_bc"})
   if not flag then
      error("Invalid name for item supplied to WallBC_WithSlip2 constructor.", 2)
   end
   o = BoundaryCondition.new(self, o)
   -- In a turbulence model sense, a slip wall is NOT a wall
   -- We mean a wall where the speed of the gas matches the speed of the wall
   o.is_wall_with_viscous_effects = false
   o.convective_flux_computed_in_bc = true
   o.ghost_cell_data_available = true
   o.preReconAction = { InternalCopyThenReflect:new() }
   o.postConvFluxAction = { UpdateEnergyWallNormalVelocity:new()}
   o.preSpatialDerivActionAtBndryFaces = { CopyCellData:new()}
   o.is_configured = true
   return o
end

-- Select the default slip-wall boundary condition flavour.
WallBC_WithSlip = WallBC_WithSlip0


WallBC_NoSlip_FixedT0 = BoundaryCondition:new()
WallBC_NoSlip_FixedT0.type = "wall_no_slip_fixed_t"
function WallBC_NoSlip_FixedT0:new(o)
   local flag = type(self)=='table' and self.type=='wall_no_slip_fixed_t'
   if not flag then
      error("Make sure that you are using WallBC_NoSlip_FixedT0:new{}"..
               " and not WallBC_NoSlip_FixedT0.new{}", 2)
   end
   o = o or {}
   flag = checkAllowedNames(o, {"Twall", "wall_function", "field_bc",
                                "catalytic_type", "wall_massf_composition",
                                "label", "group", "is_design_surface", "num_cntrl_pts",
                                "user_post_diff_flux"})
   if not flag then
      error("Invalid name for item supplied to WallBC_NoSlip_FixedT0 constructor.", 2)
   end
   o = BoundaryCondition.new(self, o)
   o.ghost_cell_data_available = true
   o.is_wall_with_viscous_effects = true
   o.preReconAction = { InternalCopyThenReflect:new() }
   o.preSpatialDerivActionAtBndryFaces = { CopyCellData:new(), ZeroSlipWallVelocity:new(),
					   FixedT:new{Twall=o.Twall}}

   if o.catalytic_type == "fixed_composition" then
      o.preSpatialDerivActionAtBndryFaces[#o.preSpatialDerivActionAtBndryFaces+1] =
         FixedComposition:new{wall_massf_composition=convertSpeciesTableToArray(o.wall_massf_composition)}
   elseif o.catalytic_type == "equilibrium" then
      o.preSpatialDerivActionAtBndryFaces[#o.preSpatialDerivActionAtBndryFaces+1] =
         EquilibriumComposition:new{}
   end
   o.preSpatialDerivActionAtBndryFaces[#o.preSpatialDerivActionAtBndryFaces+1] =
      UpdateThermoTransCoeffs:new()

   if config.turbulence_model ~= "none" then
      o.preSpatialDerivActionAtBndryFaces[#o.preSpatialDerivActionAtBndryFaces+1] = WallTurbulent:new()
      if o.wall_function and config.turbulence_model == "k_omega" then
         -- Only makes sense to add a wall function if the k-omega model is active.
         o.preSpatialDerivActionAtBndryFaces[#o.preSpatialDerivActionAtBndryFaces+1] =
            WallFunctionInterfaceEffect:new{}
         o.preSpatialDerivActionAtBndryCells = { WallFunctionCellEffect:new() }
      end
   end

   if o.user_post_diff_flux then
      o.postDiffFluxAction = {UserDefinedFlux:new{fileName=o.user_post_diff_flux}}
   end
   o.is_configured = true
   return o
end

WallBC_NoSlip_FixedT1 = BoundaryCondition:new()
WallBC_NoSlip_FixedT1.type = "wall_no_slip_fixed_t"
function WallBC_NoSlip_FixedT1:new(o)
   local flag = type(self)=='table' and self.type=='wall_no_slip_fixed_t'
   if not flag then
      error("Make sure that you are using WallBC_NoSlip_FixedT1:new{}"..
               " and not WallBC_NoSlip_FixedT1.new{}", 2)
   end
   o = o or {}
   flag = checkAllowedNames(o, {"Twall", "wall_function","field_bc",
                                "catalytic_type", "wall_massf_composition",
                                "label", "group", "is_design_surface", "num_cntrl_pts"})
   if not flag then
      error("Invalid name for item supplied to WallBC_NoSlip_FixedT1 constructor.", 2)
   end
   o = BoundaryCondition.new(self, o)
   o.ghost_cell_data_available = false
   o.is_wall_with_viscous_effects = true
   o.preReconAction = {}
   o.preSpatialDerivActionAtBndryFaces = { CopyCellData:new(), ZeroSlipWallVelocity:new(),
					   FixedT:new{Twall=o.Twall}}

   if o.catalytic_type == "fixed_composition" then
      o.preSpatialDerivActionAtBndryFaces[#o.preSpatialDerivActionAtBndryFaces+1] =
         FixedComposition:new{wall_massf_composition=convertSpeciesTableToArray(o.wall_massf_composition)}
   elseif o.catalytic_type == "equilibrium" then
      o.preSpatialDerivActionAtBndryFaces[#o.preSpatialDerivActionAtBndryFaces+1] =
         EquilibriumComposition:new{}
   end
   o.preSpatialDerivActionAtBndryFaces[#o.preSpatialDerivActionAtBndryFaces+1] =
      UpdateThermoTransCoeffs:new()

   if config.turbulence_model ~= "none" then
      o.preSpatialDerivActionAtBndryFaces[#o.preSpatialDerivActionAtBndryFaces+1] = WallTurbulent:new()
      if o.wall_function and config.turbulence_model == "k_omega" then
	 -- Only makes sense to add a wall function if the k-omega model is active.
	 o.preSpatialDerivActionAtBndryFaces[#o.preSpatialDerivActionAtBndryFaces+1] =
	    WallFunctionInterfaceEffect:new{}
	 o.preSpatialDerivActionAtBndryCells = { WallFunctionCellEffect:new() }
      end
   end
   o.is_configured = true
   return o
end

WallBC_NoSlip_FixedT = WallBC_NoSlip_FixedT0

WallBC_NoSlip_UserDefinedT = BoundaryCondition:new()
WallBC_NoSlip_UserDefinedT.type = "wall_no_slip_user_defined_t"
function WallBC_NoSlip_UserDefinedT:new(o)
   local flag = type(self)=='table' and self.type=='wall_no_slip_user_defined_t'
   if not flag then
      error("Make sure that you are using WallBC_NoSlip_UserDefinedT:new{}"..
               " and not WallBC_NoSlip_UserDefinedT.new{}", 2)
   end
   o = o or {}
   flag = checkAllowedNames(o, {"Twall", "wall_function","field_bc",
                                "catalytic_type", "wall_massf_composition",
                                "label", "group", "is_design_surface", "num_cntrl_pts",
                                "thermionic_emission"})
   if not flag then
      error("Invalid name for item supplied to WallBC_NoSlip_UserDefinedT constructor.", 2)
   end
   o = BoundaryCondition.new(self, o)
   o.ghost_cell_data_available = true
   o.is_wall_with_viscous_effects = true
   o.preReconAction = { InternalCopyThenReflect:new() }
   o.preSpatialDerivActionAtBndryFaces = { CopyCellData:new(), ZeroSlipWallVelocity:new(),
                                           UserDefinedInterface:new{fileName=o.Twall}}

   if o.catalytic_type == "fixed_composition" then
      o.preSpatialDerivActionAtBndryFaces[#o.preSpatialDerivActionAtBndryFaces+1] =
         FixedComposition:new{wall_massf_composition=convertSpeciesTableToArray(o.wall_massf_composition)}
   elseif o.catalytic_type == "equilibrium" then
      o.preSpatialDerivActionAtBndryFaces[#o.preSpatialDerivActionAtBndryFaces+1] =
         EquilibriumComposition:new{}
   end
   o.preSpatialDerivActionAtBndryFaces[#o.preSpatialDerivActionAtBndryFaces+1] =
      UpdateThermoTransCoeffs:new()

   if config.turbulence_model ~= "none" then
      o.preSpatialDerivActionAtBndryFaces[#o.preSpatialDerivActionAtBndryFaces+1] = WallTurbulent:new()
      if o.wall_function and config.turbulence_model == "k_omega" then
         -- Only makes sense to add a wall function if the k-omega model is active.
         o.preSpatialDerivActionAtBndryFaces[#o.preSpatialDerivActionAtBndryFaces+1] =
            WallFunctionInterfaceEffect:new{}
         o.preSpatialDerivActionAtBndryCells = { WallFunctionCellEffect:new() }
      end
   end

   if o.thermionic_emission == "true" then
      o.postDiffFluxAction = {ThermionicElectronFlux:new{Ar=1.20e6, phi=2.0}}
   end
   o.is_configured = true
   return o
end


WallBC_ThermionicEmission = BoundaryCondition:new()
WallBC_ThermionicEmission.type = "wall_thermionic_emission"
function WallBC_ThermionicEmission:new(o)
   local flag = type(self)=='table' and self.type=='wall_thermionic_emission'
   if not flag then
      error("Make sure that you are using WallBC_ThermionicEmission:new{}"..
               " and not WallBC_ThermionicEmission.new{}", 2)
   end
   o = o or {}
   flag = checkAllowedNames(o, {"emissivity", "Ar", "phi", "ThermionicEmissionActive",
                                "catalytic_type", "wall_massf_composition", "field_bc",
                                "label", "group", "is_design_surface", "num_cntrl_pts"})
   if not flag then
      error("Invalid name for item supplied to WallBC_ThermionicEmission constructor.", 2)
   end
   o = BoundaryCondition.new(self, o)
   o.ghost_cell_data_available = true
   o.is_wall_with_viscous_effects = true
   o.preReconAction = { InternalCopyThenReflect:new() }
   o.preSpatialDerivActionAtBndryFaces = { CopyCellData:new(), ZeroSlipWallVelocity:new(),
      ThermionicRadiativeEquilibrium:new{emissivity=o.emissivity, Ar=o.Ar, phi=o.phi,
                                  ThermionicEmissionActive=o.ThermionicEmissionActive,
                                  catalytic_type=o.catalytic_type,
                                  wall_massf_composition=o.wall_massf_composition}}

   -- ThermionicRadiativeEquilibrium handles setting the thermostate and transprops
   -- at the boundaries. It needs access to these routines anyway, so it might as well.
   o.is_configured = true
   return o
end

WallBC_NoSlip_Adiabatic0 = BoundaryCondition:new()
WallBC_NoSlip_Adiabatic0.type = "wall_no_slip_adiabatic"
function WallBC_NoSlip_Adiabatic0:new(o)
   local flag = type(self)=='table' and self.type=='wall_no_slip_adiabatic'
   if not flag then
      error("Make sure that you are using WallBC_NoSlip_Adiabatic0:new{}"..
               " and not WallBC_NoSlip_Adiabatic0.new{}", 2)
   end
   o = o or {}
   flag = checkAllowedNames(o, {"wall_function","field_bc",
                                "catalytic_type", "wall_massf_composition",
                                "label", "group", "is_design_surface", "num_cntrl_pts"})
   if not flag then
      error("Invalid name for item supplied to WallBC_NoSlip_Adiabatic0 constructor.", 2)
   end
   o = BoundaryCondition.new(self, o)
   o.ghost_cell_data_available = true
   o.is_wall_with_viscous_effects = true
   o.preReconAction = { InternalCopyThenReflect:new() }
   o.preSpatialDerivActionAtBndryFaces = { CopyCellData:new(), ZeroSlipWallVelocity:new() }

   -- For the adiabatic wall we only need an UpdateThermoTransCoeffs if catalytic effects are present
   if o.catalytic_type == "fixed_composition" then
      o.preSpatialDerivActionAtBndryFaces[#o.preSpatialDerivActionAtBndryFaces+1] =
         FixedComposition:new{wall_massf_composition=convertSpeciesTableToArray(o.wall_massf_composition)}
      o.preSpatialDerivActionAtBndryFaces[#o.preSpatialDerivActionAtBndryFaces+1] =
         UpdateThermoTransCoeffs:new()
   elseif o.catalytic_type == "equilibrium" then
      o.preSpatialDerivActionAtBndryFaces[#o.preSpatialDerivActionAtBndryFaces+1] =
         EquilibriumComposition:new{}
      o.preSpatialDerivActionAtBndryFaces[#o.preSpatialDerivActionAtBndryFaces+1] =
         UpdateThermoTransCoeffs:new()
   end

   if config.turbulence_model ~= "none" then
      o.preSpatialDerivActionAtBndryFaces[#o.preSpatialDerivActionAtBndryFaces+1] = WallTurbulent:new()
      if o.wall_function and config.turbulence_model == "k_omega" then
	 -- Only makes sense to add a wall function if the k-omega model is active.
	 o.preSpatialDerivActionAtBndryFaces[#o.preSpatialDerivActionAtBndryFaces+1] =
	    AdiabaticWallFunctionInterfaceEffect:new{}
	 o.preSpatialDerivActionAtBndryCells = { WallFunctionCellEffect:new() }
      end
   end
   o.is_configured = true
   return o
end

WallBC_NoSlip_Adiabatic1 = BoundaryCondition:new()
WallBC_NoSlip_Adiabatic1.type = "wall_no_slip_adiabatic"
function WallBC_NoSlip_Adiabatic1:new(o)
   local flag = type(self)=='table' and self.type=='wall_no_slip_adiabatic'
   if not flag then
      error("Make sure that you are using WallBC_NoSlip_Adiabatic1:new{}"..
               " and not WallBC_NoSlip_Adiabatic1.new{}", 2)
   end
   o = o or {}
   flag = checkAllowedNames(o, {"wall_function","field_bc",
                                "catalytic_type", "wall_massf_composition",
                                "label", "group", "is_design_surface", "num_cntrl_pts"})
   if not flag then
      error("Invalid name for item supplied to WallBC_NoSlip_Adiabatic1 constructor.", 2)
   end
   o = BoundaryCondition.new(self, o)
   o.ghost_cell_data_available = false
   o.is_wall_with_viscous_effects = true
   o.preReconAction = {}
   o.preSpatialDerivActionAtBndryFaces = { CopyCellData:new(), ZeroSlipWallVelocity:new() }

   -- For the adiabatic wall we only need an UpdateThermoTransCoeffs if catalytic effects are present
   if o.catalytic_type == "fixed_composition" then
      o.preSpatialDerivActionAtBndryFaces[#o.preSpatialDerivActionAtBndryFaces+1] =
         FixedComposition:new{wall_massf_composition=convertSpeciesTableToArray(o.wall_massf_composition)}
      o.preSpatialDerivActionAtBndryFaces[#o.preSpatialDerivActionAtBndryFaces+1] =
         UpdateThermoTransCoeffs:new()
   elseif o.catalytic_type == "equilibrium" then
      o.preSpatialDerivActionAtBndryFaces[#o.preSpatialDerivActionAtBndryFaces+1] =
         EquilibriumComposition:new{}
      o.preSpatialDerivActionAtBndryFaces[#o.preSpatialDerivActionAtBndryFaces+1] =
         UpdateThermoTransCoeffs:new()
   end

   if config.turbulence_model ~= "none" then
      o.preSpatialDerivActionAtBndryFaces[#o.preSpatialDerivActionAtBndryFaces+1] = WallTurbulent:new()
      if o.wall_function and config.turbulence_model == "k_omega" then
	 -- Only makes sense to add a wall function if the k-omega model is active.
	 o.preSpatialDerivActionAtBndryFaces[#o.preSpatialDerivActionAtBndryFaces+1] =
	    AdiabaticWallFunctionInterfaceEffect:new{}
	 o.preSpatialDerivActionAtBndryCells = { WallFunctionCellEffect:new() }
      end
   end
   o.is_configured = true
   return o
end

WallBC_NoSlip_Adiabatic = WallBC_NoSlip_Adiabatic0

WallBC_TranslatingSurface_FixedT = BoundaryCondition:new()
WallBC_TranslatingSurface_FixedT.type = "wall_translating_surface_fixed_t"
function WallBC_TranslatingSurface_FixedT:new(o)
   local flag = type(self)=='table' and self.type=='wall_translating_surface_fixed_t'
   if not flag then
      error("Make sure that you are using WallBC_TranslatingSurface_FixedT:new{}"..
               " and not WallBC_TranslatingSurface_FixedT.new{}", 2)
   end
   o = o or {}
   flag = checkAllowedNames(o, {"v_trans", "Twall", "label", "group", "is_design_surface",
                                "num_cntrl_pts","field_bc"})
   if not flag then
      error("Invalid name for item supplied to WallBC_TranslatingSurface_FixedT constructor.", 2)
   end
   o = BoundaryCondition.new(self, o)
   o.preReconAction = { InternalCopyThenReflect:new() }
   -- Fill in missing components for v_trans
   if type(o.v_trans) == "table" then
      flag = checkAllowedNames(o.v_trans, {"x", "y", "z"})
      if not flag then
         error("Table representing v_trans should have only named cartesian components.", 2)
      end
   end
   o.v_trans.x = o.v_trans.x or 0.0
   o.v_trans.y = o.v_trans.y or 0.0
   o.v_trans.z = o.v_trans.z or 0.0
   o.preSpatialDerivActionAtBndryFaces = { CopyCellData:new(),
					   TranslatingSurface:new{v_trans=o.v_trans},
					   FixedT:new{Twall=o.Twall},
					   UpdateThermoTransCoeffs:new() }
   if config.turbulence_model ~= "none" then
      o.preSpatialDerivActionAtBndryFaces[#o.preSpatialDerivActionAtBndryFaces+1] = WallTurbulent:new()
   end
   o.is_configured = true
   return o
end

WallBC_TranslatingSurface_Adiabatic = BoundaryCondition:new()
WallBC_TranslatingSurface_Adiabatic.type = "wall_translating_surface_adiabatic"
function WallBC_TranslatingSurface_Adiabatic:new(o)
   local flag = type(self)=='table' and self.type=='wall_translating_surface_adiabatic'
   if not flag then
      error("Make sure that you are using WallBC_TranslatingSurface_Adiabatic:new{}"..
               " and not WallBC_TranslatingSurface_Adiabatic.new{}", 2)
   end
   o = o or {}
   flag = checkAllowedNames(o, {"v_trans", "label", "group", "is_design_surface", "num_cntrl_pts","field_bc"})
   if not flag then
      error("Invalid name for item supplied to WallBC_TranslatingSurface_Adiabatic constructor.", 2)
   end
   o = BoundaryCondition.new(self, o)
   o.preReconAction = { InternalCopyThenReflect:new() }
   -- Fill in missing components for v_trans
   if type(o.v_trans) == "table" then
      flag = checkAllowedNames(o.v_trans, {"x", "y", "z"})
      if not flag then
         error("Table representing v_trans should have only named cartesian components.", 2)
      end
   end
   o.v_trans.x = o.v_trans.x or 0.0
   o.v_trans.y = o.v_trans.y or 0.0
   o.v_trans.z = o.v_trans.z or 0.0
   o.preSpatialDerivActionAtBndryFaces = { CopyCellData:new(),
					   TranslatingSurface:new{v_trans=o.v_trans} }
   if config.turbulence_model ~= "none" then
      o.preSpatialDerivActionAtBndryFaces[#o.preSpatialDerivActionAtBndryFaces+1] = WallTurbulent:new()
   end
   o.is_configured = true
   return o
end

WallBC_RotatingSurface_FixedT = BoundaryCondition:new()
WallBC_RotatingSurface_FixedT.type = "wall_rotating_surface_fixed_t"
function WallBC_RotatingSurface_FixedT:new(o)
   local flag = type(self)=='table' and self.type=='wall_rotating_surface_fixed_t'
   if not flag then
      error("Make sure that you are using WallBC_RotatingSurface_FixedT:new{}"..
               " and not WallBC_RotatingSurface_FixedT.new{}", 2)
   end
   o = o or {}
   flag = checkAllowedNames(o, {"r_omega", "centre", "Twall", "label", "group",
                                "is_design_surface", "num_cntrl_pts","field_bc"})
   if not flag then
      error("Invalid name for item supplied to WallBC_RotatingSurface_FixedT constructor.", 2)
   end
   o = BoundaryCondition.new(self, o)
   o.preReconAction = { InternalCopyThenReflect:new() }
   -- Fill in missing components for r_omega and centre
   if type(o.r_omega) == "table" then
      flag = checkAllowedNames(o.r_omega, {"x", "y", "z"})
      if not flag then
         error("Table representing r_omega should have only named cartesian components.", 2)
      end
   end
   o.r_omega.x = o.r_omega.x or 0.0
   o.r_omega.y = o.r_omega.y or 0.0
   o.r_omega.z = o.r_omega.z or 0.0
   if type(o.centre) == "table" then
      flag = checkAllowedNames(o.centre, {"x", "y", "z"})
      if not flag then
         error("Table representing centre should have only named cartesian components.", 2)
      end
   end
   o.centre.x = o.centre.x or 0.0
   o.centre.y = o.centre.y or 0.0
   o.centre.z = o.centre.z or 0.0
   o.preSpatialDerivActionAtBndryFaces = { CopyCellData:new(),
					   RotatingSurface:new{r_omega=o.r_omega, centre=o.centre},
					   FixedT:new{Twall=o.Twall},
					   UpdateThermoTransCoeffs:new() }
   if config.turbulence_model ~= "none" then
      o.preSpatialDerivActionAtBndryFaces[#o.preSpatialDerivActionAtBndryFaces+1] = WallTurbulent:new()
   end
   o.is_configured = true
   return o
end

WallBC_RotatingSurface_Adiabatic = BoundaryCondition:new()
WallBC_RotatingSurface_Adiabatic.type = "wall_rotating_surface_adiabatic"
function WallBC_RotatingSurface_Adiabatic:new(o)
   local flag = type(self)=='table' and self.type=='wall_rotating_surface_adiabatic'
   if not flag then
      error("Make sure that you are using WallBC_RotatingSurface_Adiabatic:new{}"..
               " and not WallBC_RotatingSurface_Adiabatic.new{}", 2)
   end
   o = o or {}
   flag = checkAllowedNames(o, {"r_omega", "centre", "label", "group",
                                "is_design_surface", "num_cntrl_pts","field_bc"})
   if not flag then
      error("Invalid name for item supplied to WallBC_RotatingSurface_Adiabatic constructor.", 2)
   end
   o = BoundaryCondition.new(self, o)
   o.preReconAction = { InternalCopyThenReflect:new() }
   -- Fill in missing components for r_omega and centre
   if type(o.r_omega) == "table" then
      flag = checkAllowedNames(o.r_omega, {"x", "y", "z"})
      if not flag then
         error("Table representing r_omega should have only named cartesian components.", 2)
      end
   end
   o.r_omega.x = o.r_omega.x or 0.0
   o.r_omega.y = o.r_omega.y or 0.0
   o.r_omega.z = o.r_omega.z or 0.0
   if type(o.centre) == "table" then
      flag = checkAllowedNames(o.centre, {"x", "y", "z"})
      if not flag then
         error("Table representing centre should have only named cartesian components.", 2)
      end
   end
   o.centre.x = o.centre.x or 0.0
   o.centre.y = o.centre.y or 0.0
   o.centre.z = o.centre.z or 0.0
   o.preSpatialDerivActionAtBndryFaces = { CopyCellData:new(),
					   RotatingSurface:new{r_omega=o.r_omega, centre=o.centre} }
   if config.turbulence_model ~= "none" then
      o.preSpatialDerivActionAtBndryFaces[#o.preSpatialDerivActionAtBndryFaces+1] = WallTurbulent:new()
   end
   o.is_configured = true
   return o
end

InFlowBC_Supersonic = BoundaryCondition:new()
InFlowBC_Supersonic.type = "inflow_supersonic"
function InFlowBC_Supersonic:new(o)
   local flag = type(self)=='table' and self.type=='inflow_supersonic'
   if not flag then
      error("Make sure that you are using InFlowBC_Supersonic:new{}"..
               " and not InFlowBC_Supersonic.new{}", 2)
   end
   o = o or {}
   flag = checkAllowedNames(o, {"flowState", "flowCondition", "x0", "y0", "z0", "r", "label", "group","field_bc"})
   if not flag then
      error("Invalid name for item supplied to InFlowBC_Supersonic constructor.", 2)
   end
   if o.flowState == nil then o.flowState = o.flowCondition end -- look for old name
   -- Look for conical-flow parameters.
   o.x0 = o.x0 or 0.0
   o.y0 = o.y0 or 0.0
   o.z0 = o.z0 or 0.0
   o.r = o.r or 0.0
   -- Configure parameters and actions.
   o = BoundaryCondition.new(self, o)
   o.is_wall_with_viscous_effects = false
   o.preReconAction = {
      FlowStateCopy:new{flowState=o.flowState, x0=o.x0, y0=o.y0, z0=o.z0, r=o.r}
   }
   o.preSpatialDerivActionAtBndryFaces = {
      FlowStateCopyToInterface:new{flowState=o.flowState, x0=o.x0, y0=o.y0, z0=o.z0, r=o.r}
   }
   o.is_configured = true
   return o
end

InFlowBC_StaticProfile = BoundaryCondition:new()
InFlowBC_StaticProfile.type = "inflow_static_profile"
function InFlowBC_StaticProfile:new(o)
   local flag = type(self)=='table' and self.type=='inflow_static_profile'
   if not flag then
      error("Make sure that you are using InFlowBC_StaticProfile:new{}"..
               " and not InFlowBC_StaticProfile.new{}", 2)
   end
   o = o or {}
   flag = checkAllowedNames(o, {"filename", "fileName", "match", "label", "group","field_bc"})
   if not flag then
      error("Invalid name for item supplied to InFlowBC_StaticProfile constructor.", 2)
   end
   o = BoundaryCondition.new(self, o)
   o.is_wall_with_viscous_effects = false
   o.match = o.match or "xyz-to-xyz"
   o.filename = o.filename or o.fileName
   o.preReconAction = { FlowStateCopyFromProfile:new{filename=o.filename, match=o.match} }
   o.preSpatialDerivActionAtBndryFaces = {
      FlowStateCopyFromProfileToInterface:new{filename=o.filename, match=o.match}
   }
   o.is_configured = true
   return o
end

InFlowBC_Transient = BoundaryCondition:new()
InFlowBC_Transient.type = "inflow_transient"
function InFlowBC_Transient:new(o)
   local flag = type(self)=='table' and self.type=='inflow_transient'
   if not flag then
      error("Make sure that you are using InFlowBC_Transient:new{}"..
               " and not InFlowBC_Transient.new{}", 2)
   end
   o = o or {}
   flag = checkAllowedNames(o, {"filename", "fileName", "label", "group","field_bc"})
   if not flag then
      error("Invalid name for item supplied to InFlowBC_Transient constructor.", 2)
   end
   o = BoundaryCondition.new(self, o)
   o.is_wall_with_viscous_effects = false
   o.filename = o.filename or o.fileName
   o.preReconAction = { FlowStateCopyFromHistory:new{filename=o.filename} }
   o.preSpatialDerivActionAtBndryFaces = {
      FlowStateCopyFromHistoryToInterface:new{filename=o.filename}
   }
   o.is_configured = true
   return o
end

InFlowBC_Synthetic = BoundaryCondition:new()
InFlowBC_Synthetic.type = "inflow_synthetic"
function InFlowBC_Synthetic:new(o)
   local flag = type(self)=='table' and self.type=='inflow_synthetic'
   if not flag then
      error("Make sure that you are using InFlowBC_Synthetic:new{}"..
               " and not InFlowBC_Synthetic.new{}", 2)
   end
   o = o or {}
   flag = checkAllowedNames(o, {"filename", "fileName", "label", "group","field_bc"})
   if not flag then
      error("Invalid name for item supplied to InFlowBC_Synthetic constructor.", 2)
   end
   o = BoundaryCondition.new(self, o)
   o.is_wall_with_viscous_effects = false
   o.filename = o.filename or o.fileName
   o.preReconAction = { SynthesiseFlowState:new{filename=o.filename} }
   o.preSpatialDerivActionAtBndryFaces = {
      SynthesiseFlowStateToInterface:new{filename=o.filename}
   }
   o.is_configured = true
   return o
end

InFlowBC_ConstFlux = BoundaryCondition:new()
InFlowBC_ConstFlux.type = "inflow_const_flux"
function InFlowBC_ConstFlux:new(o)
   local flag = type(self)=='table' and self.type=='inflow_const_flux'
   if not flag then
      error("Make sure that you are using InFlowBC_ConstFlux:new{}"..
               " and not InFlowBC_ConstFlux.new{}", 2)
   end
   o = o or {}
   flag = checkAllowedNames(o, {"flowState", "flowCondition", "x0", "y0", "z0", "r", "label", "group","field_bc"})
   if not flag then
      error("Invalid name for item supplied to InFlowBC_ConstFlux constructor.", 2)
   end
   if o.flowState == nil then o.flowState = o.flowCondition end -- look for old name
   -- Look for conical-flow parameters.
   o.x0 = o.x0 or 0.0
   o.y0 = o.y0 or 0.0
   o.z0 = o.z0 or 0.0
   o.r = o.r or 0.0
   -- Configure parameters and actions.
   o = BoundaryCondition.new(self, o)
   o.is_wall_with_viscous_effects = false
   o.convective_flux_computed_in_bc = true
   o.ghost_cell_data_available = false
   o.postConvFluxAction = { ConstFlux:new{flowState=o.flowState, x0=o.x0, y0=o.y0, z0=o.z0, r=o.r} }
   o.preSpatialDerivActionAtBndryFaces = { CopyCellData:new() }
   o.is_configured = true
   return o
end

InFlowBC_ShockFitting = BoundaryCondition:new()
InFlowBC_ShockFitting.type = "inflow_shock_fitting"
function InFlowBC_ShockFitting:new(o)
   local flag = type(self)=='table' and self.type=='inflow_shock_fitting'
   if not flag then
      error("Make sure that you are using InFlowBC_ShockFitting:new{}"..
               " and not InFlowBC_ShockFitting.new{}", 2)
   end
   o = o or {}
   flag = checkAllowedNames(o, {"flowState", "flowCondition", "x0", "y0", "z0", "r", "label", "group","field_bc"})
   if not flag then
      error("Invalid name for item supplied to InFlowBC_ShockFitting constructor.", 2)
   end
   if o.flowState == nil then o.flowState = o.flowCondition end -- look for old name
   -- Look for conical-flow parameters.
   o.x0 = o.x0 or 0.0
   o.y0 = o.y0 or 0.0
   o.z0 = o.z0 or 0.0
   o.r = o.r or 0.0
   -- Configure parameters and actions.
   o = BoundaryCondition.new(self, o)
   o.is_wall_with_viscous_effects = false
   o.convective_flux_computed_in_bc = true
   o.ghost_cell_data_available = false
   o.postConvFluxAction = { ConstFlux:new{flowState=o.flowState, x0=o.x0, y0=o.y0, z0=o.z0, r=o.r} }
   o.preSpatialDerivActionAtBndryFaces = { CopyCellData:new() }
   o.is_configured = true
   return o
end

InFlowBC_FromStagnation = BoundaryCondition:new()
InFlowBC_FromStagnation.type = "inflow_from_stagnation_condition"
function InFlowBC_FromStagnation:new(o)
   local flag = type(self)=='table' and self.type=='inflow_from_stagnation_condition'
   if not flag then
      error("Make sure that you are using InFlowBC_FromStagnation:new{}"..
               " and not InFlowBC_FromStagnation.new{}", 2)
   end
   o = o or {}
   flag = checkAllowedNames(o, {"stagnationState", "stagCondition", "fileName", "filename",
                                "direction_type", "direction_x", "direction_y", "direction_z",
                                "alpha", "mass_flux", "relax_factor","field_bc",
                                "label", "group"})
   if not flag then
      error("Invalid name for item supplied to InFlowBC_FromStagnation constructor.", 2)
   end
   if o.stagnationState == nil then o.stagnationState = o.stagCondition end -- look for old name
   o.fileName = o.fileName or o.filename
   o = BoundaryCondition.new(self, o)
   o.is_wall_with_viscous_effects = false
   o.preReconAction = { FromStagnation:new{stagnationState=o.stagnationState,
                                           fileName=o.fileName,
					   direction_type=o.direction_type,
					   direction_x=o.direction_x,
					   direction_y=o.direction_y,
					   direction_z=o.direction_z,
					   alpha=o.alpha, beta=o.beta,
					   mass_flux=o.mass_flux,
					   relax_factor=o.relax_factor} }
   o.preSpatialDerivActionAtBndryFaces = { CopyCellData:new() }
   o.is_configured = true
   return o
end

InOutFlowBC_Ambient = BoundaryCondition:new()
InOutFlowBC_Ambient.type = "inoutflow_ambient"
function InOutFlowBC_Ambient:new(o)
   local flag = type(self)=='table' and self.type=='inoutflow_ambient'
   if not flag then
      error("Make sure that you are using InOutFlowBC_Ambient:new{}"..
               " and not InOutFlowBC_Ambient.new{}", 2)
   end
   o = o or {}
   flag = checkAllowedNames(o, {"flowState", "flowCondition", "label", "group","field_bc"})
   if not flag then
      error("Invalid name for item supplied to InOutFlowBC_Ambient constructor.", 2)
   end
   if o.flowState == nil then o.flowState = o.flowCondition end -- look for old name
   o = BoundaryCondition.new(self, o)
   o.is_wall_with_viscous_effects = false
   o.preReconAction = { FromUpwindCopy:new{flowState=o.flowState} }
   o.preSpatialDerivActionAtBndryFaces = { CopyCellData:new() }
   o.is_configured = true
   return o
end

InOutFlowBC_DualState = BoundaryCondition:new()
InOutFlowBC_DualState.type = "inoutflow_dualstate"
function InOutFlowBC_DualState:new(o)
   local flag = type(self)=='table' and self.type=='inoutflow_dualstate'
   if not flag then
      error("Make sure that you are using InOutFlowBC_DualState:new{}"..
               " and not InOutFlowBC_DualState.new{}", 2)
   end
   o = o or {}
   flag = checkAllowedNames(o, {"flowState1", "flowState2", "p", "n", "label", "group", "field_bc"})
   if not flag then
      error("Invalid name for item supplied to InOutFlowBC_DualState constructor.", 2)
   end
   if o.flowState1 == nil then
      error("Need to supply flowState1 for InOutFlowBC_DualState:new{}", 2)
   end
   if o.flowState2 == nil then
      error("Need to supply flowState2 for InOutFlowBC_DualState:new{}", 2)
   end
   if o.p == nil then
      error("Need to supply point p as a Vector3 or simple table for InOutFlowBC_DualState:new{}", 2)
   end
   o.p.x = o.p.x or 0.0
   o.p.y = o.p.y or 0.0
   o.p.z = o.p.z or 0.0
   if o.n == nil then
      error("Need to supply direction vector n as a Vector3 or simple table for InOutFlowBC_DualState:new{}", 2)
   end
   o.n.x = o.n.x or 0.0
   o.n.y = o.n.y or 1.0
   o.n.z = o.n.z or 0.0
   o = BoundaryCondition.new(self, o)
   o.is_wall_with_viscous_effects = false
   o.preReconAction = { FromUpwindCopyDualState:new{flowState1=o.flowState1, flowState2=o.flowState2, p=o.p, n=o.n} }
   o.preSpatialDerivActionAtBndryFaces = { CopyCellData:new() }
   o.is_configured = true
   return o
end

OutFlowBC_SimpleExtrapolate = BoundaryCondition:new()
OutFlowBC_SimpleExtrapolate.type = "outflow_simple_extrapolate"
function OutFlowBC_SimpleExtrapolate:new(o)
   local flag = type(self)=='table' and self.type=='outflow_simple_extrapolate'
   if not flag then
      error("Make sure that you are using OutFlowBC_SimpleExtrapolate:new{}"..
               " and not OutFlowBC_SimpleExtrapolate.new{}", 2)
   end
   o = o or {}
   flag = checkAllowedNames(o, {"xOrder", "label", "group","field_bc"})
   if not flag then
      error("Invalid name for item supplied to OutFlowBC_SimpleExtrapolate constructor.", 2)
   end
   o = BoundaryCondition.new(self, o)
   o.is_wall_with_viscous_effects = false
   o.preReconAction = { ExtrapolateCopy:new{xOrder = o.xOrder} }
   o.preSpatialDerivActionAtBndryFaces = { CopyCellData:new() }
   o.is_configured = true
   return o
end

OutFlowBC_SimpleFlux = BoundaryCondition:new()
OutFlowBC_SimpleFlux.type = "outflow_simple_flux"
function OutFlowBC_SimpleFlux:new(o)
   local flag = type(self)=='table' and self.type=='outflow_simple_flux'
   if not flag then
      error("Make sure that you are using OutFlowBC_SimpleFlux:new{}"..
               " and not OutFlowBC_SimpleFlux.new{}", 2)
   end
   o = o or {}
   flag = checkAllowedNames(o, {"label", "group","field_bc"})
   if not flag then
      error("Invalid name for item supplied to OutFlowBC_SimpleFlux constructor.", 2)
   end
   o = BoundaryCondition.new(self, o)
   o.is_wall_with_viscous_effects = false
   o.preReconAction = { ExtrapolateCopy:new{xOrder = o.xOrder} }
   o.preSpatialDerivActionAtBndryFaces = { CopyCellData:new() }
   o.postConvFluxAction = { SimpleOutflowFlux:new() }
   o.ghost_cell_data_available = true
   o.convective_flux_computed_in_bc = true
   o.is_configured = true
   return o
end
-- Old name is retained as an alias.
OutFlowBC_Simple = OutFlowBC_SimpleFlux

OutFlowBC_FixedP = BoundaryCondition:new()
OutFlowBC_FixedP.type = "outflow_fixed_p"
function OutFlowBC_FixedP:new(o)
   local flag = type(self)=='table' and self.type=='outflow_fixed_p'
   if not flag then
      error("Make sure that you are using OutFlowBC_FixedP:new{}"..
               " and not OutFlowBC_FixedP.new{}", 2)
   end
   o = o or {}
   flag = checkAllowedNames(o, {"xOrder", "p_outside", "label", "group","field_bc"})
   if not flag then
      error("Invalid name for item supplied to OutFlowBC_FixedP constructor.", 2)
   end
   o = BoundaryCondition.new(self, o)
   o.is_wall_with_viscous_effects = false
   o.preReconAction = { ExtrapolateCopy:new{xOrder = o.xOrder},
			FixedP:new{p_outside=o.p_outside} }
   o.preSpatialDerivActionAtBndryFaces = { CopyCellData:new() }
   o.is_configured = true
   return o
end

OutFlowBC_FixedPT = BoundaryCondition:new()
OutFlowBC_FixedPT.type = "outflow_fixed_p_and_t"
function OutFlowBC_FixedPT:new(o)
   local flag = type(self)=='table' and self.type=='outflow_fixed_p_and_t'
   if not flag then
      error("Make sure that you are using OutFlowBC_FixedPT:new{}"..
               " and not OutFlowBC_FixedPT.new{}", 2)
   end
   o = o or {}
   flag = checkAllowedNames(o, {"xOrder", "p_outside", "T_outside",
                                "label", "group","field_bc"})
   if not flag then
      error("Invalid name for item supplied to OutFlowBC_FixedPT constructor.", 2)
   end
   o = BoundaryCondition.new(self, o)
   o.is_wall_with_viscous_effects = false
   o.preReconAction = { ExtrapolateCopy:new{xOrder = o.xOrder},
			FixedPT:new{p_outside=o.p_outside, T_outside=o.T_outside} }
   o.preSpatialDerivActionAtBndryFaces = { CopyCellData:new() }
   o.is_configured = true
   return o
end

ExchangeBC_FullFace = BoundaryCondition:new()
ExchangeBC_FullFace.type = "exchange_over_full_face"
function ExchangeBC_FullFace:new(o)
   local flag = type(self)=='table' and self.type=='exchange_over_full_face'
   if not flag then
      error("Make sure that you are using ExchangeBC_FullFace:new{}"..
               " and not ExchangeBC_FullFace.new{}", 2)
   end
   o = o or {}
   flag = checkAllowedNames(o, {"otherBlock", "otherFace", "orientation",
                                "reorient_vector_quantities", "Rmatrix",
                                "label", "group"})
   if not flag then
      error("Invalid name for item supplied to ExchangeBC_FullFace constructor.", 2)
   end
   o = BoundaryCondition.new(self, o)
   o.is_wall_with_viscous_effects = false
   o.preReconAction = { FullFaceCopy:new{otherBlock=o.otherBlock,
                                         otherFace=o.otherFace,
                                         orientation=o.orientation,
                                         reorient_vector_quantities=o.reorient_vector_quantities,
                                         Rmatrix=o.Rmatrix} }
   o.preSpatialDerivActionAtBndryFaces = { UpdateThermoTransCoeffs:new() }
   o.is_configured = true
   return o
end

ExchangeBC_FullFacePlusUDF = BoundaryCondition:new()
ExchangeBC_FullFacePlusUDF.type = "exchange_over_full_face_plus_udf"
function ExchangeBC_FullFacePlusUDF:new(o)
   local flag = type(self)=='table' and self.type=='exchange_over_full_face_plus_udf'
   if not flag then
      error("Make sure that you are using ExchangeBC_FullFacePlusUDF:new{}"..
               " and not ExchangeBC_FullFacePlusUDF.new{}", 2)
   end
   o = o or {}
   flag = checkAllowedNames(o, {"otherBlock", "otherFace", "orientation",
                                "reorient_vector_quantities", "Rmatrix",
                                "fileName", "filename", "label", "group"})
   if not flag then
      error("Invalid name for item supplied to ExchangeBC_FullFacePlusUDF constructor.", 2)
   end
   o = BoundaryCondition.new(self, o)
   o.is_wall_with_viscous_effects = false
   o.preReconAction = { FullFaceCopy:new{otherBlock=o.otherBlock,
                                         otherFace=o.otherFace,
                                         orientation=o.orientation,
                                         reorient_vector_quantities=o.reorient_vector_quantities,
                                         Rmatrix=o.Rmatrix},
                        UserDefinedGhostCell:new{fileName=o.fileName}
   }
   o.preSpatialDerivActionAtBndryFaces = {
      UserDefinedInterface:new{fileName=o.fileName},
      UpdateThermoTransCoeffs:new()
   }
   o.is_configured = true
   return o
end

ExchangeBC_MappedCell = BoundaryCondition:new()
ExchangeBC_MappedCell.type = "exchange_using_mapped_cells"
function ExchangeBC_MappedCell:new(o)
   local flag = type(self)=='table' and self.type=='exchange_using_mapped_cells'
   if not flag then
      error("Make sure that you are using ExchangeBC_MappedCell:new{}"..
               " and not ExchangeBC_MappedCell.new{}", 2)
   end
   o = o or {}
   flag = checkAllowedNames(o, {"cell_mapping_from_file", "fileName", "filename",
                                "transform_position", "c0", "n", "alpha", "delta",
                                "list_mapped_cells",
                                "reorient_vector_quantities", "Rmatrix",
                                "label", "group"})
   if not flag then
      error("Invalid name for item supplied to ExchangeBC_MappedCell constructor.", 2)
   end
   o = BoundaryCondition.new(self, o)
   o.is_wall_with_viscous_effects = false
   o.fileName = o.fileName or o.filename
   o.preReconAction = { MappedCellCopy:new{cell_mapping_from_file=o.cell_mapping_from_file,
                                           fileName=o.fileName,
                                           transform_position=o.transform_position,
                                           c0=o.c0, n=o.n, alpha=o.alpha, delta=o.delta,
                                           list_mapped_cells=o.list_mapped_cells,
                                           reorient_vector_quantities=o.reorient_vector_quantities,
                                           Rmatrix=o.Rmatrix} }
   o.preSpatialDerivActionAtBndryFaces = { UpdateThermoTransCoeffs:new() }
   o.is_configured = true
   return o
end

UserDefinedGhostCellBC = BoundaryCondition:new()
UserDefinedGhostCellBC.type = "user_defined_ghost_cell"
function UserDefinedGhostCellBC:new(o)
   local flag = type(self)=='table' and self.type=='user_defined_ghost_cell'
   if not flag then
      error("Make sure that you are using UserDefinedGhostCellBC:new{}"..
               " and not UserDefinedGhostCellBC.new{}", 2)
   end
   o = o or {}
   flag = checkAllowedNames(o, {"fileName", "filename", "label", "group",
                                "is_design_surface", "num_cntrl_pts","field_bc"})
   if not flag then
      error("Invalid name for item supplied to UserDefinedGhostCellBC constructor.", 2)
   end
   o = BoundaryCondition.new(self, o)
   o.fileName = o.fileName or o.filename
   o.preReconAction = { UserDefinedGhostCell:new{fileName=o.fileName} }
   o.preSpatialDerivActionAtBndryFaces = {
      UserDefinedInterface:new{fileName=o.fileName},
      UpdateThermoTransCoeffs:new()
   }
   o.is_configured = true
   return o
end
-- Keep the old name.
UserDefinedBC = UserDefinedGhostCellBC

UserDefinedFluxBC = BoundaryCondition:new()
UserDefinedFluxBC.type = "user_defined_flux"
function UserDefinedFluxBC:new(o)
   local flag = type(self)=='table' and self.type=='user_defined_flux'
   if not flag then
      error("Make sure that you are using UserDefinedFluxBC:new{}"..
               " and not UserDefinedFluxBC.new{}", 2)
   end
   o = o or {}
   flag = checkAllowedNames(o, {"fileName", "filename", "funcName", "funcname",
                                "label", "group", "is_design_surface", "num_cntrl_pts","field_bc"})
   if not flag then
      error("Invalid name for item supplied to UserDefinedFluxBC constructor.", 2)
   end
   o = BoundaryCondition.new(self, o)
   o.fileName = o.fileName or o.filename
   o.funcName = o.funcName or o.funcname
   o.funcName = o.funcName or 'convectiveFlux'
   o.postConvFluxAction = { UserDefinedFlux:new{fileName=o.fileName, funcName=o.funcName} }
   o.ghost_cell_data_available = false
   o.convective_flux_computed_in_bc = true
   o.is_configured = true
   return o
end

WallBC_AdjacentToSolid = BoundaryCondition:new()
WallBC_AdjacentToSolid.type = "wall_adjacent_to_solid"
function WallBC_AdjacentToSolid:new(o)
   local flag = type(self)=='table' and self.type=='wall_adjacent_to_solid'
   if not flag then
      error("Make sure that you are using WallBC_AdjacentToSolid:new{}"..
               " and not WallBC_AdjacentToSolid.new{}", 2)
   end
   o = o or {}
   flag = checkAllowedNames(o, {"otherBlock", "otherFace", "orientation",
                                "label", "group", "is_design_surface", "num_cntrl_pts","field_bc"})
   if not flag then
      error("Invalid name for item supplied to WallBC_AdjacentToSolid constructor.", 2)
   end
   o = BoundaryCondition.new(self, o)
   o.group = "adjacent_to_solid"
   o.is_wall_with_viscous_effects = true
   o.preReconAction = { GasSolidFullFaceCopy:new{otherBlock=o.otherBlock,
                                                 otherFace=o.otherFace,
                                                 orientation=o.orientation},
                        InternalCopyThenReflect:new() }
   o.preSpatialDerivActionAtBndryFaces = { CopyCellData:new(),
                                           ZeroSlipWallVelocity:new(),
                                           TemperatureFromGasSolidInterface:new{otherBlock=o.otherBlock,
                                                                                otherFace=o.otherFace,
                                                                                orientation=o.orientation},
                                           UpdateThermoTransCoeffs:new()}
   if config.turbulence_model ~= "none" then
      o.preSpatialDerivActionAtBndryFaces[#o.preSpatialDerivActionAtBndryFaces+1] = WallTurbulent:new()
   end
   o.postDiffFluxAction = { EnergyFluxFromAdjacentSolid:new{otherBlock=o.otherBlock,
							    otherFace=o.otherFace,
							    orientation=o.orientation }
   }
   o.is_configured = true
   return o
end

WallBC_AdjacentToSolid2 = BoundaryCondition:new()
WallBC_AdjacentToSolid2.type = "wall_adjacent_to_solid2"
function WallBC_AdjacentToSolid2:new(o)
   local flag = type(self)=='table' and self.type=='wall_adjacent_to_solid2'
   if not flag then
      error("Make sure that you are using WallBC_AdjacentToSolid2:new{}"..
               " and not WallBC_AdjacentToSolid2.new{}", 2)
   end
   o = o or {}
   flag = checkAllowedNames(o, {"otherBlock", "otherFace", "orientation", "catalytic_type", "wall_massf_composition",
                                "label", "group", "is_design_surface", "num_cntrl_pts","field_bc"})
   if not flag then
      error("Invalid name for item supplied to WallBC_AdjacentToSolid2 constructor.", 2)
   end
   o = BoundaryCondition.new(self, o)
   o.group = "adjacent_to_solid"
   o.is_wall_with_viscous_effects = true
   o.preReconAction = { GasSolidFullFaceCopy:new{otherBlock=o.otherBlock,
                                                 otherFace=o.otherFace,
                                                 orientation=o.orientation},
                        InternalCopyThenReflect:new() }
   o.preSpatialDerivActionAtBndryFaces = { CopyCellData:new(),
                                           ZeroSlipWallVelocity:new(),
                                           TemperatureFromGasSolidInterface2:new{otherBlock=o.otherBlock,
                                                                                 otherFace=o.otherFace,
                                                                                 orientation=o.orientation}}

   if o.catalytic_type == "fixed_composition" then
      o.preSpatialDerivActionAtBndryFaces[#o.preSpatialDerivActionAtBndryFaces+1] =
         FixedComposition:new{wall_massf_composition=convertSpeciesTableToArray(o.wall_massf_composition)}
   elseif o.catalytic_type == "equilibrium" then
      o.preSpatialDerivActionAtBndryFaces[#o.preSpatialDerivActionAtBndryFaces+1] = EquilibriumComposition:new{}
   end

   o.preSpatialDerivActionAtBndryFaces[#o.preSpatialDerivActionAtBndryFaces+1] = UpdateThermoTransCoeffs:new()

   if config.turbulence_model ~= "none" then
      o.preSpatialDerivActionAtBndryFaces[#o.preSpatialDerivActionAtBndryFaces+1] = WallTurbulent:new()
   end
   o.postDiffFluxAction = { }
   o.is_configured = true
   return o
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

SolidBIE_TemperatureAndFluxFromSolidGasInterface = SolidBoundaryInterfaceEffect:new{}
SolidBIE_TemperatureAndFluxFromSolidGasInterface.type = "temperature_and_flux_from_gas_solid_interface"
function SolidBIE_TemperatureAndFluxFromSolidGasInterface:tojson()
   local str = string.format('          {"type": "%s"} ', self.type)
   return str
end

SolidBIE_FixedT = SolidBoundaryInterfaceEffect:new{Twall=300.0}
SolidBIE_FixedT.type = "fixed_temperature"
function SolidBIE_FixedT:tojson()
   local str = string.format('          {"type": "%s", ', self.type)
   str = str .. string.format('"Twall": %.18e }', self.Twall)
   return str
end

SolidBIE_CopyAdjacentCellT = SolidBoundaryInterfaceEffect:new{}
SolidBIE_CopyAdjacentCellT.type = "copy_adjacent_cell_temperature"
function SolidBIE_CopyAdjacentCellT:tojson()
   return string.format('          {"type": "%s"} ', self.type)
end

SolidBIE_UserDefined = SolidBoundaryInterfaceEffect:new{fileName='user-defined-solid-bc.lua'}
SolidBIE_UserDefined.type = "user_defined"
function SolidBIE_UserDefined:tojson()
   local str = string.format('          {"type": "%s", ', self.type)
   str = str .. string.format('"filename": "%s" }', self.fileName)
   return str
end

SolidBIE_ConnectionBoundary = SolidBoundaryInterfaceEffect:new{otherBlock=-1,
							       otherFace=nil,
							       orientation=-1}
SolidBIE_ConnectionBoundary.type = "connection_boundary"
function SolidBIE_ConnectionBoundary:tojson()
   local str = string.format('          {"type" : "%s", ', self.type)
   str = str .. string.format('"otherBlock" : %d, ', self.otherBlock)
   str = str .. string.format('"otherFace" : "%s", ', self.otherFace)
   str = str .. string.format('"orientation" : %d ', self.orientation)
   str = str .. ' } '
   return str
end

SolidGhostCellEffect = {
   type = ""
}
function SolidGhostCellEffect:new(o)
   o = o or {}
   setmetatable(o, self)
   self.__index = self
   return o
end

SolidGCE_SolidGhostCellFullFaceCopy = SolidGhostCellEffect:new{otherBlock=-1,
							       otherFace=nil,
							       orientation=-1}
SolidGCE_SolidGhostCellFullFaceCopy.type = "solid_full_face_copy"
function SolidGCE_SolidGhostCellFullFaceCopy:tojson()
   local str = string.format('          {"type" : "%s", ', self.type)
   str = str .. string.format('"otherBlock" : %d, ', self.otherBlock)
   str = str .. string.format('"otherFace" : "%s", ', self.otherFace)
   str = str .. string.format('"orientation" : %d ', self.orientation)
   str = str .. ' } '
   return str
end

SolidGasFullFaceCopy = SolidGhostCellEffect:new{otherBlock=-1,
                                                otherFace=nil,
                                                orientation=-1}
SolidGasFullFaceCopy.type = "solid_gas_full_face_copy"
function SolidGasFullFaceCopy:tojson()
   local str = string.format('          {"type" : "%s", ', self.type)
   str = str .. string.format('"otherBlock" : %d, ', self.otherBlock)
   str = str .. string.format('"otherFace" : "%s", ', self.otherFace)
   str = str .. string.format('"orientation" : %d ', self.orientation)
   str = str .. ' } '
   return str
end

GasSolidFullFaceCopy = SolidGhostCellEffect:new{otherBlock=-1,
                                                otherFace=nil,
                                                orientation=-1}
GasSolidFullFaceCopy.type = "gas_solid_full_face_copy"
function GasSolidFullFaceCopy:tojson()
   local str = string.format('          {"type" : "%s", ', self.type)
   str = str .. string.format('"otherBlock" : %d, ', self.otherBlock)
   str = str .. string.format('"otherFace" : "%s", ', self.otherFace)
   str = str .. string.format('"orientation" : %d ', self.orientation)
   str = str .. ' } '
   return str
end

SolidBoundaryFluxEffect = {
   type = ""
}
function SolidBoundaryFluxEffect:new(o)
   o = o or {}
   setmetatable(o, self)
   self.__index = self
   return o
end

SolidBFE_ZeroFlux = SolidBoundaryFluxEffect:new{}
SolidBFE_ZeroFlux.type = "zero_flux"
function SolidBFE_ZeroFlux:tojson()
   local str = string.format('          {"type": "%s"} ', self.type)
   return str
end

SolidBFE_ConstantFlux = SolidBoundaryFluxEffect:new{fluxValue=0.0}
SolidBFE_ConstantFlux.type = "constant_flux"
function SolidBFE_ConstantFlux:tojson()
   local str = string.format('          {"type": "%s", ', self.type)
   str = str .. string.format('"flux_value": %.18e }', self.fluxValue)
   return str
end

SolidBFE_ConstantFluxFromSolidGasInterface = SolidBoundaryFluxEffect:new{}
SolidBFE_ConstantFluxFromSolidGasInterface.type = "constant_flux_from_solid_gas_interface"
function SolidBFE_ConstantFluxFromSolidGasInterface:tojson()
   local str = string.format('          {"type": "%s"} ', self.type)
   return str
end

SolidBFE_UserDefined = SolidBoundaryFluxEffect:new{fileName='user-defined-solid-bc.lua'}
SolidBFE_UserDefined.type = "user_defined"
function SolidBFE_UserDefined:tojson()
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
   is_solid_domain_bc = true,
   setsFluxDirectly = false,
   preSpatialDerivActionAtBndryCells = {},
   preSpatialDerivActionAtBndryFaces = {},
   postFluxAction = {}
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
   str = str .. string.format('        "sets_flux_directly": %s,\n',
                              tostring(self.setsFluxDirectly))
   str = str .. '        "pre_spatial_deriv_action_at_bndry_cells": [\n'
   for i,effect in ipairs(self.preSpatialDerivActionAtBndryCells) do
      str = str .. effect:tojson()
      -- Extra code to deal with annoying JSON trailing comma deficiency
      if i ~= #self.preSpatialDerivActionAtBndryCells then str = str .. "," end
   end
   str = str .. '\n        ],\n'
   str = str .. '        "pre_spatial_deriv_action_at_bndry_faces": [\n'
   for i,effect in ipairs(self.preSpatialDerivActionAtBndryFaces) do
      str = str .. effect:tojson()
      -- Extra code to deal with annoying JSON trailing comma deficiency
      if i ~= #self.preSpatialDerivActionAtBndryFaces then str = str .. "," end
   end
   str = str .. '\n        ],\n'
   str = str .. '        "post_flux_action": [\n'
   for i,effect in ipairs(self.postFluxAction) do
      str = str .. effect:tojson()
      -- Extra code to deal with annoying JSON trailing comma deficiency
      if i ~= #self.postFluxAction then str = str .. "," end
   end
   str = str .. '\n        ]\n    }'
   return str
end

SolidFixedTBC = SolidBoundaryCondition:new()
SolidFixedTBC.type = "SolidFixedT"
function SolidFixedTBC:new(o)
   o = SolidBoundaryCondition.new(self, o)
   o.preSpatialDerivActionAtBndryFaces = { SolidBIE_FixedT:new{Twall=o.Twall} }
   return o
end

SolidAdiabaticBC = SolidBoundaryCondition:new()
SolidAdiabaticBC.type = "SolidAdiabatic"
function SolidAdiabaticBC:new(o)
   o = SolidBoundaryCondition.new(self, o)
   o.preSpatialDerivActionAtBndryFaces = { SolidBIE_CopyAdjacentCellT:new{} }
   o.postFluxAction = { SolidBFE_ZeroFlux:new{} }
   o.setsFluxDirectly = true
   return o
end

SolidConstantFluxBC = SolidBoundaryCondition:new()
SolidConstantFluxBC.type = "SolidConstantFlux"
function SolidConstantFluxBC:new(o)
   o = SolidBoundaryCondition.new(self, o)
   o.preSpatialDerivActionAtBndryFaces = { SolidBIE_CopyAdjacentCellT:new{} }
   o.postFluxAction = { SolidBFE_ConstantFlux:new{fluxValue=o.fluxValue} }
   o.setsFluxDirectly = true
   return o
end


SolidUserDefinedBC = SolidBoundaryCondition:new()
SolidUserDefinedBC.type = "SolidUserDefined"
function SolidUserDefinedBC:new(o)
   o = SolidBoundaryCondition.new(self, o)
   if (o.specifyFlux) then
      o.postFluxAction = { SolidBFE_UserDefined:new{fileName=o.fileName} }
      o.preSpatialDerivActionAtBndryFaces = { SolidBIE_CopyAdjacentCellT:new{} }
      o.setsFluxDirectly = true
   else
      o.preSpatialDerivActionAtBndryFaces = { SolidBIE_UserDefined:new{fileName=o.fileName} }
   end
   return o
end

SolidAdjacentToGasBC = SolidBoundaryCondition:new()
SolidAdjacentToGasBC.type = "SolidAdjacentToGas"
function SolidAdjacentToGasBC:new(o)
   o = SolidBoundaryCondition.new(self, o)
   o.setsFluxDirectly = true
   o.preSpatialDerivActionAtBndryCells = { SolidGasFullFaceCopy:new{otherBlock=o.otherBlock,
                                                                    otherFace=o.otherFace,
                                                                    orientation=o.orientation}}
   o.preSpatialDerivActionAtBndryFaces = { SolidBIE_TemperatureAndFluxFromSolidGasInterface:new{}}
   return o
end

SolidAdjacentToGasBC2 = SolidBoundaryCondition:new()
SolidAdjacentToGasBC2.type = "SolidAdjacentToGas2"
function SolidAdjacentToGasBC2:new(o)
   o = SolidBoundaryCondition.new(self, o)
   o.preSpatialDerivActionAtBndryCells = { SolidGasFullFaceCopy:new{otherBlock=o.otherBlock,
                                                                    otherFace=o.otherFace,
                                                                    orientation=o.orientation}}
   o.preSpatialDerivActionAtBndryFaces = { SolidBIE_CopyAdjacentCellT:new{} }
   o.postFluxAction = { SolidBFE_ConstantFluxFromSolidGasInterface:new{} }
   o.setsFluxDirectly = true
   return o
end

SolidConnectionBoundaryBC = SolidBoundaryCondition:new()
SolidConnectionBoundaryBC.type = "SolidConnectionBoundary"
function SolidConnectionBoundaryBC:new(o)
   o = SolidBoundaryCondition.new(self, o)
   o.setsFluxDirectly = true
   o.preSpatialDerivActionAtBndryFaces = { SolidBIE_ConnectionBoundary:new{otherBlock=o.otherBlock,
                                                                           otherFace=o.otherFace,
                                                                           orientation=o.orientation} }
   return o
end

SolidFullFaceCopyBoundaryBC = SolidBoundaryCondition:new()
SolidFullFaceCopyBoundaryBC.type = "SolidFullFaceCopyBoundary"
function SolidFullFaceCopyBoundaryBC:new(o)
   o = SolidBoundaryCondition.new(self, o)
   o.setsFluxDirectly = false
   o.preSpatialDerivActionAtBndryCells = { SolidGCE_SolidGhostCellFullFaceCopy:new{otherBlock=o.otherBlock,
                                                                                   otherFace=o.otherFace,
                                                                                   orientation=o.orientation} }
   return o
end

