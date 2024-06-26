-- Lua containers for a solid thermal model.
--
-- Author: RJG
-- Date: 2024-04-26

SolidThermalModel = {
   type = ""
}
function SolidThermalModel:new(o)
   o = o or {}
   setmetatable(o, self)
   self.__index = self
   return o
end

ConstantPropertiesModel = SolidThermalModel:new{rho=nil, k=nil, Cp=nil}
ConstantPropertiesModel.type = "constant_properties"
function ConstantPropertiesModel:tojson()
   local str = '  {\n'
   str = str .. string.format('  "type": "%s",\n', self.type)
   str = str .. string.format('  "rho": %.18e,\n', self.rho)
   str = str .. string.format('  "k": %.18e,\n', self.k)
   str = str .. string.format('  "Cp": %.18e\n', self.Cp)
   str = str .. '   }'
   return str
end

LinearVariationModel = SolidThermalModel:new{e_ref=0.0, min=nil, max=nil}
LinearVariationModel.type = "linear_variation"
function LinearVariationModel:tojson()
   local str = '  {\n'
   str = str .. string.format('   "type": "%s",\n', self.type)
   str = str .. string.format('   "e_ref": %.18e,\n', self.e_ref)
   str = str .. string.format('   "min": { "T": %.18e, "rho": %.18e, "k": %.18e, "Cp": %.18e },\n',
      self.min.T, self.min.rho, self.min.k, self.min.Cp)
   str = str .. string.format('   "max": { "T": %.18e, "rho": %.18e, "k": %.18e, "Cp": %.18e }\n',
      self.max.T, self.max.rho, self.max.k, self.max.Cp)
   str = str .. '   }'
   return str
end

TabulatedPropertiesModel = SolidThermalModel:new{rho=nil, filename=''}
TabulatedPropertiesModel.type = "tabulated_properties"
function TabulatedPropertiesModel:tojson()
   local str = '  {\n'
   str = str .. string.format('   "type": "%s",\n', self.type)
   str = str .. string.format('   "rho": %.18e,\n', self.rho)
   str = str .. string.format('   "filename": "%s"\n', self.filename)
   str = str .. '   }'
   return str
end

function registerSolidModels(t)
   for k,v in pairs(t) do
      _solidModels[k] = v
   end
end

return {
   ConstantPropertiesModel = ConstantPropertiesModel,
   registerSolidModels = registerSolidModels
}
