# gas_files.cmake

set (GAS_MODEL_FILES_D
  package.d
  co2gas.d
  co2gas_sw.d
  gas_model.d
  gas_state.d
  init_gas_model.d
  ideal_gas.d
  ideal_helium.d
  cubic_gas.d
  cea_gas.d
  physical_constants.d
  sf6virial.d
  therm_perf_gas.d
  therm_perf_gas_equil.d
  very_viscous_air.d
  uniform_lut.d
  uniform_lut_plus_ideal.d
  adaptive_lut_CEA.d
  ideal_air_proxy.d
  ideal_gas_ab.d
  two_temperature_reacting_argon.d
  two_temperature_argon_plus_ideal.d
  ideal_dissociating_gas.d
  two_temperature_air.d
  two_temperature_nitrogen.d
  two_temperature_dissociating_nitrogen.d
  vib_specific_nitrogen.d
  fuel_air_mix.d
  equilibrium_gas.d
  steam.d
  pseudo_species_gas.d
  pseudo_species.d
  electronically_specific_gas.d
  electronic_species.d
  two_temperature_gasgiant.d
  )

set (GAS_MODEL_FILES_F
  ideal_air_fortran.f
  )
