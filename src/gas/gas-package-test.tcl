#!/usr/bin/tclsh
# gas-package-test.tcl
#
# Testing of the gas pacakge.
#
# RJG, 2015-10-23
#
# This script can be started by hand, but usually we will
# get the makefile to execute this when issued with a
# 'make test'. To execute by hand:
#    1. ./gas-package-test.tcl
#    2. tclsh gas-package-test.tcl

package require tcltest 2.0
namespace import ::tcltest::*
configure -verbose {pass start body error}

# Tests in top-level gas/ area
puts "-------------------------------------"
puts "   Top-level gas/ tests"
puts "-------------------------------------"

test init-gas-model-test {Testing init_gas_model.d} -body {
    exec ./init_gas_model_test
} -result {} -returnCodes {0}

test init-gas-model-complex-test {Testing init_gas_model.d complex flavour} -body {
    exec ./init_gas_model_complex_test
} -result {} -returnCodes {0}

test ideal-gas-test {Testing ideal_gas.d} -body {
    exec ./ideal_gas_test
} -result {} -returnCodes {0}

test ideal-gas-complex-test {Testing ideal_gas.d complex flavour} -body {
    exec ./ideal_gas_complex_test
} -result {} -returnCodes {0}

test ideal-helium-test {Testing ideal_helium.d} -body {
    exec ./ideal_helium_test
} -result {} -returnCodes {0}

test cubic-gas-test {Testing cubic_gas.d} -body {
    exec ./cubic_gas_test
} -result {} -returnCodes {0}

test cea-gas-test {Testing cea_gas.d} -body {
    exec ./cea_gas_test
} -result {} -returnCodes {0}

test composite-gas-test {Testing composite_gas.d} -body {
    exec ./composite_gas_test
} -result {} -returnCodes {0}

test therm-perf-gas-test {Testing therm_perf_gas.d} -body {
    exec ./therm_perf_gas_test
} -result {} -returnCodes {0}

test therm-perf-gas-complex-test {Testing therm_perf_gas.d complex flavour} -body {
    exec ./therm_perf_gas_complex_test
} -result {} -returnCodes {0}

test therm-perf-gas-equil-test {Testing therm_perf_gas_equil.d} -body {
    exec ./therm_perf_gas_equil_test > LOGFILE_THERM_PERF_EQUIL_TEST
} -result {} -returnCodes {0}

test very-viscous-air-test {Testing very_viscous_air.d} -body {
    exec ./very_viscous_air_test
} -result {} -returnCodes {0}

test uniform-lut-test {Testing uniform_lut.d} -body {
    exec ./uniform_lut_test > LOGFILE_UNIFORM_LUT_TEST
} -result {} -returnCodes {0}

test adaptive-lut-test {Testing adaptive_lut.d} -body {
    exec ./adaptive_lut_CEA_test
} -result {} -returnCodes {0}

test uniform-lut-plus-ideal-test {Testing uniform_lut_plus_ideal.d} -body {
    exec cp sample-data/uniform-lut-plus-ideal-air-gas-model.lua ./uniform-lut-plus-ideal-air-gas-model.lua
    exec cp sample-data/cea-lut-air-version-test.lua ./cea-lut-air-version-test.lua
    exec cp sample-data/ideal-air-gas-model.lua ./ideal-air-gas-model.lua
    exec ./uniform_lut_plus_ideal_test > LOGFILE_UNIFORM_LUT_PLUS_IDEAL_TEST
} -result {} -returnCodes {0}

test ideal-gas-ab-test {Testing ideal_gas_ab.d} -body {
    exec ./ideal_gas_ab_test
} -result {} -returnCodes {0}

test vib-specific-nitrogen-test {Testing vib_specific_nitrogen.d} -body {
    exec ./vib_specific_nitrogen_test
} -result {} -returnCodes {0}

test vib-specific-co-test {Testing vib-specific-co.d} -body {
    exec ./vib_specific_co_test
} -result {} -returnCodes {0}

test vib-specific-co-mixture-test {Testing vib-specific-co.d} -body {
    exec ./vib_specific_co_mixture_test
} -result {} -returnCodes {0}


test ideal-dissociating-gas-test {Testing ideal_dissociating_gas.d} -body {
    exec ./ideal_dissociating_gas_test
} -result {} -returnCodes {0}

test two-temperature-reacting-argon-test {Testing Daniel's two-T argon gas model} -body {
    exec ./two_temperature_reacting_argon_test
} -result {} -returnCodes {0}

test two-temperature-argon-plus-ideal-test {Testing two_temperature_argon_plus_ideal.d} -body {
    exec cp sample-data/two-temperature-argon-plus-ideal-air-gas-model.lua ./two-temperature-argon-plus-ideal-air-gas-model.lua
    exec cp sample-data/two-temperature-reacting-argon-model.lua ./two-temperature-reacting-argon-model.lua
    exec cp sample-data/ideal-air-gas-model.lua ./ideal-air-gas-model.lua
    exec ./two_temperature_argon_plus_ideal_test > LOGFILE_TWO_TEMPERATURE_ARGON_PLUS_IDEAL_TEST
} -result {} -returnCodes {0}

#test fuel_air_mix-test {Testing fuel_air_mix.d} -body {
#    exec ./fuel_air_mix_test
#} -result {} -returnCodes {0}

test equilibrium-gas-test {Testing equilibrium_gas.d} -body {
    exec ./equilibrium_gas_test
} -result {} -returnCodes {0}

test two-temperature-gasgiant-test {Testing Daisy's two-T H2-He gas model} -body {
    exec ./two_temperature_gasgiant_test
} -result {} -returnCodes {0}

puts "-------------------------------------"
puts "   thermo/ tests"
puts "-------------------------------------"

test cea-thermo-curves-test {Testing thermo/cea_thermo_curves.d} -body {
    exec ./cea_thermo_curves_test
} -result {} -returnCodes {0}

test cea-thermo-curves-complex-test {Testing thermo/cea_thermo_curves_complex.d} -body {
    exec ./cea_thermo_curves_complex_test
} -result {} -returnCodes {0}

test perf-gas-mix-eos-test {Testing thermo/perf_gas_mix_eos.d} -body {
    exec ./perf_gas_mix_eos_test
} -result {} -returnCodes {0}

test therm-perf-gas-mix-eos-test {Testing thermo/therm_perf_gas_mix_eos.d} -body {
    exec ./therm_perf_gas_mix_eos_test
} -result {} -returnCodes {0}

puts "-------------------------------------"
puts "   diffusion/ tests"
puts "-------------------------------------"

test cea-therm-cond-test {Testing diffusion/cea_therm_cond.d} -body {
    exec ./cea_therm_cond_test
} -result {} -returnCodes {0}

test cea-viscosity-test {Testing diffusion/cea_viscosity.d} -body {
    exec ./cea_viscosity_test
} -result {} -returnCodes {0}

test sutherland-therm-cond-test {Testing diffusion/sutherland_therm_cond.d} -body {
    exec ./sutherland_therm_cond_test
} -result {} -returnCodes {0}

test sutherland-viscosity-test {Testing diffusion/sutherland_viscosity.d} -body {
    exec ./sutherland_viscosity_test
} -result {} -returnCodes {0}

test wilke-mixing-therm-cond-test {Testing diffusion/wilke_mixing_therm_cond.d} -body {
    exec ./wilke_mixing_therm_cond_test
} -result {} -returnCodes {0}

test wilke-mixing-viscosity-test {Testing diffusion/wilke_mixing_viscosity.d} -body {
    exec ./wilke_mixing_viscosity_test
} -result {} -returnCodes {0}


puts "-----------------------------------------"
puts "   Lua wrapped functions (gas-calc) test "
puts "-----------------------------------------"

test gas-calc-test {Testing gas-calc} -body {
    exec ./gas-calc wrapped-gas-model-test.lua
} -result {} -returnCodes {0}

puts ""
puts "=====================================  SUMMARY  ====================================="
cleanupTests
puts "====================================================================================="

