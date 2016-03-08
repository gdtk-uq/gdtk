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

test gas-model-test {Testing gas_model.d} -body {
    exec ./gas_model_test
} -result {} -returnCodes {0}

test gas-model-util-test {Testing gas_model_util.d} -body {
    exec ./gas_model_util_test
} -result {} -returnCodes {0}

test ideal-gas-test {Testing ideal_gas.d} -body {
    exec ./ideal_gas_test
} -result {} -returnCodes {0}

test therm-perf-gas-test {Testing therm_perf_gas.d} -body {
    exec ./therm_perf_gas_test
} -result {} -returnCodes {0}

test very-viscous-air-test {Testing very_viscous_air.d} -body {
    exec ./very_viscous_air_test
} -result {} -returnCodes {0}

test co2gas-sw-test {Testing co2gas_sw.d} -body {
    exec ./co2gas_sw_test > LOGFILE_CO2GAS_SW_TEST
} -result {} -returnCodes {0}

test uniform-lut-test {Testing uniform_lut.d} -body {
    exec ./uniform_lut_test > LOGFILE_UNIFORM_LUT_TEST
} -result {} -returnCodes {0}

test adaptive-lut-test {Testing uniform_lut.d} -body {
    exec ./adaptive_lut_CEA_test
} -result {} -returnCodes {0}

puts "-------------------------------------"
puts "   thermo/ tests"
puts "-------------------------------------"

test cea-thermo-curves-test {Testing thermo/cea_thermo_curves.d} -body {
    exec ./cea_thermo_curves_test
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

test sutherland-therm-cond-test {Testing diffusion/sutherland_therm_cone.d} -body {
    exec ./sutherland_therm_cond_test
} -result {} -returnCodes {0}

test sutherland-viscosity-test {Testing diffusion/sutherland_viscosity.d} -body {
    exec ./sutherland_viscosity_test
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

