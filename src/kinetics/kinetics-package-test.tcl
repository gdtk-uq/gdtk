#!/usr/bin/tclsh
# Testing of the kinetics pacakge.
#
# RJG, 2015-10-23
#

package require tcltest 2.0
namespace import ::tcltest::*
configure -verbose {pass start body error}

test chemistry-update-test {Testing chemistry_update.d} -body {
    exec ./chemistry_update_test
} -result {} -returnCodes {0}

test equilibrium-update-test {Testing equilibrium_update.d} -body {
    exec ./equilibrium_update_test
} -result {} -returnCodes {0}

test rate-constant-test {Testing rate_constant.d} -body {
    exec ./rate_constant_test
} -result {} -returnCodes {0}

test reaction-test {Testing reaction.d} -body {
    exec ./reaction_test
} -result {} -returnCodes {0}

test reaction-mechanism-test {Testing reaction_mechanism.d} -body {
    exec ./reaction_mechanism_test
} -result {} -returnCodes {0}

test two-temperature-argon-kinetics-test {Testing Daniel Smith's two-temperature argon reaction mechanism.} -body {
    exec ./two_temperature_argon_kinetics_test
} -result {} -returnCodes {0}

test vib-specific-CO-kientics-test {Testing Nick Gibbons implementation of a state-specific carbon monoxide model.} -body {
    exec ./two_temperature_argon_kinetics_test
} -result {} -returnCodes {0}

test two-temperature-argon-with-ideal-gas-test {Testing Daniel Smith's two-temperature argon reaction mechanism with ideal-gas present.} -body {
    exec ./two_temperature_argon_with_ideal_gas_test
} -result {} -returnCodes {0}

test two-temperature-gasgiant-kinetics-test {Testing Daisy's two-temperature Gas-Giant reaction mechanism.} -body {
    exec ./two_temperature_gasgiant_kinetics_test
} -result {} -returnCodes {0}

puts ""
puts "==========================================  SUMMARY  =========================================="
cleanupTests
puts "==============================================================================================="

