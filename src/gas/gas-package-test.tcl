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

test gas-model-test {Testing gas_model.d} -body {
    exec ./gas_model_test
} -result {} -returnCodes {0}

test ideal-gas-test {Testing ideal_gas.d} -body {
    exec ./ideal_gas_test
} -result {} -returnCodes {0}

cleanupTests
