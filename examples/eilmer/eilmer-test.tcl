#!/usr/bin/tclsh
# eilmer-test.tcl
#
# Smoke test the Eilmer4 code.
#
# PJ, 11-Jan-2011, 12-Jul-2011 (Eilmer3 version)
#     2015-10-22 Port to Eilmer4
# 
# Presently, we work through the specified directories and explicitly invoke 
# test scripts.  Of course, there must be a better way to do this using the 
# tcltest module.
# Short tests are those which may take up to 20 minutes on my workstation.
# I don't want to wait more than an hour, or so, for the set of short tests
# to run.
#
# Usage:
# 1. ./eilmer-text.tcl
# 2. tclsh eilmer-test.tcl
# 3. tclsh eilmer-test.tcl --long-tests
# 4. ./eilmer-test.tcl --dummy-run

set long_tests 0
set dummy_run 0
for {set i 0} {$i < $argc} {incr i} {
    set arg [lindex $argv $i]
    if {[string first "long" $arg] >= 0} {
        set long_tests 1
    }
    if {[string first "dummy" $arg] >= 0} {
        set dummy_run 1
    }
}
set test_scripts [list "2D/sharp-cone-20-degrees/sg/cone20.test"]
lappend test_scripts "2D/sharp-cone-20-degrees/usg/cone20-usg.test"
lappend test_scripts "3D/sod-shock-tube/sg/sod.test"
lappend test_scripts "3D/sod-shock-tube/usg/sod.test"
# lappend test_scripts "2D/sod/N2-O2/sod.test"
lappend test_scripts "2D/sphere-sawada/fixed-grid/ss3.test"
lappend test_scripts "2D/nozzle-conical-back/back.test"
# lappend test_scripts "2D/methane-reactor/psr.test"
lappend test_scripts "2D/channel-with-bump/bump.test"
lappend test_scripts "2D/manufactured-solution/sg/smoke-tests/mms-euler.test"
lappend test_scripts "2D/manufactured-solution/sg/smoke-tests/mms-ns-div-theorem.test"
lappend test_scripts "2D/manufactured-solution/sg/smoke-tests/mms-ns-least-sq-at-vtxs.test"
lappend test_scripts "2D/manufactured-solution/sg/smoke-tests/mms-ns-least-sq-at-faces.test"
lappend test_scripts "2D/manufactured-solution/usg/mms-euler.test"
lappend test_scripts "2D/cht-manufactured-solution/spatial-verification/smoke-tests/single-thread.test"
# lappend test_scripts "2D/unstructured-manufactured-solution/mms-ns.test"
lappend test_scripts "2D/cylinder-shockfitting/cyl-sf.test"
# lappend test_scripts "2D/odw/odw.test"
lappend test_scripts "2D/duct-hydrogen-combustion/bittker.test"
# lappend test_scripts "2D/radiating-cylinder/gray-gas/MC/cyl.test"
# lappend test_scripts "2D/giordano/inf_cyl.test"
lappend test_scripts "3D/simple-ramp/sg/ramp.test"
lappend test_scripts "2D/cylinder-dlr-n90/n90.test"
#lappend test_scripts "2D/sphere-lobb/smoke-test/lobb.test"
# lappend test_scripts "2D/radiating-cylinder/Argon/MC/cyl.test"
# lappend test_scripts "2D/nenzfr-Mach4-nozzle-eq/nozzle-eq-marching.test"
# lappend test_scripts "2D/nenzfr-Mach4-nozzle-eq/nozzle-eq.test"
if {$long_tests} {
    puts "Do long tests as well as short tests..."
    # lappend test_scripts "2D/turb-flat-plate/turb_flat_plate.test"
    # lappend test_scripts "2D/nenzfr-Mach4-nozzle-noneq/nozzle-noneq.test"
    # lappend test_scripts "2D/Rutowski-hemisphere/Ms_12.70/Rutowski-short.test"
} else {
    puts "Do short tests only..."
}
set original_dir [pwd]
foreach test_script $test_scripts {
    cd [file dir $test_script]
    puts "[exec date] [pwd]"
    if { !$dummy_run } {
        source [file tail $test_script]
    }
    cd $original_dir
}
puts "[exec date] Finished tests."
