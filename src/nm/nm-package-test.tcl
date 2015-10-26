#!/usr/bin/tclsh
# Testing of the nm pacakge.
#
# RJG & PJ, 2015-10-24
#

package require tcltest 2.0
namespace import ::tcltest::*
configure -verbose {pass start body error}

test bbla-test {Testing bbla.d} -body {
    exec ./bbla_test
} -result {} -returnCodes {0}

test bracketing-test {Testing bracketing.d} -body {
    exec ./bracketing_test
} -result {} -returnCodes {0}

test gaussquad-test {Testing gaussquad.d} -body {
    exec ./gaussquad_test
} -result {} -returnCodes {0}

test linesearch-test {Testing linesearch.d} -body {
    exec ./linesearch_test
} -result {} -returnCodes {0}

test newtoncotes-test {Testing newtoncotes.d} -body {
    exec ./newtoncotes_test
} -result {} -returnCodes {0}

test ridder-test {Testing ridder.d} -body {
    exec ./ridder_test
} -result {} -returnCodes {0}

test brent-test {Testing brent.d} -body {
    exec ./brent_test
} -result {} -returnCodes {0}

test rungekutta-test {Testing rungekutta.d} -body {
    exec ./rungekutta_test
} -result {} -returnCodes {0}

puts ""
puts "=======================================  SUMMARY  ======================================="
cleanupTests
puts "========================================================================================="

