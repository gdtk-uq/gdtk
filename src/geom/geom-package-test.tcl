#!/usr/bin/tclsh
# geom-package-test.tcl
#
# Testing of the geometry pacakge.
#
# RJG & PJ, 2015-10-24
#

package require tcltest 2.0
namespace import ::tcltest::*
configure -verbose {pass start body error}

set module_names [list vector3 projection properties \
		      line arc bezier helix polyline \
		      coonspatch aopatch \
		      tfivolume sweptsurfacevolume slabvolume wedgevolume \
		      sgrid \
		      univariatefunctions svg ]

foreach name $module_names {
    test $name-test "Testing $name.d" \
	-body "exec [join [list ./ $name _test] {}]" \
	-result {} -returnCodes {0}
}

# test bbla-test {Testing bbla.d} -body {
#     exec ./bbla_test
# } -result {} -returnCodes {0}

puts ""
puts "=======================================  SUMMARY  ======================================="
cleanupTests
puts "========================================================================================="

