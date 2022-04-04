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

set module_names [list vector3 vector3_complex \
                      projection projection_complex \
                      properties properties_complex \
		      line line_complex \
                      arc arc_complex \
                      bezier bezier_complex \
                      helix helix_complex \
                      polyline polyline_complex \
                      xspline \
                      xsplinelsq \
                      svgpath \
		      coonspatch coonspatch_complex \
                      aopatch aopatch_complex \
                      gmopatch \
                      bezierpatch bezierpatch_complex \
                      controlpointpatch \
                      beziertrianglepatch \
		      spherepatch spherepatch_complex \
		      cubepatch cubepatch_complex \
		      tfivolume tfivolume_complex \
                      sweptsurfacevolume sweptsurfacevolume_complex \
                      twosurfacevolume twosurfacevolume_complex \
                      wedgevolume wedgevolume_complex \
                      slabvolume slabvolume_complex \
		      sgrid sgrid_complex \
		      usgrid usgrid_complex \
		      univariatefunctions svg ]
# [TODO] PJ 2021-12-14: get gmopatch_complex test working.

foreach name $module_names {
    test ${name}_test "Testing $name.d" \
	-body "exec [join [list ./ $name _test] {}]" \
	-result {} -returnCodes {0}
}

# test bbla_test {Testing bbla.d} -body {
#     exec ./bbla_test
# } -result {} -returnCodes {0}

puts ""
puts "=======================================  SUMMARY  ======================================="
cleanupTests
puts "========================================================================================="

