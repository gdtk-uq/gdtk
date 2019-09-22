#!/usr/bin/env ruby
# eilmer-test.rb
#
# Smoke test the Eilmer4 code.
# 
# Presently, we work through the specified directories and
# explicitly invoke the test scripts that reside there. 
# Short tests are those which may take up to 20 minutes on my workstation.
# I don't want to wait more than an hour, or so, for the set of short tests
# to run.
#
# Usage:
# 1. ./eilmer-text.rb
# 2. ruby eilmer-test.rb
#
# PJ, 2011-01-11 Eilmer3 version
#     2015-10-22 Port to Eilmer4
#     2019-05-16 Ruby version.

require 'date'

long_tests = false

test_scripts = []
test_scripts << "2D/sharp-cone-20-degrees/sg/cone20-test.rb"
test_scripts << "2D/sharp-cone-20-degrees/sg-mpi/cone20-mpi-test.rb"
test_scripts << "2D/sharp-cone-20-degrees/usg/cone20-usg-test.rb"
test_scripts << "2D/sharp-cone-20-degrees/usg-su2/cone20-usg-su2-test.rb"
gpmetis_exe = `which gpmetis`
if gpmetis_exe.length > 0 then
  puts "Found gpmetis"
  test_scripts << "2D/sharp-cone-20-degrees/usg-metis/cone20-usg-metis-test.rb"
end
test_scripts << "2D/moving-grid/piston-w-const-vel/simple/piston-test.rb"
test_scripts << "2D/moving-grid/piston-w-const-vel/simple-2/piston-test.rb"
test_scripts << "2D/moving-grid/piston-w-const-vel/simple-reverse/piston-test.rb"
test_scripts << "2D/moving-grid/piston-in-tube/piston-1-block/pit1-test.rb"
test_scripts << "2D/moving-grid/piston-in-tube/piston-2-block/pit2-test.rb"
test_scripts << "3D/sod-shock-tube/sg/sod-test.rb"
test_scripts << "3D/sod-shock-tube/usg/sod-usg-test.rb"
test_scripts << "3D/sod-shock-tube/ideal-air-eq-air/sod-test.rb"
test_scripts << "3D/connection-test/connection-shared-memory-test.rb"
test_scripts << "2D/vortex/vtx-test.rb"
test_scripts << "2D/reactor-n2/reactor-test.rb"
test_scripts << "2D/sphere-sawada/fixed-grid/ss3-test.rb"
test_scripts << "2D/nozzle-conical-back/back-test.rb"
test_scripts << "2D/channel-with-bump/bump-test.rb"
test_scripts << "2D/manufactured-solution/sg/smoke-tests/mms-euler-test.rb"
test_scripts << "2D/manufactured-solution/sg/smoke-tests/mms-ns-div-theorem-test.rb"
test_scripts << "2D/manufactured-solution/sg/smoke-tests/mms-ns-least-sq-at-vtxs-test.rb"
test_scripts << "2D/manufactured-solution/sg/smoke-tests/mms-ns-least-sq-at-faces-test.rb"
test_scripts << "2D/manufactured-solution/usg/mms-euler-test.rb"
test_scripts << "2D/cht-manufactured-solution/spatial-verification/smoke-tests/single-thread-test.rb"
test_scripts << "2D/shock-fitting/cylinder/cyl-sf-test.rb"
test_scripts << "2D/shock-fitting/cylinder-mpi/cyl-sf-mpi-test.rb"
test_scripts << "2D/oblique-detonation-wave/odw-test.rb"
test_scripts << "2D/duct-hydrogen-combustion/bittker-test.rb"
test_scripts << "2D/cylinder-giordano/two-temperature/inf_cyl-test.rb"
test_scripts << "3D/simple-ramp/sg/ramp-test.rb"
test_scripts << "2D/cylinder-dlr-n90/cpu-chem/n90-test.rb"
test_scripts << "2D/binary-diffusion/bd-test.rb"
# test_scripts << "2D/sphere-lobb/smoke-test/lobb.test"
# test_scripts << "2D/radiating-cylinder/Argon/MC/cyl.test"
test_scripts << "2D/nozzle-shock-tunnel-t4m4/t4m4-test.rb"
if long_tests then
  puts "Do long tests as well as short tests..."
  # test_scripts << "2D/turb-flat-plate/turb_flat_plate.test"
  # test_scripts << "2D/nenzfr-Mach4-nozzle-noneq/nozzle-noneq.test"
  # test_scripts << "2D/Rutowski-hemisphere/Ms_12.70/Rutowski-short.test"
else
  puts "Do short tests only..."
end

time_start = DateTime.now
original_dir = Dir.getwd
test_scripts.each do |ts|
  Dir.chdir(File.dirname(ts))
  puts "#{DateTime.now} #{Dir.getwd}"
  case File.extname(ts)
  when ".rb"
    cmd = "ruby"
  when ".test", ".tcl"
    cmd = "tclsh"
  else
    puts "Dodgy extension for test script."
    cmd = ""
  end
  if cmd.length > 0
    cmd << " " + File.basename(ts)
    puts "cmd= " + cmd.to_s
    system(cmd)
  end
  Dir.chdir(original_dir)
end
time_end = DateTime.now
puts "#{time_end} Finished tests."
delta = time_end - time_start
h = (delta*24).to_i
m = (delta*24*60).to_i - h*60
s = (delta*24*3600).to_i - h*3600 - m*60
puts "Elapsed time: #{h} hr, #{m} min, #{s} sec."
