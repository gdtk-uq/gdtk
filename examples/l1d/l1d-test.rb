#!/usr/bin/env ruby
# l1d-test.rb
#
# Smoke test the L1d4 code.
# Presently, we work through the specified directories and
# explicitly invoke the test scripts that reside there.
#
# Usage:
# 1. ./l1d-test.rb
# 2. ruby l1d-test.rb
#
# PJ, 2020-05-21, adapted from eilmer-test.rb

require 'date'

test_scripts = []
test_scripts << "sod/sod-test.rb"
test_scripts << "sod-he-air/sod-he-air-test.rb"
test_scripts << "piston/piston-test.rb"
test_scripts << "piston-with-valve/piston-with-valve-test.rb"

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
