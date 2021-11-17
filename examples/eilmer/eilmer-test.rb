#! /usr/bin/env ruby
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
# 1. ./eilmer-test.rb
# 2. ruby eilmer-test.rb
#
# PJ, 2011-01-11 Eilmer3 version
#     2015-10-22 Port to Eilmer4
#     2019-05-16 Ruby version.
#     2021-11-17 refactored to work with a file of test names.

require 'date'

gpmetis_exe = `which gpmetis`
found_gpmetis = gpmetis_exe.length > 0

f = File.new('test-names.txt', 'r')
lines = f.readlines()
f.close
test_scripts = []
lines.each do |txt|
  name = txt.chomp
  test_scripts << name if ((not name.include?("metis")) or found_gpmetis)
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
