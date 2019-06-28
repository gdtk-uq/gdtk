#! /usr/bin/env ruby
# ramp-test.rb
# Functional test for a 3D, inviscid flow.
# This exercises quite a few of the basic functions of the 3D code.
# PJ, 2015-10-26, 2019-06-24
#
require 'test/unit'
require 'open3'

class TestRamp < Test::Unit::TestCase
  def test_0_prep
    cmd = "prep-gas ideal-air.inp ideal-air-gas-model.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "e4shared --prep --job=ramp"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_1_run
    cmd = "e4shared --run --job=ramp --verbosity=1 --max-cpus=2"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    steps = 0
    lines = o.split("\n")
    lines.each do |txt|
      if txt.match('final-t=') then
        items = txt.split(' ')
        steps = items[1].to_i
      end
    end
    assert((steps - 862).abs < 3, "Failed to take correct number of steps.")
  end

  def test_2_shock_angle
    cmd = 'e4shared --custom-script --script-file=estimate_shock_angle.lua'
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    shock_angle = 0.0
    average_deviation = 1000.0 # something big
    lines = o.split("\n")
    lines.each do |txt|
      if txt.match('shock_angle_deg') then
        items = txt.split(' ')
        shock_angle = items[1].to_f
      end
      if txt.match('average_deviation_metres') then
        items = txt.split(' ')
        average_deviation = items[1].to_f
      end
    end
    assert((shock_angle - 57.0).abs < 1.0, "Failed to get correct shock angle.")
    assert((average_deviation).abs < 0.002, "Failed to find a straight shock.")
  end

  def test_3_estimate_forces
    cmd = 'e4shared --custom-script --script-file=estimate_ramp_force.lua'
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    f_x = 0.0; f_y = 0.0; f_z = 0.0
    lines = o.split("\n")
    lines.each do |txt|
      if txt.match('force=') then
        # expected line looks like
        # force=	Vector3([2217, 3.92703e-14, -12573.3])	Newtons
        str = txt.match(/\[(.+)\]/)[1]
        items = str.split(',')
        f_x = items[0].to_f; f_y = items[1].to_f; f_z = items[2].to_f
      end
    end
    assert((f_x - 2215.0)/2215.0 < 0.003, "Failed to get correct f_x.")
    assert(f_y.abs < 0.001, "Failed to get correct f_y.")
    assert((f_z + 12560.0)/12560.0 < 0.002, "Failed to get correct f_z.")
  end
end
