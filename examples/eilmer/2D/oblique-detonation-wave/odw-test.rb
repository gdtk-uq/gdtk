#! /usr/bin/env ruby
# odw-test.rb
# Functional test for the Powers and Aslam oblique detonation wave.
# This exercises quite a few of the basic functions of the 2D code
# plus the generic connection to a custom reacting gas.
# PJ, 2019-06-28
#
require 'test/unit'
require 'open3'

class TestODW < Test::Unit::TestCase
  def test_0_prep
    cmd = "e4shared --prep --job=odw"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_1_run
    cmd = "e4shared --run --job=odw"
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
    assert((steps - 2113).abs < 6, "Failed to take correct number of steps.")
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
    assert((shock_angle - 45.0).abs < 0.5, "Failed to get correct shock angle.")
    assert((average_deviation).abs < 0.002, "Failed to find a straight shock.")
  end
end
