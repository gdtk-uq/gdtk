#! /usr/bin/env ruby
# plate-test.rb
# Tests for the solid plate example.
# KAD, 2020-09-24
#
require 'test/unit'
require 'open3'

class TestPlate < Test::Unit::TestCase
  def test_0_prep
    cmd = "prep-gas ideal-air.inp ideal-air-gas-model.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "e4shared --prep --job=plate"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_1_run
    cmd = "e4shared --run --job=plate"
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
    assert((steps - 1000).abs < 3, "Failed to take correct number of steps.")
  end

  def test_2_post_temperature_points
    cmd = 'e4shared --post --job=plate --tindx-plot=last --probe=0.05,0.01,0.0'
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    lines = o.split("\n")
    temp = 0.0
    solid_data_line = false
    lines.each do |txt|
      if (solid_data_line)
        items = txt.split(' ')
        temp = items[5].to_f
        solid_data_line = false # we have passed the solid data line
      end
      if txt.match('#')
        solid_data_line = true # the next line is the solid data
      end
    end
    assert((temp - 606.56).abs < 1.0, "Failed to see correct temperature.")
  end
end
