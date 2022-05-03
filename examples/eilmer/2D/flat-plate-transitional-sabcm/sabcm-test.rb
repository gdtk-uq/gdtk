#! /usr/bin/env ruby
# sabcm-test.rb
# Tests for the transitional Spalart Allmaras Turbulence model
# NNG, 2020-05-08
#
require 'test/unit'
require 'open3'

class TestSABCM < Test::Unit::TestCase
  def test_0_prep
    cmd = "prep-gas ideal-air.inp ideal-air-gas-model.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "e4shared --prep --job=fp"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_1_run
    cmd = "e4-nk-shared --job=fp --verbosity=1"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    steps = 0
    lines = o.split("\n")
    fstep_line_num = 0
    line_count = 0
    lines.each do |txt|
      line_count += 1
      if txt.match('STOPPING') then
        fstep_line_num = line_count + 2
      end
      if line_count == fstep_line_num then
        items = txt.split(' ')
        steps = items[1].to_i
      end
    end
    assert((steps - 1131).abs < 40, "Failed to take correct number of steps.")
  end

  def test_2_drag_force
    cmd = 'python3 compute_transition_location.py fp'
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    lines = o.split("\n")
    location = 0
    lines.each do |txt|
      if txt.match('Transition Location:') then
        items = txt.split(' ')
        location = items[2].to_f
      end
    end
    assert((location - 0.2578125).abs < 1.0e-03, "Failed to compute correct transition location.")
  end
end
