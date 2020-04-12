#! /usr/bin/env ruby
# cone20-test.rb
# Tests for the mabey flat plate implicit solver example.
# KD, 2019-12-17
#
require 'test/unit'
require 'open3'

class TestMabey < Test::Unit::TestCase
  def test_0_prep
    cmd = "prep-gas ideal-air.inp ideal-air-gas-model.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "e4shared --prep --job=mabey"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_1_run
    cmd = "e4-nk-shared --job=mabey --verbosity=1"
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
    assert((steps - 839).abs < 3, "Failed to take correct number of steps.")
  end

  def test_2_drag_force
    cmd = 'python compute_drag.py'
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    lines = o.split("\n")
    force = 0
    lines.each do |txt|
      if txt.match('F_t:') then
        items = txt.split(' ')
        force = items[1].to_f
      end
    end
    assert((force - 25.0356587647).abs < 1.0e-02, "Failed to compute correct drag force.")
  end
end
