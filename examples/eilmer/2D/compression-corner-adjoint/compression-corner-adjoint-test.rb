#! /usr/bin/env ruby
# compression-corner-test.rb
# Tests for the compression corner adjoint solver example.
# KD, 2019-12-18
#
require 'test/unit'
require 'open3'

class TestCompCorner < Test::Unit::TestCase
  def test_0_prep
    cmd = "prep-gas ideal-air.inp ideal-air-gas-model.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "e4shared --prep --job=comp-corner"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_1_run
    cmd = "e4-nk-shared --job=comp-corner --verbosity=1"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "e4ssc --job=comp-corner --adjoint-method --adjoint-verification=true"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    steps = 0
    lines = o.split("\n")
    fstep_line_num = 0
    line_count = 0
    lines.each do |txt|
      if txt.match('restarts') then
        items = txt.split(' ')
        steps = items[0].to_i
      end
    end
    assert((steps - 9).abs < 1, "Failed to take correct number of restarts.")
  end

  def test_2_adjoint
    cmd = 'e4ssc --job=comp-corner --adjoint-method --adjoint-verification=true'
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    lines = o.split("\n")
    grad1 = 0
    grad2 = 0
    lines.each do |txt|
      if txt.match('gradient for variable 1:') then
        items = txt.split(':')
        grad1 = items[1].to_f
      end
      if txt.match('gradient for variable 2:') then
        items = txt.split(':')
        grad2 = items[1].to_f
      end
    end
    assert((grad1 - -1.8577754782790569e+02).abs < 1.0e-10, "Failed to compute correct gradient for variable 1.")
    assert((grad2 - -5.0900383692758169e+00).abs < 1.0e-10, "Failed to compute correct gradient for variable 2.")
  end

end
