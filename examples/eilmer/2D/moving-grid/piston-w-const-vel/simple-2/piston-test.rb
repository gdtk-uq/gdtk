#! /usr/bin/env ruby
# piston-test.rb
# Tests for the piston-with-constant-velocity example.
# PJ, 2019-05-22
#
require 'test/unit'
require 'open3'

class TestPiston < Test::Unit::TestCase
  def test_0_prep
    cmd = "prep-gas ideal-air.inp ideal-air-gas-model.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "e4shared --prep --job=piston"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_1_run
    cmd = "e4shared --run --job=piston"
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
    assert((steps - 506).abs < 10, "Failed to take correct number of steps.")
  end

  def test_2_post
    # Check post-shock conditions.
    cmd = 'e4shared --post --job=piston --probe=0.28,0,0'
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    # Sift through the output from probe, looking for our data line.
    # The data line that we want comes immediately after a comment line
    # that specifies the varible names.
    posx = nil; velx = nil; p = nil; temp=nil
    lines = o.split("\n")
    found_comment_line = false
    lines.each do |txt|
      if found_comment_line then
        # Next line is my data line and...
        items = txt.split(' ')
        posx = items[0].to_f; velx = items[5].to_f
        p = items[8].to_f; temp = items[19].to_f
        # once I have extracted values from it, get out.
        break
      end
      if txt.match('^# ') then
        found_comment_line = true
      end
    end
    assert((posx - 0.283).abs < 1.0e-2, "Failed to get correct location.")
    assert((velx - 293.0).abs < 1.0, "Failed to see correct velocity.")
    assert((p - 30.257e3).abs < 500.0, "Failed to see correct pressure.")
    assert((temp - 398.5).abs < 1.0, "Failed to see correct temperature.")
  end
end
