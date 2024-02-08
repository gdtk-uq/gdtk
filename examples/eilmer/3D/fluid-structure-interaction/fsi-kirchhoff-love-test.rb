#! /usr/bin/env ruby
# fsi-euler-bernoulli-test.rb
# Simply test that the beam tip displacement is what we expect it to be.
# Lachlan Whyborn, 08/02/2024

require 'test/unit'
require 'open3'

class TestFSIKirchhoffLove < Test::Unit::TestCase
  def test_0_prep
    cmd = "prep-gas ideal-air.inp ideal-air.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "e4shared --prep --job=kirchhoff_love"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_1_run
    cmd = "e4shared --run --job=kirchhoff_love"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_2_tip_deflection
    # Expected displacement at corners
    expected_displacement = -3.963672e-8

    # Check top and bottom corners
    cmd = "tail -n 1 FSI/hist/0010.dat"
    o, e, s = Open3.capture3(*cmd.split)
    final_vals = o.split(' ')
    computed_displacement = final_vals[1].to_f
    percent_diff = (expected_displacement - computed_displacement) / expected_displacement 
    assert((percent_diff).abs < 1e-3, "Failed to see correct bottom corner displacement.")
    
    # Check that the expected and computed tip displacements are within 0.1%
    cmd = "tail -n 1 FSI/hist/0120.dat"
    o, e, s = Open3.capture3(*cmd.split)
    final_vals = o.split(' ')
    computed_displacement = final_vals[1].to_f
    percent_diff = (expected_displacement - computed_displacement) / expected_displacement 
    assert((percent_diff).abs < 1e-3, "Failed to see correct bottom corner displacement.")
  end
end
