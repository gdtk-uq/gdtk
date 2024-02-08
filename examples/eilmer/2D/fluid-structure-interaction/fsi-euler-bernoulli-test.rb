#! /usr/bin/env ruby
# fsi-euler-bernoulli-test.rb
# Simply test that the beam tip displacement is what we expect it to be.
# Lachlan Whyborn, 08/02/2024

require 'test/unit'
require 'open3'

class TestFSIEulerBernoulli < Test::Unit::TestCase
  def test_0_prep
    cmd = "prep-gas ideal-air.inp ideal-air.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "e4shared --prep --job=euler_bernoulli"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_1_run
    cmd = "e4shared --run --job=euler_bernoulli"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_2_tip_deflection
    cmd = "tail -n 1 FSI/hist/0010.dat"
    o, e, s = Open3.capture3(*cmd.split)
    final_vals = o.split(' ')
    # Check that the expected and computed tip displacements are within 0.1%
    expected_displacement = -6e-8
    computed_displacement = final_vals[1].to_f
    percent_diff = (expected_displacement - computed_displacement) / expected_displacement
    assert((percent_diff).abs < 1e-3, "Failed to see correct tip displacement.")
  end
end
