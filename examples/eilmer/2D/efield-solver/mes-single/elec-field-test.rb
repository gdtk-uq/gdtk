#! /usr/bin/env ruby
# elec-field-test.rb
# Tests for the electric field solver, in single block mode, using an exact solution.
# NNG, 2021-11-15
#
require 'test/unit'
require 'open3'

class TestElec < Test::Unit::TestCase
  def test_0_prep
    cmd = "prep-gas ideal-air.inp ideal-air-gas-model.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "e4shared --prep --job=elec"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_1_run
    cmd = "e4shared --run --job=elec --max-cpus=1"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_2_exact_solution
    cmd = 'python3 mes.py'
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    lines = o.split("\n")
    rms = 0
    lines.each do |txt|
      if txt.match('RMS Error:') then
        items = txt.split(' ')
        rms = items[2].to_f
      end
    end
    assert((rms - 0.00026664520908412055).abs < 1.0e-06, "Failed to compute correct RMS.")
  end
end
