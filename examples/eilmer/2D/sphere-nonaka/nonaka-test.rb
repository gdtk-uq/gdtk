#! /usr/bin/env ruby
# nonaka-test.rb
# This test covers the new multi-temperature relaxation capability, using a 2T gas solved over a shock-fitted sphere.
# @author: Nick N. Gibbons (2021/05/04)

require 'test/unit'
require 'open3'

class TestNonaka < Test::Unit::TestCase
  def test_0_prep
    cmd = "prep-gas air-5sp-2T.inp air-5sp-2T.gas"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "prep-chem air-5sp-2T.gas GuptaEtAl-air-2T.lua air-5sp-6r-2T.chem"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "prep-kinetics air-5sp-2T.gas air-5sp-6r-2T.chem air-VT-energy-exchange.lua air-VT.exch"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "e4shared --job=nonaka --prep"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_1_run
    cmd = "mpirun -np 4 e4mpi --job=nonaka --run --report-residuals"
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
    assert((steps - 4713).abs < 100, "Failed to take correct number of steps.")
  end

  def test_2_shock_stand_off
    cmd = 'e4shared --custom-script --script-file=shock-shape.lua'
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = 'python3 compute-error.py'
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    shock_location = 0.0
    lines = o.split("\n")
    rms = 0
    lines.each do |txt|
      if txt.match('Shock Standoff RMS:') then
        items = txt.split(' ')
        rms = items[3].to_f
      end
    end
    # Check that we are in the same place as before relative to the experimental data
    rms_ref = 0.011668210375449803
    assert((rms - rms_ref).abs < 1.0e-4,
           "Failed to get correct shock position.")
  end
end
