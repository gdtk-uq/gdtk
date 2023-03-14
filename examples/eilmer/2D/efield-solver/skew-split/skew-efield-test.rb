#! /usr/bin/env ruby
# elec-field-test.rb
# Tests for the electric field solver, in multi-block mode, using an known answer
# NNG, 2022-02-16
#
require 'test/unit'
require 'open3'

class TestElec < Test::Unit::TestCase
  def test_0_prep
    cmd = "prep-gas ideal-air.inp ideal-air.lua"
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
    lines = o.split("\n")
    current_in = 0
    lines.each do |txt|
      if txt.match('    Current in:') then
        items = txt.split(' ')
        current_in = items[2].to_f
      end
    end
    assert((current_in - 1.011589).abs < 1.0e-06, "Failed to compute correct in-current.")

    lines = o.split("\n")
    current_out = 0
    lines.each do |txt|
      if txt.match('    Current out:') then
        items = txt.split(' ')
        current_out = items[2].to_f
      end
    end
    assert((current_out - 1.011589).abs < 1.0e-06, "Failed to compute correct out-current.")

    lines = o.split("\n")
    solves = 0
    lines.each do |txt|
      if txt.match('    Solve Complete:') then
        items = txt.split(' ')
        solves = items[3].split('=')[1].split('/')[0].to_i
      end
    end
    assert(solves==180, "Failed to compute efield in the corrent number of steps.")
  end
end
