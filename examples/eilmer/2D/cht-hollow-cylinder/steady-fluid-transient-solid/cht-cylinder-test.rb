#! /usr/bin/env ruby
# cht-cylinder-test.rb
#
# This test covers the steady fluid transient solid conjugate heat transfer solver,
# the test computes the heat soak into a hollow cylinder subjected to hypersonic flow for 5 seconds.
#
# author: Kyle A. Damm
# date:   2023-11-09

require 'test/unit'
require 'open3'

class TestCHTCyl < Test::Unit::TestCase
  def test_0_prep
    cmd = "prep-gas gm-ideal-air.inp ideal-air.gas"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "e4shared --prep --job=cylinder"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_1_run
    cmd = "mpirun -np 4 e4-cht-dist --job=cylinder --snapshot-start=0"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_2_surface_temperature
    cmd = 'python3 compute_error.py'
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    rms = 0.0
    lines = o.split("\n")
    lines.each do |txt|
      if txt.match('Surface Temperature RMS:') then
        items = txt.split(' ')
        rms = items[3].to_f
      end
    end
    # Check that we have computed the same temperature profile as before
    assert(rms.abs < 1.0e-6, "Failed to get correct surface temperature profile.")
  end
end
