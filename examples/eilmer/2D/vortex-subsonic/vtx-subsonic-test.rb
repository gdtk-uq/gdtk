#! /usr/bin/env ruby
# vtx-subsonic-test.rb
# Tests for the subsonic vortex example with Lachlan's alpha-split flux calculator.
# PJ, 2019-11-23
#
require 'test/unit'
require 'open3'

class TestVtxSubsonic < Test::Unit::TestCase
  def test_0_prep
    cmd = "prep-gas ideal-air.inp ideal-air-gas-model.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "e4shared --prep --job=vortex"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_1_run
    cmd = "mpirun -np 4 e4mpi --run --job=vortex --verbosity=1"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    steps = 0
    lines = o.split("\n")
    lines.each do |txt|
      if txt.match('final-t=') then
        items = txt.split(' ')
        steps = items[1].to_i
        break
      end
    end
    assert((steps - 10939).abs < 30, "Failed to take correct number of steps.")
  end

  def test_2_post_Lx_norms
    cmd = 'e4shared --post --job=vortex --tindx-plot=last'
    cmd = cmd + ' --ref-soln=vortex-flow-spec.lua --norms="rho"'
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    lines = o.split("\n")
    errL1 = 0.0; errL2 = 0.0
    lines.each do |txt|
      if txt.match('L2=') then
        # Found line with norms; extract value.
        items = txt.split(' ')
        errL1 = items[1].to_f
        errL2 = items[3].to_f
        break
      end
    end
    assert((errL1 - 2.328e-4).abs < 0.1e-4, "Failed to see expected L1 density error.")
    assert((errL2 - 6.274e-4).abs < 0.1e-4, "Failed to see expected L2 density error.")
  end
end
