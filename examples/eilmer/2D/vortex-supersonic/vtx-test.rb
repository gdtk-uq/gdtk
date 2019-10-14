#! /usr/bin/env ruby
# vtx-test.rb
# Tests for the supersonic vortex example.
# PJ, 2019-09-23
#
require 'test/unit'
require 'open3'

class TestVtx < Test::Unit::TestCase
  def test_0_prep
    cmd = "prep-gas ideal-air.inp ideal-air-gas-model.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "e4shared --prep --job=vtx"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_1_run
    cmd = "mpirun -np 4 e4mpi --run --job=vtx --verbosity=1"
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
    assert((steps - 2761).abs < 3, "Failed to take correct number of steps.")
  end

  def test_2_post_L2_norms
    cmd = 'e4shared --post --job=vtx --tindx-plot=last --ref-soln=udf-vortex-flow.lua --norms="p"'
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    lines = o.split("\n")
    pressure = 0.0
    lines.each do |txt|
      if txt.match('L2=') then
        # Found line with norms; extract value.
        items = txt.split(' ')
        pressure = items[3].to_f
        break
      end
    end
    assert((pressure - 800.0).abs < 100.0, "Failed to see expected pressure error.")
    #
    cmd = 'e4shared --post --job=vtx --tindx-plot=last --ref-soln=udf-vortex-flow.lua --norms="T"'
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    lines = o.split("\n")
    temperature = 0.0
    lines.each do |txt|
      if txt.match('L2=') then
        # Found line with norms; extract value.
        items = txt.split(' ')
        temperature = items[3].to_f
        break
      end
    end
    assert((temperature - 0.405).abs < 0.10, "Failed to see expected temperature error.")
  end
end
