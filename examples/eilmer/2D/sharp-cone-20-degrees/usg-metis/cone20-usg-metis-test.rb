#! /usr/bin/env ruby
# cone20-usg-metis-test.rb
# Tests for the conical-flow sharp-cone, unstructured-grid, metis partitioning.
# Replaces Tcl test script and is a bit more tolerant
# of changes to column numbers for the flow data.
# PJ, 2019-06-25
#
require 'test/unit'
require 'open3'

class TestCone20 < Test::Unit::TestCase
  def test_0_prep
    # Clean up old generated files.
    files = Dir.glob("block_*")
    files.each do |f| File.delete(f) end
    files = Dir.glob("mapped_cells")
    files.each do |f| File.delete(f) end
    # Now ready for clean start.
    cmd = "ugrid_partition cone20.su2 mapped_cells 4 2"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "prep-gas ideal-air.inp ideal-air-gas-model.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "e4shared --prep --job=cone20"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_1_run
    cmd = "e4shared --run --job=cone20 --verbosity=1 --max-cpus=4"
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
    assert((steps - 2008).abs < 3, "Failed to take correct number of steps.")
  end

  def test_2_post_free_stream_flow
    cmd = 'e4shared --post --job=cone20 --tindx-plot=last --add-vars="mach" --probe=0.4,0.5,0.0'
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    lines = o.split("\n")
    columns = {'a'=>0, 'M'=>0, 'p'=>0, 'T'=>0}
    values = {'a'=>0.0, 'M'=>0.0, 'p'=>0.0, 'T'=>0.0}
    lines.each do |txt|
      if columns['a'] > 0
        # Already found names, must be data line.
        items = txt.split(' ')
        values['a'] = items[columns['a']-1].to_f
        values['p'] = items[columns['p']-1].to_f
        values['T'] = items[columns['T']-1].to_f
        values['M'] = items[columns['M']-1].to_f
        break
      end
      if txt.match('pos.x')
        # Found variable names, extract column numbers.
        columns['a'] = txt.match(/(\d+):a\s/)[1].to_i
        columns['p'] = txt.match(/(\d+):p\s/)[1].to_i
        columns['T'] = txt.match(/(\d+):T\s/)[1].to_i
        columns['M'] = txt.match(/(\d+):M_local/)[1].to_i
      end
    end
    assert((values['a'] - 666.0).abs < 1.0, "Failed to see correct sound speed.")
    assert((values['p'] - 95.84e3).abs < 500.0, "Failed to see correct pressure.")
    assert((values['T'] - 1103.0).abs < 1.0, "Failed to see correct temperature.")
    assert((values['M'] - 1.50).abs < 0.02, "Failed to see correct Mach number.")
  end

  def test_3_post_pressure_on_cone_surface
    cmd = 'e4shared --post --job=cone20 --tindx-plot=last --add-vars="mach" --probe=0.7,0.182,0.0'
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    lines = o.split("\n")
    column_p = 0
    value_p = 0.0
    lines.each do |txt|
      if column_p > 0
        # Already found names, must be data line.
        items = txt.split(' ')
        value_p = items[column_p-1].to_f
        break
      end
      if txt.match('pos.x')
        # Found variable names, extract column number.
        column_p = txt.match(/(\d+):p\s/)[1].to_i
      end
    end
    q_inf = 0.5 * 0.3029 * 1.0e6 # dynamic pressure for free stream
    p_cone_surface = 95.84e3 + 0.387 * q_inf
    assert((value_p - p_cone_surface).abs < 3.0e3, "Failed to see correct surface pressure.")
  end

  def test_4_post_shock_angle
    cmd = 'e4shared --custom-script --script-file=estimate_shock_angle.lua'
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    shock_angle = 0.0
    average_deviation = 1000.0 # something big
    lines = o.split("\n")
    lines.each do |txt|
      if txt.match('shock_angle_deg') then
        items = txt.split(' ')
        shock_angle = items[1].to_f
      end
      if txt.match('average_deviation_metres') then
        items = txt.split(' ')
        average_deviation = items[1].to_f
      end
    end
    assert((shock_angle - 49.547).abs < 1.0, "Failed to get correct shock angle.")
    assert((average_deviation).abs < 0.004, "Failed to find a straight shock.")
  end
end
