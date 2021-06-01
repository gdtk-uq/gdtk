#! /usr/bin/env ruby
# cone20-test.rb
# Tests for the conical-flow sharp-cone example with the new flow format.
# Replaces Tcl test script and is a bit more tolerant
# of changes to column numbers for the flow data.
# Very minor changes from PJ's original (PJ, 2019-06-24).
# At this stage, the custom processing does not work with the new flow format,
# so tests 3 and 4 are dropped.
# Lachlan Whyborn, 2021-06-01
#
require 'test/unit'
require 'open3'

class TestCone20_nff < Test::Unit::TestCase
  def test_0_prep
    cmd = "prep-gas ideal-air.inp ideal-air-gas-model.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "e4shared --prep --job=cone20-nff"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_1_run
    cmd = "e4shared --run --job=cone20-nff --verbosity=1 --max-cpus=2"
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
    assert((steps - 833).abs < 3, "Failed to take correct number of steps.")
  end

  def test_2_post_free_stream_flow
    cmd = 'e4shared --post --job=cone20-nff --tindx-plot=last --add-vars="mach" --probe=0.4,0.5,0.0'
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
end
