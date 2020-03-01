#! /usr/bin/env ruby
# bump-test.rb
# Functional test for inviscid flow in a channel.
# This exercises quite a few of the basic functions of the block-marching code.
# PJ, 2019-06-26
#
require 'test/unit'
require 'open3'

class TestBump < Test::Unit::TestCase
  def test_0_prep
    cmd = "prep-gas ideal-air.inp ideal-air-gas-model.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "e4shared --prep --job=bump --verbosity=1 --only-blocks=0..<96"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "e4shared --prep --job=bump --verbosity=1 --only-blocks=96..<192"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_1_run
    cmd = "e4shared --run --job=bump --verbosity=1 --max-cpus=4"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_2_post_exit_flow
    cmd = 'e4shared --post --job=bump --tindx-plot=last --add-vars="mach" --probe=2.921,0.364,0.0'
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
    assert((values['a'] - 349.3).abs < 1.0, "Failed to see correct sound speed.")
    assert((values['p'] - 121000).abs < 1000.0, "Failed to see correct pressure.")
    assert((values['T'] - 303.5).abs < 1.0, "Failed to see correct temperature.")
    assert((values['M'] - 1.52).abs < 0.05, "Failed to see correct Mach number.")
  end

end
