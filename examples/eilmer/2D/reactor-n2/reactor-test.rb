#! /usr/bin/env ruby
# reactor-test.rb
# This exercises the finite-rate chemistry coupled to the gas-dynamics.
# PJ, 2019-06-27
#
require 'test/unit'
require 'open3'

class TestReactor < Test::Unit::TestCase
  def test_0_prep
    cmd = 'prep-gas nitrogen-2sp.inp nitrogen-2sp.lua'
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = 'prep-chem nitrogen-2sp.lua nitrogen-2sp-2r.lua chem.lua'
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_1_run
    cmd = 'e4shared --prep --job=reactor'
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = 'e4shared --run --job=reactor --verbosity=1 --max-cpus=1'
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
    assert((steps - 20001).abs < 3, "Failed to take correct number of steps.")
  end

  def test_2_reaction_products
    f = File.new("hist/reactor-blk-0-cell-0.dat.0", "r")
    lines = f.readlines
    f.close
    names = ['p', 'T', 'massf_N2', 'massf_N']
    refs = {'p'=>1.494e5, 'T'=>6381.7, 'massf_N2'=>0.8760, 'massf_N'=>0.1240}
    columns = {}
    values = {}
    # Extract column numbers from first line.
    txt = lines[0]
    columns['p'] = txt.match(/(\d+):p\s/)[1].to_i
    columns['T'] = txt.match(/(\d+):T\s/)[1].to_i
    columns['massf_N2'] = txt.match(/(\d+):massf\[0\]-N2\s/)[1].to_i
    columns['massf_N'] = txt.match(/(\d+):massf\[1\]-N\s/)[1].to_i
    # Extract the values from the last line and check.
    txt = lines[-1]
    items = txt.split(' ')
    names.each do |n|
      values[n] = items[columns[n]-1].to_f
      assert((values[n]-refs[n]).abs/refs[n] < 1.0e-3, "Failed to see correct value for #{n}")
    end
  end
end
