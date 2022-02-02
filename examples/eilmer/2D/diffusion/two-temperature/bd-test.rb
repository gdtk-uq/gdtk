#! /usr/bin/env ruby
# bd-test.rb
# Functional test for binary diffusion of two gases.
# This exercises the two temperature diffusion process
# NNG, 2022-02-02; RJG, 2017-07-02; PJ, 2019-06-24
#
require 'test/unit'
require 'open3'

class TestCone20 < Test::Unit::TestCase
  def test_0_prep
    cmd = "prep-gas n2-3sp-gm.inp n2-3sp-gm.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "prep-chem n2-3sp-gm.lua n2-blank-rr.inp n2-blank-rr.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "prep-kinetics n2-3sp-gm.lua n2-blank-rr.lua n2-blank-ee.inp n2-blank-ee.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "e4shared --prep --job=bd"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_1_run
    cmd = "e4shared --run --job=bd"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_2_composition_at_x
    cmd = 'e4shared --post --job=bd --tindx-plot=last --probe=0.5e-5,0.0,0.0'
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    lines = o.split("\n")
    names = ['T_modes[0]', 'massf_N2', 'massf_N2+', 'massf_e-']
    refs = {'T_modes[0]'=>2000.0,
            'massf_N2'=>0.9973563683756822,
            'massf_N2+'=>0.00264357985473305,
            'massf_e-'=>5.176958476351787e-08}
    tols = {'T_modes[0]'=>1.0e-3, 'massf_N2'=>2.5e-2, 'massf_N2+'=>2.5e-2, 'massf_e-'=>2.5e-2}
    columns = {}
    values = {}
    lines.each do |txt|
      if columns['T_modes[0]']
        # Already found names, must be data line.
        items = txt.split(' ')
        names.each do |n|
          values[n] = items[columns[n]-1].to_f
        end
        break
      end
      if txt.match('pos.x')
        # Found variable names, extract column numbers.
        columns['T_modes[0]'] = txt.match(/(\d+):T_modes\[0\]\s/)[1].to_i
        columns['massf_N2'] = txt.match(/(\d+):massf\[\d+\]-N2\s/)[1].to_i
        columns['massf_N2+'] = txt.match(/(\d+):massf\[\d+\]-N2\+\s/)[1].to_i
        columns['massf_e-'] = txt.match(/(\d+):massf\[\d+\]-e-+\s/)[1].to_i
      end
    end
    names.each do |n|
      assert((values[n] - refs[n]).abs/refs[n] < tols[n], "Failed to see correct #{n}.")
    end
  end

end
