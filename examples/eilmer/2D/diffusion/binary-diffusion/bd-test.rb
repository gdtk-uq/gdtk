#! /usr/bin/env ruby
# bd-test.rb
# Functional test for binary diffusion of two gases.
# This exercises the laminar diffusion model.
# RJG, 2017-07-02; PJ, 2019-06-24
#
require 'test/unit'
require 'open3'

class TestCone20 < Test::Unit::TestCase
  def test_0_prep
    cmd = "prep-gas gas-model.inp thermally-perfect-N2-O2-mix.lua"
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
    names = ['rho', 'massf_N2', 'massf_O2']
    refs = {'rho'=>1.373172040772266067,
            'massf_N2'=>1.819097532464550027e-01,
            'massf_O2'=>8.180902467535450251e-01}
    tols = {'rho'=>5.0e-3, 'massf_N2'=>2.5e-2, 'massf_O2'=>2.5e-2}
    columns = {}
    values = {}
    lines.each do |txt|
      if columns['rho']
        # Already found names, must be data line.
        items = txt.split(' ')
        names.each do |n|
          values[n] = items[columns[n]-1].to_f
        end
        break
      end
      if txt.match('pos.x')
        # Found variable names, extract column numbers.
        columns['rho'] = txt.match(/(\d+):rho\s/)[1].to_i
        columns['massf_N2'] = txt.match(/(\d+):massf\[\d+\]-N2\s/)[1].to_i
        columns['massf_O2'] = txt.match(/(\d+):massf\[\d+\]-O2\s/)[1].to_i
      end
    end
    names.each do |n|
      assert((values[n] - refs[n]).abs/refs[n] < tols[n], "Failed to see correct #{n}.")
    end
  end

end
