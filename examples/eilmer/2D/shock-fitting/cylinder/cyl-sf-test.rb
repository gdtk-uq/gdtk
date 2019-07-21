#! /usr/bin/env ruby
# cyl-sf-test.rb
# Tests the shock-fitting boundary condition and grid motion code.
# PJ, 2019-06-26
#
require 'test/unit'
require 'open3'

class TestCyl_SF < Test::Unit::TestCase
  def test_0_prep
    cmd = "prep-gas ideal-air.inp ideal-air-gas-model.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "e4shared --prep --job=cyl-sf"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_1_run
    cmd = "e4shared --run --job=cyl-sf --verbosity=1 --max-cpus=4"
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
    assert((steps - 5185).abs < 100, "Failed to take correct number of steps.")
  end

  def test_2_post_shock_condition
    # With this shock-fitting simulation, we take the first point 
    # in the grid as (being close enough to) the position of the shock.
    # For a free stream Mach number of 7, Billig's correlation gives
    # delta/R = 0.4246.
    # Normal shock jump tables give p2/p1=57, T2/T1=10.47, M2=0.3974
    ref = {'x'=>-1.425, 'p'=>5.7e6, 'T'=>3141.0, 'M'=>0.3974}
    cmd = 'e4shared --post --job=cyl-sf --tindx-plot=last --add-vars="mach" --probe=-1.425,0.0,0.0'
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    lines = o.split("\n")
    columns = {'x'=>0, 'M'=>0, 'p'=>0, 'T'=>0}
    values = {'x'=>0.0, 'M'=>0.0, 'p'=>0.0, 'T'=>0.0}
    lines.each do |txt|
      if columns['x'] > 0
        # Already found names, must be data line.
        items = txt.split(' ')
        values['x'] = items[columns['x']-1].to_f
        values['p'] = items[columns['p']-1].to_f
        values['T'] = items[columns['T']-1].to_f
        values['M'] = items[columns['M']-1].to_f
        break
      end
      if txt.match('pos.x')
        # Found variable names, extract column numbers.
        columns['x'] = txt.match(/(\d+):pos.x\s/)[1].to_i
        columns['p'] = txt.match(/(\d+):p\s/)[1].to_i
        columns['T'] = txt.match(/(\d+):T\s/)[1].to_i
        columns['M'] = txt.match(/(\d+):M_local/)[1].to_i
      end
    end
    assert((values['x'] - ref['x']).abs < 3.0e-2, "Failed to see correct shock position.")
    assert((values['p'] - ref['p']).abs/ref['p'] < 3.0e-2, "Failed to see correct pressure.")
    assert((values['T'] - ref['T']).abs/ref['T'] < 1.0e-2, "Failed to see correct temperature.")
    assert((values['M'] - ref['M']).abs/ref['M'] < 0.025, "Failed to see correct Mach number.")
  end

end
