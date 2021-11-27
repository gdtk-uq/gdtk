#! /usr/bin/env ruby
# cyl50-usg-test.rb
# This is a viscous boundary layer case with an unstructured-grid.
# PJ, 2019-11-14
#
require 'test/unit'
require 'open3'

class TestCyl50_USG < Test::Unit::TestCase
  def test_0_prep
    cmd = "prep-gas ideal-air.inp ideal-air-gas-model.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "e4shared --prep --job=cyl50"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_1_run
    cmd = "e4shared --run --job=cyl50 --verbosity=1 --max-cpus=4"
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
    # Revised 2021-08-06, 2021-11-27 PJ.
    # For the transient code with explicit updates, expect 187980 steps.
    # It is very much lower with the implicit update scheme and slightly
    # sensitive to the model of the boundary condition.
    assert((steps - 2044).abs < 50, "Failed to take correct number of steps.")
  end

  def test_2_peak_T_in_boundary_layer
    ref = {'rho'=>3.47e-3, 'p'=>257.3, 'T'=>257.0, 'velx'=>300.0}
    cmd = 'e4shared --post --job=cyl50 --tindx-plot=last --add-vars="mach" --probe=0.917,0.00775,0.0'
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    lines = o.split("\n")
    columns = {'rho'=>0, 'velx'=>0, 'p'=>0, 'T'=>0}
    values = {'rho'=>0.0, 'velx'=>0.0, 'p'=>0.0, 'T'=>0.0}
    lines.each do |txt|
      if columns['rho'] > 0
        # Already found names, must be data line.
        items = txt.split(' ')
        values['rho'] = items[columns['rho']-1].to_f
        values['p'] = items[columns['p']-1].to_f
        values['T'] = items[columns['T']-1].to_f
        values['velx'] = items[columns['velx']-1].to_f
        break
      end
      if txt.match('pos.x')
        # Found variable names, extract column numbers.
        columns['rho'] = txt.match(/(\d+):rho\s/)[1].to_i
        columns['p'] = txt.match(/(\d+):p\s/)[1].to_i
        columns['T'] = txt.match(/(\d+):T\s/)[1].to_i
        columns['velx'] = txt.match(/(\d+):vel.x/)[1].to_i
      end
    end
    assert((values['rho'] - ref['rho']).abs/ref['rho'] < 2.0e-2, "Failed to see correct density")
    assert((values['p'] - ref['p']).abs/ref['p'] < 1.0e-2, "Failed to see correct pressure.")
    assert((values['T'] - ref['T']).abs/ref['T'] < 2.0e-2, "Failed to see correct temperature.")
    assert((values['velx'] - ref['velx']).abs/ref['velx'] < 1.0e-1, "Failed to see correct x velocity.")
  end

end
