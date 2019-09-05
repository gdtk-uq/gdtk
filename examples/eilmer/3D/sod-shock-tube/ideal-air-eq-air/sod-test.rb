#! /usr/bin/env ruby
# sod-test.rb
# This exercises the LuaFnVolume class for user-defined geometry
# and the LUT-plus-ideal-gas model in a structured-grid example.
# PJ, 2019-07-11
#
require 'test/unit'
require 'open3'

class TestSod < Test::Unit::TestCase
  def test_0_prep
    cmd = "cp #{ENV['DGD_REPO']}/src/gas/sample-data/cea-lut-air-version-test.lua ."
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "prep-gas ideal-air.inp ideal-air-gas-model.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "e4shared --prep --job=sod"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_1_run
    cmd = "e4shared --run --job=sod --verbosity=1 --max-cpus=1"
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
    assert((steps - 115).abs < 3, "Failed to take correct number of steps.")
  end

  def test_2_post_shock_region
    ref = {'rho'=>0.2647, 'p'=>30.2e3, 'T'=>398.0, 'velx'=>293.0}
    cmd = 'e4shared --post --job=sod --tindx-plot=last --add-vars="mach" --probe=0.78,0.025,0.025'
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
    assert((values['rho'] - ref['rho']).abs/ref['rho'] < 1.0e-2, "Failed to see correct density.")
    assert((values['p'] - ref['p']).abs/ref['p'] < 1.0e-2, "Failed to see correct pressure.")
    assert((values['T'] - ref['T']).abs/ref['T'] < 1.0e-2, "Failed to see correct temperature.")
    assert((values['velx'] - ref['velx']).abs/ref['velx'] < 1.0e-2, "Failed to see correct x velocity.")
  end

  def test_3_expanded_driver_region
    ref = {'rho'=>0.4271, 'p'=>30.2e3, 'T'=>247.0, 'velx'=>293.0}
    cmd = 'e4shared --post --job=sod --tindx-plot=last --add-vars="mach" --probe=0.6,0.025,0.025'
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
    assert((values['rho'] - ref['rho']).abs/ref['rho'] < 1.0e-2, "Failed to see correct density.")
    assert((values['p'] - ref['p']).abs/ref['p'] < 1.0e-2, "Failed to see correct pressure.")
    assert((values['T'] - ref['T']).abs/ref['T'] < 1.0e-2, "Failed to see correct temperature.")
    assert((values['velx'] - ref['velx']).abs/ref['velx'] < 1.0e-2, "Failed to see correct x velocity.")
  end

end
