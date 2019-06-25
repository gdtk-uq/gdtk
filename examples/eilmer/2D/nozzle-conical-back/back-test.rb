#! /usr/bin/env ruby
# back-test.rb
# This exercises the SubsonicInBC that forms the upstream boundary for
# the simulation of the Back et al. nozzle.
#
# PJ, 2019-06-25 Ruby port
# Tcl versions:03-Mar-2014
#     2015-10-22 ported to Eilmer4
#     2015-10-27 relax tolerance to accommodate the less steady mass-flux solution.
#     2019-05-07 relax further to accommodate the different reconstruction
#     associated with the WallBC_WithSlip1 boundary condition.
#     We are sampling the exit flow in the most sensitive region, near the axis.
#
require 'test/unit'
require 'open3'

class TestBack < Test::Unit::TestCase
  def test_0_prep
    cmd = "prep-gas ideal-air.inp ideal-air-gas-model.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "e4shared --prep --job=back"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_1_run
    cmd = "e4shared --run --job=back --verbosity=1 --max-cpus=2"
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
    assert((steps - 5400).abs < 100, "Failed to take correct number of steps.")
  end

  def test_2_post_subsonic_region
    cmd = 'e4shared --post --job=back --tindx-plot=last --add-vars="mach" --probe=-0.075,0.0,0.0'
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    lines = o.split("\n")
    columns = {'rho'=>0, 'M'=>0, 'p'=>0, 'T'=>0}
    values = {'rho'=>0.0, 'M'=>0.0, 'p'=>0.0, 'T'=>0.0}
    lines.each do |txt|
      if columns['rho'] > 0
        # Already found names, must be data line.
        items = txt.split(' ')
        values['rho'] = items[columns['rho']-1].to_f
        values['p'] = items[columns['p']-1].to_f
        values['T'] = items[columns['T']-1].to_f
        values['M'] = items[columns['M']-1].to_f
        break
      end
      if txt.match('pos.x')
        # Found variable names, extract column numbers.
        columns['rho'] = txt.match(/(\d+):rho\s/)[1].to_i
        columns['p'] = txt.match(/(\d+):p\s/)[1].to_i
        columns['T'] = txt.match(/(\d+):T\s/)[1].to_i
        columns['M'] = txt.match(/(\d+):M_local/)[1].to_i
      end
    end
    # The following values were extracted from the solution
    # as it was computed on PJ's computer 03-Mar-2014.
    ref = {'rho'=>5.743503e+00, 'p'=>4.927381e+05, 'T'=>2.988554e+02, 'M'=>1.384033e-01}
    assert((values['rho'] - ref['rho']).abs/ref['rho'] < 0.01, "Failed to see correct density.")
    assert((values['p'] - ref['p']).abs/ref['p'] < 0.01, "Failed to see correct pressure.")
    assert((values['T'] - ref['T']).abs/ref['T'] < 0.01, "Failed to see correct temperature.")
    assert((values['M'] - ref['M']).abs/ref['M'] < 0.01, "Failed to see correct Mach number.")
  end

  def test_3_post_supersonic_region
    cmd = 'e4shared --post --job=back --tindx-plot=last --add-vars="mach" --probe=0.075,0.0,0.0'
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    lines = o.split("\n")
    columns = {'rho'=>0, 'M'=>0, 'p'=>0, 'T'=>0}
    values = {'rho'=>0.0, 'M'=>0.0, 'p'=>0.0, 'T'=>0.0}
    lines.each do |txt|
      if columns['rho'] > 0
        # Already found names, must be data line.
        items = txt.split(' ')
        values['rho'] = items[columns['rho']-1].to_f
        values['p'] = items[columns['p']-1].to_f
        values['T'] = items[columns['T']-1].to_f
        values['M'] = items[columns['M']-1].to_f
        break
      end
      if txt.match('pos.x')
        # Found variable names, extract column numbers.
        columns['rho'] = txt.match(/(\d+):rho\s/)[1].to_i
        columns['p'] = txt.match(/(\d+):p\s/)[1].to_i
        columns['T'] = txt.match(/(\d+):T\s/)[1].to_i
        columns['M'] = txt.match(/(\d+):M_local/)[1].to_i
      end
    end
    # The following values were extracted from the solution
    # as it was computed on PJ's computer 03-Mar-2014.
    ref = {'rho'=>3.904117e-01, 'p'=>1.142634e+04, 'T'=>1.019545e+02, 'M'=>3.116296e+00}
    assert((values['rho'] - ref['rho']).abs/ref['rho'] < 0.05, "Failed to see correct density.")
    assert((values['p'] - ref['p']).abs/ref['p'] < 0.05, "Failed to see correct pressure.")
    assert((values['T'] - ref['T']).abs/ref['T'] < 0.02, "Failed to see correct temperature.")
    assert((values['M'] - ref['M']).abs/ref['M'] < 0.02, "Failed to see correct Mach number.")
  end

end
