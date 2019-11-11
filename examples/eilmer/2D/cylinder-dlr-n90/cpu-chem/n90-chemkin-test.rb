#! /usr/bin/env ruby
# n90-test.rb
# This exercises the nonequilibrium thermochemistry with the flow of
# dissociating nitrogen over a cylinder -- Chemkin thermochemistry.
# PJ, 2015-11-10, 2019-11-12
#
require 'test/unit'
require 'open3'

class TestN90Chemkin < Test::Unit::TestCase
  def test_0_prep
    cmd = "prep-gas nitrogen-2sp.inp nitrogen-2sp.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "chemkin2eilmer n2-dissociation.inp nitrogen-2sp.lua n2-dissociation.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "prep-chem nitrogen-2sp.lua n2-dissociation.lua e4-chem.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "e4shared --prep --job=n90"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_1_run
    cmd = "e4shared --run --job=n90 --verbosity=1 --max-cpus=4"
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
    assert((steps - 4428).abs < 30, "Failed to take correct number of steps.")
  end

  def test_2_post_shock_state
    cmd = 'e4shared --post --job=n90 --tindx-plot=5 --probe=-0.010,0.0,0.0'
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    names = ['rho', 'vel.x', 'p', 'massf_N', 'T']
    # The following values were extracted from the solution, with a nicer grid,
    # as it was computed on PJ's computer 2019-11-12.
    refs = {'rho'=>1.870e-02, 'vel.x'=>5.924e+02, 'p'=>5.353e+04,
            'massf_N'=>2.098e-02, 'T'=>9.448e+03}
    tols = {'rho'=>0.01, 'vel.x'=>0.1, 'p'=>0.01, 'massf_N'=>0.15, 'T'=>0.01}
    columns = {}
    values = {}
    lines = o.split("\n")
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
        columns['vel.x'] = txt.match(/(\d+):vel.x\s/)[1].to_i
        columns['p'] = txt.match(/(\d+):p\s/)[1].to_i
        columns['T'] = txt.match(/(\d+):T/)[1].to_i # T is at end of line
        columns['massf_N'] = txt.match(/(\d+):massf\[\d+\]-N\s/)[1].to_i
      end
    end
    names.each do |n|
      assert((values[n] - refs[n]).abs/refs[n] < tols[n], "Failed to see correct #{n}.")
    end
  end
end
