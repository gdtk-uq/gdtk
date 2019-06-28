#! /usr/bin/env ruby
# bittker-test.rb
# Functional test for reacting inviscid flow in a pipe.
# This exercises the finite-rate chemistry coupled to the gas dynamics.
# PJ, 21-Jun-2011, adapted to Eilmer4 18-Mar-2016, 2019-06-24
#
require 'test/unit'
require 'open3'

class TestBittker < Test::Unit::TestCase
  def test_0_prep
    cmd = "prep-gas combusting-species.inp h2-o2-n2-9sp.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "prep-chem h2-o2-n2-9sp.lua Bittker-Scullin.lua h2-o2-n2-9sp-18r.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "e4shared --prep --job=bittker"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_1_run
    cmd = "e4shared --run --job=bittker --verbosity=1 --max-cpus=2"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_2_combustion_products
    # Probe the solution near the end of the duct
    # and pull a few values from the line of data
    # that is produced by the postprocessor.
    cmd = 'e4shared --post --job=bittker --tindx-plot=last --probe=0.098,0.05,0.0'
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    names = ['rho', 'vel.x', 'p', 'T', 'massf_O', 'massf_H2', 'massf_H2O']
    # The following values were extracted from the solution
    # as it was computed on PJ's computer 20-Mar-2016 for cfl_value = 0.5.
    refs = {'rho'=>1.588144e-01, 'vel.x'=>4.479598e+03, 'p'=>1.498976e+05,
            'T'=>2485.2, 'massf_O'=>1.477749e-02, 'massf_H2'=>6.661693e-03,
            'massf_H2O'=>1.591853e-01}
    tols = {'rho'=>5.0e-3, 'vel.x'=>5.0e-3, 'p'=>5.0e-3, 'T'=>1.0e-2,
            'massf_O'=>5.0e-2, 'massf_H2'=>2.5e-2, 'massf_H2O'=>1.0e-2}
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
        columns['vel.x'] = txt.match(/(\d+):vel\.x\s/)[1].to_i
        columns['p'] = txt.match(/(\d+):p\s/)[1].to_i
        # Note that we expect T at the end of the line.
        columns['T'] = txt.match(/(\d+):T/)[1].to_i
        columns['massf_O'] = txt.match(/(\d+):massf\[\d\]-O\s/)[1].to_i
        columns['massf_H2'] = txt.match(/(\d+):massf\[\d\]-H2\s/)[1].to_i
        columns['massf_H2O'] = txt.match(/(\d+):massf\[\d\]-H2O\s/)[1].to_i
      end
    end
    names.each do |n|
      assert((values[n] - refs[n]).abs/refs[n] < tols[n], "Failed to see correct #{n}.")
    end
  end
end
