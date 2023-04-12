#! /usr/bin/env ruby
# MHD-ShockTube-Test.rb
# Test the MHD capability within Eilmer,
# specifically the MHD HLLE flux method
# from Vince and Daryl. Specific conditions
# are the third problem from Dai and Woodward,
# "A Simple Finite Difference Scheme for 
# Multidimensional Magnetohydrodynamical
# Equations", JCP 142. Figure 7 and Table 2.

require 'test/unit'
require 'open3'

class TestMHDSolver < Test::Unit::TestCase
  def test_0_prep
    cmd = "prep-gas ideal-He.inp ideal-He.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "e4shared --prep --job=MHDShockTube"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_1_run
    cmd = "e4shared --run --job=MHDShockTube --verbosity=1 --max-cpus=2"
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
    assert((steps - 880).abs < 10, "Failed to take correct number of steps.")
  end

  def test_2_region_conditions
    cmd = "e4shared --post --job=MHDShockTube --tindx-plot=last --probe=0.3,0.25,0.0;0.5,0.25,0.0;0.7,0.25,0.0"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    lines = o.split("\n")
    columns = {'p'=>0, 'vel.x'=>0, 'vel.y'=>0, 'vel.z'=>0, 'B.y'=>0, 'B.z'=>0, 'rho'=>0}
    values = {'p'=>[], 'vel.x'=>[], 'vel.y'=>[], 'vel.z'=>[], 'B.y'=>[], 'B.z'=>[], 'rho'=>[]}
    lines.each do |txt|
      if columns['p'] > 0
        # Already found names, must be data line.
        items = txt.split(' ')
        values['p'].append(items[columns['p']-1].to_f)
        values['vel.x'].append(items[columns['vel.x']-1].to_f)
        values['vel.y'].append(items[columns['vel.y']-1].to_f)
        values['vel.z'].append(items[columns['vel.z']-1].to_f)
        values['B.y'].append(items[columns['B.y']-1].to_f * Math.sqrt(Math::PI * 4))
        values['B.z'].append(items[columns['B.z']-1].to_f * Math.sqrt(Math::PI * 4))
        values['rho'].append(items[columns['rho']-1].to_f)
      end
      if txt.match('pos.x')
        # Found variable names, extract column numbers.
        columns['p'] = txt.match(/(\d+):p\s/)[1].to_i
        columns['vel.x'] = txt.match(/(\d+):vel.x\s/)[1].to_i
        columns['vel.y'] = txt.match(/(\d+):vel.y\s/)[1].to_i
        columns['vel.z'] = txt.match(/(\d+):vel.z\s/)[1].to_i
        columns['B.y'] = txt.match(/(\d+):B.y\s/)[1].to_i
        columns['B.z'] = txt.match(/(\d+):B.z\s/)[1].to_i
        columns['rho'] = txt.match(/(\d+):rho\s/)[1].to_i
      end
    end

    # The values from Dai and Woodward for each region
    ref_values = {'p'=>[1.84083e2, 1.84083e2, 1.84083e2], 'vel.x'=>[1.63777e-5, 0.0, 1.63777e-5], 'vel.y'=>[-2.36600e-1, 0.0, 2.36599e-1], 'vel.z'=>[-5.71203e-1, 0.0, -5.71203e-1], 'B.y'=>[4.00346, 5.66174, 4.00346], 'B.z'=>[-4.00346, -1.23741e-7, 4.00346], 'rho'=>[3.90913, 3.90913, 3.90913]}
    ref_tol = {'p'=>1e1, 'vel.x'=>1e-1, 'vel.y'=>1e-2, 'vel.z'=>1e-2, 'B.y'=>1e-1, 'B.z'=>1e-1, 'rho'=>1e-1}

    for i in 0..2
      ref_values.each do |key, value|
        assert((values[key][i] - ref_values[key][i]).abs < ref_tol[key], "Failed to see correct #{key} in region #{i}\n. Got #{values[key][i]}, expected #{ref_values[key][i]}")
      end
    end
  end
end


