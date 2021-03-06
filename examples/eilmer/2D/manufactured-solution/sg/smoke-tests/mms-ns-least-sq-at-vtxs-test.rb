#! /usr/bin/env ruby
# mms-ns-least-sq-at-vtxs-test.rb
# Run the method of manufactured solutions test case to exercise
# various options in the code. These cases are all run on a single
# grid. Our tests involve checking the density norms against
# known good values. These known good values come from an order
# of accuracy test.
# RJG, 2015-11-19; PJ, 2019-06-24
#
require 'test/unit'
require 'open3'

class TestMMS_NS_LeastSqAtVtx < Test::Unit::TestCase
  def test_0_prep
    cmd = "cp case-ns-least-sq-at-vtxs.txt case.txt"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "python3 make_lua_files.py"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "e4shared --prep --job=mms"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_1_run
    cmd = "e4shared --run --job=mms"
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
    assert((steps - 1507).abs < 3, "Failed to take correct number of steps.")
  end

  def test_2_check_norms
    cmd = 'e4shared --job=mms --post --tindx-plot=20 --ref-soln=ref-soln.lua --norms="rho"'
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    lines = o.split("\n")
    names = ['L1', 'L2', 'Linf']
    # For symmetric weighting of vertex gradients.
    # refs = {'L1'=>1.46279480e-03, 'L2'=>1.77472463e-03, 'Linf'=>4.13957472e-03}
    # For upwind weighting of vertex gradients. PJ 2021-09-13
    refs = {'L1'=>2.35516834e-03, 'L2'=>2.77413701e-03, 'Linf'=>5.71746417e-03}
    values = {'L1'=>0.0, 'L2'=>0.0, 'Linf'=>0.0}
    lines.each do |txt|
      if txt.match('Linf')
        # Must be our data line.
        items = txt.strip().split(' ')
        values['L1'] = items[1].to_f
        values['L2'] = items[3].to_f
        values['Linf'] = items[5].to_f
        break
      end
    end
    names.each do |n|
      assert((values[n]-refs[n]).abs/refs[n] < 1.0e-3, "Failed to see correct #{n}.")
    end
  end
end
