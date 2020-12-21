#! /usr/bin/env ruby
# mms-SA-test.rb
# Run the method of manufactured solutions test case to exercise
# the Navier-Stokes equations using the unstructured (steady-state) solver.
# Our tests involve checking the density norms against
# known good values. These known good values come from an order
# of accuracy test.
# DMB, 2021-02-01;
#
require 'test/unit'
require 'open3'

class TestMMS_SA < Test::Unit::TestCase
  def test_0_prep
    cmd = "python3 make_lua_files.py"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "e4shared --prep --job=mms"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_1_run
    cmd = "e4-nk-shared --job=mms"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    # print o
    steps = 0
    lines = o.split("\n")
    lines.each do |txt|
      if txt.match('STEP=') then
        items = txt.split(' ')
        steps = items[1].to_i
      end
    end
    assert((steps - 9).abs < 3, "Failed to take correct number of steps.")
  end

  def test_2_check_norms
    cmd = 'e4shared --job=mms --post --tindx-plot=last --ref-soln=ref-soln.lua --norms="rho" > norms.txt'
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    # print o
    lines = o.split("\n")
    names = ['L1', 'L2', 'Linf']
    refs = {'L1'=>3.973558775758266950e-03, 'L2'=>4.181660063469227508e-03, 'Linf'=>9.182445521693782808e-03}
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
