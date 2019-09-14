#! /usr/bin/env ruby
# mms-euler-test.rb
# Run the method of manufactured solutions test case for an
# unstructured grid. These known good values come from an order
# of accuracy test.
#
# KAD, 2016-05-03; PJ, 2019-06-24
#
require 'test/unit'
require 'open3'

class TestMMS_Euler_USG < Test::Unit::TestCase
  def test_0_prep
    cmd = "python3 make_source_terms.py"
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
    # recalibrate 2018-10-06, 2018-11-28
    assert((steps - 2623).abs < 3, "Failed to take correct number of steps.")
  end

  def test_2_check_norms
    cmd = 'e4shared --job=mms --post --tindx-plot=20 --ref-soln=udf-bc.lua --norms="rho"'
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    lines = o.split("\n")
    names = ['L1', 'L2', 'Linf']
    # While experimenting with having only a single ghost cell,
    # we recalibrate these values so that the test does not fail. PJ 2016-07-29
    refs = {'L1'=>1.730145442812528172e-03,
            'L2'=>2.436181746051084129e-03,
            'Linf'=>6.587021426121220102e-03}
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
