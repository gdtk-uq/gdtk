#! /usr/bin/env ruby
# single-thread-test.rb
# Run the method of manufactured solutions to test the coupling
# of the gas and solid domains. This test executes using a single thread.
#
# RJG, 2016-06-14; PJ, 2019-06-24
#
require 'test/unit'
require 'open3'

class TestCHT_MMS_Spatial < Test::Unit::TestCase
  def test_0_prep
    cmd = "python3 make_lua_files.py"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "e4shared --prep --job=coupled-mms"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_1_run
    cmd = "e4shared --run --job=coupled-mms --max-cpus=1"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_2_check_norms_gas
    cmd = 'e4shared --job=coupled-mms --post --tindx-plot=last --ref-soln=ref-soln.lua --norms="T"'
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    lines = o.split("\n")
    names = ['L1', 'L2', 'Linf']
    refs = {'L1'=>1.535024800974120751e-01,
            'L2'=>2.105657110574111757e-01,
            'Linf'=>8.343791649878085082e-01}
    values = {'L1'=>0.0, 'L2'=>0.0, 'Linf'=>0.0}
    lines.each do |txt|
      if txt.match('Linf')
        # Must be our data line, in the batch of fluid blocks.
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

  def test_2_check_norms_solid
    cmd = 'e4shared --job=coupled-mms --post --tindx-plot=last --ref-soln=ref-solid-soln.lua --norms="T"'
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    lines = o.split("\n")
    names = ['L1', 'L2', 'Linf']
    refs = {'L1'=>5.601069822072819449e-01,
            'L2'=>7.483531389162372260e-01,
            'Linf'=>1.768327973065368042e+00}
    values = {'L1'=>0.0, 'L2'=>0.0, 'Linf'=>0.0}
    foundState = 0
    lines.each do |txt|
      if txt.match('Linf')
        foundState += 1
      end
      if txt.match('Linf') and foundState==2
        # Must be our data line, in the second batch of (solid) blocks.
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
