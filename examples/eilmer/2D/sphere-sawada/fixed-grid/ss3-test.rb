#! /usr/bin/env ruby
# ss3-test.rb
# 2D axisymmetric flow with look-up table model for air.
# PJ, 2019-06-27
#
require 'test/unit'
require 'open3'

class TestSS3 < Test::Unit::TestCase
  def test_0_prep
    cmd = "cp #{ENV['DGD_REPO']}/src/gas/sample-data/cea-lut-air-version-test.lua ./cea-lut-air.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "e4shared --prep --job=ss3"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_1_run
    cmd = "e4shared --run --job=ss3 --verbosity=1 --max-cpus=4"
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
    assert((steps - 11137).abs < 100, "Failed to take correct number of steps.")
  end

  def test_2_shock_stand_off
    cmd = 'e4shared --post --job=ss3 --tindx-plot=last --slice-list=0,:,0,: '+
          '--output-file=ss3_stag_line.data'
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = 'awk -f locate_shock.awk ss3_stag_line.data'
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    shock_location = 0.0
    lines = o.split("\n")
    lines.each do |txt|
      if txt.match('shock-location=') then
        items = txt.split(' ')
        shock_location = items[1].to_f
      end
    end
    # Compute shock stand-off im mm and then check that the 
    # relative error is less than 15%, which corresponds to 
    # a distance of two cells.
    stand_off = shock_location* (-1000) - 31.8
    assert((stand_off - 2.66).abs/2.66 < 0.15, "Failed to get correct shock location.")
  end
end
