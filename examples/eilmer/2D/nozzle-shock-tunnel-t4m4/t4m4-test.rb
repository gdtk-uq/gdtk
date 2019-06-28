#! /usr/bin/env ruby
# t4m4-test.rb
# Functional test for axisymmetric T4 Mach 4 nozzle flow, block-marching version.
# This is Wilson's Mach 7 flight condition at q of 50 kPa.
# PJ, 2018-04-14, 2019-06-27
#
require 'test/unit'
require 'open3'

class TestT4M4 < Test::Unit::TestCase
  def test_0_prep
    cmd = "cp #{ENV['DGD_REPO']}/src/gas/sample-data/cea-lut-air-version-test.lua ./cea-lut-air.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "cp #{ENV['DGD_REPO']}/src/gas/sample-data/cea-air13species-gas-model.lua ./cea-air13species-gas-model.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "e4shared-debug --prep --job=t4m4"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_1_run
    cmd = "e4shared --run --job=t4m4 --verbosity=1 --max-cpus=4"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_2_post
    cmd = 'e4shared --post --job=t4m4 --tindx-plot=last --add-vars="mach,pitot,total-p,total-h" '+
          '--output-file=nozzle-exit.data --extract-line="0.5114,0,0,0.5114,0.044,0,20"'
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def is_within_tolerance(v, v_ref, tol)
    (v - v_ref).abs < (tol*v_ref)
  end

  def mean(my_list)
    my_list.inject(:+) / my_list.length
  end
  
  def test_3_exit_condition
    # The nozzle exit conditions 444K, 1997m/s, 0.125kg/m**3.
    f = File.new('nozzle-exit.data', 'r')
    lines = f.readlines
    f.close
    myT_list = []; vel_list = []; density_list = []
    lines.each do |txt|
      if !txt.match(/^#/) then
        items = txt.split(' ')
        ypos = items[1].to_f
        if ypos > 0.044 then break end
        myT_list << items[19].to_f
        vel_list << items[5].to_f
        density_list << items[4].to_f
      end
    end
    assert(is_within_tolerance(mean(myT_list), 443.7, 0.01), "Failed to get correct temperature.")
    assert(is_within_tolerance(mean(vel_list), 1996.6, 0.01), "Failed to get correct velocity.")
    assert(is_within_tolerance(mean(density_list), 0.1252, 0.01), "Failed to get correct density.")
  end
end
