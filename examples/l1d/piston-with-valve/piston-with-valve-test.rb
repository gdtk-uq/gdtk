#! /usr/bin/env ruby
# piston-with-valve-test.rb
# Tests for the piston-in-tube, with valve, example for L1d4.
# PJ, 2020-04-28
#
require 'test/unit'
require 'open3'

class TestPiston < Test::Unit::TestCase
  def test_0_prep
    cmd = "prep-gas ideal-air.inp ideal-air-gas-model.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "l1d4-prep --job=piston-with-valve"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_1_run
    cmd = "l1d4 --run-simulation --job=piston-with-valve"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    sim_time = 0.0
    lines = o.split("\n")
    lines.each do |txt|
      if txt.match('Step=700') then
        items = txt.split(' ')
        sim_time_items = items[1].split('=')
        sim_time = sim_time_items[1].to_f
      end
    end
    assert((sim_time - 0.03931).abs < 0.001, "Incorrect sim_time at step 700.")
  end

  def test_2_post
    cmd = "l1d4 --piston-history --job=piston-with-valve --pindx=0"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    f = File.new("piston-0000-history.data", "r")
    txt = f.readlines[-1]
    f.close
    items = txt.split(' ')
    t = items[1].to_f
    x = items[2].to_f
    v = items[3].to_f
    assert((t - 0.040).abs < 0.0001, "Failed to reach correct final time.")
    assert((x - 6.559).abs < 0.1, "Failed to reach correct position.")
    assert((v - 276.7).abs < 1.0, "Failed to reach correct velocity.")
  end

  def test_3_energies
    f = File.new("piston/energies.data", "r")
    txt = f.readlines
    f.close
    items = txt[1].split(' ')
    e0 = items[-1].to_f
    items = txt[-1].split(' ')
    e1 = items[-1].to_f
    assert((e1 - e0).abs/e0 < 0.0001, "Failed to conserve energy.")
  end

end
