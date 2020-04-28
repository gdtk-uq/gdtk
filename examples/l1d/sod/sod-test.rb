#! /usr/bin/env ruby
# sod-test.rb
# Tests for the Sod shock-tube example for L1d4.
# PJ, 2020-04-28
#
require 'test/unit'
require 'open3'

class TestSod < Test::Unit::TestCase
  def test_0_prep
    cmd = "prep-gas ideal-air.inp ideal-air-gas-model.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "l1d4-prep --job=sod"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_1_run
    cmd = "l1d4 --run-simulation --job=sod"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    sim_time = 0.0
    lines = o.split("\n")
    lines.each do |txt|
      if txt.match('Step=200') then
        items = txt.split(' ')
        sim_time_items = items[1].split('=')
        sim_time = sim_time_items[1].to_f
      end
    end
    assert((sim_time - 5.421e-4).abs < 1.0e-5, "Inorrect sim_time at step 200.")
  end

  def test_2_post
    cmd = "l1d4 --time-slice --job=sod --tindx=60"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    f = File.new("slug-0001-tindx-0060-cells.data", "r")
    lines = f.readlines
    f.close
    columns = {'vel'=>0, 'p'=>0, 'T'=>0}
    values = {'vel'=>0.0, 'p'=>0.0, 'T'=>0.0}
    items = lines[0].split
    # Header line starts with a sharp character.
    columns['vel'] = items.index('vel') - 1
    columns['p'] = items.index('p') - 1
    columns['T'] = items.index('T') - 1
    lines.each do |txt|
      items = txt.split(' ')
      x = items[0].to_f
      if x >= 0.80 then
        values['vel'] = items[columns['vel']].to_f
        values['p'] = items[columns['p']].to_f
        values['T'] = items[columns['T']].to_f
        break
      end
    end
    assert((values['vel'] - 293.3).abs < 1.0, "Incorrect post-shock velocity.")
    assert((values['p'] - 30.31e3).abs < 100.0, "Incorrect post-shock pressure.")
    assert((values['T'] - 398.0).abs < 1.0, "Incorrect post-shock temperature.")
  end

end
