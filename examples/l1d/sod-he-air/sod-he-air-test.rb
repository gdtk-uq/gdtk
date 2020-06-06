#! /usr/bin/env ruby
# sod-test.rb
# Tests for the Sod shock-tube, with helium driving air.
# PJ, 2020-05-21
#
require 'test/unit'
require 'open3'

class TestSod < Test::Unit::TestCase
  def test_0_prep
    cmd = "prep-gas ideal-helium.inp ideal-helium-gas-model.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "cp #{ENV['DGD_REPO']}/src/gas/sample-data/cea-lut-air-version-test.lua ./cea-lut-air.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "l1d4-prep --job=sod-he-air"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_1_run
    cmd = "l1d4 --run-simulation --job=sod-he-air"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    sim_time = 0.0
    lines = o.split("\n")
    lines.each do |txt|
      if txt.match('Step=350') then
        items = txt.split(' ')
        sim_time_items = items[1].split('=')
        sim_time = sim_time_items[1].to_f
      end
    end
    assert((sim_time - 3.776e-4).abs < 1.0e-5, "Incorrect sim_time at step 350.")
  end

  def test_2_post
    cmd = "l1d4 --time-slice --job=sod-he-air --tindx=40"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    f = File.new("slug-0001-tindx-0040-cells.data", "r")
    lines = f.readlines
    f.close
    columns = {'vel'=>0, 'p'=>0, 'T'=>0, 'rho'=>0}
    values = {'vel'=>0.0, 'p'=>0.0, 'T'=>0.0, 'rho'=>0.0}
    items = lines[0].split
    # Header line starts with a sharp character.
    columns['vel'] = items.index('vel') - 1
    columns['p'] = items.index('p') - 1
    columns['T'] = items.index('T') - 1
    columns['rho'] = items.index('rho') - 1
    lines.each do |txt|
      items = txt.split(' ')
      x = items[0].to_f
      if x >= 0.75 then
        values['vel'] = items[columns['vel']].to_f
        values['p'] = items[columns['p']].to_f
        values['T'] = items[columns['T']].to_f
        values['rho'] = items[columns['rho']].to_f
        break
      end
    end
    assert((values['vel'] - 444.0).abs < 1.0, "Incorrect post-shock velocity.")
    assert((values['p'] - 48472.4).abs < 100.0, "Incorrect post-shock pressure.")
    assert((values['T'] - 486.15).abs < 1.0, "Incorrect post-shock temperature.")
    assert((values['rho'] - 0.347342).abs < 0.001, "Incorrect post-shock density.")
  end

  def test_3_energies
    f = File.new("sod-he-air/energies.data", "r")
    txt = f.readlines
    f.close
    items = txt[1].split(' ')
    e0 = items[-1].to_f
    items = txt[-1].split(' ')
    e1 = items[-1].to_f
    assert((e1 - e0).abs/e0 < 0.0001, "Failed to conserve energy.")
  end

end
