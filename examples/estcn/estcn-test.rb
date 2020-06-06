#! /usr/bin/env ruby
# estcn-test.rb
# Tests for the ESTCN program.
# PJ, 2020-06-07
#
require 'test/unit'
require 'open3'

class TestESTCN < Test::Unit::TestCase
  def test_0_prep
    cmd = "cp #{ENV['DGD_REPO']}/src/gas/sample-data/cea-air5species-gas-model.lua ."
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "cp #{ENV['DGD_REPO']}/src/gas/sample-data/cea-lut-air-version-test.lua ./cea-lut-air.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_1_shock_tunnel_lut
    cmd = "estcn  --task=stn --gas=cea-lut-air.lua --T1=300 --p1=125.0e3 --Vs=2414 --pe=34.37e6 --ar=27.0"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    enthalpy = 0.0; velocity = 0.0
    found_state_7 = false
    pressure = 0.0; temperature = 0.0
    lines = o.split("\n")
    lines.each do |txt|
      if txt.match('Enthalpy difference') then
        items = txt.split(' ')
        enthalpy = items[-2].to_f
      end
      if txt.match('V7') then
        items = txt.split(' ')
        velocity = items[1].to_f
      end
      # The last thing that we do is look for the pressure and temperature.
      if found_state_7 then
        items = txt.split(', ')
        pressure_item = items[1]
        pressure = pressure_item.split('=')[-1].to_f
        temperature_item = items[2]
        temperature = temperature_item.split('=')[-1].to_f
        break
      end
      if txt.match('State 7') then
        found_state_7 = true
      end
    end
    assert((enthalpy - 5.42908e+06).abs/5.42908e+06 < 1.0e-3, "Incorrect enthalpy.")
    assert((pressure - 93702.4).abs/94702.4 < 1.0e-3, "Incorrect pressure.")
    assert((temperature - 1283.6).abs/1283.6 < 1.0e-3, "Incorrect temperature.")
  end

end
