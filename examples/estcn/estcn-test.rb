#! /usr/bin/env ruby
# estcn-test.rb
# Tests for the ESTCN program.
# PJ, 2020-06-07
#
require 'test/unit'
require 'open3'

class TestESTCN < Test::Unit::TestCase
  def test_0_prepare_gas_models
    cmd = "cp #{ENV['DGD_REPO']}/src/gas/sample-data/ideal-air-gas-model.lua ."
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "cp #{ENV['DGD_REPO']}/src/gas/sample-data/cea-air5species-gas-model.lua ."
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "cp #{ENV['DGD_REPO']}/src/gas/sample-data/cea-air13species-gas-model.lua ."
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "cp #{ENV['DGD_REPO']}/src/gas/sample-data/cea-lut-air-version-test.lua ./cea-lut-air.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "prep-gas he-n2-o2-mix.inp he-n2-o2-mix-gas-model.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_1_shock_tunnel_lut
    cmd = "estcn --task=stn --gas=cea-lut-air.lua --T1=300 --p1=125.0e3 --Vs=2414 --pe=34.37e6 --ar=27.0"
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
      if found_state_7 && pressure == 0.0 then
        items = txt.split(' ')
        pressure = items[1].to_f
        temperature = items[4].to_f
      end
      if txt.match('State 7') then
        # Note that we want the data from the following line.
        found_state_7 = true
      end
    end
    assert((enthalpy - 5.42908e+06).abs/5.42908e+06 < 1.0e-3, "Incorrect enthalpy.")
    assert((pressure - 93702.4).abs/93702.4 < 1.0e-3, "Incorrect pressure.")
    assert((temperature - 1283.6).abs/1283.6 < 1.0e-3, "Incorrect temperature.")
    assert((velocity - 2950.34).abs/2950.34 < 1.0e-3, "Incorrect velocity.")
  end

  def test_2_shock_tunnel_cea2_gas
    cmd = "estcn --task=stn --gas=cea-air5species-gas-model.lua --T1=300 --p1=125.0e3 --Vs=2414 --pe=34.37e6 --ar=27.0"
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
      if found_state_7 && pressure == 0.0 then
        items = txt.split(' ')
        pressure = items[1].to_f
        temperature = items[4].to_f
      end
      if txt.match('State 7') then
        # Note that we want the data from the following line.
        found_state_7 = true
      end
    end
    assert((enthalpy - 5.43258e+06).abs/5.43258e+06 < 1.0e-3, "Incorrect enthalpy.")
    assert((pressure - 93940.5).abs/93940.5 < 1.0e-3, "Incorrect pressure.")
    assert((temperature - 1280.98).abs/1280.98 < 1.0e-3, "Incorrect temperature.")
    assert((velocity - 2949.81).abs/2949.81 < 1.0e-3, "Incorrect velocity.")
  end

  def test_3_incident_shock_cea2_gas
    cmd = "estcn --task=ishock --gas=cea-air13species-gas-model.lua --p1=59 --T1=283 --Vs=11296"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    velocity = 0.0
    found_state_2 = false
    pressure = 0.0; temperature = 0.0
    lines = o.split("\n")
    lines.each do |txt|
      if txt.match('Vg') then
        items = txt.split(' ')
        velocity = items[4].to_f
      end
      if found_state_2 && pressure == 0.0 then
        items = txt.split(' ')
        pressure = items[1].to_f
        temperature = items[4].to_f
      end
      if txt.match('State 2') then
        # Note that we want the data from the following line.
        found_state_2 = true
      end
    end
    assert((pressure - 86686.0).abs/86686.0 < 1.0e-3, "Incorrect pressure.")
    assert((temperature - 12034.0).abs/12034.0 < 1.0e-3, "Incorrect temperature.")
    assert((velocity - 10561.1).abs/10561.1 < 1.0e-3, "Incorrect velocity.")
  end

  def test_4_incident_shock_ideal_air_gas
    cmd = "estcn --task=ishock --gas=ideal-air-gas-model.lua --p1=59 --T1=283 --Vs=11296"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    velocity = 0.0
    found_state_2 = false
    pressure = 0.0; temperature = 0.0
    lines = o.split("\n")
    lines.each do |txt|
      if txt.match('Vg') then
        items = txt.split(' ')
        velocity = items[4].to_f
      end
      if found_state_2 && pressure == 0.0 then
        items = txt.split(' ')
        pressure = items[1].to_f
        temperature = items[4].to_f
      end
      if txt.match('State 2') then
        # Note that we want the data from the following line.
        found_state_2 = true
      end
    end
    assert((pressure - 77204.1).abs/77204.1 < 1.0e-3, "Incorrect pressure.")
    assert((temperature - 61994.8).abs/61994.8 < 1.0e-3, "Incorrect temperature.")
    assert((velocity - 9404.94).abs/9404.94 < 1.0e-3, "Incorrect velocity.")
  end

  def test_5_pitot_pressure
    cmd = "estcn --gas=cea-lut-air.lua --task=pitot --p1=93.6e3 --T1=1284 --V1=2.95e3"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    lines = o.split("\n")
    items = lines[-2].split(' ')
    pressure= items[1].to_f
    temperature = items[4].to_f
    assert((pressure - 2.1462e+06).abs/2.1462e+06 < 1.0e-3, "Incorrect pressure.")
    assert((temperature - 3875.76).abs/3875.76 < 1.0e-3, "Incorrect temperature.")
  end

  def test_6_cone_pressure
    cmd = "estcn --gas=cea-lut-air.lua --task=cone --sigma-deg=15 --p1=93.6e3 --T1=1284 --V1=2.95e3"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    lines = o.split("\n")
    items = lines[-2].split(' ')
    pressure = items[1].to_f
    temperature = items[4].to_f
    angle_deg = 0.0
    lines.each do |txt|
      if txt.match('Shock angle') then
        items = txt.split(' ')
        angle_deg = items[-2].to_f
      end
    end
    assert((pressure - 271013.0).abs/271013.0 < 1.0e-3, "Incorrect pressure.")
    assert((temperature - 1680.1).abs/1680.1 < 1.0e-3, "Incorrect temperature.")
    assert((angle_deg - 21).abs/21 < 1.0e-3, "Incorrect shock angle.")
  end

  def test_7_total_pressure
    cmd = "estcn --gas=cea-lut-air.lua --task=total --p1=93.6e3 --T1=1284 --V1=2.95e3"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    lines = o.split("\n")
    items = lines[-2].split(' ')
    pressure = items[1].to_f
    temperature = items[4].to_f
    assert((pressure - 3.4273e+07).abs/3.4273e+07 < 1.0e-3, "Incorrect pressure.")
    assert((temperature - 4160.5).abs/4160.5 < 1.0e-3, "Incorrect temperature.")
  end

  def test_8_incident_shock_thermally_perfect_air_gas
    cmd = "estcn --task=ishock --gas=he-n2-o2-mix-gas-model.lua --T1=300 --p1=125.0e3 --molef=0.0,0.79,0.21 --Vs=2414.0"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    velocity = 0.0
    found_state_2 = false
    pressure = 0.0; temperature = 0.0
    lines = o.split("\n")
    lines.each do |txt|
      if txt.match('Vg') then
        items = txt.split(' ')
        velocity = items[4].to_f
      end
      if found_state_2 && pressure == 0.0 then
        items = txt.split(' ')
        pressure = items[1].to_f
        temperature = items[4].to_f
      end
      if txt.match('State 2') then
        # Note that we want the data from the following line.
        found_state_2 = true
      end
    end
    assert((pressure - 7.8815e+06).abs/7.8815e+06 < 1.0e-3, "Incorrect pressure.")
    assert((temperature - 2792.97).abs/2792.97 < 1.0e-3, "Incorrect temperature.")
    assert((velocity - 2057.56).abs/2057.56 < 1.0e-3, "Incorrect velocity.")
  end

end
