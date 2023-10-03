#! /usr/bin/env ruby
# nenzf1d-test.rb
# Exercise the NENZF1d program for a number of cases.
# PJ, 2021-10-21
#
require 'test/unit'
require 'open3'

class TestNENZF1D < Test::Unit::TestCase
  def test_0_prep_gas_files_air
    cmd = "cp #{ENV['DGD_REPO']}/src/gas/sample-data/cea-air5species-gas-model.lua ."
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "cp #{ENV['DGD_REPO']}/src/gas/sample-data/air-5sp-eq.lua ."
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "cp #{ENV['DGD_REPO']}/examples/kinetics/air-chemistry-1T/air-5sp-1T.inp ."
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "cp #{ENV['DGD_REPO']}/examples/kinetics/air-chemistry-1T/GuptaEtAl-air-reactions.lua ."
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "prep-gas air-5sp-1T.inp air-5sp-1T.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "prep-chem air-5sp-1T.lua GuptaEtAl-air-reactions.lua air-5sp-1T-reactions.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "cp #{ENV['DGD_REPO']}/examples/kinetics/air-chemistry-2T/air-5sp-gas-model.lua ./air-5sp-2T.inp"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "cp #{ENV['DGD_REPO']}/examples/kinetics/air-chemistry-2T/GuptaEtAl-air-reactions-2T.lua ./GuptaEtAl-air-2T.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "cp #{ENV['DGD_REPO']}/examples/kinetics/air-chemistry-2T/air-energy-exchange.lua ./air-energy-exchange.inp"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "sed -i s+11-species-air+5-species-air+g GuptaEtAl-air-2T.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "prep-gas air-5sp-2T.inp air-5sp-2T.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "prep-chem air-5sp-2T.lua GuptaEtAl-air-2T.lua air-5sp-6r-2T.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "prep-kinetics air-5sp-2T.lua air-5sp-6r-2T.lua air-energy-exchange.inp air-energy-exchange.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_1_t4m8_11742
    cmd = "nenzf1d t4m8-11742.yaml"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    found_exit_condition = false
    enthalpy = 0.0; mach = 0.0; temperature = 0.0; pressure = 0.0
    massf_NO = 0.0
    lines = o.split("\n")
    lines.each do |txt|
      if txt.match('Exit condition') then
        found_exit_condition = true
      end
      if txt.match('H5s-H1') then
        items = txt.split(' ')
        enthalpy = items[1].to_f
      end
      if found_exit_condition then
        if txt.match('Mach') then
          items = txt.split(' ')
          mach = items[1].to_f
        end
        if txt.match('temperature') then
          items = txt.split(' ')
          temperature = items[1].to_f
        end
        if txt.match('pressure') then
          items = txt.split(' ')
          pressure = items[1].to_f
        end
        if txt.match('massf\[NO\]') then
          items = txt.split(' ')
          massf_NO = items[1].to_f
        end
      end
    end
    assert((enthalpy - 4.79689).abs < 0.01, "Failed to compute correct enthalpy.")
    assert(found_exit_condition, "Failed to find exit condition in output.")
    assert((mach - 7.15261).abs < 0.01, "Failed to compute correct Mach number.")
    assert((temperature - 432.096).abs < 0.1, "Failed to compute correct temperature.")
    assert((pressure - 4.69785).abs < 0.01, "Failed to compute correct pressure.")
    assert((massf_NO - 0.068366).abs < 0.001, "Failed to compute correct mass-fraction of NO.")
  end

  def test_2_t4m7_air_2T
    cmd = "nenzf1d t4m7-air-2T.yaml"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    found_exit_condition = false
    enthalpy = 0.0; mach = 0.0; temperature = 0.0; pressure = 0.0
    massf_NO = 0.0; tvib = 0.0
    lines = o.split("\n")
    lines.each do |txt|
      if txt.match('Exit condition') then
        found_exit_condition = true
      end
      if txt.match('H5s-H1') then
        items = txt.split(' ')
        enthalpy = items[1].to_f
      end
      if found_exit_condition then
        if txt.match('Mach') then
          items = txt.split(' ')
          mach = items[1].to_f
        end
        if txt.match('temperature') then
          items = txt.split(' ')
          temperature = items[1].to_f
        end
        if txt.match('pressure') then
          items = txt.split(' ')
          pressure = items[1].to_f
        end
        if txt.match('massf\[NO\]') then
          items = txt.split(' ')
          massf_NO = items[1].to_f
        end
        if txt.match('T_modes\[0\]') then
          items = txt.split(' ')
          tvib = items[1].to_f
        end
      end
    end
    assert((enthalpy - 2.47976).abs < 0.01, "Failed to compute correct enthalpy.")
    assert(found_exit_condition, "Failed to find exit condition in output.")
    assert((mach - 7.64642).abs < 0.01, "Failed to compute correct Mach number.")
    assert((temperature - 219.023).abs < 0.1, "Failed to compute correct temperature.")
    assert((pressure - 2.81165).abs < 0.01, "Failed to compute correct pressure.")
    assert((massf_NO - 0.00972722).abs < 0.001, "Failed to compute correct mass-fraction of NO.")
    assert((tvib - 1385.79).abs < 0.1, "Failed to compute correct vibrational temperature.")
  end

  def test_3_t4m8_11742_python
    cmd = "python3 py-t4m8-11742.py"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    found_exit_condition = false
    enthalpy = 0.0; mach = 0.0; temperature = 0.0; pressure = 0.0
    massf_NO = 0.0
    # puts o
    lines = o.split("\n")
    lines.each do |txt|
      if txt.match('Exit condition') then
        found_exit_condition = true
      end
      if txt.match('H5s-H1') then
        items = txt.split(' ')
        enthalpy = items[1].to_f
      end
      if found_exit_condition then
        if txt.match('Mach') then
          items = txt.split(' ')
          mach = items[1].to_f
        end
        if txt.match('temperature') then
          items = txt.split(' ')
          temperature = items[1].to_f
        end
        if txt.match('pressure') then
          items = txt.split(' ')
          pressure = items[1].to_f
        end
        if txt.match('massf\[NO\]') then
          items = txt.split(' ')
          massf_NO = items[1].to_f
        end
      end
    end
    assert((enthalpy - 4.79689).abs < 0.01, "Failed to compute correct enthalpy.")
    assert(found_exit_condition, "Failed to find exit condition in output.")
    assert((mach - 7.15261).abs < 0.01, "Failed to compute correct Mach number.")
    assert((temperature - 432.096).abs < 0.1, "Failed to compute correct temperature.")
    assert((pressure - 4.69785).abs < 0.01, "Failed to compute correct pressure.")
    assert((massf_NO - 0.068366).abs < 0.001, "Failed to compute correct mass-fraction of NO.")
  end
end
