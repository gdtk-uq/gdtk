#! /usr/bin/env ruby
# sabcm-test.rb
# Tests using the steepening wave problem.
# NNG, 2023-12-13
#
require 'test/unit'
require 'open3'

class TestSWP < Test::Unit::TestCase
  def test_0_prep
    cmd = "rm -rf ausmdv"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)

    cmd = "mkdir -p ausmdv"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_1_run32
    cmd = "cp -r blank/swp.lua ./"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)

    cmd = "cp -r blank/ideal-air.inp ./"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)

    cmd = 'sed -i s/N=128/N=32/ swp.lua'
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)

    cmd = 'sed -i s/asf/ausmdv/ swp.lua'
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)

    cmd = "prep-gas ideal-air.inp ideal-air.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)

    cmd = "e4shared --prep --job=swp"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)

    cmd = "e4shared --run --job=swp"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)

    cmd = "mkdir -p ausmdv/032"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)

    cmd = "mv config flow grid hist ideal-air.inp ideal-air.lua loads solid swp.lua ausmdv/032"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_2_run64
    cmd = "cp -r blank/swp.lua ./"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)

    cmd = "cp -r blank/ideal-air.inp ./"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)

    cmd = 'sed -i s/N=128/N=64/ swp.lua'
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)

    cmd = 'sed -i s/asf/ausmdv/ swp.lua'
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)

    cmd = "prep-gas ideal-air.inp ideal-air.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)

    cmd = "e4shared --prep --job=swp"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)

    cmd = "e4shared --run --job=swp"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)

    cmd = "mkdir -p ausmdv/064"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)

    cmd = "mv config flow grid hist ideal-air.inp ideal-air.lua loads solid swp.lua ausmdv/064"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_3_run128
    cmd = "cp -r blank/swp.lua ./"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)

    cmd = "cp -r blank/ideal-air.inp ./"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)

    cmd = 'sed -i s/N=128/N=128/ swp.lua'
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)

    cmd = 'sed -i s/asf/ausmdv/ swp.lua'
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)

    cmd = "prep-gas ideal-air.inp ideal-air.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)

    cmd = "e4shared --prep --job=swp"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)

    cmd = "e4shared --run --job=swp"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)

    cmd = "mkdir -p ausmdv/128"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)

    cmd = "mv config flow grid hist ideal-air.inp ideal-air.lua loads solid swp.lua ausmdv/128"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_4_order
    cmd = 'python3 scripts/post_process.py ausmdv/032 ausmdv/064 ausmdv/128'
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    lines = o.split("\n")

    pL0=0; pL2=0; rhoL0=0; rhoL2=0; velL0=0; velL2=0;
    lines.each do |txt|
      if txt.match('p L0') then
        items = txt.split('=')
        pL0 = items[1].to_f
      end

      if txt.match('p L2') then
        items = txt.split('=')
        pL2 = items[1].to_f
      end

      if txt.match('rho L0') then
        items = txt.split('=')
        rhoL0 = items[1].to_f
      end

      if txt.match('rho L2') then
        items = txt.split('=')
        rhoL2 = items[1].to_f
      end

      if txt.match('vel.x L0') then
        items = txt.split('=')
        velL0 = items[1].to_f
      end

      if txt.match('vel.x L2') then
        items = txt.split('=')
        velL2 = items[1].to_f
      end
    end
    assert((pL0-2.03).abs < 1.0e-03, "Incorrect pressure L0 norm.")
    assert((pL2-2.05).abs < 1.0e-03, "Incorrect pressure L2 norm.")
    assert((rhoL0-1.96).abs < 1.0e-03, "Incorrect density L0 norm.")
    assert((rhoL2-2.03).abs < 1.0e-03, "Incorrect density L2 norm.")
    assert((velL0-2.03).abs < 1.0e-03, "Incorrect velocity L0 norm.")
    assert((velL2-2.08).abs < 1.0e-03, "Incorrect velocity L2 norm.")
  end
end
