#! /usr/bin/env ruby
# lobb-test.rb
# This exercises the nonequilibrium thermochemistry plus shock-fitting
# with the flow of reacting air over a sphere.
# RJG, 2016-01-25 (Tcl test) PJ, 2019-11-12 (Ruby port)
#
require 'test/unit'
require 'open3'

class TestLobb < Test::Unit::TestCase
  def test_0_prep
    cmd = "cp #{ENV['DGD_REPO']}/examples/kinetics/air-chemistry-1T/air-5sp-1T.inp air-5sp.inp"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "cp #{ENV['DGD_REPO']}/examples/kinetics/air-chemistry-1T/GuptaEtAl-air-reactions.lua ."
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "prep-gas air-5sp.inp air-5sp.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "prep-chem air-5sp.lua GuptaEtAl-air-reactions.lua air-5sp-6r.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "e4shared --prep --job=lobb"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_1_run
    cmd = "e4shared --run --job=lobb --verbosity=1 --max-cpus=4"
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
    # steps observed for simulation with backward-euler update at 2021-11-28
    assert((steps - 15530).abs < 140, "Failed to take correct number of steps.")
  end

  def test_2_shock_stand_off
    cmd = 'e4shared --job=lobb --custom-post --script-file=shock-detachment.lua'
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    shock_location = 0.0
    lines = o.split("\n")
    lines.each do |txt|
      if txt.match('shock-detachment=') then
        items = txt.split(' ')
        shock_location = items[1].to_f
      end
    end
    shock_ref = 0.0005384 # changed by PJ 2021-11-28
    # PJ 2021-08-05 Doing a bit more work on the shock-fitting code
    # shifts the position of the vertex a little, to a position of 0.0005500m.
    # However, looking at the temperature field near the axis shows that
    # the first face (length nearly 0.00045, yes, close to the stand-off distance)
    # varies in x-position across a range of 0.00002m, or a relative change of 0.037.
    # Small (and hopefully benign) changes to the simulation code will most likey
    # lead to changes in shock location that are a reasonable fraction of 0.00002m.
    # Let's just change the tolerance on shock location to be a couple of percent.
    assert((shock_location - shock_ref).abs/shock_ref < 0.02,
           "Failed to get correct shock location.")
  end
end
