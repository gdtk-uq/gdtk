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
    # steps observed for simulation at 2021-02-17
    assert((steps - 25513).abs < 200, "Failed to take correct number of steps.")
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
    # The old value 0.0005409 is probably better.
    # We need to fix the shock-fitting that is currently producing
    # (1) a 'stiff' final bit of shock toward the axis and
    # (2) a wobbly final bit of shock far from the body.
    # PJ 2019-11-09
    # shock_ref = 0.0005426
    shock_ref = 0.0005438 # changed by PJ 2021-02-17
    assert((shock_location - shock_ref).abs < 3.0e-6,
           "Failed to get correct shock location.")
  end
end
