#! /usr/bin/env ruby
# inf_cyl-test.rb
# Functional test for 2D, inviscid flow in thermal nonequilibrium.
# This exercises quite a few of the basic functions of the 2D code
# plus Rowan's two-temperature nitrogen model.
#
# PJ, 13-Jan-2011
#     01-Aug-2017 Ported to Eilmer4
#     2019-05-07 Increase tolerances to accommodate no-ghost-cell BC. 
#     2019-06-24 Ruby flavour.
#
require 'test/unit'
require 'open3'

class TestInfCyl < Test::Unit::TestCase
  def test_0_prep
    cmd = "e4shared --prep --job=inf_cyl"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_1_run
    cmd = "e4shared --run --job=inf_cyl --verbosity=1 --max-cpus=4"
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
    assert((steps - 7124).abs < 50, "Failed to take correct number of steps.")
  end

  def test_2_stagnation_point_temperatures
    cmd = 'e4shared --job=inf_cyl --post --slice-list=0,:,0,0 '+
          '--output-file=stag-prof-50Pa-Blackman.data'
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    f = File.new('stag-prof-50Pa-Blackman.data', 'r')
    lines = f.readlines()
    f.close
    #
    columns = {}
    txt = lines[0]
    columns['T_tr'] = txt.match(/(\d+):T\s/)[1].to_i
    columns['T_v'] = txt.match(/(\d+):T_modes\[0\]\s/)[1].to_i
    #
    values = {}
    txt = lines[-1]
    items = txt.split(' ')
    values['T_tr'] = items[columns['T_tr']-1].to_f
    values['T_v'] = items[columns['T_v']-1].to_f
    #
    assert((values['T_tr'] - 2536.4).abs < 20.0, "Failed to see correct translational T.")
    assert((values['T_v'] - 2349.0).abs < 60.0, "Failed to see correct vibrational T.")
  end
end
