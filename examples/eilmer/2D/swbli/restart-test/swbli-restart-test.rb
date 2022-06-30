#! /usr/bin/env ruby
# swbli-restart-test.rb
# Test the NK solvers restart capability using the swbli case
# Note that we turn off both the preconditioner and keep a
# constant CFL, since both of these introduce hysteresis in
# the solver. Is it worth addressing this.
# NNG, 2022-06-30

require 'test/unit'
require 'open3'

class TestRestart < Test::Unit::TestCase
  def test_0_prep
    cmd = "prep-gas gm-air.inp gm-air.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "e4shared --prep --job=swbli"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_1_run
    cmd = "e4-nk-shared --max-cpus=8 --job=swbli"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    steps = 0
    lines = o.split("\n")
    fstep_line_num = 0
    line_count = 0
    lines.each do |txt|
      line_count += 1
      if txt.match('STOPPING') then
        fstep_line_num = line_count + 1
      end
      if line_count == fstep_line_num then
        items = txt.split(' ')
        steps = items[1].to_i
      end
    end
    assert((steps - 60).abs < 1, "Failed to take correct number of steps.")
  end

  def test_2_run
    cmd = "e4-nk-shared --max-cpus=8 --job=swbli --snapshot-start=1"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    steps = 0
    lines = o.split("\n")
    fstep_line_num = 0
    line_count = 0
    lines.each do |txt|
      line_count += 1
      if txt.match('STOPPING') then
        fstep_line_num = line_count + 1
      end
      if line_count == fstep_line_num then
        items = txt.split(' ')
        steps = items[1].to_i
      end
    end
    assert((steps - 60).abs < 1, "Failed to take correct number of steps.")
  end

  def test_3_post
    cmd = 'e4shared --custom-post --script-file=post.lua'
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    lines = o.split("\n")
    pdiff = 1e99
    lines.each do |txt|
      if txt.match('pdiff:') then
        items = txt.split(' ')
        pdiff = items[1].to_f
      end
    end
    assert(pdiff < 1.0e-06, "Failed to compute a clean restart.")
  end
end
