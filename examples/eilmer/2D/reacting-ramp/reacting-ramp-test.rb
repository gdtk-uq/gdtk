#! /usr/bin/env ruby
# reacting-ramp-test.rb
# Tests for the NK solver using the svan_albada limiter
# NNG, 2023-02-21
#
require 'test/unit'
require 'open3'

class TestRR < Test::Unit::TestCase
  def test_0_grid
    cmd = "e4shared --custom-script --script-file=su2-grid/su2-grid-gen.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "rm -rf su2-grid/block_0_ramp15.su2"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "rm -rf su2-grid/block_1_ramp15.su2"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "rm -rf su2-grid/block_2_ramp15.su2"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "rm -rf su2-grid/block_3_ramp15.su2"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)

    cmd = "ugrid_partition ramp15.su2 mapped_cells 4 2"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)

    cmd = "mv block_0_ramp15.su2 block_1_ramp15.su2 block_2_ramp15.su2 block_3_ramp15.su2 su2-grid"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_1_prep
    cmd = "prep-gas species.inp gas-model.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "prep-chem gas-model.lua reaction_mechanism.lua reac-file.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "e4shared --job=ramp15 --prep"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_2_run
    cmd = "e4-nk-shared --job=ramp15 --verbosity=1"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    steps = 0
    lines = o.split("\n")
    fstep_line_num = 0
    line_count = 0
    lines.each do |txt|
      line_count += 1
      if txt.match('STOPPING') then
        fstep_line_num = line_count + 2
      end
      if line_count == fstep_line_num then
        items = txt.split(' ')
        steps = items[1].to_i
      end
    end
    assert((steps - 551).abs < 20, "Failed to take correct number of steps.")
  end

end
