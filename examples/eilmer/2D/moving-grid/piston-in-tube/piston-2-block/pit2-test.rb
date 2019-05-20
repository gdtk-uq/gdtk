# pit2-test.rb
# Tests for the piston-in-tube example, 2 FluidBlock case.
# PJ, 2019-05-20
#
require 'test/unit'
require 'open3'

class TestPIT2 < Test::Unit::TestCase
  def test_0_prep
    cmd = "prep-gas ideal-air.inp ideal-air-gas-model.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "e4shared --prep --job=pit2"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_1_run
    cmd = "mpirun -np 2 e4mpi --job=pit2 --run"
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
    assert((steps - 831).abs < 3, "Failed to take correct number of steps.")
  end

  def test_2_post
    f = File.new("piston.data", "r")
    txt = f.readlines[-1]
    f.close
    items = txt.split(' ')
    t = items[0].to_f
    x = items[1].to_f
    v = items[2].to_f
    assert((t - 0.040).abs < 0.0001, "Failed to reach correct final time.")
    assert((x - 6.550).abs < 0.1, "Failed to reach correct position.")
    assert((v - 276.6).abs < 1.0, "Failed to reach correct velocity.")
  end

  def test_3_post
    cmd = 'e4shared --custom-script --script-file=balanceCheck.lua'
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    energy_error = 1.0e6
    lines = o.split("\n")
    lines.each do |txt|
      if txt.match('Energy-error =') then
        items = txt.split(' ')
        energy_error = items[2].to_f
      end
    end
    assert((energy_error - 229).abs < 100.0, "Failed to get small energy error.")
  end
end
