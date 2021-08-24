#! /usr/bin/env ruby
# flow-format-test.rb
# This exercises the file format IO using the forward facing step flow example
# This test simply checks for the successful completion of a sequence of commands
#
# DB, 2021-08-24 - Creation

require 'test/unit'
require 'open3'

class TestBack < Test::Unit::TestCase
 
  def test_1_restart
    Dir.chdir('restart_test'){
      cmds = [
        "prep-gas ideal-air.inp ideal-air-gas-model.lua",

        "e4shared --prep --job=ffs",
        "e4shared --run --job=ffs --verbosity=1 --max-cpus=3",
        "e4shared --post --job=ffs --tindx-plot=all --vtk-xml",

        "e4shared --prep --job=ffs_restart_1",
        "e4shared --run --job=ffs_restart_1 --verbosity=1 --max-cpus=3",
        "e4shared --post --job=ffs_restart_1 --tindx-plot=all --vtk-xml",

        "e4shared --prep --job=ffs_restart_2",
        "e4shared --run --job=ffs_restart_2 --verbosity=1 --max-cpus=3",
        "e4shared --post --job=ffs_restart_2 --tindx-plot=all --vtk-xml"
      ]

      cmds.each do |cmd|
        o, e, s = Open3.capture3(*cmd.split)
        assert_equal(s.success?, true)
      end

    }
  end

  def test_2_multi_output
    
    Dir.chdir('multi_output'){
      cmds = [
        "prep-gas ideal-air.inp ideal-air-gas-model.lua",

        "e4shared --prep --job=ffs",
        "e4shared --run --job=ffs --verbosity=1 --max-cpus=3",

        "e4shared --post --job=ffs --tindx-plot=all --vtk-xml",
        
        "e4shared --post --job=ffs --tindx-plot=\"1,2,3,4,5\" --vtk-xml --plotTag=average",
        "e4shared --post --job=ffs --tindx-plot=last --vtk-xml --plotTag=DFT",
        "e4shared --post --job=ffs --tindx-plot=\"1,2,3,4,5\" --vtk-xml --plotTag=gradient",
      ]

      cmds.each do |cmd|
        o, e, s = Open3.capture3(*cmd.split)
        assert_equal(s.success?, true)
      end
    }
  end

end
