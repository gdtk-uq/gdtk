#! /usr/bin/env ruby
# connection-shared-memory-test.rb
# Functional test for structured block connections in 3D.
# RJG, 2018-01-23, PJ, 2019-06-27
#
require 'test/unit'
require 'open3'

class TestConnection < Test::Unit::TestCase
  def test_0_generate_ref_result
    cmd = 'python3 generate-ref-result.py'
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_1_run_connection_tests
    cmd = 'python3 run-cases.py'
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    assert(!o.match("FAILED"), "Failed connection test")
  end
end
