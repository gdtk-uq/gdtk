# test_oblique_shock.rb
#
# $ cp ${DGD_REPO}/src/gas/sample-data/cea-air13species-gas-model.lua .
# $ ruby test_oblique_shock.rb
#
# PJ, 2019-12-01
#
$LOAD_PATH << '~/gdtkinst/lib'
require 'gdtk/gas'
require 'test/unit'
require 'open3'

class TestObliqueShock < Test::Unit::TestCase
  def test_0_prep
    cmd = "cp #{ENV['DGD_REPO']}/src/gas/sample-data/air-13sp-eq-gas-model.lua ."
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "cp #{ENV['DGD_REPO']}/src/gas/sample-data/air-13sp-1T-input.lua ."
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "prep-gas air-13sp-1T-input.lua air-13sp-gas-model.lua"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_1_normal_shock_given_speed
    m1 = 1.5
    # puts "Oblique-shock demo for m1=%g" % [m1]

    gmodel = GasModel.new('air-13sp-eq-gas-model.lua')
    state1 = GasState.new(gmodel)
    state1.p = 100.0e3 # Pa
    state1.T = 300.0 # K ideal air, not high T
    state1.update_thermo_from_pT()
    state1.update_sound_speed()
    # puts "state1: %s" % state1
    v1 = m1*state1.a
    # puts "v1=%g" % [v1]
    assert((state1.rho - 1.1566).abs < 1.0e-3, "Failed to compute correct initial density.")
    assert((v1 - 521.58).abs < 1.0e-1, "Failed to compute correct initial velocity.")

    beta = 45.0 * Math::PI/180.0
    # puts "  given beta(degrees)=%g" % [beta*180/Math::PI]
    state2 = GasState.new(gmodel)
    flow = GasFlow.new(gmodel)
    theta, v2 = flow.theta_oblique(state1, v1, beta, state2)
    # puts "  theta=%g degrees" % [theta*180/Math::PI]
    # puts "  v2=%g" % [v2]
    # puts "  state2: %s" % [state2]
    assert((state2.rho - 1.2748).abs/state2.rho < 1.0e-2, "Failed to compute correct post-shock density.")
    assert((theta*180/Math::PI - 2.784).abs < 1.0e-3, "Failed to compute correct deflection angle.")
    assert((v2 - 498.0).abs < 1.0, "Failed to compute correct post-shock velocity.")

    # puts "Oblique shock angle from deflection."
    beta2 = flow.beta_oblique(state1, v1, theta)
    # puts "  beta2(degrees)=%g" % [beta2*180/Math::PI]
    assert((beta2*180/Math::PI - 45.0).abs < 1.0e-3, "Failed to compute correct shock angle.")
  end
end
