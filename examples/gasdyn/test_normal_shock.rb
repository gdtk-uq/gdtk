# test_normal_shock.rb
#
# $ ruby test_normal_shock.rb
#
# PJ, 2019-11-28
#     2023-11-15 Change to use Nick's equilibrium model and automate tests.
#
$LOAD_PATH << '~/gdtkinst/lib'
require 'gdtk/gas'
require 'test/unit'
require 'open3'

class TestNormalShock < Test::Unit::TestCase
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
    gmodel = GasModel.new('air-13sp-eq-gas-model.lua')
    state1 = GasState.new(gmodel)
    state1.p = 125.0e3 # Pa
    state1.T = 300.0 # K
    state1.update_thermo_from_pT()
    state1.update_sound_speed()
    # puts "state1: %s" % state1
    assert((state1.rho - 1.44578).abs < 1.0e-3, "Failed to compute correct initial density.")

    # puts "normal shock, given shock speed"
    vs = 2414.0
    # puts "vs=%g" % [vs]
    state2 = GasState.new(gmodel)
    flow = GasFlow.new(gmodel)
    v2, vg = flow.normal_shock(state1, vs, state2)
    # puts "v2=%g vg=%g" % [v2, vg]
    # puts "state2: %s" % [state2]
    assert((state2.p - 7.29e6).abs/state2.p < 1.0e-2, "Failed to compute correct post-shock pressure.")
    assert((v2 - 361.16).abs < 1.0, "Failed to compute correct post-shock velocity.")
    assert((vg - 2052.8).abs < 1.0, "Failed to compute correct lab-frame velocity.")
  end

  def test_2_normal_shock_given_pressure
    gmodel = GasModel.new('air-13sp-eq-gas-model.lua')
    state1 = GasState.new(gmodel)
    state1.p = 125.0e3 # Pa
    state1.T = 300.0 # K
    state1.update_thermo_from_pT()
    state1.update_sound_speed()
    # puts "state1: %s" % state1
    assert((state1.rho - 1.44578).abs < 1.0e-3, "Failed to compute correct initial density.")

    # puts "normal shock, given pressure ratio"
    p2p1 = 58.317
    # puts "p2p1=%g" % [p2p1]
    state2 = GasState.new(gmodel)
    flow = GasFlow.new(gmodel)
    vs, v2, vg = flow.normal_shock_p2p1(state1, p2p1, state2)
    # puts "vs=%g v2=%g vg=%g" % [vs, v2, vg]
    # puts "state2: %s" % [state2]
    assert((state2.rho - 9.6635).abs/state2.rho < 1.0e-2, "Failed to compute correct post-shock density.")
    assert((vs - 2414.0).abs < 1.0, "Failed to compute correct shock speed.")
    assert((v2 - 361.16).abs < 1.0, "Failed to compute correct post-shock velocity.")
    assert((vg - 2052.8).abs < 1.0, "Failed to compute correct lab-frame velocity.")
    
    # puts "reflected shock"
    state5 = GasState.new(gmodel)
    vr_b = flow.reflected_shock(state2, vg, state5)
    # puts "vr_b=%g" % [vr_b]
    # puts "state5: %s" % [state5]
    assert((vr_b - 572.9).abs < 1.0, "Failed to compute correct reflected-shock speed.")
    assert((state5.p - 59.37e6).abs/state5.p < 1.0e-2, "Failed to compute correct reflected-shock pressure.")
  end
end
