# test_ideal_gas_flow.py
#
# To run on Linux Mint 21.1:
# $ PYTHONDONTWRITEBYTECODE=1 pytest-3
#
# PJ, 2019-12-29, Ported from the cfpylib collection.
#     2023-06-12 convert to pytest

import math
import gdtk.ideal_gas_flow as igf
import pytest

# This is used to change to local directory so that subprocess runs nicely.
@pytest.fixture(autouse=True)
def change_test_dir(request, monkeypatch):
    monkeypatch.chdir(request.fspath.dirname)


def test_0_isentropic_flow():
    print("Begin test of isentropic flow ratios...")
    M = 2.0
    print("Computed: M=%g: A/Astar=%g, T0/T=%g, p0/p=%g, r0/r=%g" %
          (M, igf.A_Astar(M), igf.T0_T(M), igf.p0_p(M), igf.r0_r(M)))
    print("Expected: M=2, A/Astar=1.687, T0/T=1.80, p0/p=7.824, r0/r=4.347")
    assert math.isclose(igf.A_Astar(M), 1.687, rel_tol=1.0e-2, abs_tol=1.0e-5), "isentropic flow, area change from Mach number"
    assert math.isclose(igf.T0_T(M), 1.80, rel_tol=1.0e-2, abs_tol=1.0e-5), "isentropic flow, temperature ratio"
    assert math.isclose(igf.p0_p(M), 7.824, rel_tol=1.0e-2, abs_tol=1.0e-5), "isentropic flow, pressure ratio"
    assert math.isclose(igf.r0_r(M), 4.347, rel_tol=1.0e-2, abs_tol=1.0e-5), "isentropic flow, density ratio"
    print("")
    return

def test_1_normal_shock():
    print("Normal shock jump...")
    M = 2.0
    print("Computed: M=%g: M2=%g, T2/T1=%g, p2/p1=%g, r2/r1=%g" %
          (M, igf.m2_shock(M), igf.T2_T1(M), igf.p2_p1(M), igf.r2_r1(M)))
    print("Expected: M1=2, M2=0.5774, T2/T1=1.687, p2/p1=4.50, r2/r1=2.667")
    assert math.isclose(igf.m2_shock(M), 0.5774, rel_tol=1.0e-2, abs_tol=1.0e-5), "normal shock jump, Mach number"
    assert math.isclose(igf.T2_T1(M), 1.687, rel_tol=1.0e-2, abs_tol=1.0e-5), "normal shock, temperature ratio"
    assert math.isclose(igf.p2_p1(M), 4.50, rel_tol=1.0e-2, abs_tol=1.0e-5), "normal shock, pressure ratio"
    assert math.isclose(igf.r2_r1(M), 2.667, rel_tol=1.0e-2, abs_tol=1.0e-5), "normal shock, density ratio"
    print("")
    return

def test_2_rayleigh_line():
    print("Rayleigh-line flow...")
    M = 2.0
    print("Computed: M=%g: T0/Tstar=%g, T/Tstar=%g, p/pstar=%g, r/rstar=%g" %
          (M, igf.T0_T0star(M), igf.T_Tstar(M), igf.p_pstar(M), igf.r_rstar(M)))
    print("Expected: M=2, T0/T0star=0.7934, T/Tstar=0.5289, p/pstar=0.3636, r/rstar=0.6875")
    assert math.isclose(igf.T0_T0star(M), 0.7934, rel_tol=1.0e-2, abs_tol=1.0e-5), "Rayleigh-line, total temperature ratio"
    assert math.isclose(igf.T_Tstar(M), 0.5289, rel_tol=1.0e-2, abs_tol=1.0e-5), "Rayleigh-line, static temperature ratio"
    assert math.isclose(igf.p_pstar(M), 0.3636, rel_tol=1.0e-2, abs_tol=1.0e-5), "Rayleigh-line, pressure ratio"
    assert math.isclose(igf.r_rstar(M), 0.6875, rel_tol=1.0e-2, abs_tol=1.0e-5), "Rayleigh-line, density ratio"
    print("Inverse calculation: T0/T0star=%g --> M=%g" %
          (igf.T0_T0star(M), igf.M_Rayleigh(igf.T0_T0star(M))))
    assert math.isclose(igf.M_Rayleigh(igf.T0_T0star(M)), M, rel_tol=1.0e-2, abs_tol=1.0e-5), "Rayleigh-line, inverse calculation"
    print("")
    return

def test_3_prandtl_meyer():
    print("Prandtl-Meyer function...")
    M = 2
    print("Computed: M=%g --> nu=%g; Inverse: M=%g <-- nu=%g" %
          (M, igf.PM1(M), igf.PM2(1.1481), 1.1481))
    print("Expected: M=2 --> nu=0.4604; Inverse: M=4 <-- nu=1.1481")
    assert math.isclose(igf.PM1(M), 0.4604, rel_tol=1.0e-2, abs_tol=1.0e-5), "Prandtl-Meyer function"
    assert math.isclose(igf.PM2(1.1481), 4, rel_tol=1.0e-2, abs_tol=1.0e-5), "Prandtl-Meyer function, inverse"
    print("")
    return

def test_4_oblique_shock():
    print("Oblique shock relations may not quite match (data is from chart)...")
    M = 2.0
    beta = math.radians(44.0); theta = math.radians(14.0); # from chart, M=2
    print("Computed: M1=%g, theta(beta=%g)=%g, beta(theta=%g)=%g" % \
          (M, beta, igf.theta_obl(M, beta), theta, igf.beta_obl(M, theta)))
    assert math.isclose(igf.theta_obl(M, beta), theta, rel_tol=1.0e-2, abs_tol=1.0e-5), "oblique shock, theta from beta"
    assert math.isclose(igf.beta_obl(M, theta), beta, rel_tol=1.0e-2, abs_tol=1.0e-5), "oblique shock, beta from theta"
    print("Conditions behind shock:")
    print("M2=%g, expected 1.482 (from chart, 14 degree deflection)" %
          igf.M2_obl(M, beta, theta))
    assert math.isclose(igf.M2_obl(M, beta, theta), 1.482, rel_tol=1.0e-2, abs_tol=1.0e-5), "oblique shock, Mach number after"
    print("Computed: T2/T1=%g, p2/p1=%g, r2/r1=%g" %
          (igf.T2_T1_obl(M, beta), igf.p2_p1_obl(M, beta), igf.r2_r1_obl(M, beta)))
    print("Expected: T2/T1=1.249, p2/p1=2.088, r2/r1=1.673 (approx. normal-shock table M=1.390)")
    assert math.isclose(igf.T2_T1_obl(M, beta), 1.249, rel_tol=1.0e-2, abs_tol=1.0e-5), "oblique shock, temperature ratio"
    assert math.isclose(igf.p2_p1_obl(M, beta), 2.088, rel_tol=1.0e-2, abs_tol=1.0e-5), "oblique shock, pressure ratio"
    assert math.isclose(igf.r2_r1_obl(M, beta), 1.673, rel_tol=1.0e-2, abs_tol=1.0e-5), "oblique shock, density ratio"
    print("v2/v1=%g, p02/p01=%g" %
          (igf.v2_v1_obl(M, beta), igf.p02_p01_obl(M, beta)))
    print("Expected: v2/v1=0.8304=sin(B)/sin(B-d)*r1/r2")
    assert math.isclose(igf.v2_v1_obl(M, beta), 0.8304, rel_tol=1.0e-2, abs_tol=1.0e-5), "oblique shock, velocity ratio"
    print("")
    return

def test_5_taylor_maccoll():
    M1 = 1.5; p1 = 100.0e3; T1 = 300.0; R = 287.1; g = 1.4; rho1 = p1/(R*T1)
    print("Taylor-Maccoll cone flow demo with M1=%g" % M1)
    print("for M1=1.5, beta=49 degrees, expect theta=20 degrees from NACA1135.")
    a1 = math.sqrt(1.4*287*T1)
    V1 = M1 * a1
    beta = math.radians(49.0)
    theta_c, V_c, p_c, T_c = igf.theta_cone(V1, p1, T1, beta)
    print("theta_c(degrees)=", math.degrees(theta_c), "expected 20 degrees, surface speed V_c=", V_c)
    print("surface pressure coefficient=", (p_c - p1)/(0.5*rho1*V1*V1), "expected 0.385")
    print("p_c: %g, T_c: %g" % (p_c, T_c))
    assert math.isclose(math.degrees(theta_c), 20.0, rel_tol=1.0e-2, abs_tol=1.0e-5), "conical shock, cone angle"
    assert math.isclose((p_c - p1)/(0.5*rho1*V1*V1), 0.385, rel_tol=1.0e-2, abs_tol=1.0e-5), "conical shock, surface pressure coefficient"
    print("")
    print("Conical shock from cone with half-angle 20 degrees in M1=", M1)
    beta = igf.beta_cone(V1, p1, T1, math.radians(20.0))
    print("sigma(degrees)=", math.degrees(beta), "expected 49 degrees")
    assert math.isclose(math.degrees(beta), 49.0, rel_tol=1.0e-2, abs_tol=1.0e-5), "conical shock, shock angle from deflection"
    print("Repeat above test, but call beta_cone2()")
    beta = igf.beta_cone2(M1, math.radians(20.0))
    print("sigma(degrees)=", math.degrees(beta), "expected 49 degrees")
    assert math.isclose(math.degrees(beta), 49.0, rel_tol=1.0e-2, abs_tol=1.0e-5), "conical shock, shock angle from deflection and Mach"
    #
    print("Done.")
    return

if __name__ == '__main__':
    test_0_isentropic_flow()
    test_1_normal_shock()
    test_2_rayleigh_line()
    test_3_prandtl_meyer()
    test_4_oblique_shock()
    test_5_taylor_maccoll()
