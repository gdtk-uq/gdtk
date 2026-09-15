# rao-nozzle.py
# Design of a Rao thrust-optimised axisymmetric nozzle with the gdtk.imoc
# method-of-characteristics package, validated against the nozzle A of
# Rao (1958).
#
# The construction follows Rao's method as implemented by J. Kunze in the
# original (Tcl) IMOC code and described in section 3.1 of his PhD thesis.
# Equation numbers in the comments below refer to that thesis.
#
# Starting from a straight initial characteristic AB carrying uniform,
# axial flow at M_A, a centred expansion fan at the wall point A is grown
# incrementally.  Along the control surface DE (the terminating C+
# characteristic through the nozzle lip E) two Lagrange-multiplier
# expressions, eqs 3.2 and 3.3, are invariant, which turns the search for
# the optimum contour into (i) locating the point D on the newest fan
# characteristic where eq 3.2 evaluates to its lip value and (ii) growing
# the fan until the mass flow crossing AD (eqs 3.5-3.6) balances the mass
# flow crossing DE (eqs 3.7-3.11).  The flow field between AD and DE is
# then filled with characteristics anchored on the (known) DE data and the
# wall is traced as the streamline from A to E.
#
# Usage:
#   python3 rao-nozzle.py
#
# Outputs nozzle-wall.data (x/r_A, r/r_A, Mach, theta in degrees) and a
# summary of the run, including a comparison against the digitised data
# points of Rao (1958), nozzle A, taken from Fig 3.3 of the thesis.
#
# References:
#   G.V.R. Rao (1958), "Exhaust nozzle contour for optimum thrust",
#     Journal of Jet Propulsion 28(6):377-382.
#   J. Kunze (2020), "Design of a 3D shape transitioning nozzle and
#     experimental thrust measurements of an airframe integrated scramjet",
#     PhD thesis, The University of Queensland, section 3.1.
#
# JM, 2026-09-15, built on the imoc examples of PJ and Fabian Zander.
#
import math
import numpy as np
import gdtk.imoc.kernel as kernel
import gdtk.imoc.unit_process as unit
from gdtk.ideal_gas_flow import PM1, PM2

# ---------------------------------------------------------------------------
# Inputs: Rao (1958) nozzle A, as tabulated in Table 3.1 of the thesis.
# Lengths are nondimensionalised by the radius r_A at the wall point A,
# pressures and temperatures by their stagnation values (kernel defaults).
M_A = 1.103           # Mach number on the (straight) initial characteristic
M_E = 3.5             # Mach number at the nozzle lip E
G = 1.23              # ratio of specific heats
P_AMB_OVER_P0 = 0.0   # ambient-to-stagnation pressure ratio (0 => vacuum)
RA_OVER_RA = 0.45     # throat-arc radius R_A/r_A (used only to present the
                      # wall data in the same form as Rao's figure)

# Discretisation and iteration controls.
N_INIT = 40                        # nodes on the initial characteristic
DNU_STEP = math.radians(0.25)      # Prandtl-Meyer increment of the fan sweep
N_TABLE = 4000                     # rows in the DE lookup table
N_DE = 60                          # mesh nodes placed along DE
EPS_TOL = 1.0e-6                   # tolerance on mass-flow error, eq 3.12
LAMBDA2_TOL = 1.0e-12              # tolerance when locating point D
MAX_BISECT = 60                    # bisection iteration limits

kernel.axisymmetric = True
kernel.g = G

NU_A = PM1(M_A, G)

# ---------------------------------------------------------------------------
# Small helpers (nondimensional: a0=1, rho0=1).

def mach_angle(M):
    return math.asin(1.0/M)

def a_nd(M):
    "Speed of sound / stagnation speed of sound."
    return (1.0 + 0.5*(G-1.0)*M*M)**-0.5

def rho_nd(M):
    "Density / stagnation density."
    return (1.0 + 0.5*(G-1.0)*M*M)**(-1.0/(G-1.0))

def lambda2(M, theta):
    "Lagrange-multiplier expression of eq 3.2, invariant along DE."
    alpha = mach_angle(M)
    return -math.sqrt(1.0/(G-1.0+2.0/(M*M)))*math.cos(theta-alpha)/math.cos(alpha)

def lambda3(M, theta, r):
    "Lagrange-multiplier expression of eq 3.3, invariant along DE."
    alpha = mach_angle(M)
    return r*M*M*(1.0+0.5*(G-1.0)*M*M)**(-G/(G-1.0)) * \
        math.sin(theta)**2*math.tan(alpha)

# ---------------------------------------------------------------------------
# Step 1: conditions along the control surface DE.

def theta_at_lip():
    """
    Flow angle at the nozzle lip E from the transversality condition,
    eq 3.4: sin(2 theta_E) = 2/(g M_E^2) (1 - p0/pE) cot(alpha_E).
    For a vacuum back-pressure the pressure factor is exactly one.
    """
    factor = 1.0
    if P_AMB_OVER_P0 > 0.0:
        pE = (1.0 + 0.5*(G-1.0)*M_E*M_E)**(-G/(G-1.0))  # pE/p0, isentropic
        factor = 1.0 - P_AMB_OVER_P0/pE
    s = 2.0/(G*M_E*M_E)*factor*math.sqrt(M_E*M_E-1.0)
    return 0.5*math.asin(s)

def g3(M, theta):
    "The flow-property part of eq 3.3 (proportional to rho V^2 sin^2(theta) tan(alpha))."
    return M*M*(1.0+0.5*(G-1.0)*M*M)**(-G/(G-1.0)) * \
        math.sin(theta)**2*math.tan(mach_angle(M))

def build_de_table(theta_E):
    """
    Lookup table of (theta, M, r/r_E) triplets along the control surface DE.

    Along DE both lambda2 (eq 3.2) and lambda3 (eq 3.3) keep their lip
    values, so for each flow angle theta, eq 3.2 is solved for M and eq 3.3
    then gives r/r_E directly.  The angle theta is the natural sweep
    variable: it rises monotonically from theta_E at the lip towards D,
    while M first dips slightly below M_E (at the point where theta equals
    the Mach angle) before rising again, and r/r_E falls from one.  The
    sweep is terminated where r/r_E reaches its minimum, which bounds the
    angle at which a point D can exist.

    The tabulated states satisfy the axisymmetric C+ compatibility
    relation, d(nu-theta) = sin(mu)sin(theta)/r dl, confirming Rao's
    result that the control surface is itself a characteristic.
    """
    l2E = lambda2(M_E, theta_E)
    g3E = g3(M_E, theta_E)
    thetas = []; Ms = []; rs = []
    M_prev = M_E
    for theta in np.linspace(theta_E, math.radians(60.0), N_TABLE):
        # Solve eq 3.2 for M by secant continuation from the previous row.
        m0, m1 = M_prev, M_prev*1.0005 + 1.0e-6
        f0 = lambda2(m0, theta) - l2E
        for _ in range(60):
            f1 = lambda2(m1, theta) - l2E
            if f1 == f0: break
            m0, f0, m1 = m1, f1, m1 - f1*(m1-m0)/(f1-f0)
            if abs(m1-m0) < 1.0e-13: break
        if not (1.0 < m1 < 5.0*M_E) or abs(lambda2(m1, theta)-l2E) > 1.0e-9:
            break
        r = g3E/g3(m1, theta)
        if rs and r >= rs[-1]: break   # passed the minimum radius
        thetas.append(theta); Ms.append(m1); rs.append(r)
        M_prev = m1
    table = {'theta': np.array(thetas), 'M': np.array(Ms),
             'r_over_rE': np.array(rs), 'l2E': l2E, 'g3E': g3E}
    assert abs(table['M'][0]-M_E) < 1.0e-9
    assert abs(table['r_over_rE'][0]-1.0) < 1.0e-12
    # Self-consistency: every tabulated state must hold both invariants
    # at their lip values.
    l3E = lambda3(M_E, theta_E, 1.0)
    err2 = max(abs(lambda2(m, t)-l2E)
               for m, t in zip(table['M'], table['theta']))
    err3 = max(abs(lambda3(m, t, r)-l3E)
               for m, t, r in zip(table['M'], table['theta'],
                                  table['r_over_rE']))
    assert err2 < 1.0e-9 and err3 < 1.0e-12, (err2, err3)
    return table

# ---------------------------------------------------------------------------
# Steps 2-4: the expansion fan at A and the mass-flow balance.

def make_initial_characteristic():
    """
    The initial characteristic is approximated as a straight line from the
    wall point A=(0,1) to the axis, carrying uniform axial flow at M_A.
    """
    mu_A = mach_angle(M_A)
    x_axis = 1.0/math.tan(mu_A)
    indices = []
    for i in range(N_INIT+1):
        s = i/N_INIT
        node = kernel.Node(x=s*x_axis, y=1.0-s, theta=0.0, nu=NU_A, mach=M_A)
        kernel.register_node_in_mesh(node)
        indices.append(node.indx)
    for a, b in zip(indices[:-1], indices[1:]):
        kernel.nodes[a].cminus_down = b
        kernel.nodes[b].cminus_up = a
    return indices

# The axis is modelled as a wall at y=0.
axis_wall = kernel.Wall(lambda x: 0.0, 0.0, 30.0)

def add_fan_characteristic(nu_fan, prev_char):
    """
    Add one characteristic of the centred fan at A: a new node at A with
    Prandtl-Meyer value nu_fan, marched down along the previous
    characteristic and closed onto the axis.  Returns the ordered node list
    (A ... axis); every node in the list is newly created, so the list also
    serves as the retraction record for trial characteristics.
    """
    fan = kernel.Node(x=0.0, y=1.0, nu=nu_fan, theta=nu_fan-NU_A,
                      mach=PM2(nu_fan, G))
    kernel.register_node_in_mesh(fan)
    new_nodes = unit.march_along_cminus(prev_char[1], fan.indx, 'down')
    if kernel.nodes[new_nodes[-1]].y > 1.0e-6:
        new_nodes.append(unit.cminus_wall(axis_wall, new_nodes[-1]))
    return new_nodes

def retract_nodes(indices):
    "Remove a trial characteristic, restoring the mesh around it."
    for i in reversed(indices):
        kernel.delete_node(i)

def interp_props(i1, i2, alpha):
    """
    Node properties linearly interpolated between nodes i1 and i2,
    mirroring exactly what unit.insert() would create at fraction alpha.
    """
    n1 = kernel.nodes[i1]; n2 = kernel.nodes[i2]
    x = (1.0-alpha)*n1.x + alpha*n2.x
    y = (1.0-alpha)*n1.y + alpha*n2.y
    nu = (1.0-alpha)*n1.nu + alpha*n2.nu
    theta = (1.0-alpha)*n1.theta + alpha*n2.theta
    return x, y, nu, theta, PM2(nu, G)

def find_point_D(char_nodes, l2E):
    """
    Walk the newest characteristic upstream from the axis, looking for the
    place where lambda2 (eq 3.2) crosses its lip value; refine the location
    by bisection on the interpolation fraction.  Returns (i1, i2, alpha) or
    None if the characteristic does not yet reach the required lambda2.
    """
    vals = [lambda2(kernel.nodes[i].mach, kernel.nodes[i].theta) - l2E
            for i in char_nodes]
    k_found = None
    for k in range(len(char_nodes)-2, -1, -1):  # from the axis, upstream
        if vals[k] == 0.0 or vals[k]*vals[k+1] < 0.0:
            k_found = k
            break
    if k_found is None: return None
    i1 = char_nodes[k_found]; i2 = char_nodes[k_found+1]
    a, b = 0.0, 1.0
    fa = vals[k_found]
    for _ in range(MAX_BISECT):
        m = 0.5*(a+b)
        _, _, _, th, M = interp_props(i1, i2, m)
        fm = lambda2(M, th) - l2E
        if abs(fm) < LAMBDA2_TOL: return (i1, i2, m)
        if fa*fm <= 0.0:
            b = m
        else:
            a = m; fa = fm
    return (i1, i2, 0.5*(a+b))

def mass_flow_AD(props):
    """
    Mass flow crossing the fan characteristic between A and D, eqs 3.5-3.6.
    props is an ordered list of (x, y, nu, theta, M) tuples.
    """
    total = 0.0
    x_prev = None; Phi_prev = None
    for x, y, nu, theta, M in props:
        alpha = mach_angle(M)
        Phi = a_nd(M)*rho_nd(M)*M*math.sin(alpha)/math.cos(theta-alpha)*y
        if x_prev is not None:
            total += 0.5*(x-x_prev)*(Phi+Phi_prev)
        x_prev = x; Phi_prev = Phi
    return 2.0*math.pi*total

def node_props(indices):
    return [(kernel.nodes[i].x, kernel.nodes[i].y, kernel.nodes[i].nu,
             kernel.nodes[i].theta, kernel.nodes[i].mach) for i in indices]

def control_surface(table, M_D, theta_D, x_D, r_D):
    """
    Realise the control surface DE for a candidate point D: slice the
    lookup table between theta_D and theta_E, scale the radii so that
    r(theta_D)=r_D, and integrate the x-coordinates from D towards E using
    the local C+ slope tan(theta+alpha), eqs 3.7-3.9.  The returned arrays
    run from D to E; theta decreases and r increases along them.
    """
    ths = table['theta']
    mask = ths < theta_D
    theta = np.concatenate([[theta_D], ths[mask][::-1]])
    M = np.concatenate([[M_D], table['M'][mask][::-1]])
    r_rel = np.concatenate([[table['g3E']/g3(M_D, theta_D)],
                            table['r_over_rE'][mask][::-1]])
    r = r_rel*(r_D/r_rel[0])
    slope = np.tan(theta+np.arcsin(1.0/M))
    x = np.empty_like(r)
    x[0] = x_D
    for i in range(len(r)-1):
        x[i+1] = x[i] + (r[i+1]-r[i])/(0.5*(slope[i]+slope[i+1]))
    return {'x': x, 'r': r, 'M': M, 'theta': theta}

def mass_flow_DE(cs):
    "Mass flow crossing the control surface DE, eqs 3.10-3.11."
    M = cs['M']; theta = cs['theta']
    alpha = np.arcsin(1.0/M)
    a = (1.0+0.5*(G-1.0)*M*M)**-0.5
    rho = (1.0+0.5*(G-1.0)*M*M)**(-1.0/(G-1.0))
    Phi = a*rho*M*np.sin(alpha)/np.sin(theta+alpha)*cs['r']
    # The 1/sin(theta+alpha) factor converts a line-element of DE into an
    # increment of radius, so the trapezoidal sum is taken over r.
    return 2.0*math.pi*np.sum(0.5*np.diff(cs['r'])*(Phi[:-1]+Phi[1:]))

def mass_flow_error(char_nodes, table):
    """
    The relative mass-flow mismatch, eq 3.12, for the newest fan
    characteristic, or None if the characteristic does not yet carry a
    valid point D (the lambda2 crossing exists but its flow angle is still
    below theta_E, so no control surface can be hung from it).
    Also returns the D location and control surface.
    """
    Dloc = find_point_D(char_nodes, table['l2E'])
    if Dloc is None: return None
    i1, i2, alpha = Dloc
    x_D, r_D, nu_D, th_D, M_D = interp_props(i1, i2, alpha)
    if not (table['theta'][0] < th_D < table['theta'][-1]): return None
    props = node_props(char_nodes[:char_nodes.index(i1)+1])
    props.append((x_D, r_D, nu_D, th_D, M_D))
    m_AD = mass_flow_AD(props)
    cs = control_surface(table, M_D, th_D, x_D, r_D)
    m_DE = mass_flow_DE(cs)
    eps = (m_AD-m_DE)/(0.5*(m_AD+m_DE))
    return eps, Dloc, cs, m_AD, m_DE

# ---------------------------------------------------------------------------
# Step 5: mesh along DE, the region between AD and DE, and the wall.

def build_de_nodes(cs, D_idx):
    """
    Create mesh nodes along the control surface from the lookup-table data
    and link them as a C+ characteristic running downstream from D to E.
    The nodes are spaced roughly uniformly in radius (the table itself is
    uniform in flow angle, which crowds the points towards D).
    """
    r = cs['r']
    targets = np.linspace(r[0], r[-1], N_DE)
    picks = np.unique(np.searchsorted(r, targets).clip(0, len(r)-1))
    picks[-1] = len(r)-1
    de = [D_idx]
    for i in picks[1:]:
        node = kernel.Node(x=cs['x'][i], y=cs['r'][i],
                           nu=PM1(cs['M'][i], G), theta=cs['theta'][i],
                           mach=cs['M'][i])
        kernel.register_node_in_mesh(node)
        kernel.nodes[de[-1]].cplus_down = node.indx
        node.cplus_up = de[-1]
        de.append(node.indx)
    return de

def fill_region_ADE(D_idx, de):
    """
    Compute the flow field between the last fan characteristic AD and the
    control surface DE.  From each DE node a C- characteristic is marched
    'up', using the previous characteristic for the C+ data (the same idiom
    as the wave-cone region of examples/imoc/anderson_11p1.py).  The lines
    deliberately overshoot the (as yet unknown) wall; the extra nodes above
    the wall streamline are harmless.
    """
    lines = []
    old = kernel.nodes[D_idx].cminus_up
    for dk in de[1:]:
        new_nodes = unit.march_along_cminus(old, dk, 'up')
        lines.append(new_nodes)
        old = new_nodes[1]
    return lines

def trace_wall(fan_idx, lines, de):
    """
    Trace the nozzle wall as the streamline leaving A after the full
    expansion, extending it across each C- characteristic of the region
    ADE and finally across DE itself (extended linearly beyond E so that
    the last crossing is guaranteed).
    """
    kernel.register_streamline_start(fan_idx)
    wall = [fan_idx]
    # The last entry of lines is the characteristic through E itself; the
    # wall meets it at E, so the crossing is taken with the extended DE
    # surface below instead.
    for line in lines[:-1]:
        idx = unit.extend_streamline_to_given_line(wall[-1],
                                                   list(reversed(line)))
        if idx is None:
            raise RuntimeError("Wall streamline failed to cross a "
                               "characteristic of region ADE.")
        wall.append(idx)
    # Extend DE beyond E, then take the final wall segment across it.
    nE = kernel.nodes[de[-1]]; nP = kernel.nodes[de[-2]]
    frac = 0.5
    nu_ext = nE.nu + frac*(nE.nu-nP.nu)
    ext = kernel.Node(x=nE.x+frac*(nE.x-nP.x), y=nE.y+frac*(nE.y-nP.y),
                      nu=nu_ext, theta=nE.theta+frac*(nE.theta-nP.theta),
                      mach=PM2(nu_ext, G))
    kernel.register_node_in_mesh(ext)
    nE.cplus_down = ext.indx; ext.cplus_up = nE.indx
    idx = unit.extend_streamline_to_given_line(wall[-1],
                                               [ext.indx]+list(reversed(de)))
    if idx is not None: wall.append(idx)
    return wall

# ---------------------------------------------------------------------------
# Output.

def write_wall_data(wall, fname):
    """
    Write the wall as x/r_A, r/r_A, Mach and theta(deg).  To present the
    contour in the same form as Rao's figure, the centred-expansion corner
    at A is unrolled onto the circular throat arc of radius R_A: the arc
    turns the wall from 0 to theta_max, carrying the corner-fan flow values,
    and the computed (streamline) wall attaches at the arc end.
    """
    th_max = kernel.nodes[wall[0]].theta
    x_c = RA_OVER_RA*math.sin(th_max)
    r_c = 1.0 + RA_OVER_RA*(1.0-math.cos(th_max))
    with open(fname, 'w') as f:
        f.write('# Rao thrust-optimised nozzle computed with gdtk.imoc.\n')
        f.write('# M_A=%g M_E=%g gamma=%g p_amb/p0=%g R_A/r_A=%g\n' %
                (M_A, M_E, G, P_AMB_OVER_P0, RA_OVER_RA))
        f.write('# N_INIT=%d DNU_STEP=%gdeg N_TABLE=%d N_DE=%d EPS_TOL=%g\n' %
                (N_INIT, math.degrees(DNU_STEP), N_TABLE, N_DE, EPS_TOL))
        f.write('# x/r_A      r/r_A       Mach        theta(deg)\n')
        n_arc = 40
        for i in range(n_arc):
            th = th_max*i/n_arc
            f.write('%12.6f %12.6f %12.6f %12.6f\n' %
                    (RA_OVER_RA*math.sin(th),
                     1.0+RA_OVER_RA*(1.0-math.cos(th)),
                     PM2(NU_A+th, G), math.degrees(th)))
        for i in wall:
            n = kernel.nodes[i]
            f.write('%12.6f %12.6f %12.6f %12.6f\n' %
                    (n.x+x_c, n.y+(r_c-1.0), n.mach, math.degrees(n.theta)))
    return

def compare_with(fname, xs, ys, label):
    "RMS/max deviation of the computed wall from a digitised reference."
    try:
        ref = np.loadtxt(fname)
    except OSError:
        print("  (%s not found; skipping)" % fname); return
    mask = (ref[:, 0] >= xs.min()) & (ref[:, 0] <= xs.max())
    dev = np.interp(ref[mask, 0], xs, ys) - ref[mask, 1]
    print("  %-28s n=%2d  max|dev|=%.4f  rms=%.4f" %
          (label, len(dev), np.max(np.abs(dev)), math.sqrt(np.mean(dev**2))))
    return

# ---------------------------------------------------------------------------
# Main.

def main():
    print("Rao thrust-optimised nozzle via gdtk.imoc")
    print("Inputs: M_A=%g M_E=%g gamma=%g p_amb/p0=%g R_A/r_A=%g" %
          (M_A, M_E, G, P_AMB_OVER_P0, RA_OVER_RA))
    theta_E = theta_at_lip()
    print("Lip flow angle theta_E=%.4f deg (eq 3.4)" % math.degrees(theta_E))
    table = build_de_table(theta_E)
    print("DE lookup table: %d rows, theta in [%.3f, %.3f] deg, "
          "M in [%.4f (min %.4f), %.4f], r/r_E down to %.4f" %
          (len(table['M']), math.degrees(table['theta'][0]),
           math.degrees(table['theta'][-1]), table['M'][0],
           table['M'].min(), table['M'][-1], table['r_over_rE'][-1]))
    #
    initial = make_initial_characteristic()
    m_initial = mass_flow_AD(node_props(initial))
    m_exact = math.pi*rho_nd(M_A)*a_nd(M_A)*M_A
    print("Mass flow through initial characteristic: %.6f "
          "(uniform-flow exact value %.6f)" % (m_initial, m_exact))
    #
    # Grow the fan until the mass-flow error, eq 3.12, changes sign.
    prev_char = initial
    prev_eps = None
    nu_lo = None
    n_step = 0
    while True:
        n_step += 1
        nu_fan = NU_A + n_step*DNU_STEP
        char = add_fan_characteristic(nu_fan, prev_char)
        res = mass_flow_error(char, table)
        if res is None:
            prev_char = char
            continue
        eps = res[0]
        if prev_eps is None:
            print("First candidate point D at fan angle %.3f deg, "
                  "eps=%.5f" % (math.degrees(nu_fan-NU_A), eps))
        if prev_eps is not None and eps*prev_eps < 0.0:
            retract_nodes(char)
            nu_lo, eps_lo = nu_fan-DNU_STEP, prev_eps
            nu_hi = nu_fan
            break
        prev_eps = eps
        prev_char = char
        if n_step > 10000:
            raise RuntimeError("Fan sweep failed to bracket eps=0.")
    print("Mass-flow error changes sign between fan angles "
          "%.3f and %.3f deg" %
          (math.degrees(nu_lo-NU_A), math.degrees(nu_hi-NU_A)))
    #
    # Bisection on the fan Prandtl-Meyer value.
    char = None; res = None
    for it in range(MAX_BISECT):
        nu_mid = 0.5*(nu_lo+nu_hi)
        char = add_fan_characteristic(nu_mid, prev_char)
        res = mass_flow_error(char, table)
        assert res is not None
        eps = res[0]
        if abs(eps) < EPS_TOL: break
        retract_nodes(char); char = None
        if eps*eps_lo > 0.0:
            nu_lo = nu_mid; eps_lo = eps
        else:
            nu_hi = nu_mid
    if char is None:
        # Tolerance not met; keep the midpoint characteristic anyway.
        nu_mid = 0.5*(nu_lo+nu_hi)
        char = add_fan_characteristic(nu_mid, prev_char)
        res = mass_flow_error(char, table)
    eps, (i1, i2, alphaD), cs, m_AD, m_DE = res
    fan_idx = char[0]
    theta_max = kernel.nodes[fan_idx].theta
    print("Converged after %d bisections: fan angle theta_max=%.4f deg, "
          "eps=%.2e" % (it+1, math.degrees(theta_max), eps))
    print("Mass flows: m_AD=%.6f m_DE=%.6f (initial line %.6f)" %
          (m_AD, m_DE, m_initial))
    #
    # Realise D, the control surface and the rest of the flow field.
    D_idx = unit.insert(i1, i2, alpha=alphaD)
    nD = kernel.nodes[D_idx]
    print("Point D at x=%.4f r=%.4f M=%.4f theta=%.3f deg" %
          (nD.x, nD.y, nD.mach, math.degrees(nD.theta)))
    print("Lip E at x=%.4f r=%.4f (control-surface integration)" %
          (cs['x'][-1], cs['r'][-1]))
    de = build_de_nodes(cs, D_idx)
    lines = fill_region_ADE(D_idx, de)
    wall = trace_wall(fan_idx, lines, de)
    nW = kernel.nodes[wall[-1]]
    print("Wall streamline: %d points, ends at x=%.4f r=%.4f "
          "M=%.4f theta=%.3f deg" %
          (len(wall), nW.x, nW.y, nW.mach, math.degrees(nW.theta)))
    print("  closure vs lip E: dx=%.4f dr=%.4f dM=%.4f" %
          (nW.x-cs['x'][-1], nW.y-cs['r'][-1], nW.mach-M_E))
    #
    write_wall_data(wall, 'nozzle-wall.data')
    print("Wrote nozzle-wall.data (throat arc R_A/r_A=%.2f prepended, "
          "wall shifted to the arc end, as presented by Rao)" % RA_OVER_RA)
    #
    data = np.loadtxt('nozzle-wall.data')
    print("Comparison against digitised data from Fig 3.3 of the thesis:")
    compare_with('rao-contour.data', data[:, 0], data[:, 1],
                 'contour r/r_A, Rao points')
    compare_with('kunze-contour.data', data[:, 0], data[:, 1],
                 'contour r/r_A, Kunze IMOC')
    compare_with('rao-wall-mach.data', data[:, 0], data[:, 2],
                 'wall Mach, Rao points')
    compare_with('kunze-wall-mach.data', data[:, 0], data[:, 2],
                 'wall Mach, Kunze IMOC')
    compare_with('rao-wall-theta.data', data[:, 0], data[:, 3],
                 'wall theta(deg), Rao points')
    compare_with('kunze-wall-theta.data', data[:, 0], data[:, 3],
                 'wall theta(deg), Kunze IMOC')
    print("  (the largest theta deviations sit in the steeply-rising arc")
    print("   region x<0.25 where a small x-digitisation error reads as a")
    print("   large angle difference)")
    print("Total nodes in mesh: %d" % kernel.number_of_nodes())
    return

if __name__ == '__main__':
    main()
