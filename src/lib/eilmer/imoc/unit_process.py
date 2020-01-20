# unit_process.py
"""
Unit processes for making new nodes in the Isentropic Method-Of-Characteristics.
This is a regrowth of the old IMOC code that was implemented in C.

Author(s):
  Peter J.
  Centre for Hypersonics,
  School of Mechanical Engineering, U of Q.

Version:
  2020-01-09: Let's start coding in Python and see how it develops...
"""

import eilmer.imoc.kernel as kernel
import eilmer.ideal_gas_flow as igf
from math import sin, cos, sqrt, asin, atan
from eilmer.zero_solvers import secant as solve

# Tolerances for convergence checks.
max_iteration = 15
position_tolerance = 1.0e-5

def interior(node1, node2, node4):
    """
    Returns an interior point computed from two initial points.

    node1: index of initial point along C- characteristic
    node2: index of initial point along C+ characteristic
    node4: index of solution point (may have a value of -1)
    If -1 is specified as the index for node4, a new node will be
    created for the solution point.
    """
    if node1 == node2:
        raise RuntimeError("Same index given for node1 and node2.")
    n1 = kernel.nodes[node1]
    n2 = kernel.nodes[node2]
    x1 = n1.x; y1 = n1.y; pm1 = n1.nu; th1 = n1.theta; m1 = n1.mach
    x2 = n2.x; y2 = n2.y; pm2 = n2.nu; th2 = n2.theta; m2 = n2.mach
    # Mach angles
    mu1 = asin(1/m1)
    mu2 = asin(1/m2)
    # Use point 2 to get the streamline direction cosines.
    xStream = cos(th2); yStream = sin(th2)
    # Guess at some of the solution point properties.
    # The position will be way off but it is used as part of
    # the convergence check a little further on.
    x4 = 0.5*(x1+x2); y4 = 0.5*(y1+y2)
    th4 = 0.5*(th1+th2)
    mu4 = 0.5*(mu1+mu2)
    # Compute the solution point position and flow properties.
    converged = False
    iteration_count =0
    while (not converged) and (iteration_count < max_iteration):
        x4_old = x4; y4_old = y4
        # Locate solution point by assuming straight-line segments.
        sinCminus = 0.5*(sin(th1-mu1) + sin(th4-mu4))
        cosCminus = 0.5*(cos(th1-mu1) + cos(th4-mu4))
        sinCplus  = 0.5*(sin(th2+mu2) + sin(th4+mu4))
        cosCplus  = 0.5*(cos(th2+mu2) + cos(th4+mu4))
        #
        numerator = (x2-x1)*sinCplus - (y2-y1)*cosCplus;
        denominator = cosCminus*sinCplus - sinCminus*cosCplus;
        if abs(denominator) <= 1.0e-12:
            raise RuntimeError("Interior: characteristics are parallel.")
        lambdaCminus = numerator/denominator
        x4 = x1 + lambdaCminus*cosCminus
        y4 = y1 + lambdaCminus*sinCminus
        dx = x4-x4_old; dy = y4-y4_old
        change_in_position = sqrt(dx*dx + dy*dy)
        #
        # Lengths of the characteristic segments.
        dx = x4-x1; dy = y4-y1
        lengthCminus = sqrt(dx*dx + dy*dy)
        dot_product = dx*xStream + dy*yStream
        if dot_product < 0.0:
            directionCminus = -1
        else:
            directionCminus = +1
        #
        dx = x4-x2; dy = y4-y2
        lengthCplus  = sqrt(dx*dx + dy*dy)
        dot_product = dx*xStream + dy*yStream
        if dot_product < 0.0:
            directionCplus = -1
        else:
            directionCplus = +1
        #
        # Update flow properties at solution point
        # First, assume 2D planar geometry then add
        # axisymmetric contributions if flag is set.
        pm4 = 0.5*(pm1+pm2) + 0.5*(th1-th2)
        th4 = 0.5*(pm1-pm2) + 0.5*(th1+th2)
        if kernel.axisymmetric:
            if y1 < 1.0e-6 and y2 < 1.0e-6:
                raise RuntimeError("Interior: both nodes are too close to axis.")
            # Axisymmetric components.
            if y1 < 1.0e-6:
                axiterm1 = sin(mu2)*sin(th2)/y2
            else:
                axiterm1 = sin(mu1)*sin(th1)/y1
            if y2 < 1.0e-6:
                axiterm2 = sin(mu1)*sin(th1)/y1
            else:
                axiterm2 = sin(mu2)*sin(th2)/y2
            #
            axiterm4 = sin(mu4)*sin(th4)/y4
            integralCminus = 0.5*directionCminus*lengthCminus*(axiterm1+axiterm4)
            integralCplus  = 0.5*directionCplus*lengthCplus*(axiterm2+axiterm4)
            # Include axisymmetric components.
            pm4 += 0.5*(integralCminus+integralCplus)
            th4 += 0.5*(integralCminus-integralCplus)
        #
        iteration_count += 1
        converged = change_in_position < position_tolerance
    # Save the solution-point properties and connect the 
    # node into the characteristic mesh.
    if node4 == -1:
        n4 = kernel.Node()
        node4 = n4.indx
    else:
        n4 = nodes[node4]
    m4 = igf.PM2(pm4, kernel.g)
    n4.x = x4; n4.y = y4; n4.nu = pm4; n4.theta = th4; n4.mach = m4
    # We assume that the principal flow direction
    # is in the positive x-direction.
    if x4 > x2:
        n4.cplus_up = node2; n2.cplus_down = node4
    else:
        n4.cplus_down = node2; n2.cplus_up = node4
    if x4 > x1:
        n4.cminus_up = node1; n1.cminus_down = node4
    else:
        n4.cminus_down = node1; n1.cminus_up = node4
    return n4


def insert(node1, node2, node4, alpha):
    """
    Returns the node inserted between two initial points.

    node1: index of initial point
    node2: index of initial point
    node4: index of inserted point (may have a value of -1)
    alpha: fraction that node4 is like node2;
           n4.value = alpha n2.value + (1-alpha) n1.value
    If -1 is specified as the index for node4, a new node will be
    created for the inserted point.
    If node1 and node2 are adjacent nodes along a characteristic line,
    node4 will be connected in between.
    """
    if node1 == node2:
        raise RuntimeError("Same index given for node1 and node2.")
    n1 = kernel.nodes[node1]
    n2 = kernel.nodes[node2]
    if node4 == -1:
        n4 = kernel.Node()
        node4 = n4.indx
    else:
        n4 = nodes[node4]
    # Enforce a 0.0..1.0 range for alpha
    alpha = max(min(alpha, 1.0), 0.0)
    # Linearly interpolate node properties.
    n4.x = (1-alpha)*n1.x + alpha*n2.x
    n4.y = (1-alpha)*n1.y + alpha*n2.y
    n4.nu = (1-alpha)*n1.nu + alpha*n2.nu
    n4.theta = (1-alpha)*n1.theta + alpha*n2.theta
    n4.mach = igf.PM2(n4.nu, kernel.g)
    # Connect into the mesh only if nodes 1 and 2 are adjacent.
    print("node1=", node1, "n1.cminus_down=", n1.cminus_down,
          "node2=", node2, "n2.cminus_up=", n2.cminus_up)
    if (n1.cplus_down == node2) and (n2.cplus_up == node1) and (n1.cplus_down != -1):
        n4.cplus_up = node1; n1.cplus_down = node4
        n2.cplus_up = node4; n4.cplus_down = node2
    elif (n1.cplus_up == node2) and (n2.cplus_down == node1) and (n1.cplus_up != -1):
        n4.cplus_up = node2; n2.cplus_down = node4
        n1.cplus_up = node4; n4.cplus_down = node1
    elif (n1.cminus_down == node2) and (n2.cminus_up == node1) and (n1.cminus_down != -1):
        n4.cminus_up = node1; n1.cminus_down = node4
        n2.cminus_up = node4; n4.cminus_down = node2
    elif (n1.cminus_up == node2) and (n2.cminus_down == node1) and (n1.cminus_up != -1):
        n4.cminus_up = node2; n2.cminus_down = node4
        n1.cminus_up = node4; n4.cminus_down = node1
    elif (n1.czero_down == node2) and (n2.czero_up == node1) and (n1.czero_down != -1):
        n4.czero_up = node1; n1.czero_down = node4
        n2.czero_up = node4; n4.czero_down = node2
    elif (n1.czero_up == node2) and (n2.czero_down == node1) and (n1.czero_up != -1):
        n4.czero_up = node2; n2.czero_down = node4
        n1.czero_up = node4; n4.czero_down  = node1
    # Assuming success...
    return n4


def wall_position(f, x0, y0, cosx, cosy):
    """
    Returns the x,y coordinates of the point on the wall y=f(x).

    f: function y=f(x) that defines the wall.
    x0, y0: starting point of intersecting line segment
    cosx, cosy: direction cosines for the line segment
    """
    def err(dl):
        x = x0 + dl*cosx; y = y0 + dl*cosy
        return y-f(x)
    dL = solve(err, 1.0, 1.1)
    return x0+dL*cosx, y0+dL*cosy


def wall_slope(f, x, dx=1.0e-6):
    """
    Returns ths slope of the point on the wall y=f(x).

    f: function y=f(x) that defines the wall.
    x: x-loacation of point on the wall
    dx: increment for the finite-difference calculation
    """
    return (f(x+dx)-f(x))/dx


def cminus_wall(fn_wall, node1, node4):
    """
    Returns a point on the wall, computed from one initial point.

    fn_wall: user-supplied function defining wall y=f(x)
    node1: index of initial point along C- characteristic
    node4: index of solution point (may have a value of -1)
    If -1 is specified as the index for node4, a new node will be
    created for the solution point.
    """
    if not callable(fn_wall):
        raise RuntimeError("cminus_wall expects a callable function for fn_wall.")
    n1 = kernel.nodes[node1]
    x1 = n1.x; y1 = n1.y; pm1 = n1.nu; th1 = n1.theta; m1 = n1.mach
    # Mach angles
    mu1 = asin(1/m1)
    # Use point 1 to get the streamline direction cosines.
    xStream = cos(th1); yStream = sin(th1)
    # Guess at the solution point properties.
    # The position may be way off but it is used as part of
    # the convergence check a little further on.
    x4 = x1; y4 = y1
    th4 = th1; mu4 = mu1
    # Compute the solution point position and flow properties.
    converged = False
    iteration_count =0
    while (not converged) and (iteration_count < max_iteration):
        x4_old = x4; y4_old = y4
        # Locate solution point by assuming straight-line segments.
        sinCminus = 0.5*(sin(th1-mu1) + sin(th4-mu4))
        cosCminus = 0.5*(cos(th1-mu1) + cos(th4-mu4))
        #
        x4, y4 = wall_position(fn_wall, x1, y1, cosCminus, sinCminus)
        dx = x4-x4_old; dy = y4-y4_old
        change_in_position = sqrt(dx*dx + dy*dy)
        #
        # Lengths of the characteristic segment.
        dx = x4-x1; dy = y4-y1
        lengthCminus = sqrt(dx*dx + dy*dy)
        dot_product = dx*xStream + dy*yStream
        if dot_product < 0.0:
            directionCminus = -1
        else:
            directionCminus = +1
        #
        # Update flow properties at solution point
        # First, assume 2D planar geometry then add
        # axisymmetric contributions if flag is set.
        th4 = atan(wall_slope(fn_wall, x4))
        pm4 = pm1 - (th4-th1)
        if kernel.axisymmetric:
            if y1 < 1.0e-6:
                raise RuntimeError("cminus_wall: initial node is too close to axis.")
            # Axisymmetric components.
            axiterm1 = sin(mu1)*sin(th1)/y1
            if y4 < 1.0e-6:
                axiterm4 = axiterm1
            else:
                axiterm4 = sin(mu4)*sin(th4)/y4
            #
            integralCminus = 0.5*directionCminus*lengthCminus*(axiterm1+axiterm4)
            # Include axisymmetric components.
            pm4 += integralCminus
        #
        iteration_count += 1
        converged = change_in_position < position_tolerance
    # Save the solution-point properties and connect the 
    # node into the characteristic mesh.
    if node4 == -1:
        n4 = kernel.Node()
        node4 = n4.indx
    else:
        n4 = nodes[node4]
    m4 = igf.PM2(pm4, kernel.g)
    n4.x = x4; n4.y = y4; n4.nu = pm4; n4.theta = th4; n4.mach = m4
    # We assume that the principal flow direction
    # is in the positive x-direction.
    if x4 > x1:
        n4.cminus_up = node1; n1.cminus_down = node4
    else:
        n4.cminus_down = node1; n1.cminus_up = node4
    return n4


def cplus_wall(fn_wall, node2, node4):
    """
    Returns a point on the wall, computed from one initial point.

    fn_wall: user-supplied function defining wall y=f(x)
    node2: index of initial point along C+ characteristic
    node4: index of solution point (may have a value of -1)
    If -1 is specified as the index for node4, a new node will be
    created for the solution point.
    """
    if not callable(fn_wall):
        raise RuntimeError("cplus_wall expects a callable function for fn_wall.")
    n2 = kernel.nodes[node2]
    x2 = n2.x; y2 = n2.y; pm2 = n2.nu; th2 = n2.theta; m2 = n2.mach
    # Mach angles
    mu2 = asin(1/m2)
    # Use point 2 to get the streamline direction cosines.
    xStream = cos(th2); yStream = sin(th2)
    # Guess at the solution point properties.
    # The position may be way off but it is used as part of
    # the convergence check a little further on.
    x4 = x2; y4 = y2
    th4 = th2; mu4 = mu2
    # Compute the solution point position and flow properties.
    converged = False
    iteration_count =0
    while (not converged) and (iteration_count < max_iteration):
        x4_old = x4; y4_old = y4
        # Locate solution point by assuming straight-line segments.
        sinCplus = 0.5*(sin(th2+mu2) + sin(th4+mu4))
        cosCplus = 0.5*(cos(th2+mu2) + cos(th4+mu4))
        #
        x4, y4 = wall_position(fn_wall, x2, y2, cosCplus, sinCplus)
        dx = x4-x4_old; dy = y4-y4_old
        change_in_position = sqrt(dx*dx + dy*dy)
        #
        # Lengths of the characteristic segment.
        dx = x4-x2; dy = y4-y2
        lengthCminus = sqrt(dx*dx + dy*dy)
        dot_product = dx*xStream + dy*yStream
        if dot_product < 0.0:
            directionCplus = -1
        else:
            directionCplus = +1
        #
        # Update flow properties at solution point
        # First, assume 2D planar geometry then add
        # axisymmetric contributions if flag is set.
        th4 = atan(wall_slope(fn_wall, x4))
        pm4 = pm2 + (th4-th2)
        if kernel.axisymmetric:
            if y4 < 1.0e-6:
                raise RuntimeError("cplus_wall: new node is too close to axis.")
            # Axisymmetric components.
            axiterm4 = sin(mu4)*sin(th4)/y4
            if y2 < 1.0e-6:
                axiterm2 = axiterm4
            else:
                axiterm2 = sin(mu2)*sin(th2)/y2
            #
            integralCplus = 0.5*directionCplus*lengthCplus*(axiterm2+axiterm4)
            # Include axisymmetric components.
            pm4 += integralCplus
        #
        iteration_count += 1
        converged = change_in_position < position_tolerance
    # Save the solution-point properties and connect the 
    # node into the characteristic mesh.
    if node4 == -1:
        n4 = kernel.Node()
        node4 = n4.indx
    else:
        n4 = nodes[node4]
    m4 = igf.PM2(pm4, kernel.g)
    n4.x = x4; n4.y = y4; n4.nu = pm4; n4.theta = th4; n4.mach = m4
    # We assume that the principal flow direction
    # is in the positive x-direction.
    if x4 > x2:
        n4.cplus_up = node2; n2.cplus_down = node4
    else:
        n4.cplus_down = node2; n2.cplus_up = node4
    return n4
