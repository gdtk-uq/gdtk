# unit_process.py
"""
Unit processes for making new nodes in the Isentropic Method-Of-Characteristics.
This is a regrowth of the old IMOC code that was implemented in C.

Author(s):
  Peter J.
    Centre for Hypersonics,
    School of Mechanical Engineering, U of Q.
  F. Zander
    University of Southern Queensland.

Version:
  2020-01-09: Let's start coding in Python and see how it develops...
  2020-04-14: Fabian Zander has joined the effort and has ported
    functions needed for computing a nozzle profile.
"""

import gdtk.imoc.kernel as kernel
import gdtk.ideal_gas_flow as igf
from math import sin, cos, sqrt, asin, atan, pow
from gdtk.numeric.zero_solvers import secant as solve
import numpy as np


def theta_over_r(ma, mb, dx, g):
    """
    Estimate part of the axisymmetric source term from the
    Mach number change along the axis, ma to mb over distance dx.

    See PJ's workbook notes, 2020-04-16, page 62.
    """
    ex = (g+1)/(g-1)/2
    num = ma*pow((1+0.5*(g-1)*(mb**2)), ex)
    den = mb*pow((1+0.5*(g-1)*(ma**2)), ex)
    return (sqrt(num/den)-1)/dx


# Tolerances for convergence checks.
max_iteration = 15
position_tolerance = 1.0e-5


def interior(node1, node2, node4=None):
    """
    Returns the index of an interior point computed from two initial points.

    node1: index of initial point along C- characteristic
    node2: index of initial point along C+ characteristic
    node4: index of solution point (may have a value of -1 or None)
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
            raise RuntimeError("interior: characteristics are parallel.")
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
            # Axisymmetric components will need alternate evaluation at the axis.
            if y1 < 1.0e-6 and y2 < 1.0e-6:
                # Both nodes are close to the axis.
                # Use the axial variation in Mach number to estimate theta/r,
                # and use that for sin(th)/y.
                th_over_r = theta_over_r(m2, m1, x1-x2, kernel.g)
                axiterm1 = sin(mu1)*th_over_r
                axiterm2 = sin(mu2)*th_over_r
            else:
                if y1 < 1.0e-6:
                    axiterm1 = sin(mu1)*sin(th2)/y2
                else:
                    axiterm1 = sin(mu1)*sin(th1)/y1
                if y2 < 1.0e-6:
                    axiterm2 = sin(mu2)*sin(th1)/y1
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
    if (node4 is None) or (node4 == -1):
        n4 = kernel.Node()
        node4 = n4.indx
    else:
        n4 = kernel.nodes[node4]
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
    kernel.char_mesh.append(n4.indx)
    return n4.indx


def insert(node1, node2, node4=None, alpha=0.5):
    """
    Returns the index of a node (4) inserted between two existing points (1,2).

    node1: index of existing point 1
    node2: index of existing point 2
    node4: index of inserted point (may have a value of -1 or None)
    alpha: fraction that node4 is like node2;
           n4.value = alpha n2.value + (1-alpha) n1.value
    If -1 (or None) is specified as the index for node4, a new node will be
    created for the inserted node.
    If node1 and node2 are adjacent nodes along a characteristic line,
    node4 will be connected in between.
    """
    if node1 == node2:
        raise RuntimeError("Same index given for node1 and node2.")
    n1 = kernel.nodes[node1]
    n2 = kernel.nodes[node2]
    if (node4 is None) or (node4 == -1):
        n4 = kernel.Node()
        node4 = n4.indx
    else:
        n4 = kernel.nodes[node4]
    # Enforce a 0.0..1.0 range for alpha
    alpha = max(min(alpha, 1.0), 0.0)
    # Linearly interpolate node properties.
    n4.x = (1-alpha)*n1.x + alpha*n2.x
    n4.y = (1-alpha)*n1.y + alpha*n2.y
    n4.nu = (1-alpha)*n1.nu + alpha*n2.nu
    n4.theta = (1-alpha)*n1.theta + alpha*n2.theta
    n4.mach = igf.PM2(n4.nu, kernel.g)
    # Connect into the mesh only if nodes 1 and 2 are adjacent.
    # print("node1=", node1, "n1.cminus_down=", n1.cminus_down,
    #       "node2=", node2, "n2.cminus_up=", n2.cminus_up)
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
    kernel.char_mesh.append(n4.indx)
    return n4.indx


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


def cminus_wall(wall, node1, node4=None):
    """
    Returns the index of a node on the wall, computed from one initial node.

    fn_wall: user-supplied Wall object defining wall y=f(x)
    node1: index of initial point along C- characteristic
    node4: index of solution point (may have a value of -1 or None)
    If -1 (or None) is specified as the index for node4, a new node will be
    created for the solution point.
    """
    if wall.__class__ != kernel.Wall:
        raise RuntimeError("cminus_wall expects a Wall object.")
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
        x4, y4 = wall_position(wall, x1, y1, cosCminus, sinCminus)
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
        th4 = atan(wall.dfdx(x4))
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
    if (node4 is None) or (node4 == -1):
        n4 = kernel.Node()
        node4 = n4.indx
    else:
        n4 = kernel.nodes[node4]
    m4 = igf.PM2(pm4, kernel.g)
    n4.x = x4; n4.y = y4; n4.nu = pm4; n4.theta = th4; n4.mach = m4
    # We assume that the principal flow direction
    # is in the positive x-direction.
    if x4 > x1:
        n4.cminus_up = node1; n1.cminus_down = node4
    else:
        n4.cminus_down = node1; n1.cminus_up = node4
    kernel.char_mesh.append(n4.indx)
    return n4.indx


def cplus_wall(wall, node2, node4=None):
    """
    Returns the index of a node on the wall, computed from one initial node.

    wall: user-supplied Wall object defining wall y=f(x)
    node2: index of initial point along C+ characteristic
    node4: index of solution point (may have a value of -1 or None)
    If -1 is specified as the index for node4, a new node will be
    created for the solution point.
    """
    if wall.__class__ != kernel.Wall:
        raise RuntimeError("cplus_wall expects a Wall object.")
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
        x4, y4 = wall_position(wall, x2, y2, cosCplus, sinCplus)
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
        th4 = atan(wall.dfdx(x4))
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
    if (node4 is None) or (node4 == -1):
        n4 = kernel.Node()
        node4 = n4.indx
    else:
        n4 = kernel.nodes[node4]
    m4 = igf.PM2(pm4, kernel.g)
    n4.x = x4; n4.y = y4; n4.nu = pm4; n4.theta = th4; n4.mach = m4
    # We assume that the principal flow direction
    # is in the positive x-direction.
    if x4 > x2:
        n4.cplus_up = node2; n2.cplus_down = node4
    else:
        n4.cplus_down = node2; n2.cplus_up = node4
    kernel.char_mesh.append(n4.indx)
    return n4.indx


def cplus_free(node0, node2, node4=None):
    """
    Returns the index of a free-boundary node computed from one streamline node
    already on the free boundary and one interior node on a C+ characteristic.

    node0: index of the point on the free-boundary streamline
    node2: index of initial point along C+ characteristic
    node4: index of solution point (may have a value of -1 or None)
    If -1 is specified as the index for node4, a new node will be
    created for the solution point.
    """
    if node0 == node2:
        raise RuntimeError("Same index given for node0 and node2.")
    n0 = kernel.nodes[node0]
    n2 = kernel.nodes[node2]
    x0 = n0.x; y0 = n0.y; pm0 = n0.nu; th0 = n0.theta; m0 = n1.mach
    x2 = n2.x; y2 = n2.y; pm2 = n2.nu; th2 = n2.theta; m2 = n2.mach
    # Mach angles
    mu0 = asin(1/m0)
    mu2 = asin(1/m2)
    # Use point 0 to get the streamline direction cosines.
    xStream = cos(th0); yStream = sin(th0)
    # Guess at some of the solution point properties.
    # The position will be way off but it is used as part of
    # the convergence check a little further on.
    x4 = 0.5*(x0+x2); y4 = 0.5*(y0+y2)
    th4 = th0
    mu4 = mu0
    # Compute the solution point position and flow properties.
    converged = False
    iteration_count =0
    while (not converged) and (iteration_count < max_iteration):
        x4_old = x4; y4_old = y4
        # Locate solution point by assuming straight-line segments.
        sinCzero = 0.5*(sin(th0) + sin(th4))
        cosCzero = 0.5*(cos(th0) + cos(th4))
        sinCplus  = 0.5*(sin(th2+mu2) + sin(th4+mu4))
        cosCplus  = 0.5*(cos(th2+mu2) + cos(th4+mu4))
        #
        numerator = (x2-x1)*sinCplus - (y2-y1)*cosCplus;
        denominator = cosCzero*sinCplus - sinCzero*cosCplus;
        if abs(denominator) <= 1.0e-12:
            raise RuntimeError("cplus_free: streamline and characteristic are parallel.")
        lambdaCzero = numerator/denominator
        x4 = x0 + lambdaCminus*cosCzero
        y4 = y0 + lambdaCzero*sinCzero
        dx = x4-x4_old; dy = y4-y4_old
        change_in_position = sqrt(dx*dx + dy*dy)
        #
        # Lengths of the characteristic segments.
        dx = x4-x0; dy = y4-y0
        lengthCzero = sqrt(dx*dx + dy*dy)
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
        pm4 = pm0
        th4 = th2 + (pm4-pm2)
        if kernel.axisymmetric:
            # Axisymmetric components will need alternate evaluation at the axis.
            if y0 < 1.0e-6 and y2 < 1.0e-6:
                # Both nodes are close to the axis.
                # Use the axial variation in Mach number to estimate theta/r,
                # and use that for sin(th)/y.
                th_over_r = theta_over_r(m0, m1, x0-x2, kernel.g)
                axiterm2 = sin(mu2)*th_over_r
            else:
                if y2 < 1.0e-6:
                    axiterm2 = sin(mu0)*sin(th0)/y0
                else:
                    axiterm2 = sin(mu2)*sin(th2)/y2
            #
            axiterm4 = sin(mu4)*sin(th4)/y4
            integralCplus  = 0.5*directionCplus*lengthCplus*(axiterm2+axiterm4)
            # Include axisymmetric components.
            th4 -= integralCplus
        #
        iteration_count += 1
        converged = change_in_position < position_tolerance
    # Save the solution-point properties and connect the
    # node into the characteristic mesh.
    if (node4 is None) or (node4 == -1):
        n4 = kernel.Node()
        node4 = n4.indx
    else:
        n4 = kernel.nodes[node4]
    m4 = igf.PM2(pm4, kernel.g)
    n4.x = x4; n4.y = y4; n4.nu = pm4; n4.theta = th4; n4.mach = m4
    # We assume that the principal flow direction
    # is in the positive x-direction.
    if x4 > x0:
        n4.czero_up = node0; n0.czero_down = node4
    else:
        n4.czero_down = node0; n0.czero_up = node4
    if x4 > x2:
        n4.cplus_up = node2; n2.cplus_down = node4
    else:
        n4.cplus_down = node2; n2.cplus_up = node4
    kernel.char_mesh.append(n4.indx)
    return n4.indx


def cminus_free(node0, node1, node4=None):
    """
    Returns the index of a free-boundary node computed from one streamline node
    already on the free boundary and one interior node on a C- characteristic.

    node0: index of the point on the free-boundary streamline
    node1: index of initial point along C- characteristic
    node4: index of solution point (may have a value of -1 or None)
    If -1 is specified as the index for node4, a new node will be
    created for the solution point.
    """
    if node0 == node1:
        raise RuntimeError("Same index given for node0 and node1.")
    n0 = kernel.nodes[node0]
    n1 = kernel.nodes[node1]
    x0 = n0.x; y0 = n0.y; pm0 = n0.nu; th0 = n0.theta; m0 = n0.mach
    x1 = n1.x; y1 = n1.y; pm1 = n1.nu; th1 = n1.theta; m1 = n1.mach
    # Mach angles
    mu0 = asin(1/m0)
    mu1 = asin(1/m1)
    # Use point 0 to get the streamline direction cosines.
    xStream = cos(th0); yStream = sin(th0)
    # Guess at some of the solution point properties.
    # The position will be way off but it is used as part of
    # the convergence check a little further on.
    x4 = 0.5*(x1+x0); y4 = 0.5*(y1+y0)
    th4 = th0
    mu4 = mu0
    # Compute the solution point position and flow properties.
    converged = False
    iteration_count =0
    while (not converged) and (iteration_count < max_iteration):
        x4_old = x4; y4_old = y4
        # Locate solution point by assuming straight-line segments.
        sinCminus = 0.5*(sin(th1-mu1) + sin(th4-mu4))
        cosCminus = 0.5*(cos(th1-mu1) + cos(th4-mu4))
        sinCzero  = 0.5*(sin(th0) + sin(th4))
        cosCzero  = 0.5*(cos(th0) + cos(th4))
        #
        numerator = (x0-x1)*sinCzero - (y0-y1)*cosCzero;
        denominator = cosCminus*sinCzero - sinCminus*cosCzero;
        if abs(denominator) <= 1.0e-12:
            raise RuntimeError("cminus_free: streamline and characteristic are parallel.")
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
        dx = x4-x0; dy = y4-y0
        lengthCzero  = sqrt(dx*dx + dy*dy)
        #
        # Update flow properties at solution point
        # First, assume 2D planar geometry then add
        # axisymmetric contributions if flag is set.
        pm4 = pm0
        th4 = th1 + (pm1-pm4)
        if kernel.axisymmetric:
            # Axisymmetric components will need alternate evaluation at the axis.
            if y1 < 1.0e-6 and y2 < 1.0e-6:
                # Both nodes are close to the axis.
                # Use the axial variation in Mach number to estimate theta/r,
                # and use that for sin(th)/y.
                th_over_r = theta_over_r(m2, m0, x1-x0, kernel.g)
                axiterm1 = sin(mu1)*th_over_r
            else:
                if y1 < 1.0e-6:
                    axiterm1 = sin(mu0)*sin(th0)/y0
                else:
                    axiterm1 = sin(mu1)*sin(th1)/y1
            #
            axiterm4 = sin(mu4)*sin(th4)/y4
            integralCminus = 0.5*directionCminus*lengthCminus*(axiterm1+axiterm4)
            # Include axisymmetric components.
            th4 += integralCminus
        #
        iteration_count += 1
        converged = change_in_position < position_tolerance
    # Save the solution-point properties and connect the
    # node into the characteristic mesh.
    if (node4 is None) or (node4 == -1):
        n4 = kernel.Node()
        node4 = n4.indx
    else:
        n4 = kernel.nodes[node4]
    m4 = igf.PM2(pm4, kernel.g)
    n4.x = x4; n4.y = y4; n4.nu = pm4; n4.theta = th4; n4.mach = m4
    # We assume that the principal flow direction
    # is in the positive x-direction.
    if x4 > x0:
        n4.czero_up = node0; n0.czero_down = node4
    else:
        n4.czero_down = node0; n0.czero_up = node4
    if x4 > x1:
        n4.cminus_up = node1; n1.cminus_down = node4
    else:
        n4.cminus_down = node1; n1.cminus_up = node4
    kernel.char_mesh.append(n4.indx)
    return n4.indx


def streamline_intersection_weights(node0, node1, node2):
    """
    Returns the weights of the intersection point of the extended streamline
    from node0 to the line between two existing points node1 and node2.

    node0: index of a point on the streamline
    node1: index of another existing point
    node2: index of a third existing point
    """
    if node0 == node1:
        raise RuntimeError("Same index given for node0 and node1.")
    if node0 == node2:
        raise RuntimeError("Same index given for node0 and node2.")
    if node1 == node2:
        raise RuntimeError("Same index given for node1 and node2.")
    n0 = kernel.nodes[node0]
    n1 = kernel.nodes[node1]
    n2 = kernel.nodes[node2]
    # Make copies of some of the node data.
    x0 = n0.x; y0 = n0.y; th0 = n0.theta
    x1 = n1.x; y1 = n1.y; th1 = n1.theta
    x2 = n2.x; y2 = n2.y; th2 = n2.theta
    # Guess at some of the intersection point properties and then iterate.
    x4  = 0.5*(x1+x2)
    y4  = 0.5*(y1+y2)
    th4 = 0.5*(th1+th2)
    converged = False
    iteration_count =0
    while (not converged) and (iteration_count < max_iteration):
        x4_old = x4; y4_old = y4
        # Locate solution point by assuming straight-line segments.
        dx12 = x2 - x1;
        dy12 = y2 - y1;
        sinCzero = 0.5*(sin(th0) + sin(th4))
        cosCzero = 0.5*(cos(th0) + cos(th4))
        numerator = (x0-x1)*sinCzero - (y0-y1)*cosCzero
        denominator = dx12*sinCzero - dy12*cosCzero
        if abs(denominator) <= 1.0e-12:
            raise RuntimeError("streamline_intersection_weights: "+
                               "streamline and n1-to-n2-line are parallel.")
        lambda12 = numerator / denominator
        x4 = x1 + lambda12 * dx12
        y4 = y1 + lambda12 * dy12
        dx = x4 - x4_old; dy = y4 - y4_old
        change_in_position = sqrt(dx*dx + dy*dy)
        # Update the sthreamline angle at the estimated intersection point.
        th4 = (1.0-lambda12)*th1 + lambda12*th2
        #
        iteration_count += 1
        converged = change_in_position < position_tolerance
    # At this point, we do not restrict the weights to the range 0.0 to 1.0
    # for an intersection that lies between points 1 and 2.
    return [1.0-lambda12, lambda12]


def add_stream_node(node0, node1, node2, node4=None):
    """
    Returns the index of a new streamline point (node4) by extending the streamline
    from node0 to the line between points 1 and 2.
    The new point will be integrated into the characteristic mesh
    if nodes 1 and 2 are connected directly.

    Raises and exception if the intersection point is not between the points 1 and 2.
    """
    alpha1, alpha2 = streamline_intersection_weights(node0, node1, node2)
    if 0.0 <= alpha2 <= 1.0:
        # The new point lies between the nodes 1 and 2
        # so let's put it in place.
        if (node4 is None) or (node4 == -1):
            n4 = kernel.Node()
            node4 = n4.indx
        else:
            n4 = kernel.nodes[node4]
        insert(node1, node2, node4, alpha2)
        # Connect it into the streamline.
        n0 = kernel.nodes[node0]
        if n4.x > n0.x:
            n4.czero_up = node0; n0.czero_down = node4
        else:
            n4.czero_down = node0; n0.czero_up = node4
    else:
        raise RuntimeError("Intersection point is not between nodes 1 and 2.")
    return n4.indx


def step_stream_node(node0, dL, node4=None, dR=0.9, kdtree=None):
    """
    This function calculates the next node along a streamline by the length dL
    INPUTS:
            node0  - index of starting node
            dL     - length to move along the streamline. A positive value will
                     step downstream while a negative value will step upstream.
            node4  - index of solution point, specify -1 or None for new node
            dR     - Ratio, search radius divided by dL.
                     Historically 0.9 was used to avoid picking up the previous
                     streamline point, however, this is no longer required.
            kdtree - provide a scipy.spatial.kdtree object to use much faster
                     searching algorithms
    OUTPUT:
            node4 - index of solution node or None, if no nearby nodes were found.
    """
    n0 = kernel.nodes[node0]
    # Start by estimating the node4 data based on node0
    x4 = n0.x + dL * cos(n0.theta)
    y4 = n0.y + dL * sin(n0.theta)
    th4 = n0.theta
    pm4 = n0.nu
    #
    R = dR * abs(dL) # Radius of influence for finding nodes
    mu = 2.0 # Smoothing parameter for Shepard interpolation
    near_nodes = kernel.find_nodes_near(x4, y4, tol=R, max_count=10, kdtree=kdtree)
    if len(near_nodes) == 0: return None
    # Using PJs format for data handling
    x = np.zeros_like(near_nodes, dtype=np.float)
    y = np.zeros_like(near_nodes, dtype=np.float)
    nu = np.zeros_like(near_nodes, dtype=np.float)
    theta = np.zeros_like(near_nodes, dtype=np.float)
    for idx, node_idx in enumerate(near_nodes):
        x[idx] = kernel.nodes[node_idx].x
        y[idx] = kernel.nodes[node_idx].y
        theta[idx] = kernel.nodes[node_idx].theta
        nu[idx] = kernel.nodes[node_idx].nu
    # Now calculate new solution point and flow properties
    iteration_count = 0
    converged = False
    r = np.zeros_like(near_nodes)
    Xi = np.zeros_like(near_nodes)
    while (not converged) and (iteration_count < max_iteration):
        x4_old = x4
        y4_old = y4
        r = np.sqrt(x**2 + y**2)
        Xi = (1 - r / R)**mu
        sum_Xi = np.sum(Xi)
        w = Xi / sum_Xi
        pm4 = np.sum(w * nu)
        th4 = np.sum(w * theta)
        sinCzero = 0.5 * ( sin(n0.theta) + sin(th4) )
        cosCzero = 0.5 * ( cos(n0.theta) + cos(th4) )
        x4 = n0.x + cosCzero * dL
        y4 = n0.y + sinCzero * dL
        change_in_position = sqrt((x4 - x4_old)**2 + (y4 - y4_old)**2)
        iteration_count += 1
        converged = change_in_position < position_tolerance
    # Save the solution-point properties and
    # connect the node into the streamline.
    if (node4 is None) or (node4 == -1):
        n4 = kernel.Node()
        node4 = n4.indx
    else:
        n4 = kernel.nodes[node4]
    m4 = igf.PM2(pm4, kernel.g)
    n4.x = x4; n4.y = y4; n4.nu = pm4; n4.theta = th4; n4.mach = m4
    if n4.x > n0.x:
        n4.czero_up = node0
        n0.czero_down = node4
    else:
        n4.czero_down = node0
        n0.czero_up = node4
    return node4

#------------------------------------------------------------------------
# Part B: some compound processes

def march_along_cminus(old_first, new_first, direction):
    """
    Input:
    old_start: index of the starting node on the old characteristic line
    new_start: index of the starting node on the new characteristic line
    direction: 'up' or 'down'

    Returns a list of new node indices (along the new characteristic).
    """
    new_node_indices = [new_first]
    node1 = new_first
    node2 = old_first
    while True:
        node4 = interior(node1, node2)
        new_node_indices.append(node4)
        node1 = node4
        n2 = kernel.nodes[node2]
        if direction == 'up':
            node2 = n2.cminus_up
        elif direction == 'down':
            node2 = n2.cminus_down
        else:
            raise RuntimeError(f"Invalid direction for marching along Cminus: {direction}")
        if node2 is None: break
    #
    return new_node_indices


def march_along_cplus(old_first, new_first, direction):
    """
    Input:
    old_start: index of the starting node on the old characteristic line
    new_start: index of the starting node on the new characteristic line
    direction: 'up' or 'down'

    Returns a list of new node indices (along the new characteristic).
    """
    new_node_indices = [new_first]
    node1 = old_first
    node2 = new_first
    while True:
        node4 = interior(node1, node2)
        new_node_indices.append(node4)
        node2 = node4
        n1 = kernel.nodes[node1]
        if direction == 'up':
            node1 = n1.cplus_up
        elif direction == 'down':
            node1 = n1.cplus_down
        else:
            raise RuntimeError(f"Invalid direction for marching along Cplus: {direction}")
        if node1 is None: break
    #
    return new_node_indices

def get_nodes_along_characteristic(node0, direction):
    """
    Input:
    node0     - index of starting node
    direction - direction of travel to return nodes
                must be one of 'cminus_down', 'cminus_up',
                'cplus_down', 'cplus_up'

    Returns a list of node indices,
    travelling along the characteristic in the desired direction
    """
    node_indices = []
    node_indices.append(node0)
    while True:
        node1 = kernel.nodes[node_indices[-1]]
        if direction == 'cminus_down':
            if node1.cminus_down is None: break
            else: node_indices.append(node1.cminus_down)
        elif direction == 'cminus_up':
            if node1.cminus_up is None: break
            else: node_indices.append(node1.cminus_up)
        elif direction == 'cplus_down':
            if node1.cplus_down is None: break
            else: node_indices.append(node1.cplus_down)
        elif direction == 'cplus_up':
            if node1.cplus_up is None: break
            else: node_indices.append(node1.cplus_up)
        else:
            raise RuntimeError(f"Invalid characteristic direction: {direction}")
    return node_indices
