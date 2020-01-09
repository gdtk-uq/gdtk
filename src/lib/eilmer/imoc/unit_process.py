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
from math import sin, cos, sqrt, asin

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
    xStream = cos(th2);
    yStream = sin(th2);
    # Guess at some of the solution point properties.
    # The position will be way off but it is used as part of
    # the convergence check a little further on.
    x4 = 0.5*(x1+x2);
    y4 = 0.5*(y1+y2);
    th4 = 0.5*(th1+th2);
    mu4 = 0.5*(mu1+mu2);
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
