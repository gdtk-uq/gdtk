"""
billig_patch.lua -- A convenience function for bluff-body simulations.

PJ 2019-05-25
   2022-06-29 introduce quadrant==3 option
   2022-07-04 Allow different x and y scales
   2023-01-02 Python3 flavour
"""

from gdtk.billig import x_from_y
from gdtk.geom.vector3 import Vector3
from gdtk.geom.path import Line, Arc, Spline, ArcLengthParameterizedPath
from gdtk.geom.surface import CoonsPatch

def make_patch(Minf, R, quadrant=2,
               xc=0.0, yc=0.0,
               scale=None, x_scale=None, y_scale=None,
               axisymmetric=False, theta=0.0):
    """
    Construct a surface patch for use in a bluff-body simulation
    using arguments found in the supplied table.
    For the quadrant==2 case the patch is above the x-axis.

                         ++d[#d]
                       +   |
                     +     | North
           West    +       |
                  +        |
                 +      ---b
                +      /
               +      /
               +     |
             d[1]----a     c  ----> x
               South

    For the quandrant==3 case the patch is below the x-axis.

               North
             d[1]----a     c  ----> x
               +     |
               +      \
                +      \
                 +      ---b
           West   +        |
                   +       | South
                     +     |
                       +   |
                         ++d[#d]

    Input:
    Minf: Free-stream Mach number.
    R: Radius of cylinder or sphere.
    quadrant: which of the quadrants (shown above)
    xc, yc: Position of centre of body.
    xscale, yscale: Scale for accommodating thermochemical variation
                    away from ideal low-Temperature air.
    Set the scales using the most specific information available,
    eventually defaulting to a scale of 1.0 if nothing is specified.
    axisymmetric: True is body is a sphere-cone
                  False if we have a blunted plate.
    theta: Angle of aft-body with respect to freestream direction, radians.
    """
    if scale:
        if not x_scale: x_scale = scale
        if not y_scale: y_scale = scale
    if not x_scale: x_scale = 1.0
    if not y_scale: y_scale = 1.0
    #
    a = Vector3(x=xc-R, y=yc)
    b = Vector3(x=xc, y=yc+R)
    if quadrant == 3: b = Vector3(x=xc, y=yc-R)
    c = Vector3(x=xc, y=yc)
    body = Arc(a=a, b=b, c=c) if quadrant == 2 else Arc(a=b, b=a, c=c)
    #
    # In order to have a grid that fits reasonably close the the shock,
    # use Billig's shock shape correlation to generate
    # a few sample points along the expected shock position.
    print("Points on Billig's correlation.")
    xys = []
    for i,y in enumerate([0.0, 0.2, 0.4, 0.6, 1.0, 1.4, 1.6, 2.0, 2.37]):
        # y is normalized for R=1
        x = x_from_y(y*R, Minf, theta, axisymmetric, R)
        xys.append((x, y*R)) # a new coordinate pair
        # print("x=", x, "y=", y*R)
    # Scale the Billig distances, depending on the expected behaviour
    # relative to the gamma=1.4 ideal gas.
    d = [] # to keep the nodes for the shock boundary
    if quadrant == 2:
        for i, xy in enumerate(xys):
            # the outer boundary should be a little further than the shock itself
            d.append(Vector3(x=-x_scale*xy[0]+xc, y=y_scale*xy[1]+yc))
        shock = ArcLengthParameterizedPath(Spline(d))
        xaxis = Line(p0=d[0], p1=a) # shock to nose of body
        outlet = Line(p0=d[-1], p1=b) # shock to last-point of body
        patch = CoonsPatch(south=xaxis, north=outlet, west=shock, east=body)
    else:
        # For quadrant 3, the shock path progresses fro the outer point to the axis.
        for i,xy in enumerate(reversed(xys)):
            # the outer boundary should be a little further than the shock itself
            d.append(Vector3(x=-x_scale*xy[0]+xc, y=-y_scale*xy[1]+yc))
        shock = ArcLengthParameterizedPath(Spline(d))
        xaxis = Line(p0=d[-1], p1=a) # shock to nose of body
        outlet = Line(p0=d[0], p1=b) # shock to last-point of body
        patch = CoonsPatch(south=outlet, north=xaxis, west=shock, east=body)
    return {'patch':patch, 'points':{'a':a, 'b':b, 'c':c, 'd':d}}
