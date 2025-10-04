/**
 * luasurface_demo.d
 * Demonstrate the wrapped Surface objects.
 *
 * Authors: Rowan G. and Peter J.
 * Version: 2015-02-24
 */

import std.stdio;
import std.conv;
import std.string;
import util.lua;
import geom.luawrap;

void main()
{
    writeln("Begin demonstration of Lua connection to Surface objects.");
    auto L = luaL_newstate();
    luaL_openlibs(L);
    registerVector3(L);
    registerPaths(L);
    registerSurfaces(L);
    string test_code = `
print("Construct from edges")
a = Vector3:new{x=0.0, y=0.0}
b = {x=0.0, y=1.0} -- should accept a table with named coordinates, also
c = Vector3:new{x=1.0, y=0.0}
d = Vector3:new{x=1.0, y=1.0}
surf = CoonsPatch:new{north=Line:new{p0=b, p1=d}, east=Line:new{p0=c, p1=d},
                      south=Line:new{p0=a, p1=c}, west=Line:new{p0=a, p1=b}}
print("CoonsPatch representation: ", surf)
print("surf(0.5,0.5)= ", surf(0.5, 0.5))
print("vector area= ", surf:area())
--
print("Try construction using corners")
surf2 = CoonsPatch:new{p00=a, p01=b, p11=c, p10=d}
p = surf2:eval(0.5, 0.5)
print("same point p= ", p)
--
print("AO patch")
p00 = Vector3:new{x=0.0, y=0.1, z=3.0}
p10 = {x=1.0, y=0.4, z=3.0} -- should accept a table with named coordinates, also
p11 = Vector3:new{x=1.0, y=1.1, z=3.0}
p01 = Vector3:new{x=0.0, y=1.1, z=3.0}
my_aopatch = AOPatch:new{p00=p00, p10=p10, p11=p11, p01=p01}
p = my_aopatch(0.1, 0.1);
print("my_aopatch(0.1, 0.1)= ", p)
--
print("LuaFnSurface")
function myLuaFunction(r, s)
   -- Simple plane
   return {x=r, y=s, z=0.0}
end
myFnSurface = LuaFnSurface:new{luaFnName="myLuaFunction"}
print("myFnSurface= ", myFnSurface)
print("myFnSurface(0.3, 0.4)= ", myFnSurface(0.3, 0.4))
--
print("SubRangedSurface")
srs = SubRangedSurface:new{underlying_psurface=my_aopatch,
                           r0=0.0, r1=0.5, s0=0.0, s1=0.5}
print("srs(0.2,0.2)=", srs(0.2,0.2))
--
print("ChannelPatch")
cA = Line:new{p0=Vector3:new{x=0.0,y=0.0}, p1=Vector3:new{x=1.0,y=0.0}}
cB = Line:new{p0=Vector3:new{x=0.0,y=0.25}, p1=Vector3:new{x=1.0,y=1.0}}
chanp = ChannelPatch:new{south=cA, north=cB}
print("chanp= ", chanp)
print("chanp(0.5,0.5)= ", chanp(0.5, 0.5))
bpath = chanp:make_bridging_path(0.0)
print("bpath=", bpath)
--
print("RuledSurface")
cA = Line:new{p0=Vector3:new{x=0.0,y=0.0}, p1=Vector3:new{x=1.0,y=0.0}}
cB = Line:new{p0=Vector3:new{x=0.0,y=0.25}, p1=Vector3:new{x=1.0,y=1.0}}
ruledp = RuledSurface:new{edge0=cA, edge1=cB, ruled_direction="s"}
print("ruledp= ", ruledp)
print("ruledp(0.5,0.5)= ", ruledp(0.5, 0.5))
--
print("SweptPathPatch demo")
cA = Line:new{p0=Vector3:new{x=0.0,y=0.0,z=0.0}, p1=Vector3:new{x=0.0,y=1.0,z=0.0}}
cB = Line:new{p0=Vector3:new{x=1.0,y=0.25,z=0.0}, p1=Vector3:new{x=2.0,y=0.25,z=0.0}}
spp = SweptPathPatch:new{west=cA, south=cB}
print("spp= ", spp)
print("spp(0.5,0.5)= ", spp(0.5, 0.5))
--
print("SpherePatch demo")
R = 1.0
east_patch = SpherePatch:new{radius=R, centre={0.0,0.0,0.0}, face_name="east"}
mid_e = east_patch(0.5, 0.5)
p5_east = east_patch(0.0, 1.0)
p6_east = east_patch(1.0, 1.0)
print("mid_e=", mid_e)
print("p5_east=", p5_east)
print("p6_east=", p6_east)
south_patch = SpherePatch:new{radius=R, centre={0.0,0.0,0.0}, face_name="south", which_part="top"}
mid_s = south_patch(0.5, 0.0)
p5_south = south_patch(1.0, 1.0)
print("mid_s=", mid_s)
print("p5_south=", p5_south)
north_patch = SpherePatch:new{radius=R, centre={0.0,0.0,0.0}, face_name="north", which_part="top-east"}
mid_n = north_patch(0.0, 0.0)
p6_north = north_patch(1.0, 1.0)
print("mid_n=", mid_n)
print("p6_north=", p6_north)
--
print("BezierPatch demo")
L = 2.0
H = 1.0
n = 3
m = 4
dx = L/n
dy = H/m
Q = {}
for i=1,n+1 do
   Q[i] = {}
   for j=1,m+1 do
      Q[i][j] = Vector3:new{x=(i-1)*dx, y=(j-1)*dy, z=0.0}
   end
end
bezPatch = BezierPatch:new{points=Q}
print("bezPatch(0.2, 0.25)= ", bezPatch(0.2, 0.25))
print("Expected result: x=0.4, y=0.25, z=0.0")
--
print("NURBSPatch demo")
P = {}
P[1] = {{x=0.0, y=0.0, z=0.0},
        {x=0.0, y=0.0, z=0.0},
        {x=0.0, y=0.0, z=0.0},
        {x=0.0, y=0.0, z=0.0},
        {x=0.0, y=0.0, z=0.0}}
P[2] = {{x=0.0, y=0.0, z=0.0},
        {x=0.0, y=0.0, z=0.0},
        {x=0.0, y=0.0, z=0.0},
        {x=0.0, y=0.0, z=0.0},
        {x=0.0, y=0.0, z=0.0}}
P[3] = {{x=0.0, y=0.0, z=0.0},
        {x=0.0, y=2.0, z=4.0},
        {x=0.0, y=3.0, z=2.0},
        {x=0.0, y=2.0, z=0.0},
        {x=0.0, y=0.0, z=0.0}}
P[4] = {{x=0.0, y=0.0, z=0.0},
        {x=2.0, y=3.0, z=4.0},
        {x=2.0, y=4.0, z=2.0},
        {x=2.0, y=3.0, z=0.0},
        {x=0.0, y=0.0, z=0.0}}
P[5] = {{x=0.0, y=0.0, z=0.0},
        {x=4.0, y=2.0, z=4.0},
        {x=4.0, y=3.0, z=2.0},
        {x=4.0, y=2.0, z=0.0},
        {x=0.0, y=0.0, z=0.0}}
P[6] = {{x=0.0, y=0.0, z=0.0},
        {x=0.0, y=0.0, z=0.0},
        {x=0.0, y=0.0, z=0.0},
        {x=0.0, y=0.0, z=0.0},
        {x=0.0, y=0.0, z=0.0}}
P[7] = {{x=0.0, y=0.0, z=0.0},
        {x=0.0, y=0.0, z=0.0},
        {x=0.0, y=0.0, z=0.0},
        {x=0.0, y=0.0, z=0.0},
        {x=0.0, y=0.0, z=0.0}}
P[8] = {{x=0.0, y=0.0, z=0.0},
        {x=0.0, y=0.0, z=0.0},
        {x=0.0, y=0.0, z=0.0},
        {x=0.0, y=0.0, z=0.0},
        {x=0.0, y=0.0, z=0.0}}
w = {}
w[1] = {0.0, 0.0, 0.0, 0.0, 0.0}
w[2] = {0.0, 0.0, 0.0, 0.0, 0.0}
w[3] = {0.0, 1.0, 2.0, 1.0, 0.0}
w[4] = {0.0, 2.0, 6.0, 2.0, 0.0}
w[5] = {0.0, 1.0, 2.0, 1.0, 0.0}
w[6] = {0.0, 0.0, 0.0, 0.0, 0.0}
w[7] = {0.0, 0.0, 0.0, 0.0, 0.0}
w[8] = {0.0, 0.0, 0.0, 0.0, 0.0}

p = 2
q = 2

U = {0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 4.0, 5.0, 5.0, 5.0}
V = {0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 3.0}

nurbsPatch = NURBSPatch:new{points=P, weights=w,
                            u_knots=U, v_knots=V,
                            u_degree=p, v_degree=q}

print("nurbsPatch(2.5, 1.0)= ", nurbsPatch(2.5, 1.0))
print("Expected result= ", 2.0, 98./27., 68./27.)

--
print("Utility functions")
print("isSurface(my_aopatch)= ", isSurface(my_aopatch))
print("isSurface(surf2)= ", isSurface(surf2));
print("isSurface(a)= ", isSurface(a));
surf3 = makePatch{north=Line:new{p0=b, p1=d},
                  east=Line:new{p0=c, p1=d},
                  south=Line:new{p0=a, p1=c},
                  west=Line:new{p0=a, p1=b},
                  gridType="ao"}
print("surf3= ", surf3)
print("Done luasurface_demo.")
    `;
    if ( luaL_dostring(L, toStringz(test_code)) != 0 ) {
        writeln("There was a problem interpreting the test code.");
        writeln(to!string(lua_tostring(L, -1)));
    }
    writeln("Done with luageom_demo.");
}
