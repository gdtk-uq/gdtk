Db = 14.0e-3
Rc = Db/2

a = Vector3:new{x=Rc, y=0.0}
b = Vector3:new{x=0.0, y=0.0}
c = Vector3:new{x=Rc, y=Rc}
d = Vector3:new{x=-0.25*Rc, y=0.0}
e = Vector3:new{x=-0.25*Rc, y=Rc}
f = Vector3:new{x=0.5*Rc,  y=2*Rc}
g = Vector3:new{x=Rc,      y=2*Rc}

bc = Arc:new{p0=b, p1=c, centre=a}
dg = Bezier:new{points={d, e, f, g}}

function writeGridProCurve(fname, curve, npts)
   f = assert(io.open(fname, "w"))
   f:write(string.format("%d\n", npts))
   dt = 1.0/(npts-1)
   for i=1,npts do
      t = (i-1)*dt
      p = curve(t)
      f:write(string.format("%.12e %.12e %.12e\n", p.x, p.y, p.z))
   end
end

writeGridProCurve("inflow-boundary.dat", dg, 200)
writeGridProCurve("body.dat", bc, 200)

