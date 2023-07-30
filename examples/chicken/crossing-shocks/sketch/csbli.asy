// csbli.asy
// Drawing of the half domain with plate and wedge.
// PJ 2023-07-29
settings.outformat = "pdf";
settings.prc = false;
settings.render = 4;

import three;
size(8cm,0);
// Work in millimetres.
real z1 = 50.0mm;
real x1 = 50.0mm;
real L = 120.0mm;
real theta = radians(10.0);
real x2 = x1 + L*cos(theta);
real dz2 = L*sin(theta);
real z2 = z1 - dz2;
real th2 = radians(30.0);
real L2 = dz2/sin(th2);
real x3 = x2 + L2*cos(th2);
real x4 = x3 + 30.0mm;
real H = 50.0mm;

currentprojection = perspective(-100mm, 200mm, -200mm, up=Y);

guide3 plate = (0,0,0)--(x4,0,0)--(x4,0,z1)--(x3,0,z1)--(x2,0,z2)--
  (x1,0,z1)--(0,0,z1)--cycle;
guide3 wedge = (x1,0,z1)--(x2,0,z2)--(x2,H,z2)--(x1,H,z1)--cycle;

void wireframe(
  triple p000, triple p100, triple p110, triple p010,
  triple p001, triple p101, triple p111, triple p011) {
  pen mypen = linewidth(0.5pt)+lightblue;
  draw(p000--p100--p110--p010--p000, p=mypen); // bottom
  draw(p001--p101--p111--p011--p001, p=mypen); // top
  draw(p100--p110--p111--p101--p100, p=mypen); // east
  draw(p000--p010--p011--p001--p000, p=mypen); // west
  draw(p000--p100--p101--p001--p000, p=mypen); // south
  draw(p010--p110--p111--p011--p010, p=mypen); // north
}

wireframe((x1,0,0), (0,0,0), (0,0,z1), (x1,0,z1),
          (x1,H,0), (0,H,0), (0,H,z1), (x1,H,z1));
wireframe((x2,0,0), (x1,0,0), (x1,0,z1), (x2,0,z2),
          (x2,H,0), (x1,H,0), (x1,H,z1), (x2,H,z2));
wireframe((x3,0,0), (x2,0,0), (x2,0,z2), (x3,0,z1),
          (x3,H,0), (x2,H,0), (x2,H,z2), (x3,H,z1));
wireframe((x4,0,0), (x3,0,0), (x3,0,z1), (x4,0,z1),
          (x4,H,0), (x3,H,0), (x3,H,z1), (x4,H,z1));

draw(wedge, p=linewidth(1pt));
draw(surface(wedge), surfacepen=gray(2));
draw(plate, p=linewidth(1pt));
draw(surface(plate), surfacepen=lightgray);

draw((x4+5mm,0,0)--(x4+20mm,0,0), arrow=Arrow3()); label("x", (x4+20mm,0,0), N);
draw((0,H+5mm,0)--(0,H+20mm,0), arrow=Arrow3()); label("y", (0,H+20mm,0), E);
draw((0,0,z1+5mm)--(0,0,z1+20mm), arrow=Arrow3()); label("z", (0,0,z1+20mm), E);

