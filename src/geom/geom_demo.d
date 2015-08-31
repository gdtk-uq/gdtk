/**
 * geom_demo.d Demonstrate some of the behaviour of the geometric primitives.
 *
 * Author: Peter J.
 * Version: 2014-06-16
 */

import std.stdio;
import geom;

void main()
{
    writeln("Begin demonstration of the geometric primitives.");
    Vector3 a = Vector3([1.0, 2.2, 3.0]);
    Vector3 b = Vector3(1.0);
    writeln("a = ", a, ", b = ", b);
    Vector3 c = a + b;
    writeln("c = a + b = ", c);
    Vector3 d = a;
    a.refy = 99.0;
    writeln("a = ", a, ", d = ", d); 
    Vector3 e = a * 2.0;
    Vector3 f = 3 * d; // int promoted to double OK
    writeln("e = ", e, " f = ", f);
    Vector3 g = d / 3.0;
    writeln("g = ", g);
    g += f;
    writeln("g += f -> ", g);
    g /= 2.0;
    writeln("g /= 2.0 -> ", g);
    Vector3 u = unit(g*2);
    writeln("unit(g*2) = ", u, " magnitude = ", abs(u));
    Vector3 x = Vector3(1.0, 0.0, 0.0);
    Vector3 y = Vector3(0.0, 1.0, 0.0);
    Vector3 z = cross(x,y);
    writeln("z = ", z);

    Vector3 n = unit(Vector3(1.0,1.0,0.0));
    Vector3 t1 = unit(Vector3(-1.0,1.0,0.0));
    Vector3 t2 = cross(n, t1);
    Vector3 h = Vector3(1.0,0.0,1.0);
    Vector3 h_ref = Vector3(h);
    writeln("original h = ", h);
    h.transform_to_local_frame(n, t1, t2);
    writeln("in local frame h = ", h);
    h.transform_to_global_frame(n, t1, t2);
    writeln("back to global frame h = ", h);

    writeln("Try out geometric functions that build on Vector3 objects.");
    a = Vector3(1.0, 0.0, 0.0); // plane through a,b,c
    b = Vector3(1.0, 1.0, 0.0);
    c = Vector3(0.5, 0.0, 0.0);
    Vector3 qr = Vector3(3.0, 3.0, -3.0); // direction
    Vector3 q = Vector3(0.0, 0.0, 1.0); // start point
    int flag =  project_onto_plane(q, qr, a, b, c);
    writeln("projected point q = ", q);

    writeln("Done.");
}
