/** vector3_test.chpl
 * Try out some of the Vector3 functions.
 *
 * Build with:
 * $ chpl vector3_test.chpl vector3.chpl
 *
 * Run with:
 * $ ./vector3_test
 *
 * Peter J. 2024-01-22
 */

use Geom;

proc main() {
  var v1 = new Vector3(1.0, 2.0, 3.0);           // One way to initialize.
  var v2 = new Vector3(x=1.1, y=2.2, z=3.3); // Another way to initialize.
  v1 = (v1+v2)*2.0 - (new Vector3(1.0, 1.0, 1.0));
  var v3 = new Vector3(x=3.2, y=7.4, z=11.6);
  if !v1.approxEquals(v3, 1.0e-9) {
    writeln("Vector3 operator+-* error v1=", v1, " v3=", v3);
  }

  v1.x = 1.0; v1.y = 2.0; v1.z = 3.0;
  v1.div(2);
  if !approxEquals(dot(v1,v1), 14.0/4) {
    writeln("Dot product error dot(v1,v1)=", dot(v1,v1));
  }

  v1.set(1.0, 2.0, 3.0);
  v1.add(v2); v1.mul(2.0); v1.sub(new Vector3(1.0, 1.0, 1.0));
  if !v1.approxEquals(v3, 1.0e-9) {
    writeln("Vector3 add sub mul error v1=", v1, " v3=", v3);
  }

  var n = new Vector3(1.0, 1.0); n.normalize();
  var t1 = new Vector3(-1.0, 1.0); t1.normalize();
  if !t1.approxEquals(new Vector3(-1.0/sqrt(2.0), 1.0/sqrt(2.0), 0.0), 1.0e-9) {
    writeln("Vector3 normalize error t1=", t1);
  }
  var t2 = cross(n, t1);
  if !t2.approxEquals(new Vector3(0.0, 0.0, 1.0), 1.0e-9) {
    writeln("Vector3 cross product error t2=", t2);
  }

  v1.set(1.0, 0.0, 0.0);
  v1.transformToLocalFrame(n, t1, t2);
  if !v1.approxEquals(new Vector3(1.0/sqrt(2.0), -1.0/sqrt(2.0), 0.0), 1.0e-9) {
    writeln("Vector3 transformToLocalFrame error v1=", v1);
  }
  v1.transformToGlobalFrame(n, t1, t2);
  if !v1.approxEquals(new Vector3(1.0, 0.0, 0.0), 1.0e-9) {
    writeln("Vector3 transformToGlobalFrame error v1=", v1);
  }
}
