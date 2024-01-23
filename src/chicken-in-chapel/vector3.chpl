/** vector3.chpl
 * A module for working with geometric vectors in three dimensions.
 *
 * Peter J. 2024-01-22 Build just enough to play and check performance.
 */


module Geom {
  use Math;

  proc approxEquals(a: real, b: real, tol: real = 1.0e-6): bool {
    return abs(a-b)/(0.5*(abs(a)+abs(b))+1.0) <= tol;
  }

  record Vector3 {
    var x, y, z: real;

    // Methods that do not allocate, I think.

    proc ref set(x: real, y: real, z: real = 0.0) {
      this.x = x; this.y = y; this.z = z;
    }

    proc ref set(other: Vector3) {
      x = other.x; y = other.y; z = other.z;
    }

    proc ref setAsAverage(a: Vector3, b: Vector3) {
      x = 0.5*(a.x+b.x); y = 0.5*(a.x+b.x); z = 0.5*(a.x+b.x);
    }

    proc ref approxEquals(other: Vector3, tol: real = 1.0e-6): bool {
      return Geom.approxEquals(x, other.x, tol) &&
        Geom.approxEquals(y, other.y, tol) &&
        Geom.approxEquals(z, other.z, tol);
    }

    proc ref add(other: Vector3) {
      x += other.x; y += other.y; z += other.z;
    }

    proc ref sub(other: Vector3) {
      x -= other.x; y -= other.y; z -= other.z;
    }

    proc ref mul(other: real) {
      x *= other; y *= other; z *= other;
    }

    proc ref div(other: real) {
      x /= other; y /= other; z /= other;
    }

    proc ref dot(other: Vector3): real {
      return x*other.x + y*other.y + z*other.z;
    }

    // Scales the vector to unit magnitude.
    proc ref normalize() {
      var magnitude = sqrt(x*x + y*y + z*z);
      if magnitude > 0.0 {
        x /= magnitude; y /= magnitude; z /= magnitude;
      } else {
        // Clean up, in case dot product underflows.
        x = 0.0; y = 0.0; z = 0.0;
      }
      // Flush small components to zero.
      const small = 1.0e-30;
      if abs(x) < small { x = 0.0; }
      if abs(y) < small { y = 0.0; }
      if abs(z) < small { z = 0.0; }
    }

    // Transform functions used to reorient vector values in the CFD codes.

    // Rotate v from the global xyz coordinate system into the local frame
    // defined by the orthogonal unit vectors n,t1,t2.
    //
    // We assume, without checking, that these vectors do nicely define
    // such a local system.
    proc ref transformToLocalFrame(n: Vector3, t1: Vector3, t2: Vector3) {
      var v_x = this.dot(n); // normal component
      var v_y = this.dot(t1); // tangential component 1
      var v_z = this.dot(t2); // tangential component 2
      x = v_x; y = v_y; z = v_z;
    }

    // Rotate v back into the global (xyz) coordinate system.
    proc ref transformToGlobalFrame(n: Vector3, t1: Vector3, t2: Vector3) {
      var v_x = x*n.x + y*t1.x + z*t2.x; // global-x
      var v_y = x*n.y + y*t1.y + z*t2.y; // global-y
      var v_z = x*n.z + y*t1.z + z*t2.z; // global-z
      x = v_x; y = v_y; z = v_z;
    }

    // 2D flavour for change of coordinate system functions.

    proc ref transformToLocalFrame(n: Vector3, t1: Vector3) {
      var v_x = x*n.x + y*n.y;   // normal component
      var v_y = x*t1.x + y*t1.y; // tangential component 1
      x = v_x; y = v_y; z = 0.0;
    }

    // Rotate v back into the global (xy) coordinate system.
    proc ref transformToGlobalFrame(n: Vector3, t1: Vector3) {
      var v_x = x*n.x + y*t1.x; // global-x
      var v_y = x*n.y + y*t1.y; // global-y
      x = v_x; y = v_y; z = 0.0;
    }

    // Vector cross product for use in a single statement that will not make temporaries.
    proc ref cross(v1: Vector3, v2: Vector3) {
      x = v1.y * v2.z - v2.y * v1.z;
      y = v2.x * v1.z - v1.x * v2.z;
      z = v1.x * v2.y - v2.x * v1.y;
    }

  } // end record Vector3

  operator+(a: Vector3, b: Vector3): Vector3 {
    return new Vector3(a.x+b.x, a.y+b.y, a.z+b.z);
  }

  operator-(a: Vector3, b: Vector3): Vector3 {
    return new Vector3(a.x-b.x, a.y-b.y, a.z-b.z);
  }

  operator*(a: Vector3, b: real): Vector3 {
    return new Vector3(a.x*b, a.y*b, a.z*b);
  }

  operator/(a: Vector3, b: real): Vector3 {
    return new Vector3(a.x/b, a.y/b, a.z/b);
  }

  proc dot(a: Vector3, b: Vector3): real {
    return a.x*b.x + a.y*b.y + a.z*b.z;
  }

  proc cross(v1: Vector3, v2: Vector3): Vector3 {
    var v3 = new Vector3(0.0, 0.0, 0.0);
    v3.cross(v1, v2);
    return v3;
  }
}

