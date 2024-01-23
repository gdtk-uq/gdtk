/** vector3.chpl
 * A module for working with geometric vectors in three dimensions.
 *
 * Peter J. 2024-01-22 Build just enough to play and check performance.
 *          2024-01-23 Reimplemented all of the chicken functions.
 */


module Geom {
  use Math;

  proc approxEquals(a: real, b: real, tol: real = 1.0e-6): bool {
    return abs(a-b)/(0.5*(abs(a)+abs(b))+1.0) <= tol;
  }

  const zero: real = 0.0;
  const one: real = 1.0;
  const half: real = 0.5;
  const one_quarter: real = 0.25;
  const one_fifth: real = 0.2;
  const one_eighth: real = 0.125;
  const six: real = 6.0;

  record Vector3 {
    var x, y, z: real;

    // Methods that do not allocate, I think.

    proc ref set(x: real, y: real, z: real = zero) {
      this.x = x; this.y = y; this.z = z;
    }

    proc ref set(other: Vector3) {
      x = other.x; y = other.y; z = other.z;
    }

    proc ref setAsAverage(a: Vector3, b: Vector3) {
      x = half*(a.x+b.x); y = half*(a.x+b.x); z = half*(a.x+b.x);
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
        x = zero; y = zero; z = zero;
      }
      // Flush small components to zero.
      const small = 1.0e-30;
      if abs(x) < small { x = zero; }
      if abs(y) < small { y = zero; }
      if abs(z) < small { z = zero; }
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
      x = v_x; y = v_y; z = zero;
    }

    // Rotate v back into the global (xy) coordinate system.
    proc ref transformToGlobalFrame(n: Vector3, t1: Vector3) {
      var v_x = x*n.x + y*t1.x; // global-x
      var v_y = x*n.y + y*t1.y; // global-y
      x = v_x; y = v_y; z = zero;
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
    var v3 = new Vector3(zero, zero, zero);
    v3.cross(v1, v2);
    return v3;
  }

  // Quadrilateral properties of centroid, associated unit normals and area.
  //   p3-----p2
  //   |      |
  //   |      |
  //   p0-----p1
  // Resultant normal vector is up, toward you.
  // Assume that all points are in the one plane.
  proc quadProperties(p0: Vector3, p1: Vector3, p2: Vector3, p3: Vector3,
                      ref centroid: Vector3,
                      ref n: Vector3, ref t1: Vector3, ref t2: Vector3,
                      ref area: real,
                      tol: real = 1.0e-12, area_tol: real = 1.0e-20) {
    centroid.set(p0); centroid.add(p1); centroid.add(p2); centroid.add(p3);
    centroid.mul(one_quarter);
    // Compute areas via the cross products.
    // Vector3 vector_area = one_quarter * cross(p1-p0+p2-p3, p3-p0+p2-p1);
    var p01, p03, vectorArea: Vector3;
    p01.set(p1); p01.sub(p0); p01.add(p2); p01.sub(p3);
    p03.set(p3); p03.sub(p0); p03.add(p2); p03.sub(p1);
    vectorArea.cross(p01, p03); vectorArea.mul(one_quarter);
    // unit-normal and area
    // area = abs(vector_area);
    area = sqrt(vectorArea.dot(vectorArea));
    if (area > area_tol) {
      // n = unit(vectorArea);
      n.set(vectorArea); n.div(area);
      // Tangent unit-vectors:
      // t1 is parallel to side01 and side32,
      // t2 is normal to n and t1
      // t1 = unit((p1-p0)+(p2-p3)); // Works even if one edge has zero length.
      t1.set(p01);
      t1.normalize();
      // t2 = unit(cross(n, t1)); // Calling unit() to tighten up the magnitude.
      t2.cross(n, t1);
      t2.normalize();
    } else {
      // We have nothing meaningful to put into the unit vectors,
      // so, put in an arbitrary but orthogonal set.
      n.set(one,zero,zero); t1.set(zero,one,zero); t2.set(zero,zero,one);
    }
  } // end quad_properties()

  // Because of the way we lose precision when reading and writing files,
  // it may be that the vertices are not quite in their ideal position.
  // We need a couple of finite, but small, tolerances to deal with
  // collapsed volumes.
  const smallButSignificantVolume: real = 1.0e-12;
  const verySmallVolume: real = 1.0e-20;

  // For the tetrahedron geometry, we consider p0,p1,p2 the base.
  // Looking from p3 back toward the base, a counter-clockwise cycle
  // p0->p1->p2->p0 gives a positive volume.

  proc tetrahedronVolume(p0: Vector3, p1: Vector3, p2: Vector3, p3: Vector3): real {
    // Return dot(p3-p0, cross(p1-p0, p2-p0)) / 6.0;
    var d01, d02, d03, c: Vector3;
    d01.set(p1); d01.sub(p0);
    d02.set(p2); d02.sub(p0);
    c.cross(d01, d02);
    d03.set(p3); d03.sub(p0);
    return d03.dot(c) / six;
  } // end tetrahedronVolume()

  proc tetrahedronProperties(p0: Vector3, p1: Vector3, p2: Vector3, p3: Vector3,
                             ref centroid: Vector3, ref volume: real) {
    // Variant without L_min calculation.
    centroid.set(p0); centroid.add(p1); centroid.add(p2); centroid.add(p3);
    centroid.mul(one_quarter);
    volume = tetrahedronVolume(p0, p1, p2, p3);
  } // end tetrahedronProperties()

  // J. Grandy (1997) Efficient Computation of Volume of Hexahedral Cells UCRL-ID-128886.
  // Base of each dipyramid is specified clockwise from the outside.
  // pc is apex
  // pb is barycentre of base quad.
  // A base quad cycle p0->p1->p2->p3->p0 that is counterclockwise when looking
  // towards it from pc will result in a positive volume.
  // A negative volume indicates that the cycle is clockwise when looking from pc.
  proc tetragonalDipyramidVolume(p0: Vector3, p1: Vector3, p2: Vector3, p3: Vector3,
                                 pb: Vector3, pc: Vector3): real
  {
    // volume = dot(pc-pb, cross(p1-p0+p2-p3, p3-p0+p2-p1)) / 12.0;
    var bc, p01, p03, vectorArea: Vector3;
    p01.set(p1); p01.sub(p0); p01.add(p2); p01.sub(p3); p01.mul(half);
    p03.set(p3); p03.sub(p0); p03.add(p2); p03.sub(p3); p03.mul(half);
    vectorArea.cross(p01, p03);
    bc.set(pc); bc.sub(pb);
    return bc.dot(vectorArea) / 3.0;
  } // end tetragonalDipyramidVolume()

  proc pyramidProperties(p0: Vector3, p1: Vector3, p2: Vector3, p3: Vector3,
                         p4: Vector3, trueCentroid: bool,
                         ref centroid: Vector3, ref volume: real)
  {
    // p0-p1-p2-p3 is the quadrilateral base, p4 is the peak.
    // cycle of base vertices is counterclockwise, viewed from p4.
    //
    // Split into 4 tetrahedra and sum contributions to volume and moment.
    var pmB: Vector3; // Mid-point of quadrilateral base.
    pmB.set(p0); pmB.add(p1); pmB.add(p2); pmB.add(p3); pmB.mul(one_quarter);
    //
    volume = zero;
    var moment: Vector3; moment.set(zero, zero, zero);
    var tetVolume: real; var tetCentroid: Vector3;
    tetrahedronProperties(p0, p1, pmB, p4, tetCentroid, tetVolume);
    volume += tetVolume; tetCentroid.mul(tetVolume); moment.add(tetCentroid);
    tetrahedronProperties(p1, p2, pmB, p4, tetCentroid, tetVolume);
    volume += tetVolume; tetCentroid.mul(tetVolume); moment.add(tetCentroid);
    tetrahedronProperties(p2, p3, pmB, p4, tetCentroid, tetVolume);
    volume += tetVolume; tetCentroid.mul(tetVolume); moment.add(tetCentroid);
    tetrahedronProperties(p3, p0, pmB, p4, tetCentroid, tetVolume);
    volume += tetVolume; tetCentroid.mul(tetVolume); moment.add(tetCentroid);
    //
    if abs(volume) > zero { moment.div(volume); } // to get overall centroid
    if trueCentroid {
      centroid = moment;
    } else {
      // approximating the centroid via a simple averaging of vertex positions
      // has shown to be more robust when importing an unstructurd grid and
      // also appears to provide a better distribution of points for the
      // least-squares gradient estimation [KAD 2022-08-03].
      centroid.set(p0); centroid.add(p1); centroid.add(p2);
      centroid.add(p3); centroid.add(p4); centroid.mul(one_fifth);
    }
  } // end pyramidProperties()

  proc hexCellProperties(p0: Vector3, p1: Vector3, p2: Vector3, p3: Vector3,
                         p4: Vector3, p5: Vector3, p6: Vector3, p7: Vector3,
                         trueCentroid: bool,
                         ref centroid: Vector3, ref volume: real,
                         ref iLen: real, ref jLen: real, ref kLen: real)
  {
    // PJ 10-Sep-2012
    // When computing the volume of Rolf's thin, warped cells, we have to do
    // something better than splitting our cell into six tetrahedra, so we do that
    // by dividing the hex cell into six tetragonal dipyramids with the original
    // faces as the pyramid bases.
    //
    // Estimate the centroid so that we can use it as the peak
    // of each of the pyramid sub-volumes.
    centroid.set(p0); centroid.add(p1); centroid.add(p2); centroid.add(p3);
    centroid.add(p4); centroid.add(p5); centroid.add(p6); centroid.add(p7);
    centroid.mul(one_eighth);
    // Mid-points of faces.
    var pmN, pmE, pmS, pmW, pmT, pmB: Vector3;
    pmN.set(p3); pmN.add(p2); pmN.add(p6); pmN.add(p7); pmN.mul(one_quarter);
    pmE.set(p1); pmE.add(p2); pmE.add(p6); pmE.add(p5); pmE.mul(one_quarter);
    pmS.set(p0); pmS.add(p1); pmS.add(p5); pmS.add(p4); pmS.mul(one_quarter);
    pmW.set(p0); pmW.add(p3); pmW.add(p7); pmW.add(p4); pmW.mul(one_quarter);
    pmT.set(p4); pmT.add(p5); pmT.add(p6); pmT.add(p7); pmT.mul(one_quarter);
    pmB.set(p0); pmB.add(p1); pmB.add(p2); pmB.add(p3); pmB.mul(one_quarter);
    // Lengths between mid-points of faces.
    // Note that we are assuming that the hexahedron is not very skewed
    // when we later use these values as the widths of the hex cell.
    var d: Vector3;
    d.set(pmE); d.sub(pmW); iLen = sqrt(d.dot(d));
    d.set(pmN); d.sub(pmS); jLen = sqrt(d.dot(d));
    d.set(pmT); d.sub(pmB); kLen = sqrt(d.dot(d));
    // writeln("Single hexahedron divided into six tetragonal dipyramids.");
    // J. Grandy (1997) Efficient Computation of Volume of Hexahedral Cells UCRL-ID-128886.
    // Base of each dipyramid is specified clockwise from the outside.
    var subVolume: real; var subCentroid, moment: Vector3;
    volume = zero; moment.set(zero, zero, zero);
    pyramidProperties(p6, p7, p3, p2, centroid, true, subCentroid, subVolume);
    volume += subVolume; subCentroid.mul(subVolume); moment.add(subCentroid);
    pyramidProperties(p5, p6, p2, p1, centroid, true, subCentroid, subVolume);
    volume += subVolume; subCentroid.mul(subVolume); moment.add(subCentroid);
    pyramidProperties(p4, p5, p1, p0, centroid, true, subCentroid, subVolume);
    volume += subVolume; subCentroid.mul(subVolume); moment.add(subCentroid);
    pyramidProperties(p7, p4, p0, p3, centroid, true, subCentroid, subVolume);
    volume += subVolume; subCentroid.mul(subVolume); moment.add(subCentroid);
    pyramidProperties(p7, p6, p5, p4, centroid, true, subCentroid, subVolume);
    volume += subVolume; subCentroid.mul(subVolume); moment.add(subCentroid);
    pyramidProperties(p0, p1, p2, p3, centroid, true, subCentroid, subVolume);
    volume += subVolume; subCentroid.mul(subVolume); moment.add(subCentroid);
    //
    if (volume < zero && abs(volume) < smallButSignificantVolume) ||
      (volume >= zero && volume < verySmallVolume) {
        // We assume that we have a collapsed hex cell;
        // no real problem here but it may be a problem for the client code.
        // That code should test the value of volume, on return.
        volume = zero;
      }
    //
    if abs(volume) > zero { moment.div(volume); } // to get overall centroid
    if trueCentroid { centroid = moment; }
  } // end hexCellProperties()

} // end module Geom

