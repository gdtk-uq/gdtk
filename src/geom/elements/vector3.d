/**
 * vector3.d  Vector primitives for our 3D world.
 *
 * Author: Peter J.
 * Version: 2014-06-16 first cut.
 *          2015-02-18 Attempt to reduce the number of redundant object copies
 *          by declaring some of the "in Vector3" parameters as "ref const(Vector3)".
 *          This has been somewhat successful, however, a few "in" parameters remain
 *          so that the vector arithmetic is a lot cleaner in the code.
 *          Effectively this hides the creation of Vector3 temporaries
 *          that would otherwise have to appear explicitly in the code.
 *          2017-11-26 Repackage to make better use of the file system with smaller files.
 *          2018-05-29 Complex numbers accommodated.
 *          2022-01-12 Change internal storage away from an array.
 */
module geom.elements.vector3;

import std.conv;
import std.stdio;
import std.math;
import std.string;
import std.json;
import ntypes.complex;
import nm.number;

struct Vector3 {
    public number x, y, z;

    @nogc this(in number[] p)
    {
        switch ( p.length ) {
        case 0: x = y = z = 0.0; break;
        case 1: x = p[0]; y = z = 0.0; break;
        case 2: x = p[0]; y = p[1]; z = 0.0; break;
        default: x = p[0]; y = p[1]; z = p[2]; break;
        }
    }

    @nogc this(in number x, in number y=0.0, in number z=to!number(0.0))
    {
        this.x = x; this.y = y; this.z = z;
    }

    version(complex_numbers) {
        @nogc this(in double[] p)
            {
                switch ( p.length ) {
                case 0: x = y = z = 0.0; break;
                case 1: x = p[0]; y = z = 0.0; break;
                case 2: x = p[0]; y = p[1]; z = 0.0; break;
                default: x = p[0]; y = p[1]; z = p[2]; break;
                }
            }

        @nogc this(in double x, in double y=0.0, in double z=0.0)
            {
                this.x = x; this.y = y; this.z = z;
            }
    } // end version complex_numbers

    @nogc this(const(Vector3) other)
    {
        x = other.x; y = other.y; z = other.z;
    }

    @property Vector3 dup() const
    {
        return Vector3(this);
    }

    @nogc ref Vector3 set(ref const(Vector3) other) return
    // Convenience function for setting the components of an existing object.
    {
        x = other.x; y = other.y; z = other.z;
        return this;
    }

    @nogc ref Vector3 set(Vector3* other) return
    // Convenience function for setting the components of an existing object.
    {
        x = other.x; y = other.y; z = other.z;
        return this;
    }

    @nogc ref Vector3 set(number x, number y, number z=to!number(0.0)) return
    // Convenience function for setting the components of an existing object.
    // Note that we may supply just the x,y coordinates.
    {
        this.x = x; this.y = y; this.z = z;
        return this;
    }

    @nogc ref Vector3 set(number[] xyz) return
    {
        this.x = xyz[0]; this.y = xyz[1]; this.z = xyz[2];
        return this;
    }

    version(complex_numbers) {
        // We want to retain the flavour with double numbers.
        @nogc ref Vector3 set(double x, double y, double z=0.0) return
        // Convenience function for setting the components of an existing object.
        // Note that we may supply just the x,y coordinates.
        {
            this.x = x; this.y = y; this.z = z;
            return this;
        }

        @nogc ref Vector3 set(double[] xyz) return
        {
            this.x = xyz[0]; this.y = xyz[1]; this.z = xyz[2];
            return this;
        }

        @nogc ref Vector3 clear_imaginary_components() return
        // Convenience function for setting-to-zero the imaginary components of an existing object.
        {
            x.im = 0.0; y.im = 0.0; z.im = 0.0;
            return this;
        }
    } // end version complex_numbers

    @nogc ref Vector3 clear() return
    // Convenience function for setting-to-zero the components of an existing object.
    {
        x = 0.0; y = 0.0; z = 0.0;
        return this;
    }

    @nogc ref Vector3 add(number x, number y, number z, number factor) return
    // Convenience function for adding the components of an existing object.
    // This avoids the temporary associated with += (below)
    {
        x += factor*x; y += factor*y; z += factor*z;
        return this;
    }

    @nogc ref Vector3 add(ref const(Vector3) other) return
    // Convenience function for adding the components of an existing object.
    // This avoids the temporary associated with += (below)
    {
        x += other.x; y += other.y; z += other.z;
        return this;
    }

    @nogc ref Vector3 add(Vector3* other) return
    // Convenience function for adding the components of an existing object.
    // This avoids the temporary associated with += (below)
    {
        x += other.x; y += other.y; z += other.z;
        return this;
    }

    @nogc ref Vector3 add(ref const(Vector3) other, number factor) return
    // Convenience function for adding the components of an existing object, scaled.
    // This avoids the temporary associated with += (below)
    {
        x += other.x*factor; y += other.y*factor; z += other.z*factor;
        return this;
    }

    @nogc ref Vector3 add(Vector3* other, number factor) return
    // Convenience function for adding the components of an existing object, scaled.
    // This avoids the temporary associated with += (below)
    {
        x += other.x*factor; y += other.y*factor; z += other.z*factor;
        return this;
    }

    @nogc ref Vector3 scale(number factor) return
    // Convenience function for scaling the components of an existing object.
    // This avoids the temporary associated with *= (below)
    {
        x *= factor; y *= factor; z *= factor;
        return this;
    }

    version(complex_numbers) {
        // We want to retain the flavour with double numbers.

        @nogc ref Vector3 add(double x, double y, double z, double factor) return
        // Convenience function for adding the components of an existing object, scaled.
        // This avoids the temporary associated with += (below)
        {
            x += factor*x; y += factor*y; z += factor*z;
            return this;
        }

        @nogc ref Vector3 add(ref const(Vector3) other, double factor) return
        // Convenience function for adding the components of an existing object, scaled.
        // This avoids the temporary associated with += (below)
        {
            x += other.x*factor; y += other.y*factor; z += other.z*factor;
            return this;
        }

        @nogc ref Vector3 add(Vector3* other, double factor) return
        // Convenience function for adding the components of an existing object, scaled.
        // This avoids the temporary associated with += (below)
        {
            x += other.x*factor; y += other.y*factor; z += other.z*factor;
            return this;
        }

        @nogc ref Vector3 scale(double factor) return
        // Convenience function for scaling the components of an existing object.
        // This avoids the temporary associated with *= (below)
        {
            x *= factor; y *= factor; z *= factor;
            return this;
        }
    } // end version complex_numbers

    string toString() const
    {
        return format("Vector3(%s, %s, %s)", to!string(x), to!string(y), to!string(z));
    }

    // Some operators, at least those that make sense.
    Vector3 opUnary(string op)()
        if (op == "+")
    {
        return Vector3(x, y, z);
    }

    Vector3 opUnary(string op)()
        if (op == "-")
    {
        return Vector3(-x, -y, -z);
    }

    Vector3 opBinary(string op)(in Vector3 rhs) const
        if (op == "+")
    {
        return Vector3(x+rhs.x, y+rhs.y, z+rhs.z);
    }

    Vector3 opBinary(string op)(in Vector3 rhs) const
        if (op == "-")
    {
        return Vector3(x-rhs.x, y-rhs.y, z-rhs.z);
    }

    Vector3 opBinary(string op)(in number rhs) const
        if (op == "*")
    {
        return Vector3(x*rhs, y*rhs, z*rhs);
    }

    version(complex_numbers) {
        // Retain the double version.
        Vector3 opBinary(string op)(in double rhs) const
            if (op == "*")
        {
            return Vector3(x*rhs, y*rhs, z*rhs);
        }
    } // end version complex_numbers

    Vector3 opBinaryRight(string op)(in number lhs) const
        if (op == "*")
    {
        return Vector3(x*lhs, y*lhs, z*lhs);
    }

    version(complex_numbers) {
        // Retain the double version.
        Vector3 opBinaryRight(string op)(in double lhs) const
            if (op == "*")
        {
            return Vector3(x*lhs, y*lhs, z*lhs);
        }
    } // end version complex_numbers

    Vector3 opBinary(string op)(in number rhs) const
        if (op == "/")
    {
        return Vector3(x/rhs, y/rhs, z/rhs);
    }

    version(complex_numbers) {
        // Retain the double version.
        Vector3 opBinary(string op)(in double rhs) const
            if (op == "/")
        {
            return Vector3(x/rhs, y/rhs, z/rhs);
        }
    } // end version complex_numbers

    // Assignment operators. (Alexandrescu Section 7.1.5.1)
    @nogc void opAssign(in Vector3 rhs)
    {
        x = rhs.x; y = rhs.y; z = rhs.z;
    }

    // Combined assignment operators do change the original object.
    @nogc void opOpAssign(string op)(in Vector3 rhs)
        if (op == "+")
    {
        x += rhs.x; y += rhs.y; z += rhs.z;
    }

    @nogc void opOpAssign(string op)(in Vector3 rhs)
        if (op == "-")
    {
        x -= rhs.x; y -= rhs.y; z -= rhs.z;
    }

    @nogc void opOpAssign(string op)(in number rhs)
        if (op == "*")
    {
        x *= rhs; y *= rhs; z *= rhs;
    }

    @nogc void opOpAssign(string op)(in number rhs)
        if (op == "/")
    {
        x /= rhs; y /= rhs; z /= rhs;
    }

    version(complex_numbers) {
        // Retain the double version.
        @nogc void opOpAssign(string op)(in double rhs)
            if (op == "*")
        {
            x *= rhs; y *= rhs; z *= rhs;
        }

        @nogc void opOpAssign(string op)(in double rhs)
            if (op == "/")
        {
            x /= rhs; y /= rhs; z /= rhs;
        }
    } // end version complex_numbers

    // Other vector-specific operations.

    /**
     * Scales the vector to unit magnitude.
     */
    @nogc void normalize()
    {
        number magnitude = sqrt(this.dot(this));
        if (magnitude > 0.0) {
            x /= magnitude; y /= magnitude; z /= magnitude;
        } else {
            // Clean up, in case dot() underflows.
            x = y = z = 0.0;
        }
        // Flush small components to zero.
        const double small = 1.0e-30;
        version(complex_numbers) {
            if (fabs(x.re) < small && fabs(x.im) < small) { x = 0.0; }
            if (fabs(y.re) < small && fabs(y.im) < small) { y = 0.0; }
            if (fabs(z.re) < small && fabs(z.im) < small) { z = 0.0; }
        } else {
            if (fabs(x) < small) { x = 0.0; }
            if (fabs(y) < small) { y = 0.0; }
            if (fabs(z) < small) { z = 0.0; }
        }
    }

    @nogc number dot(ref const(Vector3) other) const
    {
        return x*other.x + y*other.y + z*other.z;
    }

    // Transform functions used to reorient vector values in the CFD codes.

    /**
     * Rotate v from the global xyz coordinate system into the local frame
     * defined by the orthogonal unit vectors n,t1,t2.
     *
     * We assume, without checking, that these vectors do nicely define
     * such a local system.
     */
    @nogc void transform_to_local_frame(ref const(Vector3) n,
                                        ref const(Vector3) t1,
                                        ref const(Vector3) t2)
    {
        number v_x = this.dot(n); // normal component
        number v_y = this.dot(t1); // tangential component 1
        number v_z = this.dot(t2); // tangential component 2
        x = v_x; y = v_y; z = v_z;
    }

    /**
     * Rotate v back into the global (xyz) coordinate system.
     */
    @nogc void transform_to_global_frame(ref const(Vector3) n,
                                         ref const(Vector3) t1,
                                         ref const(Vector3) t2)
    {
        number v_x = x*n.x + y*t1.x + z*t2.x; // global-x
        number v_y = x*n.y + y*t1.y + z*t2.y; // global-y
        number v_z = x*n.z + y*t1.z + z*t2.z; // global-z
        x = v_x; y = v_y; z = v_z;
    }

    // 2D flavour for change of coordinate system functions.

    @nogc void transform_to_local_frame(ref const(Vector3) n,
                                        ref const(Vector3) t1)
    {
        number v_x = x*n.x + y*n.y;   // normal component
        number v_y = x*t1.x + y*t1.y; // tangential component 1
        x = v_x; y = v_y; z = 0.0;
    }

    /**
     * Rotate v back into the global (xy) coordinate system.
     */
    @nogc void transform_to_global_frame(ref const(Vector3) n,
                                         ref const(Vector3) t1)
    {
        number v_x = x*n.x + y*t1.x; // global-x
        number v_y = x*n.y + y*t1.y; // global-y
        x = v_x; y = v_y; z = 0.0;
    }

    // Change of coordinate system; rotation with translation.

    // Transform coordinates from global frame to local (dash) frame.
    // Local frame is defined by unit vectors (n, t1 and t2) at location c.
    @nogc void transform_to_local_frame(ref const(Vector3) n,
                                        ref const(Vector3) t1,
                                        ref const(Vector3) t2,
                                        ref const(Vector3) c)
    {
        x -= c.x; y -= c.y; z -= c.z; // shift to local origin
        number v_x = this.dot(n); // normal component
        number v_y = this.dot(t1); // tangential component 1
        number v_z = this.dot(t2); // tangential component 2
        x = v_x; y = v_y; z = v_z;
    }

    /**
     * Rotate v back into the global (xyz) coordinate system.
     */
    @nogc void transform_to_global_frame(ref const(Vector3) n,
                                         ref const(Vector3) t1,
                                         ref const(Vector3) t2,
                                         ref const(Vector3) c)
    {
        number v_x = x*n.x + y*t1.x + z*t2.x + c.x; // global-x
        number v_y = x*n.y + y*t1.y + z*t2.y + c.y; // global-y
        number v_z = x*n.z + y*t1.z + z*t2.z + c.z; // global-z
        x = v_x; y = v_y; z = v_z;
    }

    /**
     * General matrix transformation (used when rotating flowstate vectors).
     */
    @nogc void apply_matrix_transform(ref const(number[]) Rmatrix)
    {
        // Write out the matrix multiplication, long-hand.
        number old_x = x; number old_y = y; number old_z = z;
        x = Rmatrix[0]*old_x + Rmatrix[1]*old_y + Rmatrix[2]*old_z;
        y = Rmatrix[3]*old_x + Rmatrix[4]*old_y + Rmatrix[5]*old_z;
        z = Rmatrix[6]*old_x + Rmatrix[7]*old_y + Rmatrix[8]*old_z;
    }
    version(complex_numbers) {
        // Retain the flavour with double numbers in the matrix.
        @nogc void apply_matrix_transform(ref const(double[]) Rmatrix)
        {
            // Write out the matrix multiplication, long-hand.
            number old_x = x; number old_y = y; number old_z = z;
            x = Rmatrix[0]*old_x + Rmatrix[1]*old_y + Rmatrix[2]*old_z;
            y = Rmatrix[3]*old_x + Rmatrix[4]*old_y + Rmatrix[5]*old_z;
            z = Rmatrix[6]*old_x + Rmatrix[7]*old_y + Rmatrix[8]*old_z;
        }
    }

    /**
     * Compute mirror-image location for plane defined by point and normal.
     */
    @nogc ref Vector3 mirror_image(ref const(Vector3) point,
                                   ref const(Vector3) normal) return
    {
        Vector3 n = Vector3(normal.x, normal.y, normal.z); n.normalize();
        // Construct tangents to the plane.
        Vector3 different = n + Vector3(1.0, 1.0, 1.0);
        Vector3 t1; cross(t1, n, different); t1.normalize();
        Vector3 t2; cross(t2, n, t1); t2.normalize();
        // Mirror image the vector in a frame local to the plane.
        transform_to_local_frame(n, t1, t2, point);
        x = -x;
        transform_to_global_frame(n, t1, t2, point);
        return this;
    }

    /**
     * Rotate point about the z-axis by angle dtheta, in radians.
     */
    @nogc ref Vector3 rotate_about_zaxis(double dtheta) return
    {
        double theta = atan2(y.re, x.re) + dtheta;
        double r = sqrt(x.re*x.re + y.re*y.re);
        x = r * cos(theta);
        y = r * sin(theta);
        return this;
    }

    /**
     * Alternative implementation for rotation in (x,y)-plane.
     */
    @nogc ref Vector3 rotate2d(double dtheta) return
    {
        double sn = sin(dtheta);
        double cs = cos(dtheta);
        number _x = x*cs - y*sn;
        number _y = y*cs + x*sn;
        x=_x; y=_y;
        return this;
    }
} // end class Vector3

/**
 * Returns the distance between two points.
 */
@nogc
double distance_between(ref const(Vector3) v1, ref const(Vector3) v2)
{
    number d = sqrt((v1.x-v2.x)^^2 + (v1.y-v2.y)^^2 + (v1.z-v2.z)^^2);
    return d.re;
}

/**
 * Returns the scalar dot product of two vectors.
 */
@nogc
number dot(ref const(Vector3) v1, ref const(Vector3) v2)
{
    number result = 0.0;
    // Maybe we should be careful with underflow and overflow...
    result = v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
    return result;
}

/**
 * Returns magnitude of the vector.
 */
@nogc
number abs(ref const(Vector3) v)
{
    return sqrt(v.dot(v));
}

/**
 * Returns a unit vector in the same direction as v.
 */
Vector3 unit(ref const(Vector3) v)
{
    Vector3 v2 = Vector3(v);
    v2.normalize();
    return v2;
}

/**
 * Vector cross product for use in a single statement that will not make temporaries.
 */
@nogc
void cross(ref Vector3 v3, ref const(Vector3) v1, ref const(Vector3) v2)
{
    v3.x = v1.y * v2.z - v2.y * v1.z;
    v3.y = v2.x * v1.z - v1.x * v2.z;
    v3.z = v1.x * v2.y - v2.x * v1.y;
}

/**
 * Vector cross product for use in Vector3 expressions.
 * We need to keep the "in" qualifiers.
 */
@nogc
Vector3 cross(in Vector3 v1, in Vector3 v2)
{
    Vector3 v3;
    cross(v3, v1, v2);
    return v3;
}

/**
 * Component forms
 */

@nogc
number dot_product(number ax, number ay, number az, number bx, number by, number bz)
{
    return ax*bx + ay*by + az*bz;
}

@nogc
void cross_product(number ax, number ay, number az, number bx, number by, number bz,
                   ref number cx, ref number cy, ref number cz)
{
    cx = ay*bz - az*by;
    cy = az*bx - ax*bz;
    cz = ax*by - ay*bx;
    return;
}

@nogc number dot(number ax, number ay, number az, ref const(Vector3) other)
{
    return ax*other.x + ay*other.y + az*other.z;
}

// Transform functions used to reorient vector values in the CFD codes.

/**
 * Rotate v from the global xyz coordinate system into the local frame
 * defined by the orthogonal unit vectors n,t1,t2.
 *
 * We assume, without checking, that these vectors do nicely define
 * such a local system.
 */
@nogc void transform_to_local_frame(ref number ax, ref number ay, ref number az,
                                    ref const(Vector3) n,
                                    ref const(Vector3) t1,
                                    ref const(Vector3) t2)
{
    number v_x = dot(ax, ay, az, n); // normal component
    number v_y = dot(ax, ay, az, t1); // tangential component 1
    number v_z = dot(ax, ay, az, t2); // tangential component 2
    ax = v_x; ay = v_y; az = v_z;
}

/**
 * Rotate v back into the global (xyz) coordinate system.
 */
@nogc void transform_to_global_frame(ref number ax, ref number ay, ref number az,
                                     ref const(Vector3) n,
                                     ref const(Vector3) t1,
                                     ref const(Vector3) t2)
{
    number v_x = ax*n.x + ay*t1.x + az*t2.x; // global-x
    number v_y = ax*n.y + ay*t1.y + az*t2.y; // global-y
    number v_z = ax*n.z + ay*t1.z + az*t2.z; // global-z
    ax = v_x; ay = v_y; az = v_z;
}

/**
 * Returns true if all of the components of two vectors are approximately equal.
 */
@nogc
bool approxEqualVectors(in Vector3 v1, in Vector3 v2,
                        double maxRelDiff=1.0e-2, double maxAbsDiff=1.0e-5)
{
    return (approxEqualNumbers(v1.x, v2.x, maxRelDiff, maxAbsDiff) &&
            approxEqualNumbers(v1.y, v2.y, maxRelDiff, maxAbsDiff) &&
            approxEqualNumbers(v1.z, v2.z, maxRelDiff, maxAbsDiff));
}

version(vector3_test) {
    import util.msg_service;
    int main() {
        // Check that we have separate data with the correct values.
        Vector3 a = Vector3([1.0, 2.2, 3.0]);
        Vector3 b = Vector3(1.0);
        assert(a.x == 1.0, failedUnitTest());
        assert(a.y == 2.2, failedUnitTest());
        assert(a.z == 3.0, failedUnitTest());
        assert(a.x == b.x, failedUnitTest());
        assert(b.y == 0.0, failedUnitTest());
        assert(b.z == 0.0, failedUnitTest());
        b.set(a);
        assert(a.x == b.x && a.y == b.y && a.z == b.z, failedUnitTest());
        b.set(1.0, 0.0, 0.0);
        assert(b.x == 1.0 && b.y == 0.0 && b.z == 0.0, failedUnitTest());

        // Check operators
        b = -a;
        assert(b.x == -a.x && b.y == -a.y && b.z == -a.z, failedUnitTest());

        b = Vector3(1.0);
        Vector3 c = a + b;
        assert(c.y == a.y+b.y, failedUnitTest());
        c = a - b;
        assert(c.y == a.y-b.y, failedUnitTest());
        Vector3 d = a.dup;
        a.y = 99.0;
        assert(a.y == 99.0 && d.y == 2.2, failedUnitTest());
        Vector3 d2 = a;
        a.y = 3.3;
        assert(a.y == 3.3 && d2.y == 99.0, failedUnitTest());

        Vector3 e = a * 2.0;
        Vector3 f = 3.0 * d;
        assert(e.z == 6.0 && f.z == 9.0, failedUnitTest());
        Vector3 g = d / 3.0;
        assert(g.z == 1.0, failedUnitTest());

        g += f;
        assert(g.z == 10.0, failedUnitTest());
        g /= 2.0;
        assert(g.z == 5.0, failedUnitTest());

        a = Vector3(1.0, 0.0, 0.0);
        a.rotate_about_zaxis(PI/4);
        assert(approxEqualVectors(a, Vector3(0.7071, 0.7071, 0)), failedUnitTest());

        a = Vector3(1.0, 0.0, 0.0);
        Vector3 point = Vector3(0.0, 1.0, 0.0);
        Vector3 normal = Vector3(0.0, 1.0, 0.0);
        a.mirror_image(point, normal);
        assert(approxEqualVectors(a, Vector3(1.0, 2.0, 0)), failedUnitTest());

        Vector3 u = unit(g);
        assert(approxEqualNumbers(abs(u), to!number(1.0)), failedUnitTest());

        Vector3 x = Vector3(1.0, 0.0, 0.0);
        Vector3 y = Vector3(0.0, 1.0, 0.0);
        Vector3 z = cross(x,y);
        Vector3 zref = Vector3(0.0,0.0,1.0);
        assert(approxEqualVectors(z, zref), failedUnitTest());

        Vector3 n = Vector3(1.0,1.0,0.0); n = unit(n);
        Vector3 t1 = Vector3(-1.0,1.0,0.0); t1 = unit(t1);
        Vector3 t2 = cross(n, t1);
        Vector3 h = Vector3(1.0,0.0,1.0);
        Vector3 h_ref = Vector3(h);
        h.transform_to_local_frame(n, t1, t2);
        assert(approxEqualVectors(h, Vector3(sqrt(1.0/2.0), -sqrt(1.0/2.0), 1.0)),
               failedUnitTest());
        h.transform_to_global_frame(n, t1, t2);
        assert(approxEqualVectors(h, h_ref), failedUnitTest());
        //
        // 2D variants, starting with fresh values.
        n.set(1.0,1.0); n = unit(n);
        t1.set(-1.0,1.0); t1 = unit(t1);
        h_ref.set(1.0,0.0);
        h.set(h_ref);
        h.transform_to_local_frame(n, t1);
        assert(approxEqualVectors(h, Vector3(sqrt(1.0/2.0), -sqrt(1.0/2.0))),
               failedUnitTest());
        h.transform_to_global_frame(n, t1);
        assert(approxEqualVectors(h, h_ref), failedUnitTest());

        Vector3 a45 = Vector3(cos(to!number(PI)/4),sin(to!number(PI)/4));
        Vector3 a60 = Vector3(cos(to!number(PI)/3),sin(to!number(PI)/3));
        assert(approxEqualVectors(a45.rotate2d(15.0*PI/180), a60), failedUnitTest());
        Vector3 a30 = Vector3(cos(to!number(PI)/6),sin(to!number(PI)/6));
        assert(approxEqualVectors(a30.rotate2d(30.0*PI/180), a60), failedUnitTest());

        return 0;
    }
} // end vector3_test

Vector3 getJSONVector3(JSONValue jsonData, string key, Vector3 defaultValue)
// Read a Vector3 value as an array of 3 floating-point values.
{
    Vector3 value;
    try {
        auto json_values = jsonData[key].array;
        foreach (i, json_val; json_values) {
            switch (i) {
            case 0: value.x = to!double(json_val.floating); break;
            case 1: value.y = to!double(json_val.floating); break;
            case 2: value.z = to!double(json_val.floating); break;
            default:
            }
        }
    } catch (Exception e) {
        value = defaultValue;
    }
    return value;
} // end getJSONVector3()
