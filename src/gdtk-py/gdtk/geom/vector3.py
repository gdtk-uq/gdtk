# vector3.py
"""
A vector class for working in 2 or 3 dimensions.

Peter Jacobs and Rowan Gollan.

Versions:
  2020-02-05 Build enough to make Points and Paths
  2020-07-30 Integrate minimal_geometry from cfcfd3 project.
  2022-10-07 Arrayification by NNG, (Newcastle Interchange Line, NSW)
"""

import numbers
import numpy as np

class Vector3(object):
    """
    A vector with 3 Cartesian coordinates.
    """
    __slots__ = ['x', 'y', 'z']

    def __init__(self, x, y=None, z=None):
        """
        Accept the coordinates as a list of numbers,
        a dictionary of named numbers, or as individual numbers.

        For example:
        >>> from gdtk.geom.vector3 import Vector3
        >>> p0 = Vector3(x=1.0, y=2.0, z=3.0)
        >>> p1 = Vector3(1.0, 2.0, 3.0)
        >>> p2 = Vector3([1.0, 2.0, 3.0])
        >>> p3 = Vector3({'x':1.0, 'y':2.0, 'z':3.0})
        >>> p4 = Vector3(p3)
        """
        if isinstance(x, list) or isinstance(x, tuple):
            self.x = x[0]
            if len(x) > 1:
                self.y = x[1]
                if len(x) > 2:
                    self.z = x[2]
                else:
                    self.z = 0.0*self.y
            else:
                raise Exception("Need to supply at least 2 coordinates.")
        elif isinstance(x, dict):
            self.x = x["x"] if 'x' in list(x.keys()) else 0.0
            self.y = x["y"] if 'y' in list(x.keys()) else 0.0
            self.z = x["z"] if 'z' in list(x.keys()) else 0.0
        elif isinstance(x, Vector3):
            self.x = x.x
            self.y = x.y
            self.z = x.z
        elif isinstance(x, np.ndarray):
            self.x = x
            if isinstance(y, numbers.Real):
                self.y = x*0.0 + y
            else:
                self.y = y
            if isinstance(z, numbers.Real):
                self.z = x*0.0 + z
            else:
                self.z = z
        else:
            # Presume that we have been given at least 2 numbers.
            self.x = x
            if y is None: raise Exception("Expected a number for y coordinate.")
            self.y = y
            if z is None:
                self.z = 0.0*self.y
            else:
                self.z = z
        args_are_not_real = not(isinstance(self.x, numbers.Real)
                                and isinstance(self.y, numbers.Real)
                                and isinstance(self.z, numbers.Real))
        args_are_not_array= not(isinstance(self.x, np.ndarray)
                                and isinstance(self.y, np.ndarray)
                                and isinstance(self.z, np.ndarray))
        if args_are_not_real and args_are_not_array:
            print(self.x, self.y, self.z, type(self.x), type( self.y), type(self.z))
            raise Exception("Elements should be Real numbers or arrays.")
        return

    def __repr__(self):
        "Returns string representation."
        return f"Vector3(x={self.x}, y={self.y}, z={self.z})"

    def __abs__(self):
        "Returns magnitude."
        return np.sqrt(self.x**2 + self.y**2 + self.z**2)

    def __pos__(self):
        "Returns copy of self."
        return Vector3(self)

    def __neg__(self):
        "Returns negative copy of self."
        return Vector3(-self.x, -self.y, -self.z)

    def __add__(self, other):
        "Returns new Vector3 for self+other"
        return Vector3(self.x+other.x, self.y+other.y, self.z+other.z)

    def __iadd__(self, other):
        "Returns self += other Vector3."
        self.x += other.x; self.y += other.y; self.z += other.z
        return self

    def __sub__(self, other):
        "Returns new Vector3 self-other."
        return Vector3(self.x-other.x, self.y-other.y, self.z-other.z)

    def __isub__(self, other):
        "self -= other Vector3"
        self.x -= other.x; self.y -= other.y; self.z -= other.z
        return self

    def __mul__(self, other):
        "Returns new Vector3 self*number"
        if isinstance(other, numbers.Real) or isinstance(other, np.ndarray):
            return Vector3(self.x*other, self.y*other, self.z*other)
        else:
            return NotImplemented

    def __rmul__(self, other):
        "Returns new Vector3 number*self."
        return self * other

    def __imul__(self, other):
        "Returns self *= other number."
        if isinstance(other, numbers.Real) or isinstance(other, np.ndarray):
            self.x *= other; self.y *= other; self.z *= other
            return self
        else:
            return NotImplemented

    def __truediv__(self, other):
        "Returns new Vector3 self/number"
        if isinstance(other, numbers.Real) or isinstance(other, np.ndarray):
            return Vector3(self.x/other, self.y/other, self.z/other)
        else:
            return NotImplemented

    def __itruediv__(self, other):
        "Returns self /= other number."
        if isinstance(other, numbers.Real) or isinstance(other, np.ndarray):
            self.x /= other; self.y /= other; self.z /= other
            return self
        else:
            return NotImplemented

    def unit(self):
        "Returns new unit vector."
        mag = abs(self)
        if np.any(mag == 0.0):
            raise ValueError("Zero magnitude vector has no defined direction.")
        return Vector3(self.x/mag, self.y/mag, self.z/mag)

    def normalize(self):
        "Scales self to unit magnitude."
        mag = abs(self)
        self /= mag
        return

    def dot(self, other):
        "Returns dot product of self with other Vector3."
        if isinstance(other, Vector3):
            return self.x*other.x + self.y*other.y + self.z*other.z
        else:
            raise Exception("dot() not implemented for {}".format(type(other)))

    def transform_to_local_frame(self, n, t1, t2, c=None):
        """
        Change coordinates into the local right-handed (RH) system at point c.

        We trust that n, t1 and t2 are normalized and for a correct RH system.
        """
        if c is not None: self -= c
        newX = self.dot(n); newY = self.dot(t1); newZ = self.dot(t2)
        self.x = newX; self.y = newY; self.z = newZ
        return self

    def transform_to_global_frame(self, n, t1, t2, c=None):
        """
        Change coordinates out of the local right-handed (RH) system at point c.

        We trust that n, t1 and t2 are normalized and for a correct RH system.
        """
        newX = self.x*n.x + self.y*t1.x + self.z*t2.x
        newY = self.x*n.y + self.y*t1.y + self.y*t2.y
        newZ = self.x*n.z + self.y*t1.z + self.z*t2.z
        if c is not None:
            newX += c.x; newY += c.y; newZ += c.z
        self.x = newX; self.y = newY; self.z = newZ
        return self

    def rotate_about_zaxis(self, dtheta):
        """
        Rebuild of a method in the d-side geometry library.
        """
        theta = np.arctan2(self.y, self.x) + dtheta;
        r = np.sqrt(self.x*self.x + self.y*self.y);
        self.x = r * np.cos(theta);
        self.y = r * np.sin(theta);
        return self;


    # ------- end class Vector3 ---------

VERY_SMALL_MAGNITUDE = 1.0e-200

def approxEqualVectors(a, b, rel_tol=1.0e-2, abs_tol=1.0e-5):
    """
    Test that all components are close.
    """
    return np.all([np.isclose(a.x, b.x, rtol=rel_tol, atol=abs_tol),
                np.isclose(a.y, b.y, rtol=rel_tol, atol=abs_tol),
                np.isclose(a.z, b.z, rtol=rel_tol, atol=abs_tol)])

def cross(a, b):
    """
    Returns Vector3 cross product.
    """
    ab_x = a.y*b.z - b.y*a.z
    ab_y = b.x*a.z - a.x*b.z
    ab_z = a.x*b.y - b.x*a.y
    return Vector3(ab_x, ab_y, ab_z)

def dot(a, b):
    "Returns dot product."
    return a.dot(b)

def unit(a):
    "Returns a new unit vector."
    return a.unit()

def quad_properties(p0, p1, p2, p3):
    "Returns centroid, quadrilateral-defining unit vectors, and area."
    vector_area = 0.25 * cross(p1-p0+p2-p3, p3-p0+p2-p1)
    n = vector_area.unit()
    area = abs(vector_area)
    t1 = ((p1-p0)+(p2-p3)).unit() # Works even if one edge has zero length.
    t2 = cross(n, t1).unit() # Using unit() to tighten up on magnitude.
    centroid = 0.25 * (p0 + p1 + p2 + p3)
    return centroid, n, t1, t2, area

def quad_centroid(p0, p1, p2, p3):
    "Returns centroid of quadrilateral."
    centroid, n, t1, t2, area = quad_properties(p0, p1, p2, p3)
    return centroid

def quad_area(p0, p1, p2, p3):
    "Returns area for quadrilateral."
    centroid, n, t1, t2, area = quad_properties(p0, p1, p2, p3)
    return area

def quad_normal(p0, p1, p2, p3):
    "Returns unit normal for quadrilateral."
    centroid, n, t1, t2, area = quad_properties(p0, p1, p2, p3)
    return n

def tetrahedron_properties(p0, p1, p2, p3):
    "Returns centroid and volume of tetrahedron."
    volume = dot(p3-p0, cross(p1-p0, p2-p0))/6.0
    centroid = 0.25 * (p0 + p1 + p2 + p3)
    return centroid, volume

def wedge_properties(p0, p1, p2, p3, p4, p5):
    "Returns centroid and volume for wedge."
    c1, v1 = tetrahedron_properties(p0, p4, p5, p3)
    c2, v2 = tetrahedron_properties(p0, p5, p4, p1)
    c3, v3 = tetrahedron_properties(p0, p1, p2, p5)
    volume = v1 + v2 + v3
    # TODO: It's hard to make this array agnostic... (NNG)
    #if volume < VERY_SMALL_MAGNITUDE:
    #    #print "Warning wedge_properties():"
    #    #print "Very small or negative volume: ", volume
    #    #print "Setting volume to zero."
    #    volume = 0.0
    #    centroid = (c1 + c2 + c3)/3.0
    #    return centroid, volume
    #
    centroid = (c1*v1 + c2*v2 + c3*v3)/volume
    return centroid, volume

def hexahedron_properties(p0, p1, p2, p3, p4, p5, p6, p7):
    "Returns centroid and volume for hexahedron."
    c1, v1 = wedge_properties(p0, p1, p2, p4, p5, p6)
    c2, v2 = wedge_properties(p0, p2, p3, p4, p6, p7)
    volume = v1 + v2
    # TODO: It's hard to make this array agnostic... (NNG)
    #small = volume < VERY_SMALL_MAGNITUDE
    centroid = (c1*v1 + c2*v2)/volume
    return centroid, volume

def hexahedron_volume(p0, p1, p2, p3, p4, p5, p6, p7):
    "Returns volume for hexahedron."
    c, v = hexahedron_properties(p0, p1, p2, p3, p4, p5, p6, p7)
    return v

if __name__=='__main__':
    p0 = Vector3(x=1.0, y=2.0, z=3.0)
    p1 = Vector3(1.0, 2.0, 3.0)
    p2 = Vector3([1.0, 2.0, 3.0])
    p3 = Vector3({'x':1.0, 'y':2.0, 'z':3.0})
    p4 = Vector3(p3)

    p5 = Vector3(x=np.array([1.0, 1.0, 1.0, 1.0]),
                 y=np.array([2.0, 2.0, 2.0, 2.0]),
                 z=np.array([3.0, 3.0, 3.0, 3.0]))

    assert(np.isclose(abs(p0), abs(p5)[0]))
    assert(np.isclose((p0+p0).x, (p5+p5).x[0]))
    assert(np.isclose((p0+p0).y, (p5+p5).y[0]))
    assert(np.isclose((p0+p0).z, (p5+p5).z[0]))

    assert(np.isclose((p0+p0).unit().x, (p5+p5).unit().x[0]))
    assert(np.isclose((p0+p0).unit().y, (p5+p5).unit().y[0]))
    assert(np.isclose((p0+p0).unit().z, (p5+p5).unit().z[0]))

    assert(np.isclose(p0.dot(p0), p5.dot(p5)[0]))

    assert(approxEqualVectors(p5,p5))

    b0 = Vector3(x=3.0, y=2.0, z=1.0)
    b5 = Vector3(x=np.array([3.0, 3.0, 3.0, 3.0]),
                 y=np.array([2.0, 2.0, 2.0, 2.0]),
                 z=np.array([1.0, 1.0, 1.0, 1.0]))
    pb0 = cross(p0, b0)
    pb5 = cross(p5, b5)
    assert(np.isclose((pb0).x, (pb5).x[0]))
    assert(np.isclose((pb0).y, (pb5).y[0]))
    assert(np.isclose((pb0).z, (pb5).z[0]))

    p0 = Vector3(x=0.0, y=0.0)
    p1 = Vector3(x=1.0, y=0.0)
    p2 = Vector3(x=1.0, y=1.0)
    p3 = Vector3(x=0.0, y=1.0)
    centroid, n, t1, t2, area = quad_properties(p0, p1, p2, p3)

    q0 = Vector3(x=np.array([0.0, 0.0]), y=np.array([0.0, 0.0]))
    q1 = Vector3(x=np.array([1.0, 1.0]), y=np.array([0.0, 0.0]))
    q2 = Vector3(x=np.array([1.0, 1.0]), y=np.array([1.0, 1.0]))
    q3 = Vector3(x=np.array([0.0, 0.0]), y=np.array([1.0, 1.0]))
    qcentroid, qn, qt1, qt2, qarea = quad_properties(q0, q1, q2, q3)
    assert(np.isclose((centroid).x, (qcentroid).x[0]))
    assert(np.isclose((centroid).y, (qcentroid).y[0]))
    assert(np.isclose((centroid).z, (qcentroid).z[0]))
    assert(np.isclose(area, qarea[0]))
    assert(np.isclose(area, qarea[0]))
    assert(np.isclose(area, qarea[0]))

    p0 = Vector3(x=0.0, y=0.0, z=0.0)
    p1 = Vector3(x=1.0, y=0.0, z=0.0)
    p2 = Vector3(x=1.0, y=1.0, z=0.0)
    p3 = Vector3(x=0.0, y=0.0, z=1.0)
    p4 = Vector3(x=1.0, y=0.0, z=1.0)
    p5 = Vector3(x=1.0, y=1.0, z=1.0)
    centroid, volume = wedge_properties(p0, p1, p2, p3, p4, p5)

    q0 = Vector3(x=np.array([0.0,0.0]), y=np.array([0.0,0.0]), z=np.array([0.0,0.0])) 
    q1 = Vector3(x=np.array([1.0,1.0]), y=np.array([0.0,0.0]), z=np.array([0.0,0.0]))
    q2 = Vector3(x=np.array([1.0,1.0]), y=np.array([1.0,1.0]), z=np.array([0.0,0.0]))
    q3 = Vector3(x=np.array([0.0,0.0]), y=np.array([0.0,0.0]), z=np.array([1.0,1.0]))
    q4 = Vector3(x=np.array([1.0,1.0]), y=np.array([0.0,0.0]), z=np.array([1.0,1.0]))
    q5 = Vector3(x=np.array([1.0,1.0]), y=np.array([1.0,1.0]), z=np.array([1.0,1.0]))
    qcentroid, qvolume = wedge_properties(q0, q1, q2, q3, q4, q5)
    assert(np.isclose((centroid).x, (qcentroid).x[0]))
    assert(np.isclose((centroid).y, (qcentroid).y[0]))
    assert(np.isclose((centroid).z, (qcentroid).z[0]))

    assert(np.isclose(volume, qvolume[0]))

    p0 = Vector3(x=0.0, y=0.0, z=0.0)
    p1 = Vector3(x=1.0, y=0.0, z=0.0)
    p2 = Vector3(x=1.0, y=1.0, z=0.0)
    p3 = Vector3(x=0.0, y=1.0, z=0.0)
    p4 = Vector3(x=0.0, y=0.0, z=1.0)
    p5 = Vector3(x=1.0, y=0.0, z=1.0)
    p6 = Vector3(x=1.0, y=1.0, z=1.0)
    p7 = Vector3(x=0.0, y=1.0, z=1.0)
    centroid, volume = hexahedron_properties(p0, p1, p2, p3, p4, p5, p6, p7)

    q0 = Vector3(x=np.array([0.0,0.0]), y=np.array([0.0,0.0]), z=np.array([0.0,0.0])) 
    q1 = Vector3(x=np.array([1.0,1.0]), y=np.array([0.0,0.0]), z=np.array([0.0,0.0]))
    q2 = Vector3(x=np.array([1.0,1.0]), y=np.array([1.0,1.0]), z=np.array([0.0,0.0]))
    q3 = Vector3(x=np.array([0.0,0.0]), y=np.array([1.0,1.0]), z=np.array([0.0,0.0]))
    q4 = Vector3(x=np.array([0.0,0.0]), y=np.array([0.0,0.0]), z=np.array([1.0,1.0]))
    q5 = Vector3(x=np.array([1.0,1.0]), y=np.array([0.0,0.0]), z=np.array([1.0,1.0]))
    q6 = Vector3(x=np.array([1.0,1.0]), y=np.array([1.0,1.0]), z=np.array([1.0,1.0]))
    q7 = Vector3(x=np.array([0.0,0.0]), y=np.array([1.0,1.0]), z=np.array([1.0,1.0]))
    qcentroid, qvolume = hexahedron_properties(q0, q1, q2, q3, q4, q5, q6, q7)
    assert(np.isclose((centroid).x, (qcentroid).x[0]))
    assert(np.isclose((centroid).y, (qcentroid).y[0]))
    assert(np.isclose((centroid).z, (qcentroid).z[0]))

