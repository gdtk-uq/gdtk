# vector3.py
"""
A vector class for working in 2 or 3 dimensions.

PJ, 2020-02-05
"""

import math
import numbers

class Vector3(object):
    """
    A vector with 3 Cartesian coordinates.
    """
    __slots__ = ['x', 'y', 'z']

    def __init__(self, x, y=None, z=None):
        """
        Accept the coordinates as a list of numbers,
        a dictionary of named numbers, or as individual numbers.
        """
        if isinstance(x, list) or isinstance(x, tuple):
            self.x = x[0]
            if len(x) > 1:
                self.y = x[1]
                if len(x) > 2:
                    self.z = x[2]
                else:
                    self.z = 0.0
            else:
                raise Exception("Need to supply at least 2 coordinates.")
        elif isinstance(x, dict):
            self.x = x["x"] if 'x' in x.keys() else 0.0
            self.y = x["y"] if 'y' in x.keys() else 0.0
            self.z = x["z"] if 'z' in x.keys() else 0.0
        elif isinstance(x, Vector3):
            self.x = x.x
            self.y = x.y
            self.z = x.z
        else:
            # Presume that we have been given at least 2 numbers.
            self.x = x
            if y == None: raise Exception("Expected a number for y coordinate.")
            self.y = y
            if z == None:
                self.z = 0.0
            else:
                self.z = z
        if not(isinstance(self.x, numbers.Real) and
               isinstance(self.y, numbers.Real) and
               isinstance(self.z, numbers.Real)):
            raise Exception("Elements should be Real numbers.")
        return

    def __repr__(self):
        return "Vector3(x={}, y={}, z={})".format(self.x, self.y, self.z)

    def __abs__(self):
        return math.sqrt(self.x**2 + self.y**2 + self.z**2)

    def __pos__(self):
        return self

    def __neg__(self):
        self.x = -self.x; self.y = -self.y; self.z = -self.z
        return self

    def __add__(self, other):
        return Vector3(self.x+other.x, self.y+other.y, self.z+other.z)

    def __iadd__(self, other):
        self.x += other.x; self.y += other.y; self.z += other.z
        return self

    def __sub__(self, other):
        return Vector3(self.x-other.x, self.y-other.y, self.z-other.z)

    def __isub__(self, other):
        self.x -= other.x; self.y -= other.y; self.z -= other.z
        return self

    def __mul__(self, other):
        if isinstance(other, numbers.Real):
            return Vector3(self.x*other, self.y*other, self.z*other)
        else:
            return NotImplemented

    def __rmul__(self, other):
        return self * other

    def __imul__(self, other):
        if isinstance(other, numbers.Real):
            self.x *= other; self.y *= other; self.z *= other
            return self
        else:
            return NotImplemented

    def __truediv__(self, other):
        if isinstance(other, numbers.Real):
            return Vector3(self.x/other, self.y/other, self.z/other)
        else:
            return NotImplemented

    def __itruediv__(self, other):
        if isinstance(other, numbers.Real):
            self.x /= other; self.y /= other; self.z /= other
            return self
        else:
            return NotImplemented

    def normalize(self):
        mag = abs(self)
        self /= mag
        return

    def dot(self, other):
        if isinstance(other, Vector3):
            return self.x*other.x + self.y*other.y + self.z*other.z
        else:
            raise Exception("dot() not implemented for {}".format(type(other)))

    def transform_to_local_frame(self, n, t1, t2, c=None):
        """
        Change coodinates into the local right-handed (RH) system at point c.

        We trust that n, t1 and t2 are normalized and for a correct RH system.
        """
        if c: self -= c
        newX = self.dot(n); newY = self.dot(t1); newZ = self.dot(t2)
        self.x = newX; self.y = newY; self.z = newZ
        return self

    def transform_to_global_frame(self, n, t1, t2, c=None):
        """
        Change coodinates out of the local right-handed (RH) system at point c.

        We trust that n, t1 and t2 are normalized and for a correct RH system.
        """
        newX = self.x*n.x + self.y*t1.x + self.z*t2.x
        newY = self.x*n.y + self.y*t1.y + self.y*t2.y
        newZ = self.x*n.z + self.y*t1.z + self.z*t2.z
        if c:
            newX += c.x; newY += c.y; newZ += c.z
        self.x = newX; self.y = newY; self.z = newZ
        return self

    # ------- end class Vector3 ---------


def cross(a, b):
    """
    Vector3 cross product.
    """
    ab_x = a.y*b.z - b.y*a.z
    ab_y = b.x*a.z - a.x*b.z
    ab_z = a.x*b.y - b.x*a.y
    return Vector3(ab_x, ab_y, ab_z)

