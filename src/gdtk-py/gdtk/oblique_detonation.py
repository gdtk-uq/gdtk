#!/usr/bin/env python
# oblique_detonation.py
#
# This Python script contains a class
# which encapsulates the analytical
# solution for an oblique detonation wave.
#
# The analytical solution was originally published
# by Powers and Stewart (1992) and then re-presented
# as a verfication test case by Powers and Aslam (2006).
# The form of the solution is easier to interpret
# in the 2006 paper.
#
# References:
#
# 1. Powers, J.M. and Stewart, D.S. (1992)
#    Approximate solutions for oblique detonations
#    in the hypersonic limit.
#    AIAA Journal, 30:3 pp. 726--736
#
# 2. Powers, J.M and Aslam, T.D. (2006)
#    Exact solution for multidimensional compressible
#    reactive flow for verifying numerical algorithms
#    AIAA Journal, 44:2 pp. 337--344
#
# This Python script was created by...
# Rowan J Gollan
# 23-Jul-2006
#
# Adjusted for the modern Python3 world of 2022
# and Nick's arrayified geometry bits.
# Peter J. (2022-11-15)
#
from math import cos, sin, sqrt, pow, log, fabs
import numpy as np
from gdtk.numeric.zero_solvers import secant
from gdtk.geom.vector3 import Vector3
from gdtk.geom.path import Spline

class ObliqueDetonation:

    def __init__(self, beta, T1, M1, rho1,
                 R=287.0, alpha=1000.0, gamma=6.0/5.0,
                 q=300000.0):

        self.beta = beta
        self.T1 = T1
        self.M1 = M1
        self.rho1 = rho1
        self.R = R
        self.alpha = 1000.0
        self.gamma = gamma
        self.q = q
        self.p1 = rho1*R*T1
        self.a1 = sqrt(gamma * R * T1)
        self.u1 = self.M1 * self.a1
        self.v1 = 0.0
        self.V = self.u1 * cos(self.beta)

    def get_V(self):
        return self.V

    def calculate_X(self, lmbda ):
        MsinBeta2 = (self.M1 * sin(self.beta))**2
        a1 = (1.0/ ((self.gamma + 1.0) * self.M1 * sin(self.beta))) * (self.a1 / self.alpha)
        a2 = 1.0 + self.gamma * MsinBeta2
        a3 = MsinBeta2 - 1.0
        a4 = ((2.0 * MsinBeta2) / (MsinBeta2 - 1)**2) * ((self.gamma**2 - 1.0) / self.gamma) \
            * ( self.q / (self.R*self.T1))
        #
        OneMinusA4L = 1.0 - a4*lmbda
        OneMinusA4 = 1.0 - a4
        t1 = 2.0*a3*(sqrt(OneMinusA4L) - 1.0)
        t2 = pow( (1.0/(1.0 - lmbda)), a2)
        t3 = 1.0 - sqrt((OneMinusA4L)/(OneMinusA4))
        t4 = 1.0 + sqrt( 1.0 / OneMinusA4 )
        t5 = 1.0 + sqrt((OneMinusA4L)/(OneMinusA4))
        t6 = 1.0 - sqrt( 1.0 / OneMinusA4 )
        #
        X = a1 * ( t1 + log( t2 * pow( (t3*t4) / (t5*t6) , a3*sqrt(OneMinusA4) ) ) )
        return X

    def calculate_rho(self, lmbda ):
        MsinBeta2 = (self.M1 * sin(self.beta))**2
        t1 = self.rho1 * (self.gamma + 1.0 ) * MsinBeta2
        t2 = 1.0 + self.gamma * MsinBeta2
        t3 = t2*t2
        t4 = (self.gamma + 1.0)*MsinBeta2
        t5 = ((self.gamma - 1.0)/self.gamma) * (2.0*lmbda*self.q / (self.R*self.T1))
        t6 = (self.gamma - 1.0)*MsinBeta2
        rho = t1 / ( t2 - sqrt( t3 - t4 * (2.0 + t5 + t6) ) )
        return rho

    def calculate_U(self, lmbda, rho ):
        U = self.rho1 * self.u1 * sin(self.beta) / rho
        return U

    def calculate_T(self, lmbda, rho ):
        t1 = self.p1 / (rho*self.R)
        t2 = (self.rho1*self.u1*sin(self.beta))**2 / (rho*self.R)
        t3 = 1.0/self.rho1 - 1.0/rho
        T = t1 + t2*t3
        return T

    def calculate_p(self, lmbda, rho ):
        t2 = (self.rho1*self.u1*sin(self.beta))**2
        t3 = 1.0/self.rho1 - 1.0/rho
        p = self.p1 + t2 * t3
        return p

    def calculate_Yw(self, lmbda ):
        Yw = ( self.u1*cos(self.beta) / self.alpha ) * log( 1.0 / (1.0 - lmbda) )
        return Yw

    def transform_xy_2_XY(self, x, y):
        X = x * sin(self.beta) - y * cos(self.beta)
        Y = x * cos(self.beta) + y * sin(self.beta)
        return (X, Y)

    def transform_XY_2_xy(self, X, Y):
        x = X * sin(self.beta) + Y * cos(self.beta)
        y = Y * sin(self.beta) - X * cos(self.beta)
        return (x, y)

    def transform_UV_2_uv(self, U, V):
        u = U * sin(self.beta) + V * cos(self.beta)
        v = V * sin(self.beta) - U * cos(self.beta)
        return (u, v)

    def find_XYw_from_x(self, x):
        #
        def f( lmbda ):
            X = self.calculate_X(lmbda)
            Yw = self.calculate_Yw(lmbda)
            (xg, yg) = self.transform_XY_2_xy(X, Yw)
            return (x - xg)
        #
        lmbda = secant(f, 0.0, 0.999, limits=[0.0, 0.999])
        #
        X = self.calculate_X(lmbda)
        Yw = self.calculate_Yw(lmbda)
        #
        return (X, Yw)

    def create_wall_function(self, xmin, xmax):
        def wall(t):
            # Map t --> x, where t may be a scalar or ndarray.
            x = t*(xmax - xmin)
            if (type(x) is float):
                (X, Yw) = self.find_XYw_from_x(x)
                (x, y) = self.transform_XY_2_xy(X, Yw)
            elif type(x) is np.ndarray:
                # Work through the array, element by element.
                shape = x.shape
                xflat = x.flatten()
                yflat = xflat.copy()
                for i, elementx in enumerate(xflat):
                    (X, Yw) = self.find_XYw_from_x(elementx)
                    (elementx, elementy) = self.transform_XY_2_xy(X, Yw)
                    xflat[i] = elementx; yflat[i] = elementy
                x = np.reshape(xflat, shape)
                y = np.reshape(yflat, shape)
            else:
                raise ValueError("Unexpected type for t.")
            return Vector3(x, y)
        return wall

    def solution( self, x, y):
        (X, Y) = self.transform_xy_2_XY(x, y)
        #
        if( X < 0.0 ):
            rho = self.rho1
            p = self.p1
            T = self.T1
            f = [1.0, 0.0]
            u = self.u1
            v = self.v1
        else:
            def f( lmbda ):
                X = self.calculate_X(lmbda)
                (xg, yg) = self.transform_XY_2_xy(X, Y)
                return (x - xg)
            #
            lmbda = secant(f, 0.0, 0.999, limits=[0.0, 0.999])
            #
            rho = self.calculate_rho( lmbda )
            p = self.calculate_p( lmbda, rho )
            T = self.calculate_T( lmbda, rho )
            U = self.calculate_U( lmbda, rho )
            V = self.V
            (u, v) = self.transform_UV_2_uv( U, V)
            f = [ 1.0 - lmbda, lmbda ]
        #
        return (x, y, rho, p, T, f, u, v, X, Y)


if __name__ == '__main__':
    from math import pi
    obl = ObliqueDetonation( pi/4.0, 300.0, 3.0, 1.0)
    #
    X = obl.calculate_X(0.1)
    Y = obl.calculate_Yw(0.1)
    (x, y) = obl.transform_XY_2_xy(X, Y)
    rho = obl.calculate_rho(0.1)
    p = obl.calculate_p(0.1, rho)
    T = obl.calculate_T(0.1, rho)
    U = obl.calculate_U(0.1, rho)
    #
    print("X(lmbda=0.1)= ", X)
    print("Yw(lmbda=0.1)= ", Y)
    print("x(lmbda=0.1)= ", x)
    print("y(lmbda=0.1)= ", y)
    print("rho(lmbda=0.1)= ", rho)
    print("p(lmbda=0.1)= ", p)
    print("T(lmbda=0.1)= ", T)
    print("U(lmbda=0.1)= ", U)
    #
    (X, Yw) = obl.find_XYw_from_x(x)
    print("X from x: ", X)
    print("Yw from x: ", Yw)
    #
    #spline = obl.create_wall_spline(0.0, 1.75, 1.0e-5)
    #
    print("Solution at x=0.066116, y=0.035483...")
    print(obl.solution( 0.066116, 0.035483 ))
    print("Done.")
