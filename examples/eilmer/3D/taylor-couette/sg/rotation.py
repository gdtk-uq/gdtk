import numpy
from libprep3 import *
from math import sin, cos, pi

theta0 = 0.0;
theta3 = 5.0*pi/180.0
nA=[-sin(theta0),cos(theta0),0.0]
t1A=[cos(theta0),sin(theta0), 0.0]
nB=[-sin(theta3),cos(theta3),0.0]
t1B=[cos(theta3),sin(theta3), 0.0]

if type(nA) is list: nA = Vector(nA[0], nA[1], nA[2])
if type(t1A) is list: t1A = Vector(t1A[0], t1A[1], t1A[2])
if type(nB) is list: nB = Vector(nB[0], nB[1], nB[2])
if type(t1B) is list: t1B = Vector(t1B[0], t1B[1], t1B[2])

t2A = cross(nA, t1A) # make the second tangent orthogonal to the original plane
t2A.norm()
t1A = cross(t2A, nA) # make sure that the first tangent is orthogonal
nB.norm()
t2B = cross(nB, t1A)
t2B.norm()
t1B = cross(t2B, nB)
# The rotation matrix transforms the local triplet B into the local triplet A.
mat = numpy.array([
        [nB.x,  nB.y,  nB.z,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0],
        [0.0,   0.0,   0.0,   nB.x,  nB.y,  nB.z,  0.0,   0.0,   0.0],
        [0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   nB.x, nB.y,  nB.z],
        [t1B.x, t1B.y, t1B.z, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0],
        [0.0,   0.0,   0.0,   t1B.x, t1B.y, t1B.z, 0.0,   0.0,   0.0],
        [0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   t1B.x, t1B.y, t1B.z],
        [t2B.x, t2B.y, t2B.z, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0],
        [0.0,   0.0,   0.0,   t2B.x, t2B.y, t2B.z, 0.0,   0.0,   0.0],
        [0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   t2B.x, t2B.y, t2B.z]])
rhs = numpy.array([nA.x, nA.y, nA.z, t1A.x, t1A.y, t1A.z, t2A.x, t2A.y, t2A.z])
RmatrixBtoA = numpy.linalg.solve(mat, rhs)
#print RmatrixBtoA
fp = open("Rmatrix.dat", "w")
fp.write("RmatrixBtoA = %20.20e, %20.20e, %20.20e, %20.20e, %20.20e, %20.20e, %20.20e, %20.20e, %20.20e\n" %
        (RmatrixBtoA[0], RmatrixBtoA[1], RmatrixBtoA[2],
         RmatrixBtoA[3], RmatrixBtoA[4], RmatrixBtoA[5],
         RmatrixBtoA[6], RmatrixBtoA[7], RmatrixBtoA[8]))
R_B_A = numpy.array([
        [RmatrixBtoA[0],  RmatrixBtoA[1],  RmatrixBtoA[2]],
        [RmatrixBtoA[3],  RmatrixBtoA[4],  RmatrixBtoA[5]],
        [RmatrixBtoA[6],  RmatrixBtoA[7],  RmatrixBtoA[8]]])
R_A_B = numpy.linalg.inv(R_B_A)
RmatrixAtoB = numpy.array([R_A_B[0][0], R_A_B[0][1], R_A_B[0][2], 
                           R_A_B[1][0], R_A_B[1][1], R_A_B[1][2], 
                           R_A_B[2][0], R_A_B[2][1], R_A_B[2][2]])
#print RmatrixBtoA
fp.write("RmatrixAtoB = %20.20e, %20.20e, %20.20e, %20.20e, %20.20e, %20.20e, %20.20e, %20.20e, %20.20e\n" %
        (RmatrixAtoB[0], RmatrixAtoB[1], RmatrixAtoB[2],
         RmatrixAtoB[3], RmatrixAtoB[4], RmatrixAtoB[5],
         RmatrixAtoB[6], RmatrixAtoB[7], RmatrixAtoB[8]))
