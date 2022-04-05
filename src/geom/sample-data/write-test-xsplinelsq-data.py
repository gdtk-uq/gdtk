# write-test-xsplinelsq-data.py
import math
import numpy as np
N = 101
xd = np.linspace(0.0, 2.0*math.pi, N)
yd = np.sin(xd)
with open('test-xsplinelsq.dat', 'w') as f:
    f.write('# test data for CubicSplineLsq, XSplineLsq2\n')
    for i in range(N):
        f.write("%g %g 1.0\n" % (xd[i], yd[i]))
