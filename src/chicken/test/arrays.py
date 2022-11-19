# arrays.py
# Experiment with writing the numpy arrays.
# PJ 2022-11-04

from zipfile import ZipFile
import numpy as np

a = np.zeros((2,3),dtype=float)
b = np.array([[1.0, 2.0],[3.0, 4.0]])
print("a=", a)
print("b=", b)

with ZipFile('test.zip', mode='w') as zf:
    with zf.open('a', mode='w') as fp:
        np.savetxt(fp, a.flatten())
    with zf.open('b', mode='w') as fp:
        np.savetxt(fp, b.flatten())
