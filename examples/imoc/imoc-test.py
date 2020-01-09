# imoc-test.py
# Try out the imoc kernel and unit processes.
#
# Usage:
# python3 imoc-test.py
#
# PJ, 2019-12-28, first step
#     2020-01-09, unit process for interior node
#
import eilmer.imoc.kernel as kernel
import eilmer.imoc.unit_process as unit
from eilmer.ideal_gas_flow import PM1

print("Begin imoc test...")
kernel.axisymmetric = True

print("Define a couple on nodes.")
a = kernel.Node(x=1, y=1, nu=PM1(2.1), mach=2.1, theta=0.1)
b = kernel.Node(x=1, y=0, nu=PM1(2.0), mach=2.0, theta=0.0)
print("a=", a)
print("b=", b)

print("Make a new, interior node.")
c = unit.interior(a.indx, b.indx, -1)
print("a=", a)
print("b=", b)
print("c=", c)

print("Done.")
