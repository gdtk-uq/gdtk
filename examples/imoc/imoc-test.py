# imoc-test.py
# Try out the imoc kernel and unit processes.
#
# Usage:
# python3 imoc-test.py
#
# PJ, 2019-12-28
#
import eilmer.imoc.kernel as imock

print("Begin imoc test...")

a = imock.Node()
b = imock.Node()
print("a=", a, "b=", b)

print("Done.")
