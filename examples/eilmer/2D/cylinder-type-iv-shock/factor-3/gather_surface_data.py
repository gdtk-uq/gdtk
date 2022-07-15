# gather_surface_data.py
# PJ 2022-07-02

from math import degrees, atan2

ths = []; ps = []; qs = []
tindx = 9
for blk in range(6):
    fileName = "loads/t%04d/b%04d.t%04d.loads.dat" % (tindx,blk,tindx)
    with open(fileName, 'r') as f:
        lines = f.readlines()
        for txt in lines:
            if txt.startswith('#'): continue
            items = txt.strip().split(' ')
            x = 0.008 - float(items[0])
            y = float(items[1])
            ths.append(degrees(atan2(y, x)))
            ps.append(float(items[9]))
            qs.append(-float(items[23]))

print("p_max = ", max(ps), "Pa")
print("q_max = ", max(qs)/10000., "W/cm**2")
with open("surface.data", 'w') as f:
    f.write("# angle(degrees) pressure(Pascals) heat_flux(W/cm**2)\n")
    for i in range(len(ths)):
        f.write("%f %f %f\n" % (ths[i], ps[i], qs[i]/10000.))

print("Done.")
