# Auto run for wall properties
print(' Run the ETC ON condition post.py')

print('== Extracting surface quantities')
# Compute s/R, Twall, q-wall and y+, the latter is read from the loads file and write a data
# file for plotting via gnuplot.

# 1:pos.x 2:pos.y 3:pos.z 4:area 5:q 6:tau 7:l_tau 8:m_tau 9:n_tau 
# 10:sigma 11:n.x 12:n.y 13:n.z 14:T 15:Re_wall 16:y+ 17:cellWidthNormalToSurface

xl = []
yl = []
qwl = []
Twl = []
y_plusl = []
s_on_rnl = []

import glob, os
load_files = [d for d in os.listdir('../loads/') if os.path.isdir(os.path.join('../loads/', d))]
load_files.sort()
print load_files
loads_tindx_last = load_files[-1]

for file in glob.glob("../loads/%s/*" % (loads_tindx_last)):
    print file
    with open(file) as fp:
        for line in fp.readlines():
            if line[0] == '#':
                continue
            else:
                line = line.split()
                # Append current values to lists
                xl.append(float(line[0]))
                yl.append(float(line[1]))
                qwl.append(float(line[4]))
                Twl.append(float(line[15]))
                y_plusl.append(float(line[17]))

# Sort the lists
xLL = [x for y,x in sorted(zip(xl,xl))]
yLL = [x for y,x in sorted(zip(xl,yl))]
yp = [x for y,x in sorted(zip(xl,y_plusl))]
qw = [x for y,x in sorted(zip(xl,qwl))]
Tw = [x for y,x in sorted(zip(xl,Twl))]

# Nose radius for normalisation
Rn = 0.010  # m

for i in range(len(xLL)):
    tx = xLL[i]
    ty = yLL[i]

    # Compute difference between current and last
    if i == 0:
        xdiff = tx
        ydiff = ty
    else:
        xdiff = tx - xLL[i-1]
        ydiff = ty - yLL[i-1]

    sdiff = (xdiff**2 + ydiff**2)**0.5
    # Add new s_on_rn point to list
    if i == 0:
        s_on_rnl.append(sdiff/Rn)
    else:
        s_on_rnl.append(s_on_rnl[-1] + sdiff/Rn)

# Write nicely formatted list for gnuplot
with open('surface-properties.dat', 'w') as fp:
    fp.write('# s/Rn qw T y+\n')
    for i in range(len(xl)):
        fp.write('{:.6e} {:.6e} {:.6e} {:.6e}\n'.format(
            s_on_rnl[i], qw[i], Tw[i], yp[i]))

# -----------------------------------------------
# -----------------------------------------------

## Run gnu scripts
os.system('gnuplot wall_temperature.gnuplot')
os.system('gnuplot heat_transfer.gnuplot')
os.system('gnuplot y_plus.gnuplot')