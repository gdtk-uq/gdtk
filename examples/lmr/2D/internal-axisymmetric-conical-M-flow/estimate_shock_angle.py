# Author: RJG
# Date: 2025-06-25
#
# This script can be used to estimate the shock angle
# from the Eilmer flow solution for m-flow, after building
# VTK files with Mach number added:
#
# > make post
#
# Then this script can be called:
#
# > python estimate_shock_angle.py

from gdtk.lmr import LmrConfig, SimInfo
import scipy.stats.mstats as mstats
import math

M_inf = 5.0
M_threshold = 0.995*M_inf

sim = SimInfo(LmrConfig())
ff_pv = sim.load_pvd_into_pyvista(as_point_data=True, merged=True)
shock = ff_pv.contour(isosurfaces=[M_threshold], scalars="M_local")
result = mstats.linregress(shock.points[:,0], shock.points[:,1])
print(f"shock_angle= {math.degrees(math.atan(result.slope)):.3f} degrees")
