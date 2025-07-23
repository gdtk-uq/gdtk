"""
Example script for plotting residuals
"""

from gdtk import slf
import matplotlib.pyplot as plt
from sys import argv

plt.rcParams['font.family'] = 'serif'
plt.rcParams['svg.fonttype'] = 'none'

filename = argv[1]
residuals = slf.read_log_file(filename)

fig = plt.figure(figsize=(9,4.2))
axes,axes1 = fig.subplots(1,2)

axes.set_xlabel('iters')
axes.set_title('Global Relative Residual')
axes.set_ylabel('GRR')
axes.semilogy(residuals['iter'], residuals['GRR'], 'r-')
axes.grid()

axes1.set_xlabel('iters')
axes1.set_title('Timestep')
axes1.set_ylabel('dt (sec)')
axes1.semilogy(residuals['iter'], residuals['dt'], 'b-')
axes1.grid()

plt.tight_layout()
plt.show()

