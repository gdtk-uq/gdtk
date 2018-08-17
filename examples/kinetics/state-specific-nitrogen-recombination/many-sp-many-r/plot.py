#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 17 11:38:23 2018

@author: pmariotto
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

filename = "state-specific-chem.data"
data = pd.read_csv(filename, sep=',')

t = data["t(s)"].values
sp_N2 = ["N2_eX1SIGGPlus_v{}".format(i) for i in range(0,48)]
sp_list = ['N_e2p3Minus4S'] + sp_N2

# =============================================================================
#                       Evolution of massf[sp]
# =============================================================================

plt.figure()
for sp in sp_list:
    plt.loglog(t, data[sp], label=sp)
plt.legend()
plt.xlabel("$t \quad [s]$")
plt.ylabel("$massf$")
plt.show()

# =============================================================================
#                       Temporal evolution of the Vdf
# =============================================================================
n = 10
#idx = range(0,t.shape[0],n)
idx = [0,1,2,3,50]
v = np.arange(0,48)
N2X_v = ["N2_eX1SIGGPlus_v{}".format(i) for i in v]
plt.figure()
x=v
for i in idx:
    y = [data[sp][i] for sp in N2X_v]
    plt.semilogy(x, y, "-.", label="t={:.2e}".format(t[i]))
y = [data[sp].values[-1] for sp in N2X_v]
plt.semilogy(x, y, "-.", label="t={:.2e}".format(t[-1]))

plt.xlabel("Vibrational quantum number $v$")
plt.ylabel("Vfd $massf(N_2X(v))$")
plt.legend()
plt.show()

# =============================================================================
#               Temporal evolution of integrated massf
# =============================================================================

massfN =  data['N_e2p3Minus4S']
massfN2 = [np.sum([data[sp][i] for sp in sp_N2]) for i,t in enumerate(t)]

plt.figure()
plt.plot(t, massfN, "-.", label="N")
plt.plot(t, massfN2, "-.", label="N2")
plt.xlabel("$t \quad [s]$")
plt.ylabel("$massf$")
plt.show()

# =============================================================================
#               Temporal evolution of T
# =============================================================================

plt.figure()
plt.plot(t, data['T(K)'], "-.", label="T")
plt.xlabel("$t \quad [s]$")
plt.ylabel("$T \quad [K]$")
plt.show()