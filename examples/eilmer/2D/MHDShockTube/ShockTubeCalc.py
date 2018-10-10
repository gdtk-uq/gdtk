import numpy as np
from scipy import *
from scipy.optimize import *
#Define sound speed function
def sound_speed(p, rho):
	return (p * gamma / rho) ** 0.5

#Shock tube equation- needs to iterated to find p4 = p3
def shock_tube(p4, p1, p5, rho1, rho5):
	z = (p4 / p5 - 1.)

	c1 = sound_speed(p1, rho1)
	c5 = sound_speed(p5, rho5)

	gm1 = gamma - 1.
	gp1 = gamma + 1.
	g2 = 2. * gamma

	fact1 = gm1 / g2 * (c5 / c1) * z / (1. + gp1 / g2 * z) ** 0.5
	fact2 = (1. - fact1) ** (g2 / gm1)

	return p1 * fact2 - p4


#initialise the left and right states
p1 = 1e4
rho1 = 1
u1 = 0

p5 = 0.1e4
rho5 = 0.125
u5 = 0

#initialise config of the problem
gamma = 5./3.
L = 1.0
t = 2e-3

#Solve the shock tube equation
num_guesses = 100
for pguess in linspace(p5, p1, num_guesses):
	output = fsolve(shock_tube, pguess, args=(p1, p5, rho1, rho5), full_output = True)
	p4, info, ier, mesg = output
	if ier == 1: break
	if not ier == 1:
		print ("Could not find solution for the shock tube")

#Ensure output is a number
if type(p4) is ndarray:
	p4 = p4[0]

#Compute other quantities
z = (p4 / p5 - 1.)
c5 = sound_speed(p5, rho5)

gm1 = gamma - 1
gp1 = gamma + 1
gm1fac = 0.5 * gm1 / gamma
gp1fac = 0.5 * gp1 / gamma

fact = (1 + gp1fac * z) ** 0.5

u4 = c5 * z / (gamma * fact)
rho4 = rho5 * ( 1. + gp1fac * z) / (1. + gm1fac * z)

w = c5 * fact

p3 = p4
u3 = u4
rho3 = rho1 * (p3 / p1) ** (1 / gamma)

c1 = sound_speed(p1, rho1)
c3 = sound_speed(p3, rho3)

#Calculate position of each point of interest i.e. start and finish of rarefaction wave, contact disconinuity and shock
x_rare_start = L / 2. - c1 * t
x_rare_end = L / 2. + (u3 - c3) * t
x_contact_disc = L / 2. + u3 * t
x_shock = L / 2. + w * t

N = 200

x = np.linspace(0 + float(L / N) / 2, L - float(L / N) / 2, N)
rho = np.linspace(0, L, N)
p = np.linspace(0, L, N)
u = np.linspace(0, L, N)
data = np.zeros((N, 4))
print(rho3, rho4)
f = open("ShockTubeAnalytical.dat", "wrb+")
for i in range(N):
	if x[i] < x_rare_start:
		data[i, 0] = x[i]
		data[i, 1] = u1
		data[i, 2] = p1
		data[i, 3] = rho1
	elif x[i] < x_rare_end:
		data[i, 0] = x[i]
		data[i, 1] = 2. / gp1 * (c1 + (x[i] - L / 2.) / t)
		fact = 1. - 0.5 * gm1 * data[i, 1] / c1
		data[i, 2] = p1 * fact ** (2. * gamma / gm1)
		data[i, 3] = rho1 * fact ** (2. / gm1)
	elif x[i] < x_contact_disc:
		data[i, 0] = x[i]
		data[i, 1] = u3
		data[i, 2] = p3
		data[i, 3] = rho3
	elif x[i] < x_shock:
		data[i, 0] = x[i]
		data[i, 1] = u4
		data[i, 2] = p4
		data[i, 3] = rho4
	else:
		data[i, 0] = x[i]
		data[i, 1] = u5
		data[i, 2] = p5
		data[i, 3] = rho5
np.savetxt(f, data)

#Write wave speeds to a file
waves = np.array([-c1, (u3 - c3), u3, w])
f1 = open("ShockTubeWaveSpeed.dat", "wrb+")
np.savetxt(f1, waves)



