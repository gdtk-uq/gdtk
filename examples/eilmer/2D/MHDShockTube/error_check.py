#Script used to check the accuracy of the solver
from scipy import *
from numpy import *

#Load the results from the eilmer version and analytical version
exact = loadtxt("ShockTubeAnalytical.dat")
solver = loadtxt("ShockTubeNumerical.dat")
waves = loadtxt("ShockTubeWaveSpeed.dat")

(a, b) = shape(exact)

#Start with a basic measure of error- difference between exact and numerical across the tube (RMS error)
error_u = 0.
error_p = 0.
error_rho = 0.

output = ndarray(shape=(1,4))
output[0,0] = a
for i in range(a):
	output[0, 1] += (exact[i, 1] - solver[i, 1]) ** 2
	output[0, 2] += (exact[i, 2] - solver[i, 2]) ** 2
	output[0, 3] += (exact[i, 3] - solver[i, 3]) ** 2

output[0, 1] = (output[0, 1] / a) ** 0.5
output[0, 2] = (output[0, 2] / a) ** 0.5
output[0, 3] = (output[0, 3] / a) ** 0.5

f = open("error_measure.dat", 'a+')

print(shape(output))
savetxt(f, output)

#Start by finding the flow characteristics in each zone of the tube by taking approximate measures of where the wave phenomena are
tol = 1
waveloc = ndarray(shape=(2, 4))
waveloc[0, 0] = 0
waveloc[1, 3] = a
n = 0
trigger = False

for i in range(2, a-2):
	#By taking average rate of change over 5 cells we remove some sensitivity to oscillations- we dont need the location to be that accurate
	deriv = 0.5 * (solver[i+2, 1] - solver[i-2, 1]) / (solver[i, 0] - solver[i-1, 0])
	if fabs(deriv) > tol and trigger == False:
		#wave start has been detected
		print("wave detected at")
		print(solver[i, 0])
		waveloc[1, n] = i
		trigger = True

	if fabs(deriv) < tol and trigger == True:
		#wave has been resolved
		print("wave resolved at")
		print(solver[i, 0])
		waveloc[0, n+1] = i
		trigger = False
		n += 1

#Write the velocity, pressure and density in each zone to an array- use the value half between where each zone was defined
flowstates = ndarray(shape=(4, 3))

for i in range(4):
	index = int(round((waveloc[0, i] + waveloc[1, i]) / 2))
	for j in range(3):
		flowstates[i, j] = solver[index, j+1]
	if i == 1:
		[u3_exact, p3_exact, rho3_exact] = exact[index, 1:4]
	if i == 2:
		[u4_exact, p4_exact, rho4_exact] = exact[index, 1:4]

#Put results in easy to read state- no state 2 as it is part of the rarefaction wave
[rho1, p1, u1] = flowstates[0, :]
[rho3, p3, u3] = flowstates[1, :]
[rho4, p4, u4] = flowstates[2, :]
[rho5, p5, u5] = flowstates[3, :]

#Use these results to strictly define the locations of the contact discontinuity and the shock- we say the wave occurs when the density is halfway between each state (i.e. rho3.5, rho4.5)
rho_CD = 0.5 * (rho3 + rho4)
rho_SS = 0.5 * (rho4 + rho5)

for i in range(a-1):
	#Find out where the wave occurs and linearly extrapolate
	if solver[i, 1] > rho_CD and solver[i+1, 1] < rho_CD:
		percent = (solver[i, 1] - rho_CD) / (solver[i, 1] - solver[i+1, 1])
		CD_loc = solver[i, 0] - percent * (solver[i, 0] - solver[i+1, 0])
	if solver[i, 1] > rho_SS and solver[i+1, 1] < rho_SS:
		percent = (solver[i, 1] - rho_SS) / (solver[i, 1] - solver[i+1, 1])
		SS_loc = solver[i, 0] - percent * (solver[i, 0] - solver[i+1, 0])

#Convert location to a speed
CD_speed = (CD_loc - 0.5) / 2e-3
SS_speed = (SS_loc - 0.5) / 2e-3

#Calculate errors
CD_speed_error = fabs(CD_speed - waves[2]) / waves[2]
SS_speed_error = fabs(SS_speed - waves[3]) / waves[3]

CD_rho_error = fabs(rho3 - rho3_exact) / rho3_exact
SS_rho_error = fabs(rho4 - rho4_exact) / rho4_exact

#Write to a text file
with open("error_tests.dat", "a+") as f:
	f.write(str(a) + " " + str(CD_speed_error) + " " + str(SS_speed_error) + " " + str(CD_rho_error) + " " + str(SS_rho_error) + " " + str(CD_loc) + " " + str(SS_loc) + "\n") 

"""
f = open("error_tests.dat", "a+")
error_data = ndarray(shape=(1,4))
error_data = [CD_speed_error, SS_speed_error, CD_rho_error, SS_rho_error]

savetxt(f, error_data)
"""
		
