from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
from math import pi

# Specify Number of Array Antenna Elements
N = 8
M = 8

# Normalized Wave Number
k = 2*pi

# Array Antenna Element Spacing in wavelengths (consider to be the same in both planes)
dX = 1/2
dY = 1/2

# Array Factor initiated variable for a two-dimensional planar array antenna
AF = [[0 for x in range(180)] for y in range(180)]

AngleTheta = range(180)
AnglePhi = range(180)
AngleThetaRadians = np.zeros(180)
AnglePhiRadians = np.zeros(180)

# Visible region defined as cos(theta) sin(phi)
# VisibleRegion = [[0 for x in range(180)] for y in range(180)]

# Calculate Array Factor for a two-dimensional array antenna
for m in range(M):

    for n in range(N):

        for theta in range(180):
            AngleThetaRadians = theta * pi / 180

            for phi in range(180):
                partTheta = k*np.sin(theta*pi/180)
                partPhi = m * dX * np.cos(phi * pi / 180) + n * dY * np.sin(phi * pi / 180)
                AF[theta][phi] += np.real(np.exp(1j*partTheta*partPhi))
                AnglePhiRadians[phi] = phi*pi/180


# Normalize the Array Factor and Calculate Power Pattern
NormalizedAF = [[0 for x in range(180)] for y in range(180)]
NormalizedPowerPattern = [[0 for x in range(180)] for y in range(180)]
maxAF = np.max(AF)
for theta in range(180):
    for phi in range(180):
        NormalizedAF[theta][phi] = np.abs(AF[theta][phi])/maxAF
        NormalizedPowerPattern[theta][phi] = 20*np.log10(NormalizedAF[theta][phi])

# Calculating the gain of the antenna array (assuming isotropic elements)
AveragePower =(np.average(AF))**2
MaximumPower = maxAF**2
GainDB = 10*np.log10(MaximumPower/AveragePower)

# Plotting the calculated array pattern (factor)
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot3D(AngleThetaRadians, AnglePhiRadians, NormalizedAF, 'gray')