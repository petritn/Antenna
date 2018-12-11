import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from math import pi

# Define visibility region
DIM = 181

# Specify Number of Array Antenna Elements
N = 8
M = 8

# Normalized Wave Number
k = 2*pi

# Array Antenna Element Spacing in wavelengths (consider to be the same in both planes)
dX = 1/2
dY = 1/2

# Array Factor initiated variable for a two-dimensional planar array antenna
AF = [[0 for x in range(DIM)] for y in range(DIM)]
Sxm = [[0 for x in range(DIM)] for y in range(DIM)]
Sym = [[0 for x in range(DIM)] for y in range(DIM)]

AngleTheta = np.arange(0, DIM, 1)
AnglePhi = np.arange(0, DIM, 1)
AngleThetaRadians = AnglePhiRadians = np.zeros(DIM)
VisibleRegionX = VisibleRegionY = np.zeros(DIM)

for i in range(DIM):
    AngleThetaRadians[i] = AnglePhiRadians[i] = i*pi/180
    VisibleRegionX[i] = np.cos(i*pi/180)
    VisibleRegionY[i] = np.cos(i*pi/180 - pi)
'''
# Calculate Array Factor for a two-dimensional array antenna
for m in range(M):
    for n in range(N):
        for theta in range(DIM):
            partTheta = k * np.sin(theta * pi / 180 + pi / 2)
            for phi in range(DIM):
                partPhi = m * dX * np.cos(phi * pi / 180) + n * dY * np.sin(phi * pi / 180)
                AF[theta][phi] += np.real(np.exp(1j*partTheta*partPhi))
'''

# Calculating the Array Factor of a Planar Array as a product of two linear arrays
for m in range(M):
    for theta in range(DIM):
        for phi in range(DIM):
            Sxm[theta][phi] += np.real(np.exp(1j*m*k*dX*np.sin(theta*pi/180)*np.cos(phi*pi/180)+pi/2))

for n in range(N):
    for theta in range(DIM):
        for phi in range(DIM):
            Sym[theta][phi] += np.real(np.exp(1j*n*k*dY*np.sin(theta*pi/180)*np.sin(phi*pi/180)))

for theta in range(DIM):
    for phi in range(DIM):
        AF[theta][phi] = Sxm[theta][phi] * Sym[theta][phi]

# Normalize the Array Factor and Calculate Power Pattern
NormalizedAF = [[0 for x in range(DIM)] for y in range(DIM)]
NormalizedPowerPattern = [[0 for x in range(DIM)] for y in range(DIM)]
maxAF = np.max(AF)

# Take a slice
SliceTheta = np.zeros(DIM)
SlicePhi = np.zeros(DIM)

for theta in range(DIM):
    SliceTheta[theta] = AF[theta][1]
    SlicePhi[theta] = AF[1][theta]
    for phi in range(DIM):
        NormalizedAF[theta][phi] = np.abs(AF[theta][phi])/maxAF
        NormalizedPowerPattern[theta][phi] = 20*np.log10(NormalizedAF[theta][phi])

# Calculating the gain of the antenna array (assuming isotropic elements)
AveragePower = np.average(AF)**2
MaximumPower = maxAF**2
GainDB = 10*np.log10(MaximumPower/AveragePower)

# Plotting the calculated array pattern (factor)
#plt.plot(SliceTheta)
plt.polar3d(AngleThetaRadians, AnglePhiRadians, Sxm)
plt.show()

fig = plt.figure()
ax = Axes3D(fig)
X, Y = np.meshgrid(VisibleRegionX, VisibleRegionY)
ax.plot_wireframe(X, Y, Sym, rstride=1, cstride=1)
plt.xlabel("sin($\Theta$)")
plt.ylabel("cos($\Phi$)")
plt.show()