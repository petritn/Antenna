import numpy as np
import matplotlib.pyplot as plt
import math
from mpl_toolkits.mplot3d import Axes3D
from General import sph2cart, planarAF, planarAF2
from math import pi

# Define visibility region
DIM = 181

# Specify Number of Array Antenna Elements
N = 8
M = 8

# Complex Element Excitations
Emn = np.ones((M, N), dtype=complex)

Emn[0, 0] = 2*np.exp(1j*pi)
print(Emn[0, 0])

# Normalized Wave Number
k = 2*pi

# Array Antenna Element Spacing in wavelengths (consider to be the same in both planes)
dX = 1/2
dY = 1/2

# Array Factor initiated variable for a two-dimensional planar array antenna
AF = np.zeros((DIM, DIM))
Sxm = np.zeros((DIM, DIM))
Sym = np.zeros((DIM, DIM))

AngleTheta = np.arange(0, DIM, 1)
AnglePhi = np.arange(0, DIM, 1)
AngleThetaRadians = AnglePhiRadians = np.zeros(DIM)
VisibleRegionX = VisibleRegionY = np.zeros(DIM)

for i in range(DIM):
    AngleThetaRadians[i] = AnglePhiRadians[i] = i*pi/180
    VisibleRegionX[i] = np.cos(i*pi/180)
    VisibleRegionY[i] = np.sin(i*pi/180)


# Calculate Array Factor for a two-dimensional array antenna
AF = planarAF2(M, N, dX, dY)

# Normalize the Array Factor and Calculate Power Pattern
NormalizedAF = [[0 for x in range(DIM)] for y in range(DIM)]
NormalizedPowerPattern = [[0 for x in range(DIM)] for y in range(DIM)]
maxAF = np.max(AF)

for theta in range(DIM):
    for phi in range(DIM):
        NormalizedAF[theta][phi] = np.abs(AF[theta, phi])/maxAF
        NormalizedPowerPattern[theta][phi] = 20*np.log10(NormalizedAF[theta][phi])

# Calculating the gain of the antenna array (assuming isotropic elements)
AveragePower = np.average(AF)**2
MaximumPower = maxAF**2
GainDB = 10*np.log10(MaximumPower/AveragePower)

# Plotting the calculated array pattern (factor)
fig = plt.figure()
ax = Axes3D(fig)

X = np.ones((DIM, DIM))
Y = np.ones((DIM, DIM))
Z = np.ones((DIM, DIM))

for phi in range(DIM):
    for theta in range(DIM):
        e = AF[phi][theta]

        xe, ye, ze = sph2cart(math.radians(theta), math.radians(phi), e)

        X[phi, theta] = xe
        Y[phi, theta] = ye
        Z[phi, theta] = ze

#ax.plot_wireframe(X, Y, Z, rstride=2, cstride=2)
ax.plot_surface(X, Y, Z, cmap=plt.cm.YlGnBu_r)
plt.xlabel("$\Theta$")
plt.ylabel("$\Phi$")
plt.grid()
plt.show()

