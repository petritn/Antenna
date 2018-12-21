import numpy as np
import matplotlib.pyplot as plt
import math
from mpl_toolkits.mplot3d import Axes3D
from General import sph2cart
from math import pi
from CArray import PlanarArray

# Define visibility region
DIM = 181

# Specify Number of Array Antenna Elements
N = 6
M = 6

# Complex Element Excitations
Emn = np.ones((M, N), dtype=complex)

# Emn[0, 0] = 2*np.exp(1j*pi)
# print(Emn[0, 0])

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
pa = PlanarArray(M, N, dX, dY)

AF = pa.array_factor2

# Normalize the Array Factor and Calculate Power Pattern
NormalizedAF = pa.normalized_array_factor2
NormalizedPowerPattern = pa.normalized_power_pattern2

print(pa.gain_db)

# Plotting the calculated array pattern (factor)
fig = plt.figure()
ax = Axes3D(fig)

X = np.ones((DIM, DIM))
Y = np.ones((DIM, DIM))
Z = np.ones((DIM, DIM))

for phi in range(DIM):
    for theta in range(DIM):
        e = NormalizedAF[phi][theta]

        xe, ye, ze = sph2cart(math.radians(theta), math.radians(phi), e)

        X[phi, theta] = xe
        Y[phi, theta] = ye
        Z[phi, theta] = ze

ax.plot_wireframe(X, Y, Z, rstride=2, cstride=2)
#ax.plot_surface(X, Y, Z, color='red')
plt.xlabel("$\Theta$")
plt.ylabel("$\Phi$")
plt.grid()
plt.show()

