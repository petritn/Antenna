import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from math import pi

# Define Visibility Region
DIM = 181

# Number of Array Antenna Elements
N = 8

# Amplitude of Element Excitation
Amplitude = np.ones(N)

# Phase of Element Excitation
Phase = np.zeros(N)

# Normalized Wave Number
k = 2*pi

# Array Antenna Element Spacing in wavelengths
ElementSpacing = 1/2

# Array Factor initiated variable
AF = np.zeros(DIM)

Angle = np.arange(0, DIM, 1)
AngleRadians = np.zeros(DIM)

# Visible region defined as cos(theta)
VisibleRegion = np.zeros(DIM)

# Calculate Array Factor of Linear Array
for n in range(N):
    for theta in range(DIM):
        Gamma = k*ElementSpacing*np.cos(theta*pi/180)
        AF[theta] += np.real(np.exp(1j*n*Gamma))
        AngleRadians[theta] = theta*pi/180
        VisibleRegion[theta] = np.cos(theta*pi/180)

# Normalize the Array Factor
NormalizedAF = np.zeros(DIM)
NormalizedPowerPattern = np.zeros(DIM)
maxAF = np.max(AF)
for theta in range(DIM):
    NormalizedAF[theta] = np.abs(AF[theta])/maxAF
    NormalizedPowerPattern[theta] = 20*np.log10(NormalizedAF[theta])

# Calculating the gain of the antenna array (assuming isotropic elements)
AveragePower =(np.average(AF))**2
MaximumPower = maxAF**2
GainDB = 10*np.log10(MaximumPower/AveragePower)
print(GainDB)


# Plotting the calculated array pattern (factor)
plt.figure(1)
plt.plot(VisibleRegion, NormalizedPowerPattern)
plt.grid()
plt.xlabel("cos($\Theta$)")
plt.ylabel("Normalized Power Pattern (dB)")
plt.ylim(-40, 0.0)
plt.show()

plt.figure(2)
plt.polar(AngleRadians, NormalizedAF)
plt.title("Normalized Array Factor")
plt.xlabel("$\Theta$, degrees")
plt.show()

fig = plt.figure()
ax = Axes3D(fig)

theta = np.arange(0, DIM, 1)
phi = np.arange(0, DIM, 1)

x = np.outer(np.cos(phi), np.sin(theta))
y = np.outer(np.sin(phi), np.sin(theta))

ax.plot_wireframe(x, y, NormalizedAF)
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
plt.show()

