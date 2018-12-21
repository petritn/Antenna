import numpy as np
import matplotlib.pyplot as plt
from math import pi
from CArray import LinearArray

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
#AF = np.zeros(DIM)

Angle = np.arange(0, DIM, 1)
AngleRadians = np.zeros(DIM)

# Visible region defined as cos(theta)
VisibleRegion = np.zeros(DIM)

LA = LinearArray(8, 0.5)
# Calculate Array Factor of Linear Array
AF = LA.normalized_array_factor()
NormalizedPowerPattern = LA.normalized_power_pattern()

for theta in range(DIM):
    AngleRadians[theta] = theta*pi/180
    VisibleRegion[theta] = np.cos(theta*pi/180)


# Calculating the gain of the antenna array (assuming isotropic elements)
gain = LA.gain_db()
print(gain)

# Plotting the calculated array pattern (factor)
plt.figure(1)
plt.plot(VisibleRegion, NormalizedPowerPattern)
plt.grid()
plt.xlabel("cos($\Theta$)")
plt.ylabel("Normalized Power Pattern (dB)")
plt.ylim(-40, 0.0)
plt.show()

plt.figure(1)
plt.polar(AngleRadians, AF)
plt.title("Normalized Array Factor")
plt.xlabel("$\Theta$, degrees")
plt.show()

