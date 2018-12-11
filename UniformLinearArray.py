import numpy as np
import matplotlib.pyplot as plt
from math import pi

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
AF = np.zeros(180)

Angle = range(180)
AngleRadians = np.zeros(180)
VisibleRegion = np.zeros(180) # Visible region defined as cos(theta)

# Calculate Array Factor
for n in range(N):
    for theta in range(180):
        Gamma = k*ElementSpacing*np.cos(theta*pi/180)
        AF[theta] += np.real(np.exp(1j*n*Gamma))
        AngleRadians[theta] = theta*pi/180
        VisibleRegion[theta] = np.cos(theta*pi/180)

# Normalize the Array Factor
NormalizedAF = np.zeros(180)
NormalizedPowerPattern = np.zeros(180)
maxAF = np.max(AF)
for theta in range(180):
    NormalizedAF[theta] = np.abs(AF[theta])/maxAF
    NormalizedPowerPattern[theta] = 20*np.log10(NormalizedAF[theta])

# Calculating the gain of the antenna array (assuming isotropic elements)
AveragePower =(np.average(AF))**2
MaximumPower = maxAF**2
GainDB = 10*np.log10(MaximumPower/AveragePower)


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