from CArray import LinearArray, PlanarArray
from LibraryOfFunctions import plot_pattern
import matplotlib.pyplot as plt
import numpy as np
from numpy import pi

# Create e linear array...
n = 8
d = 0.5
e = np.zeros(n, dtype=complex)

for i in range(n):
    e[i] = np.exp(-1j*i*pi)

angle = np.arange(0, 181, 1)
angle_radians = np.radians(angle)

la_ = LinearArray(n, d)

la_we = LinearArray(n, d, e)

plt.polar(angle_radians, la_.normalized_array_factor)
plt.polar(angle_radians, la_we.normalized_array_factor)
plt.show()

# Create a planar array...
pa = PlanarArray(4, 4, 0.5, 0.5)

plaf = pa.normalized_array_factor

plot_pattern(plaf)

# End

