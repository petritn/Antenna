import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Cartesian to Spherical Coordinate Transformation
def cart2sph(x, y, z):
    azimuth = np.arctan2(y,x)
    elevation = np.arctan2(z,np.sqrt(x**2 + y**2))
    r = np.sqrt(x**2 + y**2 + z**2)
    return azimuth, elevation, r

# Spherical to Cartesian Coordinate Transformation
def sph2cart(azimuth, elevation, r):
    x = r * np.cos(elevation) * np.cos(azimuth)
    y = r * np.cos(elevation) * np.sin(azimuth)
    z = r * np.sin(elevation)
    return x, y, z

# Plotting 3d antenna pattern
def plot_pattern(antenna_pattern):
    fig = plt.figure()
    ax = Axes3D(fig)

    x_ = np.ones((181, 181))
    y_ = np.ones((181, 181))
    z_ = np.ones((181, 181))

    for phi in range(181):
        for theta in range(181):
            p = antenna_pattern[phi][theta]

            xp, yp, zp = sph2cart(np.radians(theta), np.radians(phi), p)

            x_[phi, theta] = xp
            y_[phi, theta] = yp
            z_[phi, theta] = zp

    xx = (np.min(x_), np.max(x_))
    yx = (np.min(y_), np.max(y_))
    ax.plot_wireframe(x_, y_, z_, rstride=5, cstride=5)
    plt.xlim(xx[0], xx[1])
    plt.ylim(yx[0], yx[1])
    plt.xlabel("$\Theta$")
    plt.ylabel("$\Phi$")
    plt.grid()
    plt.show()

# End
