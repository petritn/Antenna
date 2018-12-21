import numpy as np
from numpy import pi

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


# Define Planar Array Factor calculation function
def planarAF(M, N, dX, dY, Emn=None):
    if Emn is None:
        arrayFactor = np.zeros((181, 181))
        for m in range(M):
            for n in range(N):
                for theta in range(181):
                    partTheta = 2 * pi * np.sin((theta-90) * pi / 180)
                    for phi in range(181):
                        partPhi = m * dX * np.cos(phi * pi / 180) + n * dY * np.sin(phi * pi / 180)
                        arrayFactor[theta, phi] += np.real(np.exp(1j * (partTheta * partPhi)))
        return arrayFactor
    else:
        arrayFactor = np.zeros((181, 181))
        for m in range(M):
            for n in range(N):
                for theta in range(181):
                    partTheta = 2 * pi * np.sin((theta) * pi / 180)
                    for phi in range(181):
                        partPhi = m * dX * np.cos(phi * pi / 180) + n * dY * np.sin(phi * pi / 180)
                        arrayFactor[theta, phi] += np.real(Emn(m,n)*np.exp(1j * (partTheta * partPhi)))
        return arrayFactor

# Define the Planar Array Factor calculation function as a product of two Linear Arrays
def planarAF2(M, N, dX, dY, Emn=None) -> np.ndarray:
    arrayFactor = np.zeros((181, 181))
    lineararrayXm = np.zeros((181, 181))
    lineararrayYn = np.zeros((181, 181))

    if Emn is None:
        for m in range(M):
            for theta in range(181):
                for phi in range(181):
                    lineararrayXm[theta, phi] += np.real(
                        np.exp(1j * m * 2 * pi * dX * np.sin((theta - 90) * pi / 180) * np.cos(phi * pi / 180) + pi / 2))

        for n in range(N):
            for theta in range(181):
                for phi in range(181):
                    lineararrayYn[theta, phi] += np.real(
                        np.exp(1j * n * 2 * pi * dY * np.sin((theta - 90) * pi / 180) * np.sin(phi * pi / 180)))

    else:
        for m in range(M):
            for theta in range(181):
                for phi in range(181):
                    lineararrayXm[theta, phi] += np.real(Emn(m, 0) *
                        np.exp(1j * m * 2 * pi * dX * np.sin((theta - 90) * pi / 180)* np.cos(phi * pi / 180) + pi / 2))

        for n in range(N):
            for theta in range(181):
                for phi in range(181):
                    lineararrayYn[theta, phi] += np.real(Emn(0, n) *
                        np.exp(1j * n * 2 * pi * dY * np.sin((theta - 90) * pi / 180) * np.sin(phi * pi / 180)))

    for theta in range(181):
        for phi in range(181):
            arrayFactor[theta, phi] = lineararrayXm[theta, phi] * lineararrayYn[theta, phi]

    return arrayFactor