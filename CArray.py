import numpy as np
from numpy import pi


class LinearArray(object):
    """Linear Antenna Array

    Parameters
    ------------
    N : int
        Number of antenna elements in the antenna array
    dX: float
        Antenna element spacing in the array.
    E: complex
        Complex antenna element excitations

    Attributes
    -----------
    Array Factor : nd-array
        Array factor of the linear antenna array.
    Gain DB : float
        Gain of the antenna array
    """

    def __init__(self, n, dx, e=None):
        if E is None:
            self.n = n
            self.dx = dx
        else:
            self.n = n
            self.dx = dx
            self.e = e

    def array_factor(self):
        """ Calculate the array factor of the antenna linear array

        Returns
        -------
        The calculated array factor of the antenna array, for the linear antenna array object

        """
        af = np.zeros(181)
        for n in range(self.N):
            for theta in range(181):
                gamma = 2 * pi * self.dX * np.cos(theta * pi / 180)
                af[theta] += np.real(np.exp(1j * n * gamma))
        return af

    def normalized_array_factor(self):
        """ Calculate the normalized array factor of the linear antenna array

        Returns
        -------
        Normalized antenna array factor

        """
        af = np.abs(self.array_factor())
        af_max = np.max(af)
        naf = af/af_max
        return naf

    def normalized_power_pattern(self):
        """ Calculate the normalized power pattern

        Returns
        -------
        Normalized power pattern of the linear array

        """
        pattern = self.normalized_array_factor()
        pp = 20 * np.log10(pattern)
        return pp

    def gain_db(self):
        """ Calculate the gain of the antenna array in dB

        Returns
        -------
        Array Antenna Gain in dB.

        """
        af = self.array_factor()
        average_power = (np.average(af)) ** 2
        max_power = np.max(af)
        gain = 10 * np.log10(max_power / average_power)
        return gain





