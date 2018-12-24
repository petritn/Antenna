import numpy as np
from numpy import pi

###############################
# Linear Array Class Definition


class LinearArray(object):
    """Linear Antenna Array

    Parameters
    ------------
    n : int
        Number of antenna elements in the antenna array
    dx: float
        Antenna element spacing in the array.
    e: list
        Complex antenna element excitations, dimension on n

    Attributes
    -----------
    Array Factor : nd-array
        Array factor of the linear antenna array.
    Gain DB : float
        Gain of the antenna array
    """

    def __init__(self, n, dx, e=None):
        if e is None:
            self.n = n
            self.dx = dx
        else:
            self.n = n
            self.dx = dx
            self.e = e

    @property
    def array_factor(self):
        """ Calculate the array factor of the antenna linear array

        :return:
        The calculated array factor of the antenna array, for the linear antenna array object

        """
        if hasattr(self, 'e'):
            af = np.zeros(181)
            for n in range(self.n):
                for theta in range(181):
                    gamma = 2 * pi * self.dx * np.cos(theta * pi / 180)
                    af[theta] += np.real(self.e[n]*np.exp(1j * n * gamma))

        else:
            af = np.zeros(181)
            for n in range(self.n):
                for theta in range(181):
                    gamma = 2 * pi * self.dx * np.cos(theta * pi / 180)
                    af[theta] += np.real(np.exp(1j * n * gamma))

        return af

    @property
    def power_pattern(self):
        return self.array_factor**2

    @property
    def normalized_array_factor(self):
        """ Calculate the normalized array factor of the linear antenna array

        :return:
        Normalized antenna array factor

        """
        af = np.abs(self.array_factor)
        af_max = np.max(af)
        naf = af / af_max
        return naf

    @property
    def normalized_power_pattern(self):
        """ Calculate the normalized power pattern

        :return:
        Normalized power pattern of the linear array in dB

        """
        pattern = self.normalized_array_factor
        pp = 20 * np.log10(pattern)
        return pp

    @property
    def gain_db(self):
        """ Calculate the gain of the antenna array in dB

        :return:
        Array Antenna Gain in dB.

        """
        af = self.array_factor
        average_power = (np.average(af)) ** 2
        max_power = np.max(af)
        gain = 10 * np.log10(max_power / average_power)
        return gain

###############################
# Planar Array Class Definition


class PlanarArray(object):
    """Planar Antenna Array

        Parameters
        ------------
        m: int
            Number of antenna elements in the X direction of antenna array
        n: int
            Number of antenna elements in the Y direction of the antenna array
        dx: float
            Antenna element spacing in the X direction of the array.
        dy: float
            Antenna element spacing in the Y direction of the array.
        e: np.array
            Complex antenna element excitations, dimension of (m, n)

        Attributes
        -----------
        Array Factor : nd-array
            Array factor of the planar antenna array.
        Gain DB : float
            Gain of the antenna array
        """

    def __init__(self, m, n, dx, dy, e=None):
        if e is None:
            self.m = m
            self.n = n
            self.dx = dx
            self.dy = dy
        else:
            self.m = m
            self.n = n
            self.dx = dx
            self.dy = dy
            self.e = e

    @property
    def array_factor(self):
        """ Calculate the array factor of the planar antenna array

        :return:

        The calculated array factor of the antenna array, for the planar antenna array object

        """
        if hasattr(self, 'e'):
            af = np.zeros((181, 181))
            for m in range(self.m):
                for n in range(self.n):
                    for theta in range(181):
                        part_theta = 2 * pi * np.sin((theta - 90) * pi / 180)
                        for phi in range(181):
                            part_phi = m * self.dx * np.cos(phi * pi / 180) + n * self.dy * np.sin(phi * pi / 180)
                            af[theta, phi] += np.real(self.e[m][n] * np.exp(1j * (part_theta * part_phi)))
            return af
        else:
            af = np.zeros((181, 181))
            for m in range(self.m):
                for n in range(self.n):
                    for theta in range(181):
                        part_theta = 2 * pi * np.sin((theta - 90) * pi / 180)
                        for phi in range(181):
                            part_phi = m * self.dx * np.cos(phi * pi / 180) + n * self.dy * np.sin(phi * pi / 180)
                            af[theta, phi] += np.real(np.exp(1j * (part_theta * part_phi)))
            return af

    @property
    def array_factor2(self):
        """ Calculate the planar array antenna factor as a product of two linear arrays

        :return:

        Calculated antenna array factor as a product of two linear arrays

        """

        af = np.zeros((181, 181))
        linear_array_x = np.zeros((181, 181))
        linear_array_y = np.zeros((181, 181))

        if hasattr(self, 'e'):
            for m in range(self.m):
                for theta in range(181):
                    for phi in range(181):
                        linear_array_x[theta, phi] += np.real(self.e[m, 0] * np.exp(1j * m * 2 * pi * self.dx *
                                                                                    np.sin((theta - 90) * pi / 180) *
                                                                                    np.cos(phi * pi / 180) + pi / 2))

            for n in range(self.n):
                for theta in range(181):
                    for phi in range(181):
                        linear_array_y[theta, phi] += np.real(self.e[0, n] * np.exp(1j * n * 2 * pi * self.dy *
                                                                                    np.sin((theta - 90) * pi / 180) *
                                                                                    np.sin(phi * pi / 180)))
        else:
            for m in range(self.m):
                for theta in range(181):
                    for phi in range(181):
                        linear_array_x[theta, phi] += np.real(np.exp(1j * m * 2 * pi * self.dx *
                                                                     np.sin((theta-90) * pi / 180) *
                                                                     np.cos(phi * pi / 180) + pi / 2))

            for n in range(self.n):
                for theta in range(181):
                    for phi in range(181):
                        linear_array_y[theta, phi] += np.real(np.exp(1j * n * 2 * pi * self.dy *
                                                                     np.sin((theta-90) * pi / 180) *
                                                                     np.sin(phi * pi / 180)))
        for theta in range(181):
            for phi in range(181):
                af[theta, phi] = linear_array_x[theta, phi] * linear_array_y[theta, phi]

        return af

    @property
    def power_pattern(self):
        """ Calculate power pattern

        :return:
        Power pattern of planar array

        """
        pp = self.array_factor**2
        return pp

    @property
    def normalized_array_factor(self):
        """ Calculate the normalized array factor

        :return:
        Normalized array factor of the planar array

        """
        af = self.array_factor
        max_af = np.max(af)
        normalized_af = np.abs(af)/max_af
        return normalized_af

    @property
    def normalized_array_factor2(self):
        """ Calculate the normalized array factor as product of two linear arrays

        :return:
        Normalized array factor of the planar array as product of two linear arrays

        """
        af = self.array_factor2
        max_af = np.max(af)
        normalized_af = np.abs(af) / max_af
        return normalized_af

    @property
    def normalized_power_pattern(self):
        """ Calculate the normalized power pattern

        :return:
        Normalized power pattern of the linear array

        """
        naf = self.normalized_array_factor
        npp = 20 * np.log10(naf)
        return npp

    @property
    def normalized_power_pattern2(self):
        """ Calculate the normalized power pattern as product of two linear arrays

        :return:
        Normalized power pattern of the linear array as product of two linear arrays

        """
        naf = self.normalized_array_factor2
        npp = 20 * np.log10(naf)
        return npp

    @property
    def gain_db(self):
        """ Calculate gain of the planar array

        :return:
        Gain of the planar array in dB

        """
        pp = self.power_pattern
        average_power = np.average(pp)
        max_power = np.max(pp)
        gain_db = 10 * np.log10(max_power/average_power)
        return gain_db

# End
