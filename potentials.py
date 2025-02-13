import numpy as np

G = 1


class SphericalPotential:
    """
    define density for all spheres using TwoPowerSphere
    Let r be spherical radius for now

    """

    def __init__(self, amp=1.0):
        self._amp = amp

    def __call__(self, r):
        return self._amp * self._evaluate(r)

    def den(self, r):
        return self._amp * self._den(r)

    def mass(self, r):
        return self._amp * self._mass(r)

    def vCirc(self, r):
        return np.sqrt(G * self.mass(r) / r)

    def radius_dis(self, r):
        return 4 * np.pi * self.den(r) * r ** 2


class Isothermal(SphericalPotential):

    def __init__(self):
        super().__init__()
        self.rho0 = 1 / (4 * np.pi)

    def _evaluate(self, r):
        return self.vCirc(r) ** 2 * np.log(r)

    def _den(self, r):
        return self.rho0 * (1/r)**2

    def _mass(self, r):
        return 4 * np.pi * self.rho0 * r


class NFW(SphericalPotential):

    def __init__(self, Rvir=1.0, conc=15):
        super().__init__()
        self.Rs = Rvir / conc
        self.conc = conc
        self.rho0 = (conc ** 3 / (4 * np.pi)) / (np.log(1 + conc) - conc / (1.0 + conc))

    def _evaluate(self, r):
        return -4 * np.pi * self.rho0 * np.log(1 + self.conc * r) / (self.conc ** 3 * r)

    def _den(self, r):
        return self.rho0 / ((r / self.Rs) * (1 + r / self.Rs)**2)

    def _mass(self, r):
        const = 4 * np.pi * self.rho0 * self.Rs ** 3
        return const * (np.log(1 + r / self.Rs) - (r/self.Rs) / (1 + r / self.Rs))
