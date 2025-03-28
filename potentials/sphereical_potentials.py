import numpy as np
G = 1

class SphericalPotential:

    def __init__(self):
        pass

    def __call__(self, r):
        return self._evaluate(r)

    def den(self, r):
        return self._den(r)

    def mass(self, r):
        return self._mass(r)

    def vCirc(self, r):
        return np.sqrt(G * self.mass(r) / r)

    def radius_dis(self, r):
        return 4 * np.pi * self.den(r) * r ** 2


class isothermal(SphericalPotential):

    def __init__(self):
        super().__init__()
        self.rho0 = 1 / (4 * np.pi)

    def _evaluate(self, r):
        return self.vCirc(r) ** 2 * np.log(r)

    def _den(self, r):
        return self.rho0 * (1 / r) ** 2

    def _mass(self, r):
        return 4 * np.pi * self.rho0 * r


class nfw(SphericalPotential):
    def __init__(self, conc=15):
        super().__init__()
        self.conc = conc
        self.rho0 = (conc ** 3 / (4 * np.pi)) * (np.log(1 + conc) - conc / (1.0 + conc)) ** -1

    def _evaluate(self, r):
        # -4 * np.pi * self.rho0 * self.Rs ** 3 * np.log(1 + r / self.Rs) / r
        return -4 * np.pi * self.rho0 * np.log(1 + r * self.conc) / (r * self.conc ** 3)

    def _den(self, r):
        return self.rho0 / (self.conc * r * (1 + self.conc * r) ** 2)

    def _mass(self, r):
        const = 4 * np.pi * self.rho0 / self.conc ** 3
        return const * (np.log(1 + r * self.conc) - ((r * self.conc) / (1 + r * self.conc)))