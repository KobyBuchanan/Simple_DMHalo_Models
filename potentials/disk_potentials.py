import numpy as np
from scipy import special
from potentials.sphereical_potentials import nfw

G = 1

class DiskPotential:
    def __init__(self):
        self.mass_disk = 1
        self.Rd = 0.2
        self.I0 = lambda x: special.i0(x / (2 * self.Rd))
        self.I1 = lambda x: special.i1(x / (2 * self.Rd))
        self.K0 = lambda x: special.k0(x / (2 * self.Rd))
        self.K1 = lambda x: special.k1(x / (2 * self.Rd))

    def __call__(self, r):
        return self._evaluate(r)

    def den(self, r):
        return self._suf_den(r)

    def mass(self, r):
        return self._mass(r)

    def vCirc(self, r):
        const = self.mass_disk * r ** 2 / (2 * self.Rd ** 3)
        bessles = self.I0(r) * self.K0(r) - self.I1(r) * self.K1(r)
        return np.sqrt(const * bessles)

    def radius_dis(self, r):
        return 2 * np.pi * r * self.den(r)


class EDNoHalo(DiskPotential):
    def __init__(self):
        super().__init__()
        self.Sig0 = self.mass_disk / (2 * np.pi * self.Rd ** 2)

    def _evaluate(self, r):
        const = np.pi * G * self.Sig0
        return const * r * (self.I1(r) * self.K0(r) - self.I0(r) * self.K1(r))

    def _suf_den(self, r):
        return self.Sig0 * np.exp(-r / self.Rd)

    def _mass(self, r):
        return 2 * np.pi * self.Sig0 * self.Rd ** 2 * (1 - np.exp(-r / self.Rd) * ((r / self.Rd) + 1))


class EDHaloOnly(DiskPotential):
    def __init__(self):
        super().__init__()
        self.halo = nfw(conc=100)
        self.Sig0 = self.mass_disk / (2 * np.pi * self.Rd ** 2)

    def _evaluate(self, r):
        return self.halo._evaluate(r)

    def _suf_den(self, r):
        return self.Sig0 * np.exp(-r / self.Rd)

    def _mass(self, r):
        return 2 * np.pi * self.Sig0 * self.Rd ** 2 * (1 - np.exp(-r / self.Rd) * ((r / self.Rd) + 1))

    def vCirc(self, r):
        return self.halo.vCirc(r)


class DiskAndHalo(DiskPotential):
    def __init__(self):
        super().__init__()
        self.halo = nfw(conc=15)
        self.Sig0 = self.mass_disk / (2 * np.pi * self.Rd ** 2)

    def _evaluate(self, r):
        const = np.pi * G * self.Sig0
        disk_pot = const * r * (self.I1(r) * self.K0(r) - self.I0(r) * self.K1(r))
        halo_pot = self.halo._evaluate(r)
        return disk_pot + halo_pot

    def _suf_den(self, r):
        return self.Sig0 * np.exp(-r / self.Rd)

    def _mass(self, r):
        return 2 * np.pi * self.Sig0 * self.Rd ** 2 * (1 - np.exp(-r / self.Rd) * ((r / self.Rd) + 1))

    def halo_vCirc(self, r):
        return self.halo.vCirc(r)