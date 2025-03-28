import numpy as np


class Sphere:
    def __init__(self):
        self.sample_space = np.linspace(0.001, 1)
        self.type = "is_sphere"

    def sample(self, n=100_000):
        r = self._sample_r(n=n)
        theta, phi = self._sample_position_angles(n=n)
        return r, theta, phi

    def _sample_r(self, n):
        M_samps = self._cum_mass_func(self.sample_space)
        random_mass_values = np.random.uniform(0, 1, n)
        M_interp = np.interp(random_mass_values, M_samps, self.sample_space)
        return M_interp

    def _cum_mass_func(self, r):
        return self._pot.mass(r) / self._pot.mass(1)

    def _sample_position_angles(self, n):
        theta = np.arcsin(np.random.uniform(-1, 1, size=n))
        phi = np.random.uniform(0, 2 * np.pi, size=n)
        return theta, phi


class Disk:
    def __init__(self):
        self.sample_space = np.linspace(0.001, 1)
        self.type = "is_cylindrical"

    def sample(self, n=100_000):
        r = self._sample_r(n=n)
        theta, z = self._sample_position(n=n)
        return r, theta, z

    def _sample_r(self, n):
        M_samps = self._cum_mass_func(self.sample_space)
        random_values = np.random.uniform(0, 1, n)
        r_samples = np.interp(random_values, M_samps, self.sample_space)
        return r_samples

    def _cum_mass_func(self, r):
        return self._pot.mass(r) / self._pot.mass(1)

    def _sample_position(self, n):
        theta = np.random.uniform(0, 2 * np.pi, size=n)
        z = np.random.normal(0, 0.3, size=n)
        return theta, z
