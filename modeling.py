import numpy as np
from scipy.integrate import quad


class Model:
    def __init__(self, pot=None, Rmax=None):
        if pot is None:
            raise OSError("Potential Must be set")
        if Rmax is None:
            self.Rmax = 1
        self._pot = pot
        self.Rmax = Rmax
        self.sample_space = np.linspace(0.001, self.Rmax)
        self.radius_dist = pot.radius_dis

    def sample(self, n=100_000):
        r = self._sample_r(n=n)
        theta, phi = self._sample_position_angles(n=n)
        #----------Velocities-----------
        v_x, v_y, v_z, sigma = self._sampleVel(r, n=n)
        return (r, theta, phi), (v_x, v_y, v_z), sigma

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

    def sigmar(self, r):
        integrand = lambda x: self._pot.den(x) * self._pot.mass(x) / x ** 2
        integral = quad(integrand, a=r, b=np.inf, limit=100)[0]
        return np.sqrt(integral / self._pot.den(r))

    def sigma_r_interp(self, r):
        sigma_values = np.array([self.sigmar(radius) for radius in self.sample_space])
        return np.interp(r, self.sample_space, sigma_values)

    def _sampleVel(self, r, n):
        sigma = self.sigma_r_interp(r)
        v_x = np.random.normal(0, sigma, n)
        v_y = np.random.normal(0, sigma, n)
        v_z = np.random.normal(0, sigma, n)
        return v_x, v_y, v_z,sigma
