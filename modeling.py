import numpy as np
from scipy.integrate import trapezoid


class Model:
    def __init__(self, pot=None, Rmax=None):
        if pot is None:
            raise OSError("Potential Must be set")
        if Rmax is None:
            self.Rmax = 1
        self._pot = pot
        self.Rmax = Rmax
        self.sample_space = np.linspace(0.01, self.Rmax, 1000)
        self.radius_dist = pot.radius_dis

    def sample(self, n=100):
        r = self._sample_r(n=n)
        theta, phi = self._sample_position_angles(n=n)
        return r, theta, phi

    def _sample_r(self, n, batch=1000):
        samples = np.array([])
        while len(samples) < n:
            x_rand = np.random.uniform(low=0, high=1, size=batch)
            y_rand = np.random.uniform(low=0, high=np.max(self.radius_dist(self.sample_space)), size=batch)
            accepted = x_rand[y_rand < self.radius_dist(x_rand)]
            samples = np.concatenate((samples,accepted))
        return samples[:n]

    def _sample_position_angles(self, n):
        theta = np.arcsin(np.random.uniform(-1, 1, size=n))
        phi = np.random.uniform(0, 2 * np.pi, size=n)
        return theta, phi

    def sigmar(self, r):
        integrand = lambda x: self._pot.den(x) * self._pot.mass(x) / x ** 2
        integrals = [integrand(x) for x in r]
        integral = trapezoid(integrals, r)
        return np.sqrt((1 / self._pot.den(r) * integral))
