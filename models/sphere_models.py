import numpy as np
from models.parent_models import Sphere
from potentials.sphereical_potentials import *

data = np.loadtxt('./NFW_dispersion.dat')
data_r, data_sigma = data[:, 0], data[:, 1]


class Isothermal(Sphere):
    def __init__(self):
        super().__init__()
        pot = isothermal()
        self._pot = pot
        self.radius_dist = pot.radius_dis

    def sample_velocities(self, r=None, n=100_000):
        sigma = 1 / np.sqrt(3)
        v_x = np.random.normal(0, sigma, n)
        v_y = np.random.normal(0, sigma, n)
        v_z = np.random.normal(0, sigma, n)
        return v_x, v_y, v_z


class NFW(Sphere):
    def __init__(self, c=15):
        super().__init__()
        pot = nfw(conc=c)
        self._pot = pot
        self.radius_dist = pot.radius_dis

    def sample_velocities(self, r=None, n=100_000):
        if r is None:
            raise OSError("NFW model requires sampled radius as an input")
        sigma = np.interp(r, data_r, data_sigma)
        v_x = np.random.normal(0, sigma, n)
        v_y = np.random.normal(0, sigma, n)
        v_z = np.random.normal(0, sigma, n)
        return v_x, v_y, v_z
