import numpy as np
from models.parent_models import Disk
from potentials.disk_potentials import *


class IsolatedExpDisk(Disk):
    def __init__(self):
        super().__init__()
        pot = EDNoHalo()
        self._pot = pot
        self.radius_dist = pot.radius_dis

    def sample_velocities(self, r=None, n=100_000):
        sigma_r_theta = 10 / 120
        sigma_z = 5 / 120
        v_r = np.random.normal(0, sigma_r_theta, n)
        v_theta = np.random.normal(0, sigma_r_theta, n)
        v_z = np.random.normal(0, sigma_z, n)
        v_theta += self._pot.vCirc(r)
        return v_r, v_theta, v_z


class MasslessDiskInHalo(Disk):
    def __init__(self):
        super().__init__()
        pot = EDHaloOnly()
        self._pot = pot
        self.radius_dist = pot.radius_dis

    def sample_velocities(self, r=None, n=100_000):
        sigma_r_theta = 10 / 120
        sigma_z = 5 / 120
        v_r = np.random.normal(0, sigma_r_theta, n)
        v_theta = np.random.normal(0, sigma_r_theta, n)
        v_z = np.random.normal(0, sigma_z, n)
        v_theta += self._pot.vCirc(r)
        return v_r, v_theta, v_z


class ED(Disk):
    def __init__(self):
        super().__init__()
        pot = DiskAndHalo()
        self._pot = pot
        self.radius_dist = pot.radius_dis

    def sample_velocities(self, r=None, n=100_000):
        sigma_r_theta = 10 / 120
        sigma_z = 5 / 120
        v_r = np.random.normal(0, sigma_r_theta, n)
        v_theta = np.random.normal(0, sigma_r_theta, n)
        v_z = np.random.normal(0, sigma_z, n)
        bulk_v_circ = self._pot.vCirc(r) + self._pot.halo_vCirc(r)
        v_theta += bulk_v_circ
        return v_r, v_theta, v_z
