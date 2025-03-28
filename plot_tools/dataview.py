import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from matplotlib.colors import LogNorm
import numpy as np
from plot_tools.plot_helpers import convert_cartesian

star_masses = 1 / 100_000


class DataView:
    def __init__(self, model_=None):
        if model_ is None:
            raise OSError("Model Must be set")
        model = model_
        self._pot = model._pot

        # model data
        self._type = model.type
        self._radius_dist = self._pot.radius_dis

        #Init Galaxy Phase_space
        self._coords = model.sample()  # defaults to spherical
        self._vels = model.sample_velocities(r=self._coords[0])  # defaults to cartesian

    def plot(self):
        x, y, z = convert_cartesian(self._coords, self._type)
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(16, 8))
        ax1.hist2d(x, y, bins=(50, 50), norm=LogNorm(), cmap='jet')
        ax1.set_aspect('equal')
        ax1.set_xlabel('x (kpc)')
        ax1.set_ylabel('y (kpc)')
        ax2.hist2d(x, z, bins=(50, 50), norm=LogNorm(), cmap='jet')
        ax2.set_xlabel('x (kpc)')
        ax2.set_ylabel('z (kpc)')
        ax2.set_aspect('equal')
        ax3.hist2d(y, z, bins=(50, 50), norm=LogNorm(), cmap='jet')
        ax3.set_xlabel('y (kpc)')
        ax3.set_ylabel('z (kpc)')
        ax3.set_aspect('equal')
        plt.show()

    def radius_distribution(self):
        if self._type == 'is_sphere':
            y_label = '$4\\pi r^2 \\rho (r)$'
        if self._type == 'is_cylindrical':
            y_label = '$2\\pi r \\Sigma (r)$'

        radii = self._coords[0]
        radiusSpace = np.linspace(0.001, 1)
        plt.hist(radii, bins=100, density=True)
        plt.plot(radiusSpace, self._radius_dist(radiusSpace))
        plt.xlabel('r')
        plt.ylabel(y_label)
        plt.show()

    def energy_ratio(self):
        v1, v2, v3 = self._vels
        radii = self._coords[0]
        velocity_vector_squared = v1 ** 2 + v2 ** 2 + v3 ** 2
        T = np.sum(0.5 * star_masses * velocity_vector_squared)
        W = np.sum(star_masses * self._pot(radii))
        print("T = ", 2 * T)
        print("|W| = ", np.abs(W))
        print("ratio = ", 2 * T / np.abs(W))

    def plot_R_den(self):
        pass

    def plot_vCirc(self, Rmax=1):
        radiusSpace = np.linspace(0.001, Rmax)
        radius = self._coords[0]
        vcirc_radius_space = self._pot.vCirc(radiusSpace)
        vcirc = self._pot.vCirc(radius)
        plt.plot(radiusSpace, vcirc_radius_space, label='Model')
        plt.scatter(radius, vcirc, alpha=0.3, c='orange', label='Sampled')
        plt.xlabel('r (kpc)')
        plt.ylabel('$v \\ (kms^{-1})$')
        plt.xlim(left=0, right=Rmax)
        plt.legend()
        plt.show()
