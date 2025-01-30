import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from matplotlib.colors import LogNorm
import numpy as np
from units import radius_space


class DataView:
    def __init__(self, table):
        self.table = table
        self.model = self.table.meta['model']
        #Position data
        #Spherical coordinates
        self.azimuth = self.table['Azimuth']
        self.elevation = self.table['Elevation']
        self.radius = self.table['Radius']
        #Cylindrical coordinates
        #Cartesian coordinates
        self.x = self.radius * np.cos(self.azimuth) * np.cos(self.elevation)
        self.y = self.radius * np.sin(self.azimuth) * np.cos(self.elevation)
        self.z = self.radius * np.sin(self.elevation)
        #velocity data
        self.vx = self.table['v_x']
        self.vy = self.table['v_y']
        self.vz = self.table['v_z']
        self.V = self.table['Velocity']
        # Energy data
        self.T = self.table['T']
        self.W = self.table['W']
        #Shell density data
        self.mid = self.table['mid'][0:1000]
        self.ShDen = self.table['ShDen'][0:1000]

    def position_plot(self, coord: str = 'Cartesian'):
        if coord == 'Cartesian':
            fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
            ax1.hist2d(self.x, self.y, bins=(50, 50), norm=LogNorm(), cmap='jet')
            ax1.set_aspect('equal')
            ax1.set_xlabel('x (kpc)')
            ax1.set_ylabel('y (kpc)')
            ax2.hist2d(self.x, self.z, bins=(50, 50), norm=LogNorm(), cmap='jet')
            ax2.set_xlabel('x (kpc)')
            ax2.set_ylabel('z (kpc)')
            ax2.set_aspect('equal')
            ax3.hist2d(self.z, self.y, bins=(50, 50), norm=LogNorm(), cmap='jet')
            ax3.set_xlabel('z (kpc)')
            ax3.set_ylabel('y (kpc)')
            ax3.set_aspect('equal')
            fig.suptitle('Position Phasespace')
            plt.show()
        if coord == 'Spherical':
            pass
        if coord == 'Cylindrical':
            pass

    def velocity_plot(self, coord: str = 'Cartesian'):
        if coord == 'Cartesian':
            plt.figure()
            plt.hist2d(self.radius, self.vx, bins=(50, 50))
            plt.xlabel('R (kpc)')
            plt.ylabel('Vx (km/s)')
            plt.figure()
            plt.hist2d(self.radius, self.vy, bins=(50, 50))
            plt.xlabel('R (kpc)')
            plt.ylabel('Vy (km/s)')
            plt.figure()
            plt.hist2d(self.radius, self.vz, bins=(50, 50))
            plt.xlabel('R (kpc)')
            plt.ylabel('Vz (km/s)')
            plt.figure()
            plt.hist2d(self.radius, self.V, bins=(50, 50))
            plt.xlabel('R (kpc)')
            plt.ylabel('V (km/s)')
            plt.show()
        if coord == 'Spherical':
            pass
        if coord == 'Cylindrical':
            pass

    def equilibrium_test(self):
        T = np.sum(self.T)
        W = np.sum(self.W)
        print(2 * T / np.abs(W))

    def density_plot(self):
        density = self.table.meta['DenFunc']
        plt.loglog(self.mid, self.ShDen, label='model')
        plt.loglog(radius_space, density(radius_space), label='true')
        plt.legend()
        plt.title(f'{self.model} model')
        plt.ylabel('Density')
        plt.xlabel('Radius (kpc)')
        plt.show()

    def acceptance_rejection_plot(self):
        density_distribution = self.table.meta['SampFunc']
        plt.plot(radius_space, density_distribution(radius_space), label="Function", color='black')
        plt.hist(self.radius, bins=50, density=True, alpha=0.2, label="Sampled Distribution")
        plt.xlim(0, 1)
        plt.ylim(0, 2.5)
        plt.xlabel("Radius")
        plt.legend()
        plt.show()


