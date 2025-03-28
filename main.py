from models.sphere_models import *
from models.disk_models import *
from plot_tools.dataview import DataView

#want 6-d phase space for each particle in model
# m[i], x[i], y[i], z[i], vx[i], vy[i], vz[i]
#assume identical mass m[i] = M_total / Number of Particles
#postions are selected for each particle to be random from a uniform distribution
#velocities are added to each particle


def main():
    #Sphere gals
    gal_iso = Isothermal()
    gal_nfw = NFW()
    plots=DataView(gal_nfw)
    plots.plot()
    plots.radius_distribution()
    plots.energy_ratio()
    plots.plot_vCirc()

    #Disk gals
    gal_disk_inHalo = MasslessDiskInHalo()
    gal_disk = IsolatedExpDisk()
    gal = ED()
    plots = DataView(gal)
    plots.plot()
    plots.radius_distribution()
    plots.plot_vCirc()
    plots.energy_ratio()
if __name__ == "__main__":
    main()