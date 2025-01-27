from units import N
from model_tools import generate_halo
from plots import DataView

#want 6-d phase space for each particle in model
# m[i], x[i], y[i], z[i], vx[i], vy[i], vz[i]
#assume identical mass m[i] = M_total / Number of Particles
#postions are selected for each particle to be random from a uniform distribution
#velocities are added to each particle


def main():
    #generates a galaxy for a given model, values are tabulated (will implement astropy tables soon)
    table = generate_halo(N, "Isothermal")

    #class object that stores data viewing tools, like phase space plots
    dv = DataView(table)

    #An example tool, this prints the equilibrium of the system
    #here equilibrium is defined as 2T / |W| = 1
    dv.equilibrium_test()

    #An example plot
    dv.position_plot()


if __name__ == "__main__":
    main()