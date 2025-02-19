import numpy as np
import matplotlib.pyplot as plt
import time

from potentials import Isothermal, NFW
from modeling import Model

#want 6-d phase space for each particle in model
# m[i], x[i], y[i], z[i], vx[i], vy[i], vz[i]
#assume identical mass m[i] = M_total / Number of Particles
#postions are selected for each particle to be random from a uniform distribution
#velocities are added to each particle


def main():
    Iso = Isothermal()
    NFW = NFW()

    start = time.time()
    model = Model(Iso, 1)

    position_data, velocity_data, sigma = model.sample()
    r = position_data[0]
    v_x = velocity_data[0]
    v_y = velocity_data[1]
    v_z = velocity_data[2]
    velocity_vectors = v_x**2 + v_y**2 + v_z**2

    T = np.sum(0.5 * velocity_vectors)
    W = -np.sum(NFW.mass(r) / r)
    print("ratio = ", 2 * T / np.abs(W))
    end = time.time()
    print(end - start)

    sorted_index = np.argsort(r)
    r_sorted = r[sorted_index]
    sigma_sorted = sigma[sorted_index]
    plt.plot(r_sorted, sigma_sorted, '--')
    plt.ticklabel_format(useOffset=False)
    plt.show()

    start = time.time()
    model = Model(NFW, 1)

    position_data, velocity_data, sigma = model.sample()
    r = position_data[0]
    v_x = velocity_data[0]
    v_y = velocity_data[1]
    v_z = velocity_data[2]
    velocity_vectors = v_x**2 + v_y**2 + v_z**2

    T = np.sum(0.5 * velocity_vectors)
    W = -np.sum(NFW.mass(r) / r)
    print("ratio = ", 2 * T / np.abs(W))
    end = time.time()
    print(end - start)

    sorted_index = np.argsort(r)
    r_sorted = r[sorted_index]
    sigma_sorted = sigma[sorted_index]
    plt.plot(r_sorted, sigma_sorted, '--')
    plt.ticklabel_format(useOffset=False)
    plt.show()


if __name__ == "__main__":
    main()