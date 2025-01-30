import sympy
import time
from astropy.table import Table
from scipy.integrate import quad
from scipy.interpolate import interp1d
from sympy import integrate, lambdify
from units import *


def radii_sampling(function, num_samples, batch=1000):
    samples = []
    while len(samples) < num_samples:
        x_rand = np.random.uniform(low=0, high=1, size=batch)
        y_rand = np.random.uniform(low=0, high=max(function(radius_space)), size=batch)
        samples += x_rand[y_rand < function(x_rand)].tolist()
    return samples[:num_samples]


def generate_coordinates(n, func):
    Azimuth = np.random.uniform(0, 2 * np.pi, size=n)
    Elevation = np.arcsin(np.random.uniform(-1, 1, size=n))
    Radial_distance = radii_sampling(func, n)
    return Azimuth, Elevation, Radial_distance


def mass_enclosed(rho, model, subs=True):
    if model == "Isothermal":
        cden = rho0_iso
    if model == "NFW":
        cden = rho0_NFW
    func = const * rho(r_prime) * r_prime ** 2
    integral = integrate(func, (r_prime, 0, r_prime))
    if subs:
        integral = integral.subs(C, c).subs(const, 4 * np.pi).subs(rho0, cden)
    lambda_function = sympy.lambdify(r_prime, integral, 'numpy')
    return sympy.simplify(integral), lambda_function


def potential_function(model):
    if model == 'Isothermal':
        rho = rho_iso
        upper_bound = R_vir
    if model == 'NFW':
        rho = rho_NFW
        upper_bound = np.inf
    mass = mass_enclosed(rho, model, subs=False)[0]
    if model == "Isothermal":
        cden = rho0_iso
    if model == "NFW":
        cden = rho0_NFW
    func = -G * mass / r_prime ** 2
    integral = integrate(func, (r_prime, r_prime, upper_bound))
    integral = integral.subs(C, c).subs(const, 4 * np.pi).subs(rho0, cden)
    lambda_function = lambdify(r_prime, integral, 'numpy')
    return sympy.simplify(integral), lambda_function


def sigma_r_calc(r, den, mass):
    integrand = lambda x: den(x) * mass(x) / x ** 2
    integral = quad(integrand, r, np.inf)
    sig_square = (1 / den(r)) * integral[0]
    return np.sqrt(sig_square)


def jeans(rho, r):
    mass = mass_enclosed(rho, 'NFW')[1]
    den = rho(r_prime).subs(C, c).subs(rho0, rho0_NFW)
    density = lambdify(r_prime, den, 'numpy')
    sigma_values = [sigma_r_calc(i, density, mass) for i in r]
    return sigma_values


def generate_velocity(n, rho,radius, model):
    mass = mass_enclosed(rho, model)[0]
    circular_velocity = sympy.sqrt(G * mass / r_prime)
    v_1 = []
    v_2 = []
    v_3 = []
    V = []
    if model == "Isothermal":
        velocity_dispersion = np.sqrt(int(circular_velocity) / 3)
        velocity_vectors = np.random.normal(loc=0, scale=velocity_dispersion, size=(n, 3))
        for i in range(len(velocity_vectors)):
            v_1.append(velocity_vectors[i][0])
            v_2.append(velocity_vectors[i][1])
            v_3.append(velocity_vectors[i][2])
            V.append(np.sqrt(v_1[i] ** 2 + v_2[i] ** 2 + v_3[i] ** 2))

    if model == "NFW":
        sigmar = jeans(rho, radius_logspace)
        sigmar_interp = interp1d(radius_logspace,sigmar)
        velocity_vectors = np.random.normal(loc=0, scale=1, size=(n, 3))
        for i in range(len(velocity_vectors)):
            v_1.append(velocity_vectors[i][0] * sigmar_interp(radius[i]))
            v_2.append(velocity_vectors[i][1] * sigmar_interp(radius[i]))
            v_3.append(velocity_vectors[i][2] * sigmar_interp(radius[i]))

    return v_1, v_2, v_3, V


def velocity_transform(coords: tuple, vels: tuple):
    c1 = coords[0]
    c2 = coords[1]
    c3 = coords[2]
    v_1_prime = []
    v_2_prime = []
    v_3_prime = []
    V_prime = []
    for i in range(len(c3)):
        theta, phi = c1[i], c2[i]
        vr, vphi, vt = vels[0][i], vels[1][i], vels[2][i]
        v_1_prime.append(vr * np.sin(theta) * np.cos(phi) - vphi * np.sin(phi) + vt * np.cos(theta) * np.cos(phi))
        v_2_prime.append(vr * np.sin(theta) * np.sin(phi) + vphi * np.cos(phi) + vt * np.cos(theta) * np.sin(phi))
        v_3_prime.append(vr * np.cos(theta) - vt * np.sin(theta))
        V_prime.append(np.sqrt(v_1_prime[i] ** 2 + v_2_prime[i] ** 2 + v_3_prime[i] ** 2))
    return v_1_prime, v_2_prime, v_3_prime, V_prime


def generate_halo(n: int, model: str):
    if model == 'Isothermal':
        rho = rho_iso
        rho_lamb = lambdify(r_prime, rho(r_prime).subs(C, c).subs(rho0, rho0_iso), 'numpy')

    if model == 'NFW':
        rho = rho_NFW
        rho_lamb = lambdify(r_prime, rho(r_prime).subs(C, c).subs(rho0, rho0_NFW), 'numpy')

    distribution_func = lambda r: 4 * np.pi * rho_lamb(r) * r ** 2
    # generate N random coordinates and velocities for each particle
    Azimuth, Elevation, Radial_distance = generate_coordinates(n, distribution_func)
    # generate velocities for N particles, see generate_velocity for details
    v_1, v_2, v_3, V = generate_velocity(n, rho, Radial_distance, model)
    if model == 'NFW':
        v_1, v_2, v_3, V = velocity_transform((Elevation, Azimuth, Radial_distance), (v_1, v_2, v_3))

    # Gravitational potential function for current model
    potential = potential_function(model)[1]

    metadata = {
        "CTime": time.ctime(time.time()),
        "model": f'{model}',
        "SampFunc": distribution_func
    }
    data = Table(meta=metadata)
    data['Index'] = np.arange(0, n, 1)
    data['Mass'] = np.ones(n) * 1 / 10
    data['Azimuth'] = Azimuth
    data['Elevation'] = Elevation
    data['Radius'] = Radial_distance
    data['Velocity'] = V
    data['v_x'] = v_1
    data['v_y'] = v_2
    data['v_z'] = v_3
    data['W'] = data['Mass'] * potential(data['Radius'])
    data['T'] = 0.5 * data['Mass'] * data['Velocity'] ** 2

    return data

#------------------Testing--------------------
