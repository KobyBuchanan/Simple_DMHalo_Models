import sympy
from astropy.table import Table
from scipy.integrate import quad
from sympy import integrate, lambdify
from units import *


def generate_coordinates(n):
    Azimuth = np.random.uniform(0, 2 * np.pi, size=n)
    Elevation = np.arcsin(np.random.uniform(-1, 1, size=n))
    Radial_distance = np.random.uniform(0, 1, size=n)
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


def generate_velocity(n, rho, model):
    mass = mass_enclosed(rho, model)[0]
    circular_velocity = float(sympy.sqrt(G * mass / r_prime))
    v_x = []
    v_y = []
    v_z = []
    V = []
    if model == "Isothermal":
        velocity_dispersion = np.sqrt(circular_velocity / 3)
        velocity_vectors = np.random.normal(loc=0, scale=velocity_dispersion, size=(n, 3))
        for i in range(len(velocity_vectors)):
            v_x.append(velocity_vectors[i][0])
            v_y.append(velocity_vectors[i][1])
            v_z.append(velocity_vectors[i][2])
            V.append(np.sqrt(v_x[i] ** 2 + v_y[i] ** 2 + v_z[i] ** 2))

    #NFW velocity not implemented
    if model == "NFW":
        velocity_vectors = np.random.normal(loc=0, scale=1, size=(n, 3))
        v_r = []
        v_phi = []
        v_theta = []
        for i in range(len(velocity_vectors)):
            sigma = 'Implement'
            v_r.append(velocity_vectors[i][0] * sigma)
            v_phi.append(velocity_vectors[i][1] * sigma)
            v_theta.append(velocity_vectors[i][2] * sigma)

    return v_x, v_y, v_z, V


def generate_halo(n: int, model: str):
    if model == 'Isothermal':
        rho = rho_iso
    if model == 'NFW':
        rho = rho_NFW

    data = Table()
    #generate N random coordinates and velocities for each particle
    Azimuth, Elevation, Radial_distance = generate_coordinates(n)
    #generate velocities for N particles, see generate_velocity for details
    v_x, v_y, v_z, V = generate_velocity(n, rho, model)
    #Gravitational potential function for current model
    potential = potential_function(model)[1]
    #----------------------------table-columns--------------------------------------------
    data['Index'] = np.arange(0, n, 1)
    data['Mass'] = np.ones(n) * 1 / 10
    data['Azimuth'] = Azimuth
    data['Elevation'] = Elevation
    data['Radius'] = Radial_distance
    data['Velocity'] = V
    data['v_x'] = v_x
    data['v_y'] = v_y
    data['v_z'] = v_z
    data['W'] = data['Mass'] * potential(data['Radius'])
    data['T'] = 0.5 * data['Mass'] * data['Velocity'] ** 2
    return data

#------------------Testing--------------------
