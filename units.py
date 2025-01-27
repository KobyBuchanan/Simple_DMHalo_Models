import numpy as np
from sympy import Symbol

#system units
G = 1
R_vir = 1
N = 100_000
radius_space = np.linspace(1 / N, 1, N)
radius_logspace = np.logspace(-5, 0, 100)
c = 15  # Milkyway concentration

#Sympy units
r_prime = Symbol('r')
C = Symbol('c')
const = Symbol('4*pi', real=True, positive=True)
rho0 = Symbol('rho0', real=True, positive=True)

#Characteristic Densities
rho0_iso = 1 / (4 * np.pi)
rho0_NFW = (c ** 3 / (4 * np.pi)) * (np.log(1 + c) - c / (1 + c)) ** -1


#Density Functions
def rho_iso(r):
    return rho0 * (1 / r) ** 2


def rho_NFW(r):
    return rho0 / (C * r * (1 + C * r) ** 2)