import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from matplotlib.colors import LogNorm

from potentials import Isothermal, NFW
from modeling import Model

Iso = Isothermal()
NFW = NFW(conc=15, mvir=1)

model = Model(NFW, 1)
R = model.sample_space
r, _, __ = model.sample(100)
sigmas = model.sigmar(R)


plt.scatter(R, sigmas)
plt.xlabel('r')
plt.ylabel('sigmar')
plt.show()