# imports
import matplotlib.pyplot as plt
from scipy.integrate import quad
from math import floor, log, sqrt
from random import uniform, gauss
import numpy as np

# dist taken from LSST DDF in figure 12 of arXiv:1409.8562
dist = [64, 258, 480, 758, 1049, 1369, 1683, 2009, 1130]
def f(z):
    if z < 0.1 or z >= 1:
        return 0
    return dist[floor(z*10)-1]

# Hubble function
def H(z):
    # provide the value for H0 in units of (km/s)/Mpc and convert it so s⁻¹
    H0 = 67.36
    H0 = H0*3.240779289*10**(-20)

    # value for Ωm
    Ωm = 0.314

    return H0*(Ωm*(1+z)**3 + 1-Ωm)**0.5

# luminosity distance
def dL(z, H):
    c = 9.715611890800001e-18  # speed of light [Gpc/s]
    return (1+z) * c * quad(lambda Z: 1/H(Z), 0, z)[0]

# get N randomly generated events from a given distribution, using rejection methods
def GetRandom(distribution, x_min, x_max, y_min, y_max, N=1):
    counter = 0
    events = []

    while counter < N:
        x = uniform(x_min, x_max)
        y = uniform(y_min, y_max)

        if y < distribution(x):
            events.append(x)
            counter += 1

    return events

# compute redshift, luminosity distances and μ's
redshifts = GetRandom(f, 0.1, 1, 0, max(dist), N=8800)
distances = [dL(i, H)*1000 for i in redshifts]  # convert to Mpc
mu = [5*log(i, 10) + 25 for i in distances]

# eq (A.1) of arXiv:2007.14335
muerrors = [sqrt((gauss(0, 0.01)*i)**2 + 0.01**2 + 0.025**2 + 0.12**2) for i in redshifts]

# distribute the errors around the most likely value using a gaussian distribution
mu = [gauss(i, j) for i, j in zip(mu, muerrors)]

# theoretical line
redshift_theoretical = np.linspace(0.1, 1, num=250)
mu_theoretical = [5*log(dL(i, H)*1000, 10) + 25 for i in redshift_theoretical]

# plot and show
plt.plot(redshift_theoretical, mu_theoretical, color="grey", zorder=3.5)
plt.errorbar(redshifts, mu, yerr=muerrors, color="blue", fmt='.', ecolor="lightblue", elinewidth=1, capsize=0, markersize=1)
plt.xlim([0.1, 1])
plt.grid(alpha=0.5)
plt.xlabel("$z$")
plt.ylabel("$\mu(z)$")
plt.show()
