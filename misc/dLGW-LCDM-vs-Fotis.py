## info:
# - compares the best-fit for the gravitational wave luminosity distance for ΛCDM and fotis-noOmegar

# imports
from scipy.special import lambertw
from scipy.optimize import fsolve
from scipy.integrate import quad
from math import exp
import gwcatalog as gwc
import matplotlib.pyplot as plt
import numpy as np

## fotis
def f(E, z, Om):
    l  = 0.5 + lambertw(-Om/(2*exp(0.5))).real   # compute lambda twice :(

    return (E**2 - 2*l)*exp(l/E**2) - Om*(1+z)**3

def E_fotis(z, Om):
    x0 = 10*z + 0.1                     # random inital value
    E = fsolve(f, x0, args=(z, Om))[0]  # solve for E(z) numerically

    return E

def dL_fotis(z, h, Om):
    H0 = 299792458*h/(2.9979*3.085678*10**25)    # H0 = 100h [km s⁻¹ Mpc⁻¹] = 100h*(km to Mpc conversion) [s⁻¹]
    l  = 0.5 + lambertw(-Om/(2*exp(0.5))).real   # lambda
    c  = 9.715611890800001e-18                   # speed of light [Gpc/s]

    # standard luminosity distance
    distance = (1+z) * (c/H0) * quad(lambda Z: 1/E_fotis(Z, Om), 0, z)[0]

    # correction
    E = E_fotis(z, Om)
    correction = exp((l/2)*(1-1/E**2)) * ((1-l)/(1-l/E**2))**0.5

    return correction*distance

## ΛCDM
def H_LCDM(z, h, Ωm):
    H0 = 299792458*h/(2.9979*3.085678*10**25)
    return H0*(Ωm*(1+z)**3 + 1-Ωm)**0.5

def dL_LCDM(z, h, Ωm):
    c = 9.715611890800001e-18
    return (1+z) * c * quad(lambda Z: 1/H_LCDM(Z, h, Ωm), 0, z)[0]


## plot
redshifts = np.linspace(0, 9, 100)
LCDM = []
fotis = []
fotis_pantheon = []

for z in redshifts:
    LCDM.append(dL_LCDM(z, 0.6995, 0.290))
    fotis.append(dL_fotis(z, 0.6989, 0.257))

plt.plot(redshifts, fotis, label="$f\,(Q)$ as Dark Energy (ET best-fit)")
plt.plot(redshifts, LCDM, label="$\Lambda$CDM (ET best-fit)")

eventredshifts, distances, errors = gwc.load("data/ET-4.csv")
plt.errorbar(eventredshifts, distances, yerr=errors, fmt=".", markersize=7.5, color="purple", ecolor="purple", elinewidth=1, capsize=2, label="ET", zorder=1.5)

plt.grid()
plt.legend()
plt.xlabel("redshift")
plt.ylabel("luminosity distance [Gpc]")
plt.show()
