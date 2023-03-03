## fQ-M=2-EM.py
# GW luminosity distance for f(Q) model with ΛCDM background: f(Q) = Q + MQ^½
# assumes h = 0.68, Ωm = 0.353 and M = 2.0311 (best fit from arXiv:2004.07867)

# imports
from scipy.optimize import fsolve
from scipy.integrate import quad
from math import exp

# description
description = "GW wave luminosity distance for f(Q) = Q + MQ^½ with h = 0.674, Ωm = 0.315 and M = 2.0331 (best fit from arXiv:2004.07867)"

# constants
c  = 9.715611890800001e-18        # speed of light [Gpc/s]

# model parameters
h  = 0.674                        # provide the value for h
Ωm = 0.315                        # matter density
M = 2.0331                        # value of M

# conver h to H0
H0 = 100*h*3.240779289*10**(-20)  # H0 = 100h [km s⁻¹ Mpc⁻¹] = 100h*(km to Mpc conversion) [s⁻¹]

# Hubble function
def H(z):
    return H0*(Ωm*(1+z)**3 + 1-Ωm)**0.5

# luminosity distance
def dL(z, H):
    return (1+z) * (c) * quad(lambda Z: 1/H(Z), 0, z)[0]
