# GW luminosity distance for a cosmological model from 2104.15123: f(Q) = Qe^(λQ₀/Q)
# uses best fit values from arXiv:2104.15123 with SNIa+CC+BAOs: Ωm = 0.353, Ωr = 0 and h = 0.68

# imports
from scipy.special import lambertw
from scipy.optimize import fsolve
from scipy.integrate import quad
from math import exp

# description
description = "GW wave luminosity distance for f(Q) = Qe^(λQ₀/Q) with Ωm = 0.353, Ωr = 0 and h = 0.68 (SNIa+CC+BAOs best-fit, see arXiv:2104.15123)"

# model parameters and constants
h  = 0.68                                    # provide the value for h
Om = 0.353                                   # matter density
H0 = 299792458*h/(2.9979*3.085678*10**25)    # H0 = 100h [km s⁻¹ Mpc⁻¹] = 100h*(km to Mpc conversion) [s⁻¹]
l  = 0.5 + lambertw(-Om/(2*exp(0.5))).real   # lambda
c  = 9.715611890800001e-18                   # speed of light [Gpc/s]

# auxiliary function to obtain E(z) for a given z numerically
def f(E, z):
    return (E**2 - 2*l)*exp(l/E**2) - Om*(1+z)**3

# Hubble function
def H(z):
    x0 = 10*z + 0.1                 # random inital value
    E = fsolve(f, x0, args=(z))[0]  # solve for E(z) numerically

    return H0*E

# luminosity distance
def dL(z, H):
    # standard luminosity distance
    distance = (1+z) * c * quad(lambda Z: 1/H(Z), 0, z)[0]

    # correction
    E = H(z)/H0
    correction = exp((l/2)*(1-1/E**2)) * ((1-l)/(1-l/E**2))**0.5

    return correction*distance
