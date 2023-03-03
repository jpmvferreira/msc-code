## fotis.py
# EM luminosity distance for a cosmological model from 2104.15123: f(Q) = Qe^(λQ₀/Q)
# uses best fit values from the same paper, using real data (SNIa+CC+BAOs): Ωm = 0.353, λ = 0.371, h = 0.68 and Ωr = (1-2λ)e^λ - Ωm ≈ 0.0209

# TO-DO: implementar distancia luminosa eletromagnetica para o modelo do fotis! fotis disse no email que era a standard, por isso imagino que seja so isto tudo sem a correcao

# imports
from scipy.optimize import fsolve
from scipy.integrate import quad
from math import exp

# description
description = "EM wave luminosity distance for f(Q) = Qe^(λQ₀/Q) with Ωm = 0.353, λ = 0.371, h = 0.68 and Ωr = (1-2λ)e^λ - Ωm ≈ 0.0209 (SNIa+CC+BAOs best-fit, see arXiv:2104.15123)"

# model parameters and constants
c  = 9.715611890800001e-18        # speed of light [Gpc/s]
l  = 0.371                        # lambda
h  = 0.68                         # provide the value for h
H0 = 100*h*3.240779289*10**(-20)  # H0 = 100h [km s⁻¹ Mpc⁻¹] = 100h*(km to Mpc conversion) [s⁻¹]
Om = 0.353                        # matter density
Or = (1-2*l)*exp(l) - Om          # compute Ωr (≈ 0.0209) since it was not provided directly in the paper

# auxiliary function to obtain E(z) for a given z numerically
def f(E, z):
    return (E**2 - 2*l)*exp(l/E**2) - Om*(1+z)**3 - Or*(1+z)**4

# Hubble function
def H(z):
    x0 = 1                          # random inital value
    E = fsolve(f, x0, args=(z))[0]  # solve E(z) using scipy

    return H0*E

# luminosity distance
def dL(z, H):
    return (1+z) * c * quad(lambda Z: 1/H(Z), 0, z)[0]
