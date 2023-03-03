## LCDM-SS.py

# imports
from scipy.integrate import quad

# description
description = "GW wave luminosity distance for ΛCDM with Ωm = 0.271 h = 0.7063 (SS catalog LISA-N30-1 best fit)"

# Hubble function
def H(z):
    # define Hubble's constant using c/H0 = (2.9979 Gpc)/h and providing the value for h
    h = 0.7063
    H0 = 299792458*h/(2.9979*3.085678*10**25)

    # value for Ωm
    Ωm = 0.271

    return H0*(Ωm*(1+z)**3 + 1-Ωm)**0.5

# luminosity distance
def dL(z, H):
    c = 9.715611890800001e-18  # speed of light [Gpc/s]
    return (1+z) * c * quad(lambda Z: 1/H(Z), 0, z)[0]
