## info:
# - tentativa de encontrar valor de h em LCDM que elimine a tensao no modelo do fotis
# - faco isso assumindo que LCDM concorda com SnIa (Om = 0.284), vou mandado valores de h aleatorios, gero a curva da distancia luminosa, e encontro o valor de {h, Om} no modelo do fotis que gera a curva da distancia luminosa de GW com o melhor da fit

import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt
from scipy.optimize import minimize, fsolve
from scipy.special import lambertw
from math import exp

# auxiliary function to obtain E(z) for a given z numerically
def f(E, z, Om, l):
    return (E**2 - 2*l)*exp(l/E**2) - Om*(1+z)**3

# Hubble function
def E(z, Om, l):
    x0 = 10*z + 0.1                     # random inital value
    Eval = fsolve(f, x0, args=(z, Om, l))[0]  # solve for E(z) numerically
    return Eval

# luminosity distance for fotis
def fotis_dL(z, h, Om):
    c  = 9.715611890800001e-18                   # speed of light [Gpc/s]
    H0 = 299792458*h/(2.9979*3.085678*10**25)    # H0 = 100h [km s⁻¹ Mpc⁻¹] = 100h*(km to Mpc conversion) [s⁻¹]
    l  = 0.5 + lambertw(-Om/(2*exp(0.5))).real   # lambda

    # standard luminosity distance
    distance = (1+z) * c * (1/H0) * quad(lambda Z: 1/E(Z, Om, l), 0, z)[0]

    # correction
    Eval = E(z, Om, l)
    correction = exp((l/2)*(1-1/Eval**2)) * ((1-l)/(1-l/Eval**2))**0.5

    return correction*distance

# luminosity distance for LCDM
def LCDM_dL(z, h, Om):
    H0 = 299792458*h/(2.9979*3.085678*10**25)
    c = 9.715611890800001e-18
    return (1+z) * c * quad(lambda Z: 1/(H0*(Om*(1+Z)**3 + 1-Om)**0.5), 0, z)[0]


print("LCDM:")
print(((0.284*(1+1)**3 + 1-0.284)**0.5)**2)

print("fotis:")
print(E(1, 0.335, 0.386156)**2)
crash.me

redshifts = np.linspace(0, 8, 100)
# fotis = [fotis_dL(i, 0.641703347982603, 0.335) for i in redshifts]
# LCDM = [LCDM_dL(i, 0.7, 0.285) for i in redshifts]
# plt.plot(redshifts, fotis, label="fotis")
# plt.plot(redshifts, LCDM, label="LCDM")
# plt.legend()
# plt.show()
# crash.me

# to-do: implementar o metodo de minimos quadrados algures
def residual(h, LCDM):
    global redshifts
    fotis = [fotis_dL(i, h, 0.335) for i in redshifts]

    S = 0
    for i in range(len(LCDM)):
        S += (LCDM[i] - fotis[i])**2

    return S

hLCDM = []
hfotis = []
Om = []
Omdiff = []
S = []
for h in np.linspace(0.62, 0.72, 6):
    LCDM = [LCDM_dL(i, h, 0.285) for i in redshifts]

    # minimize
    result = minimize(residual, [h], args=(LCDM))

    # guardar dados
    hLCDM.append(h)
    hfotis.append(float(result.x))
    S.append(residual(float(result.x), LCDM))

for l in hLCDM, hfotis, S:
    print(", ".join([str(i) for i in l]))
