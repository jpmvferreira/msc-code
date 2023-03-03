# imports
from scipy.special import lambertw
from scipy.optimize import fsolve
from scipy.integrate import quad
from math import exp, log
import numpy as np
import matplotlib.pyplot as plt
import pandas

## Fotis model
# auxiliary function to obtain E(z) for a given z numerically
def f(E, z, Om):
    l  = 0.5 + lambertw(-Om/(2*exp(0.5))).real   # lambda
    return (E**2 - 2*l)*exp(l/E**2) - Om*(1+z)**3

# Hubble function
def Efotis(z, Om):
    x0 = 10*z + 0.1                 # random inital value
    E = fsolve(f, x0, args=(z, Om))[0]  # solve for E(z) numerically

    return E

# luminosity distance
def dLfotis(z, Om):
    # standard luminosity distance
    distance = (1+z) * quad(lambda Z: 1/Efotis(Z, Om), 0, z)[0]

    return distance


## ΛCDM
def ELCDM(z, Ωm):
    return (Ωm*(1+z)**3 + 1-Ωm)**0.5

def dLLCDM(z, Ωm):
    return (1+z) * quad(lambda Z: 1/ELCDM(Z, Ωm), 0, z)[0]


## get SNIa data
with open("data/pantheon-binned.csv", "r") as file:
    columns = pandas.read_csv(file, comment="#")
    zcmb = columns["zcmb"].tolist()
    mb = columns["mb"].tolist()
    dmb = columns["dmb"].tolist()

plt.errorbar(zcmb, mb, yerr=dmb, fmt=".", markersize=7.5, color="purple", ecolor="purple", elinewidth=1, capsize=2, label="Pantheon (binned)", zorder=3.5)

## get Mscript best-fit for each model
def getMscript(dL, Om):
    data = pandas.read_csv("data/pantheon-binned.csv", comment="#")
    B = 0
    C = 0
    for i in range(0, len(data["zcmb"])):
        z = data["zcmb"][i]
        mobs = data["mb"][i]
        error = data["dmb"][i]

        delta = mobs - 5*log(dL(z, Om), 10)
        B += delta/error**2
        C += 1/error**2

    return B/C

MscriptLCDM = getMscript(dLLCDM, 0.284)
print(f"Mscript for LCDM = {MscriptLCDM}")
Mscriptfotis = getMscript(dLfotis, 0.3354)
print(f"Mscript for fotis = {Mscriptfotis}")

## get 5log(D_L(z))
zlist = np.linspace(0.014, 1.7, 1000).tolist()
distancesfotis = []
distancesLCDM = []
for z in zlist:
    distancesfotis.append(5*log(dLfotis(z, 0.3354), 10) + Mscriptfotis)
    distancesLCDM.append(5*log(dLLCDM(z, 0.284), 10) + MscriptLCDM)

plt.plot(zlist, distancesfotis, label="f(Q) as Dark Energy w/ Pantheon (binned)")
plt.plot(zlist, distancesLCDM, label="$\Lambda$CDM w/ Pantheon (binned)")

plt.grid()
plt.legend()
plt.xlabel("z")
plt.ylabel("magnitude")
plt.show()
