import pandas
from scipy.integrate import quad
import matplotlib.pyplot as plt

# Hubble function
def H(z):
    # provide the value for H0 in units of (km/s)/Mpc and convert it so s⁻¹
    H0 = 67.4
    H0 = H0*3.240779289*10**(-20)

    # value for Ωm
    Ωm = 0.315

    return H0*(Ωm*(1+z)**3 + 1-Ωm)**0.5

# luminosity distance
def dL(z, H):
    c = 9.715611890800001e-18  # speed of light [Gpc/s]
    return (1+z) * c * quad(lambda Z: 1/H(Z), 0, z)[0]

for file, label in zip(["data/ET-4-lowk.csv", "data/ET-4-highk.csv"], ["ET (low $\hat{k}$)", "ET (high $\hat{k}$)"]):
    header = pandas.read_csv(file, comment="#", nrows=0).columns.tolist()
    columns = pandas.read_csv(file, comment="#")

    redshifts = []
    sigmas = []
    for i in range(0, len(columns["redshift"])):
        redshift = columns["redshift"][i]
        dLobs = columns["luminosity_distance"][i]
        error = columns["error"][i]

        redshifts.append(redshift)
        sigmas.append(abs(dLobs - dL(redshift, H))/error)

    plt.scatter(redshifts, sigmas, label=label, zorder=3.5)

plt.legend()
plt.grid(alpha=0.5, zorder=0.5)
plt.xlabel("$z$")
plt.ylabel("$| d^{(obs)}_L(z) - d_L(z) | /\sigma^{(obs)}$")
plt.show()
