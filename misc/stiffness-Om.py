from scipy.integrate import odeint
from scipy.special import lambertw
import matplotlib.pyplot as plt
from math import exp
import numpy as np

# u'(t)
def dudt(u, t, Omega_m, l):
    f1 = exp(-l*u**2) / (l*u*(u**(-2)-2*l) - u**(-3))
    f2 = 1.5*(Omega_m)*(1+t)**2

    return f1*f2

# u(t) using ODE solver
def ut(t, Omega_m, l):
    # initial conditions plus desired time to get the value of u(t)
    t = [0, t]
    u0 = 1

    sol = odeint(dudt, u0, t, args=(Omega_m, l))

    return float(sol[1])

def main():
    Omlist = np.linspace(0.1, 1, 100)

    redshift = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
    u_obs = np.array([1.0, 0.59, 0.34, 0.23, 0.16, 0.12, 0.10, 0.08, 0.07, 0.06])

    sums = []
    for Om in Omlist:
        sum = 0
        l = 0.5 + lambertw(-(Om)/(2*exp(0.5))).real
        for i in range(0, len(redshift)):
            sum += (ut(redshift[i], Om, l) - u_obs[i])**2
        sums.append(sum)

    mini = min(sums)
    for i in range(0, len(sums)):
        if sums[i] == mini:
            j = i

    plt.plot(Omlist, sums)
    plt.scatter(Omlist[j], mini, color="red")
    plt.show()

    return

if __name__ == "__main__":
    main()
