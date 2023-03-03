## ppd-linear-regression.py
# an attempt to get a posterior predictive distribution working
# i dont really know how it should looks like -> get integrating u retard
# lets do it for a linear regression y = ax + b

## notas:
# - na grid ele considera as entradas da list como o ponto do meio

# imports
from scipy.integrate import dblquad
import matplotlib.pyplot as plt
from numpy import linspace
from math import pi, exp
import numpy as np

# gauss auxiliary function
# cuz nobody seems to have coded it yet?!
def gauss(x, mu, sigma):
    return 1/((2*pi)**0.5 * sigma) * exp(-0.5*((x-mu)/sigma)**2)

# model
def model(x, a, b):
    return a*x + b

# likelihood for a single event
def likelihood(yobs, a, b, xobs, error):
    ytheo = model(xobs, a, b)
    return gauss(yobs, ytheo, error)

# made up posterior distribution
# a ~ N(2, 0.1), b ~ N(1, 0.2)
def posterior(a, b):
    pa = gauss(a, 2, 0.1)
    pb = gauss(b, 1, 0.2)
    return pa*pb

# bounds for a and b because i can't integrate to infinity
# 3 sigmas are enough
amin = 2 - 3*0.1
amax = 2 + 3*0.1
bmin = 1 - 3*0.2
bmax = 1 + 3*0.2

# grid for observations
N = 128
xmin = -10
xmax = 10
ymin = -20
ymax = 20
xgrid = linspace(xmin, xmax, N).tolist()
ygrid = linspace(ymin, ymax, N).tolist()

pgrid = np.zeros((N, N)).tolist()
for i in range(0, N):
    for j in range(0, N):
        xobs = xgrid[i]
        yobs = ygrid[j]
        error = 0.5 # erro da observacao -> nao me consigo livrar disto! claro, pq nao existe tal coisa como uma observacao sem erro! se eu tiver uma observacao com um erro gigante, eh mais provavel encontra-la super longe!

        # "swap" x and y -> indices here work the other way around!
        pgrid[j][i] = dblquad(lambda b,a: likelihood(yobs, a, b, xobs, error)*posterior(a, b), amin, amax, bmin, bmax)[0]

fig, ax = plt.subplots()
pcm = ax.pcolormesh(xgrid, ygrid, pgrid)
plt.colorbar(pcm, ax=ax)
plt.plot([xmin, xmax], [model(xmin, 2, 1), model(xmax, 2, 1)], label="y = 2*x + 1", color="red")
plt.xlabel("x")
plt.ylabel("y")
plt.xlim((xmin, xmax))
plt.ylim((ymin, ymax))
plt.grid(alpha=0.25)
plt.legend()
plt.show()
