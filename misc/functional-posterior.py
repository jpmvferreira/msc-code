## info:
# - ref do fgivenx: https://fgivenx.readthedocs.io/en/latest/fgivenx.html#fgivenx.drivers.plot_contours
# - ref sobre ppd: https://en.wikipedia.org/wiki/Posterior_predictive_distribution ou https://vasishth.github.io/bayescogsci/book/sec-ppd.html

## to-do:
# - ver qual a diferenca entre o que o fgivenx faz e a posterior predictive distribution (ppd)

# imports
import numpy
import matplotlib.pyplot as plt
from fgivenx import plot_contours, plot_lines, plot_dkl

# model
def f(x, theta):
    m, b = theta
    return m*x + b

# posterior samples
nsamples = 1000
ms = numpy.random.normal(loc=1, scale=0.25, size=nsamples)
bs = numpy.random.normal(loc=0, scale=0.25, size=nsamples)
samples = numpy.array([(m, c) for m, c in zip(ms, bs)]).copy()

# posterior samples for a "second dataset"
ms2 = numpy.random.normal(loc=2, scale=0.25, size=nsamples)
bs2 = numpy.random.normal(loc=0.5, scale=0.25, size=nsamples)
samples2 = numpy.array([(m, c) for m, c in zip(ms2, bs2)]).copy()

# set the x range to plot on
xmin, xmax = -2, 2
nx = 100
x = numpy.linspace(xmin, xmax, nx)

# set the cache
cache = 'build/fgivenx-cache/test'

# plot
fig, axes = plt.subplots()
axes.set_ylabel(r'$P(y|x)$')
axes.set_xlabel(r'$x$')
cbar = plot_contours(f, x, samples, axes, cache=cache, contour_color_levels=[0,1,2], lines=False)
#cbar2 = plot_contours(f, x, samples2, axes, cache=cache, contour_color_levels=[0,1,2], lines=False)
#cbar = plt.colorbar(cbar, ticks=[0,1,2])
#cbar.set_ticklabels(['',r'$1\sigma$',r'$2\sigma$'])

plt.plot(x, f(x, (1,0)), label="Best-fit", color="black")

plt.grid(alpha=0.5)
axes.set_xlim(left=xmin, right=xmax)
axes.set_axisbelow(True)
plt.legend()
plt.show()
