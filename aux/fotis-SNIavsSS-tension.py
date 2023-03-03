# imports
import matplotlib.patches as mpatches
from getdist import plots, MCSamples
import matplotlib.pyplot as plt
import numpy as np
import yaml
import h5py

# flatten total steps (i.e. remove chain information) and remove warmup to a numpy array of size [steps, ndim]
def getflatsamples(samples, warmup, chains, ndim, totalsteps):
    flatsamples = np.empty([samples*chains, ndim])
    for i in range(ndim):
        start = 0
        for j in range(chains):
            flatsamples[start::chains, i] = totalsteps[warmup:, j, i]
            start += 1

    return flatsamples

# get all samples
# always leave pantheon binned to last because of the awful workaround
folders = ["fotis-noOmegar_LISA-15", "fotis-noOmegar_LISA-5", "fotis-noOmegar_ET-4", "fotis-noOmegar-SNIa_pantheon-binned"]
legends = ["LISA (median)", "LISA (best)", "ET", "Pantheon (binned)"]

mcsamples = []
for folder in folders:
    # get config
    with open(f"output/{folder}/archive/config.yml", "r") as file:
        yml_loaded = yaml.full_load(file)
        names = yml_loaded["names"]
        labels = yml_loaded["labels"]
        samples = yml_loaded["samples"]
        warmup = yml_loaded["warmup"]
        chains = yml_loaded["chains"]

    # open chain.h5
    h5file = h5py.File(f"output/{folder}/chain.h5", "r")

    # get the samples, flatten them (drop chain information), and convert it to an MCSamples object
    flatsamples = getflatsamples(samples, warmup, chains, len(names), h5file["chain"])
    mcsamples.append(MCSamples(samples=flatsamples, names=names, labels=labels))

# init and config plot
markers = None
contour_alpha = 0.5
colors = ["#006FED", "#E03424", "#008000", "#9c5500", "#9224e0", "#ed00e6", "#f2e400", "#00f2e4", "#6fd95f"]
contour_colors = ["#006FED", "#E03424", "#008000", "#9c5500", "#9224e0", "#ed00e6", "#f2e400", "#00f2e4", "#6fd95f"]
g = plots.get_subplot_plotter()
g.settings.alpha_factor_contour_lines=contour_alpha
g.settings.line_labels = False

# plot samples
g.triangle_plot(mcsamples, filled=True, markers=markers, contour_colors=contour_colors, colors=colors, param_limits={"h": [0.67, 0.77], "Omega_m": [0.13, 0.37]})

# add horizontal region to plot
band_color = colors[len(folders)-1]
g.add_y_bands(0.3353, 0.013, color=band_color, ax=(1,0), zorder=0.5, alpha1=0.7, alpha2=0.35)

# awful workaround to create custom legend without being attached to axis
# required because the legend for the 1D SNIa (or the horizontal region) isn't showing
class AnyObject:
    pass

class AnyObjectHandler:
    def __init__(self, color):
        self.color = color

    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        patch = mpatches.Rectangle([x0, y0], width, height, facecolor=self.color)
        handlebox.add_artist(patch)
        return patch

l = [AnyObject() for i in range(0, len(folders))]

d = {}
for i in range(0, len(folders)-1):
    d[l[i]] = AnyObjectHandler(colors[i])
d[l[i+1]] = AnyObjectHandler(band_color)

g.fig.legend(l, legends, handler_map=d)

# show corner plot
plt.show()
