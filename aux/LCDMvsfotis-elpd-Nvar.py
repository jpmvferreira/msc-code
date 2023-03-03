import matplotlib.pyplot as plt
from statistics import mean
from math import sqrt
import numpy as np
import sys

#plt.rcParams.update({'font.size': 12})  # for smaller plots

def get_elpd(model, dataset):
    with open(f"output/Nvar/{model}_{dataset}/criteria/PSIS-LOO-CV.log", "r") as f:
        lines = f.readlines()

    for line in lines:
        line = line.split()
        if "elpd_loo" in line:
            elpd = float(line[1])
            SE = float(line[2])
        elif "(bad)" in line:
            badcount = int(line[3])
        elif "(very" in line:
            verybadcount = int(line[4])

    if badcount or verybadcount:
        print(f"{model}_{dataset}: bad = {badcount}, very bad = {verybadcount}", file=sys.stderr)

    return elpd, SE

# pick dataset and number of events per catalog, and the index of the catalogs, here!
Nlist = [5, 10, 15, 20, 25, 30]
nlist = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
dataset="LISA"

# Nlist = [500, 750, 1000, 1250, 1500, 1750]
# nlist = [1, 2, 3, 4, 5]
# dataset="ET"

models = ["LCDM", "fotis-noOmegar"]
colors = ["#006FED", "#E03424"]
ecolors = ["#76b6ff", "#ef9991"]
labels = ["$\Lambda$CDM", "f(Q) as Dark Energy"]

# (!) crappy hardcore for two models
fig, ax = plt.subplots()
x = np.arange(len(Nlist))  # the label locations
width = len(Nlist)/35     # the width of the bars

# plot average elpd for all n (catalog ID) with same N (number of events per catalog)
for i, model in enumerate(models):
    y = []
    yerr = []
    for N in Nlist:
        elpds = []
        SEs = []
        for n in nlist:
            elpd, SE = get_elpd(model, f"{dataset}-N{N}-{n}")
            elpds.append(elpd)
            SEs.append(SE)

        y.append( abs(mean(elpds)) )
        yerr.append( sqrt(sum([i**2 for i in SEs]))/len(nlist) )

    if model == "LCDM":
        rects1 = ax.bar(x - width/1.8, y, width, align="center", yerr=yerr, color=colors[i], ecolor="black", capsize=3, label=labels[i], zorder=3.5)
    else:
        rects2 = ax.bar(x + width/1.8, y, width, align="center", yerr=yerr, color=colors[i], ecolor="black", capsize=3, label=labels[i], zorder=3.5)

plt.grid(alpha=0.5, zorder=0.5)
ax.set_xlabel("N")
ax.set_ylabel(r"|<elpd>|")
ax.set_xticks(x)
ax.set_xticklabels(Nlist)
ax.legend()

fig.tight_layout()
plt.show()

# plot elpd as a function of n (catalog ID) for each N (number of events per catalog)
for N in Nlist:
    # (!) crappy hardcore for two models
    fig, ax = plt.subplots()
    x = np.arange(len(nlist))  # the label locations
    width = len(nlist)/55     # the width of the bars

    for i, model in enumerate(models):
        elpds = []
        SEs = []

        for n in nlist:
            elpd, SE = get_elpd(model, f"{dataset}-N{N}-{n}")
            elpds.append(abs(elpd))
            SEs.append(SE)

        if model == "LCDM":
            rects1 = ax.bar(x - width/1.9, elpds, width, align="center", yerr=SEs, color=colors[i], ecolor="black", capsize=1.5, label=labels[i], zorder=3.5)
        else:
            rects2 = ax.bar(x + width/1.9, elpds, width, align="center", yerr=SEs, color=colors[i], ecolor="black", capsize=1.5, label=labels[i], zorder=3.5)

    plt.grid(alpha=0.5, zorder=0.5)
    ax.set_xlabel("n")
    ax.set_ylabel(r"|elpd|")
    ax.set_xticks(x)
    ax.set_xticklabels(nlist)
    ax.legend()

    fig.tight_layout()
    plt.show()
