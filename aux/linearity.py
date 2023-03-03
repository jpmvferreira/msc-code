## nome!
# descricao!

# imports
from scipy.stats import linregress
import matplotlib.pyplot as plt
import argparse
import os

# temp
from sigmas import readCIs

# main
def main(args):
    # get arguments
    input = args.input
    names = args.names.replace(" ", "").replace("['", "").replace("']", "").split("','")
    labels = eval(args.labels) if args.labels else [f"Batch {i}" for i in range(0, len(input))]
    regression = args.regression
    batch_regression = args.batch_regression
    print(names)

    colors = ["b", "g", "r", "c", "m", "y", "orange", "pink", "purple"]

    data1 = []
    data2 = []

    # main loop
    for i in range(0, len(input)):
        dict = readCIs(input[i])
        plt.scatter(dict[f"Δ{names[0]}"], dict[f"Δ{names[1]}"], zorder=3.5, color=colors[i], label=labels[i])

        if batch_regression:
            res = linregress(dict[f"Δ{names[0]}"], y=dict[f"Δ{names[1]}"])
            print(f"Linear regression for batch {labels[i]}")
            print(f"R² = {res.rvalue**2:.6f}")
            print(f"m = {round(res.slope, 4)}, b = {round(res.intercept, 4)}")
            print()
            plt.axline((0, res.intercept), slope=res.slope, linestyle="-", alpha = 0.15, color=colors[i])

        if regression:
            data1.extend(dict[f"Δ{names[0]}"])
            data2.extend(dict[f"Δ{names[1]}"])

    if regression:
        res = linregress(data1, y=data2)
        print(f"Linear regression for all data")
        print(f"R² = {res.rvalue**2:.6f}")
        print(f"m = {round(res.slope, 4)}, b = {round(res.intercept, 4)}")
        plt.axline((0, res.intercept), slope=res.slope, linestyle="-", alpha = 0.15, color="black")

    # fancy up and show
    plt.xlabel("$\sigma_{" + names[0] + "}$")
    plt.ylabel("$\sigma_{" + names[1] + "}$")
    plt.grid(alpha=0.5)
    plt.legend()
    plt.show()

# run if called
if __name__ == "__main__":
    # create the parser and add arguments
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", action="append", nargs="+", help="A batch of files to plot with the same color. Can be used multiple times to provide another batch of data. Each batch of data will share the same color on plot.", required=True)
    parser.add_argument("-n", "--names", type=str, help="A string with a Python style list with the LaTeX label of the parameters.", required=True)

    # custom labels
    parser.add_argument("-l", "--labels", type=str, help="A string with a Python style list with the label for each batch of dataset.")

    # regressoes
    parser.add_argument("-r", "--regression", action="store_true", help="Compute a linear regression with all provided data.")
    parser.add_argument("-br", "--batch-regression", action="store_true", help="Compute a linear regression for each batch of data.")

    # get arguments
    args = parser.parse_args()

    # call main with args
    main(args)
