# imports
import matplotlib.pyplot as plt
import argparse

# mian
def main(args):
    # get arguments
    input = args.input
    labels = eval(args.labels) if args.labels else [f"Batch {i}" for i in range(0, len(input))]
    PSIS_LOO_CV = args.PSIS_LOO_CV
    WAIC = args.WAIC
    AIC = args.AIC
    BIC = args.BIC
    DIC = args.DIC

    # pre-determined colors
    colors = ["b", "g", "r", "c", "m", "y", "orange", "pink", "purple"]

    # loop over groups
    n = 0  # nº of datasets already plotted
    for gid in range(0, len(input)):
        values = []
        errors = []

        # loop over folder of each group
        for folder in input[gid]:
            with open(f"{folder}/criteria.log", "r") as file:
                # header
                file.readline()
                file.readline()
                file.readline()

                # get right field
                if PSIS_LOO_CV:
                    text = file.readline().split(" ")
                file.readline()
                if WAIC:
                    text = file.readline().split(" ")
                file.readline()
                if AIC:
                    text = file.readline().split(" ")
                file.readline()
                if BIC:
                    text = file.readline().split(" ")
                file.readline()
                if DIC:
                    text = file.readline().split(" ")

            # append value and error
            values.append(float(text[1]))
            errors.append(float(text[3]))

        # plot each dataset for the current group
        plt.errorbar([j for j in range(n, n+len(input[gid]))], values, yerr=errors, fmt='.', color=colors[gid], ecolor=colors[gid], elinewidth=1.25, capsize=1.25, label=labels[gid])

        n += len(input[gid])  # increase nº of plotted datasets

    # fancy up and show
    plt.legend()
    plt.grid()
    plt.xticks([])
    plt.ylabel("PSIS-LOO-CV" if PSIS_LOO_CV else "WAIC" if WAIC else "AIC" if AIC else "BIC" if BIC else "DIC")
    plt.show()

# run if called
if __name__ == "__main__":
    # create the parser and add arguments
    parser = argparse.ArgumentParser()

    # (future group) required arguments
    parser.add_argument("-i", "--input", action="append", nargs="+", help="A batch of files to plot with the same color. Can be used multiple times to provide another batch of data. Each batch of data will share the same color on plot.", required=True)

    # (future group) custom labels
    parser.add_argument("-l", "--labels", type=str, help="A string with a Python style list with the label for each batch of dataset.")

    # (future grupo) criteria -> one must be required (do pre-run checks)
    parser.add_argument("--PSIS-LOO-CV", action="store_true", help="Plot the Pareto-smoothed importance sampling leave-one-out cross-validation")
    parser.add_argument("--WAIC", action="store_true", help="Plot the Watanabe-Akaike information criterion")
    parser.add_argument("--AIC", action="store_true", help="Plot the Akaike information criterion")
    parser.add_argument("--BIC", action="store_true", help="Plot the Bayesian information criterion")
    parser.add_argument("--DIC", action="store_true", help="Plot the deviance information criterion")

    # get arguments
    args = parser.parse_args()

    # call main with args
    main(args)
