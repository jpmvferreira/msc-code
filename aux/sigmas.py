## sigmas.py
# computes the best, median and worst constrains based on how small are the sigma regions for each parameter

# to-do:
# - permitir escolha dos parametros e mandar so um aviso caso um seja encontrado no texto mas nao pedido
# - imprimir best/median/worst em coluna vertical pq na horizontal nao eh mto percetivel se forem mtos
# - output em .csv? posso dar uma flag -o e assim apresento os resultados bonitos a mesma

# imports
import matplotlib.pyplot as plt
from statistics import median
import argparse
import pandas

# read from CIs.tex
def readCIs(input):
    # initialize the dictionary
    d = {}
    d["ΔTotal"] = []

    # iterate over all input files
    for i in range(0, len(input)):
        folder = input[i]
        sigmas = []

        # filter out CIs.tex
        with open(f"{folder}/CIs.tex", "r") as file:
            content = file.read()
            fileheader, sigma1, sigma2 = content.split("\n\n")

            sigma1 = sigma1.split("\n")[3:-2]  # 1σ values
            sigma2 = sigma2.split("\n")[3:-3]  # 2σ values

            for j in range(0, len(sigma1)):
                processed = sigma1[j].replace(" ", "").split("$")
                parameter = processed[1]
                value = processed[3]

                if f"Δ{parameter}" not in d.keys():
                    d[f"Δ{parameter}"] = []

                # symmetric region
                if "\pm" in value:
                    value, area = value.split("\pm")
                    value = float(value)
                    area = float(area)

                # asymmetric region
                else:
                    value, sigmalow, sigmahigh = value.replace("^{+", " ").replace("}_{-", " ").replace("}","").split()
                    value = float(value)
                    sigmalow = float(sigmalow)
                    sigmahigh = float(sigmahigh)
                    area = sigmalow + sigmahigh

                sigmas.append(area)
                d[f"Δ{parameter}"].append(area)

            # compute the total area
            sigmatot = 1
            for k in sigmas:
                sigmatot *= k

        d[f"ΔTotal"].append(sigmatot)

    return d

# main
def main(args):
    # get arguments
    input = args.input
    custom_names = eval(args.custom_names) if args.custom_names else None
    folder_names = args.folder_names
    sliced_names = args.sliced_names

    # run invalid arguments check
    if custom_names and (len(custom_names) != len (input)):
        raise Exception("Length between custom file names and input folders missmatch.")

    if bool(custom_names) + folder_names + sliced_names > 1:
        raise Exception("Several custom names provided, pick one.")

    # if folder names are requested
    if folder_names:
        custom_names = [folder.split("/")[-1] for folder in input]

    # if sliced names is requested
    if sliced_names:
        custom_names = [folder.split("/")[-1].split("_")[-1] for folder in input]

    # get data from CIs for each input folder
    dict = readCIs(input)

    # convert to panda dataframe for convenience
    df = pandas.DataFrame(data=dict, index=custom_names)
    #df.round()

    # start indexing in 1, if custom name is not provided
    if not custom_names:
        df.index = [i for i in range(1, len(df) + 1)]

    # print all Δ's
    print(df.to_string())
    print()

    # compute min/median/max
    d2 = {}
    dmin = df[df == df.min()]
    dmedian = df[df == df.median()]
    dmax = df[df == df.max()]

    for parameter in dict.keys():
        mini = " ".join( map(str, dmin[parameter].dropna().index.tolist()) )
        mediani = " ".join( map(str, dmedian[parameter].dropna().index.tolist()) )
        maxi = " ".join( map(str, dmax[parameter].dropna().index.tolist()) )

        d2[parameter] = [mini, mediani, maxi]

    mini = " ".join( map(str, dmin["ΔTotal"].dropna().index.tolist()) )
    mediani = " ".join( map(str, dmedian["ΔTotal"].dropna().index.tolist()) )
    maxi = " ".join( map(str, dmax["ΔTotal"].dropna().index.tolist()) )
    d2["ΔTotal"] = [mini, mediani, maxi]

    # pretty print min/median/max
    df2 = pandas.DataFrame(data=d2, index=["best", "median", "worst"])
    print(df2.to_string())

# run if called
if __name__ == "__main__":
    # create the parser and add arguments
    parser = argparse.ArgumentParser()

    # required arguments
    parser.add_argument("-i", "--input", nargs="*", help="Input folders", required=True)

    # custom labels
    parser.add_argument("-cn", "--custom-names", type=str, help="A string with a Python style list with the custom name for each input folder.")
    parser.add_argument("-fn", "--folder-names", action="store_true", help="Use the folder name instead of showing ordered indices. Warning: folder names must be unique.")
    parser.add_argument("-sn", "--sliced-names", action="store_true", help="Use a sliced version of the folder name (all characters after the last '_') instead of showing ordered indices.")

    # get arguments
    args = parser.parse_args()

    # call main with args
    main(args)
