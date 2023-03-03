with open('data/ET-fotis-4.csv', 'r') as t1, open('output/fotis-noOmegar_ET-Fotis-4/archive/data/ET-Fotis-4.csv', 'r') as t2:
    fileone = t1.readlines()
    filetwo = t2.readlines()

for line in filetwo:
    if line not in fileone:
        print(line)
