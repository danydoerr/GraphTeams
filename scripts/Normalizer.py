#!/usr/bin/env python

import sys
import numpy as np
import argparse as args

parser = args.ArgumentParser(description="Creates distance matrices out of Hi-C maps.")
parser.add_argument('M', metavar='Hi-C_Map', type=str, nargs='+', help="Path of all Hi-C maps to consider for the "
                                                                       + "normalization")

argus = parser.parse_args()

mapList = [np.genfromtxt(i, dtype=str) for i in argus.M]

print("Loading done")

# compute the average of all maximum entries of all matrices
maxInMats = []

for mAp in mapList:
    entries = []
    diagEntrs = []

    for i in range(len(mAp)):
        for j in range(len(mAp[i])):
            if mAp[i][j] != "NULL":
                entries.append(float(mAp[i][j]))
                if j == i + 1:
                    diagEntrs.append(float(mAp[i][j]))

    matMax = max(entries)

    # Compute the average value on the diagonal above the main diagonal of every chromosome
    avgAbDiag = 0.0

    for entr in diagEntrs:
        avgAbDiag += entr

    avgAbDiag /= float(len(diagEntrs))

    maxInMats.append(matMax)

    # change entries to distances
    for i in range(len(mAp)):
        for j in range(len(mAp[i])):
            if mAp[i][j] != "NULL":
                mAp[i][j] = str(matMax + 1.0 - float(mAp[i][j]))
            # for the main diagonal of interchromosomal maps we do not need to insert artificial distances
            if i == j and len(mAp) == len(mAp[i]):
                mAp[i][j] = str(matMax + 1.0 - avgAbDiag)

# Calculate the average maximum over all chromosomes
maxSum = 0.0

for maxi in maxInMats:
    maxSum += maxi

avgMax = maxSum / float(len(maxInMats))

overallMax = max(maxInMats)

print("Begin normalization")

# Normalize each map
for i in range(len(mapList)):
    # get normalization constant NOTE: Taking the average here leads to smaller distances in the matrices. I personally
    # prefer this
    c = avgMax / maxInMats[i]

    # Normalization
    for k in range(len(mapList[i])):
        for l in range(len(mapList[i][k])):
            # we do not consider null entries
            if mapList[i][k][l] != "NULL":
                # store new value in map
                mapList[i][k][l] = str(c * float(mapList[i][k][l]))

print("Normalisation done")

# saving the matrices
for i in range(len(mapList)):
    mapfile = open(argus.M[i] + ".dmat", "w")

    for k in range(len(mapList[i])):
        for l in range(len(mapList[i][0])):
            mapfile.write(mapList[i][k][l])

            # check if we have read the end of the line
            if l + 1 != len(mapList[i][0]):
                mapfile.write("\t")
        # check if we have reached the last line
        if k + 1 != len(mapList[i]):
            mapfile.write("\n")
    mapfile.close()
