#!/usr/bin/env python2

import numpy as np
import argparse as args

if __name__ == '__main__':
    # Read in Hi-C map

    parser = args.ArgumentParser(description="Changes the resolution of a Hi-C map.")
    parser.add_argument('N', metavar='New_Resolution', type=int, help="Resolution of the output Hi-C map.")
    parser.add_argument('M', metavar='Hi-C_Map', type=str,
                        help="Path of all Hi-C maps to consider for the normalization")

    argus = parser.parse_args()

    # We need to know weather we are dealing with an intra- or an inter-chromosomal Hi-C map. Intra-chromosomal maps do
    # have only one chromosome identifier.
    if len(argus.M.split("chr")) < 2:
        print("Chromosome identifier not set correctly.\n Are we dealing with inter- oder intra-chromosomal matrices?")
        exit(-1)
    elif len(argus.M.split("chr")) < 3:
        isInter = False
    else:
        chrs = argus.M.split("chr")

        fchr = chrs[1].split('.')[0]
        schr = chrs[2].split('.')[1]

        if fchr == schr:
            isInter = False
        else:
            isInter = True

    oldMap = np.genfromtxt(argus.M, dtype=str)

    print("Loading done\nResize Map")

    # Compute a resized matrix

    # Calculate the size ratio between old and new resolution. If the ratio is not an integer we have a problem.
    old_res = int(argus.M.split('.')[-2:][0])

    if old_res % argus.N != 0:
        print("Ratio between old and new resolution is not an integer.")
        exit(-1)
    else:
        ratio = old_res / argus.N

    # If the ratio == 1 there is nothing to resize
    if ratio == 1:
        print("Nothing to resize")
        exit(0)

    # Initialize the new map
    newMap = np.empty([len(oldMap) * ratio, len(oldMap[0]) * ratio], dtype='|S18')

    # Walk through the old map
    for i in range(len(oldMap)):
        for j in range(len(oldMap[i])):
            # Walk through the new map
            for k in range(ratio):
                for l in range(ratio):
                    # Intra-chromosomal maps should not have a value on the main diagonal
                    if not isInter and i == j and k == l:
                        newMap[i * ratio + k][j * ratio + l] = "NULL"
                    else:
                        newMap[i * ratio + k][j * ratio + l] = oldMap[i][j]

    print("Done\nSaving")

    # Save the map
    mapfile = open(argus.M + ".resized", "w")

    for i in range(len(newMap)):
        for j in range(len(newMap[i])):
            mapfile.write(newMap[i][j])

            # Check if we have read the end of the line
            if j + 1 != len(newMap[i]):
                mapfile.write("\t")
        # Check if we have reached the last line
        if i + 1 != len(newMap):
            mapfile.write("\n")
