#!/usr/bin/env python2

import sys
import argparse as args

#Reads in a cluster file and parses the relevant information into a dictionary
def readCFile(filename):
    cFile = open(filename, "r")

    #Dictionary for the cluster information
    cDict = {}

    for cluster in cFile:
        cols = cluster.split('\t')

        cSize = len(cols[1].split(';'))
        cMems = cols[0].split(';')

        #We do not want to evaluate the header file
        if len(cMems) > 1:
            #Check if a cluster of the same length as the current one was found before
            if cSize in cDict:
                cDict[cSize].append(cMems)
            else:
                cDict[cSize] = [cMems]

    return cDict

#Counts the number clusters in a dictionary
def countCs(cDict):
    num = 0

    for size in cDict.values():
        num += len(size)

    return num

#Calculates the average size of all clusters in the cluster dictionary
def calAvgCSize(cDict):
    allmems = 0.0
    num = 0.0

    for skey, cList in cDict.items():
        num += len(cList)
        allmems += float(skey) * len(cList)

    return allmems / num

#Setting up the argument parser

parser = args.ArgumentParser(description="Evaluates the found clusters of two result files.")
parser.add_argument('T', metavar='3D-Results', type=str, help="Path to the file that contains the results for the tree-dimensional data. Should be in csv format")
parser.add_argument('O', metavar='1D-Results', type=str, help="Path to the file that contains the results for the one-dimensional data. Should be in csv format")
parser.add_argument('R', metavar='Run Times', type=str, help="Path to the file that contains the run times of the algorithm.")
parser.add_argument('-d', '--delta', type=float, help="Delta which was chosen to produce the given results.")

arguments = parser.parse_args()

#Read in both files and store the clusters in dictionaries according to their size
threeDClusters = readCFile(arguments.T)
oneDClusters = readCFile(arguments.O)

#Read in run times
rTFile = open(arguments.R, "r")

# ignore header 
rTFile.readline()
# ignore memory, read runtime in HH:MM:SS format
rTRaw = rTFile.readline().rsplit('\t', 1)[1].strip()
HH, MM, SS = map(float, rTRaw.split(':'))
rT = 60 * HH + MM + SS/60

#Count the number of clusters
NumThreeD = countCs(threeDClusters)
NumOneD = countCs(oneDClusters)

#Compute the average size of the clusters
avgSizeThreeD = calAvgCSize(threeDClusters)
avgSizeOneD = calAvgCSize(oneDClusters)

#Count the number of 1D clusters that are subsets of bigger clusters in 3D
numSubsetCs = 0
avgKnowGain = 0.0

for cSize, clusterList in oneDClusters.items():
    currSize = int(cSize)

    for cluster in clusterList:
        notMem = False

        for threeDCSize, cList in threeDClusters.items():
            if int(threeDCSize) >= currSize:
                for clust in cList:
                    for mem in cluster:
                        if not mem in clust:
                            notMem =  True
                            break
                    if not notMem:
                        numSubsetCs += 1

                        avgKnowGain += float(threeDCSize) - float(currSize) 
                    else:
                        notMem = False

avgKnowGain /= float(numSubsetCs)

#Calculate the percentage of all subset clusters
perSubsetCs = float(numSubsetCs) / (float(NumOneD) / 100.0)

#print >> sys.stderr, numSubsetCs

#output all statistics
print arguments.delta, '\t', countCs(threeDClusters), '\t', countCs(oneDClusters), '\t', calAvgCSize(threeDClusters), '\t', calAvgCSize(oneDClusters), '\t', rT, '\t', avgKnowGain
