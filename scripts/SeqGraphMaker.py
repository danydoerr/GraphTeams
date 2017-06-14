#!/usr/bin/env python

import sys
import argparse as args
import numpy as np

#Function to create a graphml element that represents a weighted edge between two nodes
def cEdge(iD, orig, dest, w):
	return "    edge [\n        source " + str(orig) + "\n        target " + str(dest) + "\n        weight " + str(w) + "\n    ]\n"

#Setting up the argument parser

parser = args.ArgumentParser(description="Creates a file in graphml format for the gene teams algorithm")
parser.add_argument('S', metavar='Binsize', default=0, type=int, help="Size of bin of the Hi-C map, i.e. its resolution.")
parser.add_argument('A', metavar='AnnotationFile', type=str, help="Path to the file that contains all gene annotations. Should be in gff3 format")
parser.add_argument('H', metavar='HomologyTable', type=str, help="Path to the file that contains names of homologous genes of a gene family. Should be a tab separated file in which each line is a gene family and each column contains the names of homologous genes of a species.")
parser.add_argument('T', metavar='TargetFile', type=str, help="Path to the file in which the graph sould be written to.")
parser.add_argument('M', metavar='Hi-C_Map', type=str, nargs='+', help="Path to the file that contains a Hi-C map. Should be a tab separated file consisting of three colums containing the pair of bins and their count.")

arguments = parser.parse_args()

#Check if bin size is set
if arguments.S <= 0:
	print >> sys.stderr, "Bin sizes not or not correctly set."
	exit(-1)

#Write the header
oFile = open(arguments.T, "w")
oFile.write("graph [\n    directed 0\n")

#Reading in homology table
homologies = {}

hFile = open(arguments.H, "r")

#For each gene family...
for fam in hFile:
	#get rid of the newline character
	fam = fam.replace("\n", "")

	#split the line by column
	homs = fam.split("\t")

	#store the gene family in a dictionary
	homologies[homs[0]] = homs

hFile.close()

#Read in distances from HiC-Maps
#Dictionary that stores the different maps
hiCMaps = {}

for mAp in arguments.M:
	#find out which chromosome we are dealing with
	mapName = mAp.split('/')[-1]

	chrId = mapName.split('chr')[1].split('.')[0]

	hiCMaps[chrId] = np.genfromtxt(mAp, dtype=str)

#Counter for unique node and edge ids
nC = 0
eC = 0

#Read in annotation file

aFile = open(arguments.A, "r")

#dictionary that stores all genes according to their chromosome
geneDict = {}

#For each line in the annotation file...
for ann in aFile:
	ann = ann.replace("\n", "")

	#split the line by its columns
	infos = ann.split("\t")

	fName = infos[3]

	name = None

	#check if the gene is homologous
	for fam in homologies.values():
		if fName in fam:
			name = fam[0]

			#create node for gene
			oFile.write("    node [\n        id " + str(nC) + "\n        label " + str(nC) + '\n        class "' + name + '"\n    ]\n')

			#start of gene
			gstart = int(infos[1])
			gend = int(infos[2])

			#store genes in a dictionary according to their chromosome if we have data for that chromosome
			if infos[0] in hiCMaps:
				if infos[0] in geneDict:
					#a gene is stored as a tupel of node id, start and end
					geneDict[infos[0]].append((nC, gstart, gend))
				else:
					geneDict[infos[0]] = [(nC, gstart, gend)]

			nC += 1

#for each chromosome add edges between adjacent genes
for chrom in geneDict:
	geneList = geneDict[chrom]
	#Store the last gene processed
	lGene = None

	for gene in sorted(geneList, key=lambda x: x[1]):
		if lGene:
			#calculate middle of the genes 
			midlGene = lGene[1] + ((lGene[2] - lGene[1]) / 2)
			midGene = gene[1] + ((gene[2] - gene[1]) / 2)

			#could it happen that a gene appears in reversed orientation?
			if midlGene < 1 or midGene < 1:
				print "Wrong orientation"

			#calculate to which bins the genes belong to
			lGeneBin = lGene[1] / int(arguments.S)
			geneBin = gene[1] / int(arguments.S)

			#calculate the weight of the new edge and add it to the graph
			deltaBp = abs(midGene - midlGene)

			#find out where we have to look up the distance in the matrices
			if lGeneBin != geneBin and geneBin < len(hiCMaps[chrom]) and hiCMaps[chrom][lGeneBin][geneBin] != "NULL" and float(hiCMaps[chrom][lGeneBin][geneBin]) > 0:
				#w = float(deltaBp) * float(hiCMaps[chrom][lGeneBin][geneBin]) / float(arguments.S)
				w = float(hiCMaps[chrom][lGeneBin][geneBin])

				oFile.write(cEdge(eC, lGene[0], gene[0], w))
				eC += 1
			elif geneBin < len(hiCMaps[chrom]):
				w = float(deltaBp) * float(hiCMaps[chrom][lGeneBin][lGeneBin]) / float(arguments.S)

				oFile.write(cEdge(eC, lGene[0], gene[0], w))
				eC += 1

			lGene = gene
		else:
			lGene = gene

#Finish the output file
oFile.write("]")
oFile.close()
