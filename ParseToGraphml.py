#!/usr/bin/env python

import sys
import argparse as args

#Function to create a graphml element that represents a weighted edge between two nodes
def cEdge(iD, orig, dest, w):
	return "    edge [\n        source " + str(orig) + "\n        target " + str(dest) + "\n        weight " + str(w) + "\n    ]\n"

#Setting up the argument parser

parser = args.ArgumentParser(description="Creates a file in graphml format for the gene teams algorithm")
parser.add_argument('A', metavar='AnnotationFile', type=str, help="Path to the file that contains all gene annotations. Should be in gff3 format")
parser.add_argument('H', metavar='HomologyTable', type=str, help="Path to the file that contains names of homologous genes of a gene family. Should be a tab separated file in which each line is a gene family and each column contains the names of homologous genes of a species.")
parser.add_argument('T', metavar='TargetFile', type=str, help="Path to the file in which the graph sould be written to.")
parser.add_argument('M', metavar='Hi-C_Map', type=str, nargs='+', help="Path to the file that contains a Hi-C map. Should be a tab separated file consisting of three colums containing the pair of bins and their count.")
parser.add_argument('-s', '--binSize', default=0, type=int, help="Size of bin of the Hi-C map, i.e. its resolution.")
parser.add_argument('-o', '--mapOrder', type=str, nargs='+', help="Order of Hi-C maps which are given as parameters.")

arguments = parser.parse_args()

#Check if bin size is set
if arguments.binSize <= 0:
	print >> sys.stderr, "Bin sizes not or not correctly set."
	exit(-1)

#Check if the complete order of Hi-C maps is given
if len(arguments.mapOrder) != 2 * len(arguments.M):
	print >> sys.stderr, "Order of Hi-C maps not stated stated correctly."
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

#Reading in the annotation file
geneBins = {}

#Counter for unique node and edge ids
nC = 0
eC = 0

aFile = open(arguments.A, "r")

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
			#NOTE: A gene is assigned to the bin according to its start.
			#find out which bin corresponds to the current gene
			binNum = str((int(infos[1]) / int(arguments.binSize)) + 1)

			#calculate middle of the gene 
			mid = (int(infos[2]) - int(infos[1])) / 2

			#save to which bin the gene corresponds
			if binNum in geneBins:
				geneBins[binNum].append((infos[0], nC, mid))
			else:
				geneBins[binNum] = [(infos[0], nC, mid)]

			#increment node counter
			nC += 1

aFile.close()

#Read in Hi-C map
maps = arguments.M

#Dictionary to store amount of contacts for each region
contCounts = {}
for i in range(len(maps)):
	mFile = open(maps[i], "r")

	#a counter for save the number of lines already read in (because this is the bin number)
	l = 0

	#For each pair auf regions...
	for pair in mFile:
		l += 1

		#get rid of the newline
		pair = pair.replace('\n', '')

		#spilt by column
		cell = pair.split("\t")

		#each line has a tab at the end therefore we need to delete it first
		cell.pop()

		#first chromosome the current map belongs to
		fChr = arguments.mapOrder[2*i]

		#second chromosome the current map belongs to
		sChr = arguments.mapOrder[2*i + 1]
 
		#sum up all entries of this line
		overallContacts = 0
		for j in range(len(cell)):
			#check if entry is "NULL"
			if cell[j] == "NULL":
				overallContacts += 0
			else:
				#print "l:", l, "j:", j, "cell[j + 1]:", cell[j]
				overallContacts += float(cell[j])

			#Only consider the bins of each column if the map is from different 
			#chromosomes
			if fChr != sChr:
				curChr = sChr
				#check if there was another pair with this entry before
				if curChr in contCounts:
					if str(j) in contCounts[curChr]:
						#add the contacts to the previously observed ones
						contCounts[curChr][str(j)] += float(cell[j])
					else:
						#add a new entry for the current bin
						contCounts[curChr][str(j)] = float(cell[j])
				else:
					#add a chromosome with the current bin
					contCounts[curChr] = {str(j): float(cell[j])}

		curChr = fChr
		
		#check if there was another map for this chromosome before (s.t. there is already an entry)
		if curChr in contCounts:
			if str(l) in contCounts[curChr]:
				#add the contacts to the previously observed ones
				contCounts[curChr][str(l)] += overallContacts
			else:
				#add a new entry for the current bin
				contCounts[curChr][str(l)] = overallContacts
		else:
			#add a chromosome with the current bin
			contCounts[curChr] = {str(l): overallContacts}

	mFile.close()

#Add edges to the graph
for i in range(len(arguments.M)):
	mFile = open(arguments.M[i], "r")

	#counter for lines that are already read in to determine the current bin
	lCounter = 0

	for line in mFile:
		lCounter += 1

		line = line.replace("\n", "")

		cols = line.split('\t')

		for j in range(len(cols)):
			#we first consider the entries on the main diagonal
			if lCounter  - 1 == j and str(j + 1) in geneBins:
				genes = geneBins[str(j + 1)]

				for k in range(len(genes)):
					if genes[k][0] == arguments.mapOrder[2 * i]:
						for l in range(len(genes)):
							if k < l and genes[l][0] == arguments.mapOrder[(2 * i + 1)]:

								#calculate distance of 1bp
								dbp = float(cols[j]) / float(arguments.binSize)
								w = dbp * (abs(genes[k][2] - genes[l][2]))

								oFile.write(cEdge(eC, genes[k][1], genes[l][1], w))

								eC += 1

			#now we only want to consider all entries above the main diagonal since  
			#entries below the main diagonal are redundant
			if lCounter - 1 < j and str(j + 1) in geneBins:
				for gene in geneBins[str(j + 1)]:
					if gene[0] == arguments.mapOrder[2 * i]:
						if str(lCounter) in geneBins:
							for ogene in geneBins[str(lCounter)]:
								if ogene[0] == arguments.mapOrder[(2 * i + 1)] and j < len(cols) and cols[j] != "NULL" and float(cols[j]) > 0:
									#if we use a precalculated distance matrix we do not need to compute anything here
									w = float(cols[j])

									oFile.write(cEdge(eC, gene[1], ogene[1], w))

									eC += 1

#TODO: Verfahren fuer interchromosomale maps einbauen!

#Finish the output file
oFile.write("]")
#oFile.write("</graph>\n</graphml>")
oFile.close()
