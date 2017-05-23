#!/usr/bin/env python

import sys
import argparse as args

#Function to create a graphml element that represents a weighted edge between two nodes
def cEdge(iD, orig, dest, w):
	return "    edge [\n        source " + str(orig) + "\n        target " + str(dest) + "\n        weight " + str(w) + "\n    ]\n"
#	return '<edge id="e' + str(iD) + '" source="n' + str(orig) + '" target="n' + str(dest) + '">\n  <data key="k1">' + str(w) + '</data>\n</edge>\n'

#Setting up the argument parser

parser = args.ArgumentParser(description="Creates a file in graphml format for the gene teams algorithm")
parser.add_argument('A', metavar='AnnotationFile', type=str, help="Path to the file that contains all gene annotations. Should be in gff3 format")
parser.add_argument('H', metavar='HomologyTable', type=str, help="Path to the file that contains names of homologous genes of a gene family. Should be a tab separated file in which each line is a gene family and each column contains the names of homologous genes of a species.")
parser.add_argument('T', metavar='TargetFile', type=str, help="Path to the file in which the graph sould be written to.")
parser.add_argument('M', metavar='Hi-C_Map', type=str, nargs='+', help="Path to the file that contains a Hi-C map. Should be a tab separated file consisting of three colums containing the pair of bins and their count.")
parser.add_argument('-s', '--binSize', default=0, type=int, help="Size of bin of the Hi-C map, i.e. its resolution.")
parser.add_argument('-o', '--mapOrder', type=str, nargs='+', help="Order of Hi-C maps which are given as parameters.")
#Might be not needed anymore
#parser.add_argument('-w', '--weightConstant', default=1.0, type=float, help="Constant to adapt three-dimensional distances. Default is 1.")

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
#oFile.write('<?xml version="1.0" encoding="UTF-8"?>\n<graphml xmlns="http://graphml.graphdrawing.org/xmlns" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd">\n<key id="k0" for="node" attr.name="class" attr.type="string"/>\n<key id="k1" for="edge" attr.name="weight" attr.type="int"/>\n<graph id="G" edgedefault="undirected">\n')

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
	#TODO: Check if this is necessary
	homologies[homs[0]] = homs

hFile.close()

#Reading in the annotation file
geneBins = {}
#chrShifts = {}
#lHomGene = []

#curShift = 0

#Counter for unique node and edge ids
nC = 0
eC = 0

aFile = open(arguments.A, "r")

#For each line in the annotation file...
for ann in aFile:
	ann = ann.replace("\n", "")

	#split the line by its columns
	infos = ann.split("\t")

#	aType = infos[2]

#	#check if the current entry is a "chromosome entry"
#	if fName == "chromosome":
#		#store the shift that has to be added to the position of each gene of this
#		#chromosome
#		chrShifts[infos[0]] = curShift

#		#add the length of the chromosome for the next shift
#		curShift = int(infos[4])

	#check if the current entry is a gene
#	if aType == "gene":
	#TODO: From annotation file to annotation file it may vary how the gene name is 
	#stored. Therefore, it has to be adapted here!
	fName = infos[3]

	name = None

	#check if the gene is homologous
	for fam in homologies.values():
		if fName in fam:
			name = fam[0]

			#create node for gene
			oFile.write("    node [\n        id " + str(nC) + "\n        label " + str(nC) + '\n        class "' + name + '"\n    ]\n')
#				oFile.write('<node id=n' + str(nC) + '">\n  <data key="k0">' + name + '</data>\n</node>\n')

#				#check if the last homologous gene lied on the same chromosome
#				if lHomGene == info[0]:
#					#calculate the distance between the consecutive genes
#					dist = int(info[3]) - int(lHomGene[4])

#					if dist <= 0:
#						print >> sys.stderr, "Distanz <= 0 gefunden."

#					#create edge between both genes
#					oFile.write('<edge id="e' + str(eC) + '" source="n' + str(nC - 1) + '" target="n' + str(nC - 2) + '">\n  <data key="k1">' + str(dist) + '/<data>\n</edge>\n')

#					#increment edge counter
#					eC += 1

#				#update last gene information
#				lHomGene = info

			#TODO: A gene is assigned to the bin according to its start. Is that the right way?
			#find out which bin corresponds to the current gene
			binNum = str((int(infos[1]) / int(arguments.binSize)) + 1)

#			#calculate middle of the gene 
			mid = int(infos[1]) + ((int(infos[2]) - int(infos[1])) / 2)
#			ref = int(infos[1]) + mid
#
#			if int(binNum) != (int(infos[1]) + mid) / int(arguments.binSize):
#				ref = (int(binNum) + 1) * arguments.binSize

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

		#each map file seems to have a tab at the end of each line therefore we have to delete it first TODO: This should not be the case anymore...
		#cols.pop()

		for j in range(len(cols)):
			#we first consider the entries on the main diagonal
			if lCounter  - 1 == j and str(j + 1) in geneBins:
				genes = geneBins[str(j + 1)]

				for k in range(len(genes)):
					if genes[k][0] == arguments.mapOrder[2 * i]:
						for l in range(len(genes)):
							if k < l and genes[l][0] == arguments.mapOrder[(2 * i + 1)]:
#								#New approach to resolve gene which lie in the same bin
#								binSize = int(arguments.binSize)
#								#default in case we cannot resolve better
#								w = 1.0
#								if genes[l][2] % binSize <= binSize / 2 and j > 0:
#									share = genes[l][2] % binSize
#
#									ratio = (float(share) / binSize)
#									w = ratio * 1.0 + (1.0 - ratio) * float(cols[j - 1])
#
#								if genes[l][2] % binSize > binSize / 2 and j + 1 < len(cols):
#									share = genes[l][2] % binSize
#
#									ratio = (float(share) / binSize)
#
#									w = ratio * 1.0 + (1.0 - ratio) * float(cols[j + 1])

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
								#TODO: Mal gucken, ob das nachher in der Praxis funktioniert
								#TODO: "...and j < len(cols)..." WTF?!
								if ogene[0] == arguments.mapOrder[(2 * i + 1)] and j < len(cols) and cols[j] != "NULL" and float(cols[j]) > 0:

									w = float(cols[j])

									oFile.write(cEdge(eC, gene[1], ogene[1], w))

									eC += 1

								if ogene[0] == arguments.mapOrder[(2 * i + 1)] and j < len(cols) and cols[j] == "NULL" and lCounter == j:
									dbp = float(cols[j-1]) / float(arguments.binSize)
									w = dbp * (abs(gene[2] - ogene[2]))

									oFile.write(cEdge(eC, gene[1], ogene[1], w))

									eC += 1
								if ogene[0] == arguments.mapOrder[(2 * i + 1)] and j < len(cols) and cols[j] != "NULL" and float(cols[j]) == 0 and lCounter == j:
									dbp = float(cols[j-1]) / float(arguments.binSize)
									w = dbp * (abs(gene[2] - ogene[2]))

									oFile.write(cEdge(eC, gene[1], ogene[1], w))

									eC += 1

#TODO: Verfahren fuer interchromosomale maps einbauen!
'''
#For each bin...
for b, genes in geneBins.items():
	k = 0

	#for each common gene in that bin
	for gene in genes:
		#create edges between genes in the same bin
		for l in range(len(genes)):
			#check if genes lie on the same chromosome and are different
			if gene[0] == genes[l][0] and l > k:
				cEdge(eC, gene[1], genes[l][1], 1.0)

		k += 1

		#get order of Hi-C maps
		mOrd = arguments.mapOrder

		#find all Hi-C maps that harbor information about this gene	
		for i in range(len(mOrd)):
			if mOrd[i] == gene[0]:
				if i % 2 == 0:
					mFile = open(arguments.M[i/2], "r")

					#find all pairs with that bin
					for pair in mFile:
						pair = pair.replace("\n", "")

						cells = pair.split('\t')

						if cells[0] == b:
							#find all genes that are located at the
							#second bin
							partners = geneBins[cells[1]]

							#create an edge if the second bin is on the 
							#right chromosome
							for partner in [j for j in partners if j[0] == mOrd[i+1]]:
								if int(cells[2]) > 0:
									#calculate the edge weight
									w = contCounts[gene[0]][cells[1]] / float(cells[2])

									#create new edge
									oFile.write(cEdge(eC, gene[1], j[1], w))

									#increment edge counter
									eC += 1

					mFile.close()
				else:
					mFile = open(arguments.M[(i/2)-1], "r")

					#find all pairs with that bin
					for pair in mFile:
						pair = pair.replace("\n", "")

						cells = pair.split('\t')

						if cells[1] == b:
							#find all genes that are located at the
							#second bin
							if cells[0] in geneBins:
								partners = geneBins[cells[0]]

							#create an edge if the second bin is on the 
							#right chromosome
							for partner in [j for j in partners if j[0] == mOrd[i]]:
								if int(cells[2]) > 0:
									#calculate the edge weight
									w = contCounts[gene[0]][cells[0]] / float(cells[2])

									#create new edge
									oFile.write(cEdge(eC, gene[1], j[1], w))

									#increment edge counter
									eC += 1
					mFile.close()
'''
#Finish the output file
oFile.write("]")
#oFile.write("</graph>\n</graphml>")
oFile.close()
