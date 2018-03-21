#!/usr/bin/env python2

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF
import sys
import csv


def readHomologyFromEnsemble(data):
    res = list()

    isHeader = True
    for line in csv.reader(data, delimiter=','):
        if isHeader:
            isHeader = False
            continue
        res.append((line[0].strip(), line[3].strip()))
    return res

def mapTo(target, **source):

    sMap = [dict(map(reversed, s)) for s in source

#    for gene_id, 
    for ent in homList:
        hFile.write(humanGeneDict[ent[0]] + '\t' + ent[1]+ "\n")


if __name__ == '__main__':
    
    # setting up the argument parser
    parser = ArgumentParser(description='Reads annotations in Ensemble ' + \
            'format and produces a homology table', formatter_class=ADHF)

    parser.add_argument('annotation_file',  type=str, nargs='+',
            help='Annoation file')

    args = parser.parse_args()


    # output annotations
    out = stdout 
    for ann in annotations:
        print >> out, '\t'.join(map(str, ann))

