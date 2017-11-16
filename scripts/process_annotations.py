#!/usr/bin/env python2

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF
from sys import stderr, stdout, exit
import logging
import csv
import re

PAT_CHR = re.compile('^(?:.*;)?chromosome=([^;]+)(?:;.*)?$', re.I)
PAT_NAME = re.compile('^(?:.*;)?Name=([^;]+)(?:;.*)?$', re.I)


LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)


# supported annotation formats
ANN_FORMAT = ['ENSEMBLE', 'GFF3']


def readEnsembleData(data):
    ''' parser for ensemble annotation file '''

    res = set()

    isHeader = True
    for line in csv.reader(data, delimiter=','):
        if isHeader:
            isHeader = False
            continue
       
        # chr_id, start, end, gene_id
        res.append((line[5].strip(), int(line[6]), int(line[7]), line[3].strip()))

    return sorted(res)


def readGFF3Data(data):
    ''' parser for GFF3 annotation file '''

    chrMap = dict()
    anns = set() 
    for line in csv.reader(data, delimiter='\t'):
        if line[0].startswith('#'):
            continue
        
        if line[2] == 'region':
            m = PAT_CHR.match(line[8])
            if m:
                chrMap[line[0]] = m.group(1)
            else:
                m = PAT_NAME.match(line[8])
            if m and not chrMap.has_key(line[0]):
                chrMap[line[0]] = m.group(1)

        elif line[2] == 'gene':

            name = line[0]
            m = PAT_NAME.match(line[8])
            if m:
                name = m.group(1)
            anns.add((line[0], int(line[3]), int(line[4]), name))

    anns = list(anns)
    anns.sort()
    i = 0
    while i < len(anns):
        anns[i] = (chrMap.get(anns[i][0], anns[i][0]), ) + anns[i][1:]
        # check for overlaps and remove smaller gene annotation
        if i > 0 and anns[i-1][0] == anns[i][0] and anns[i-1][2] >= anns[i][2]:
            del anns[i]
        else:
            i += 1
            
    # inefficient, but hey, data is not that large
    anns.sort()

    return anns
           
if __name__ == '__main__':
    
    # setting up the argument parser
    parser = ArgumentParser(description='Reads annotations in Ensemble or ' + \
            'GFF3 format and produces an annotation file that will be used' + \
            ' in constructing Hi-C graphs', formatter_class=ADHF)
    parser.add_argument('-f', '--format', default=ANN_FORMAT[0], type=str,
            choices=ANN_FORMAT,
            help='Supported formats of input annotation file ')
    parser.add_argument('annotation_file',  type=str, help='Annoation file')

    args = parser.parse_args()

    # setup logging
    ch = logging.StreamHandler(stderr)
    ch.setLevel(logging.ERROR)
    ch.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))
    LOG.addHandler(ch)
   
    
    annotations = None
    if args.format == 'ENSEMBLE':
        annotations = readEnsembleData(open(args.annotation_file))
    elif args.format == 'GFF3':
        annotations = readGFF3Data(open(args.annotation_file))
    else:
        LOG.fatal('Unknown annotation file format. Exiting')
        exit(1)

    # output annotations
    out = stdout 
    for ann in annotations:
        print >> out, '\t'.join(map(str, ann))

