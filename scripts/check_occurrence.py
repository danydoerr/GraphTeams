#!/usr/bin/env python2

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF
from sys import stderr, stdout, exit, maxsize
from os.path import basename
from itertools import chain
import logging
import csv

# allow for very large fields
__l = maxsize
while True:
    try:
        csv.field_size_limit(__l)
        break
    except OverflowError:
        __l /= 10


LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)


def readRealGeneClusters(data):
    res = list()
    for line in csv.reader(data, delimiter='\t'):
        if len(line) > 1:
            res.append((line[0], line[1].split(';')))
    return res


def readGeneTeams(data):
    res = list()

    isHeader = True
    for line in csv.reader(data, delimiter='\t'):
        if isHeader:
            isHeader = False
            continue
        res.append(line[0].split(';'))
    return res


def checkOccurrence(realGeneClusters, predictedGeneClusters, threshold):
    
    res = [False] * len(realGeneClusters)

    for i, rgc in enumerate(realGeneClusters):
        if not len(rgc):
            continue
        for pgc in predictedGeneClusters:
            res[i] = len(set(rgc).intersection(pgc))/float(len(rgc)) >= \
                    threshold
            if res[i]:
                break

    return res


if __name__ == '__main__':

    # setting up the argument parser
    parser = ArgumentParser(description='Checks occurrence of real gene ' + \
            'clusters in a set of predicted gene clusters',
            formatter_class=ADHF)
    parser.add_argument('-m', '--min_occurrence', default=1, type=float, 
            help='Percent of cluster members that must at least be covered ' + \
                    'in a gene cluster occurrence')
    parser.add_argument('real_gene_clusters', type=str, 
            help='File containing biological gene clusters')
    parser.add_argument('gene_teams', nargs='+',  type=str, 
            help='Files containing predicted gene clusters')

    args = parser.parse_args()

    if args.min_occurrence < 0 or args.min_occurrence > 1:
        LOG.fatal('Parameter min_occurrence must be a value between 0 and 1')
        exit(1)


    # setup logging
    ch = logging.StreamHandler(stderr)
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(logging.Formatter('%(asctime)s %(levelname)s: %(message)s'))
    LOG.addHandler(ch)

    realGCs = readRealGeneClusters(open(args.real_gene_clusters))

    realGCnames, realGCgenes = zip(*realGCs)

    res = list()
    for f in args.gene_teams:
        predGCs = readGeneTeams(open(f))
        res.append(checkOccurrence(realGCgenes, predGCs, args.min_occurrence))
    

    out = stdout

    print >> out, '\t'.join(chain(('gene cluster', ), map(basename,
        args.gene_teams)))
    for j in xrange(len(realGCs)):
        print >> out, '\t'.join(chain((realGCnames[j], ), map(lambda x: x[j] and
            '1' or '0', res)))
         

