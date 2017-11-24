#!/usr/bin/env python2

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF
from sys import stdout, stderr, exit
from itertools import chain
from os.path import basename
from random import sample
import logging

from nearest_neighbor_go_scores import readAssociations, readGO, readGenes, \
        getAllRootPaths, constructGOTree, constructLevelMap, \
        nearestNeighborDist

LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)


def sample_nn(t, links, levels, nn_genome, genesChr, cluster_size, no_of_samples):
    """ sample GO nearest-neighbor distances  """

    all_genes = [g for g in chain(*genesChr.values()) if links.has_key(g)]
    res = list()
    for _ in xrange(no_of_samples):
        sg = sample(all_genes, cluster_size)
        sgc = nearestNeighborDist(t, links, levels, sg)
        res.append(sum(sgc[g]-nn_genome[g] for g in sg))
    res.sort()
    return res
        

if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('-c', '--compress', action='store_true', 
            help='compress output')
    parser.add_argument('go_obo_file', type=str, 
            help='Gene Ontology hierarchy in OBO format')
    parser.add_argument('go_associations', type=str,
            help='Gene Ontology assocations file that links GO terms ' + \
                    'with genes')
    parser.add_argument('genome_annotations', type=str, 
            help='genome annotation file')
    parser.add_argument('cluster_size', type=int,
            help='size of the gene cluster (measured in annotated genes)')
    parser.add_argument('no_of_samples', type=int,
            help='number of draws from the sample pool')

    args = parser.parse_args()

    # setup logging
    ch = logging.StreamHandler(stderr)
    ch.setLevel(logging.INFO)
    ch.setFormatter(logging.Formatter('%(levelname)s\t%(asctime)s\t%(message)s'))
    LOG.addHandler(ch)


    # main 
    LOG.info('loading GO OBO file')
    goData = readGO(open(args.go_obo_file))
    LOG.info('loading GO associations')
    assoData = readAssociations(open(args.go_associations))
    LOG.info('laoding gene annotations')
    genesChr = readGenes(open(args.genome_annotations))
    genes = set(chain(*genesChr.values()))
    
    pathDict = {gene: list(chain(*[map(tuple, getAllRootPaths(rep, goData)) for
        rep in reps])) for gene, reps in assoData.items() if gene in genes}

    LOG.info('found %s out of %s in GO hierarchy' %(len(pathDict), len(genes)))

    LOG.info('construcing GO hierarchy')
    t, links = constructGOTree(pathDict)
    LOG.info('top-down assignment of levels')
    levels = constructLevelMap(t)

    LOG.info('sampling %s clusters of size %s ' %(args.no_of_samples,
        args.cluster_size))
    nearest_gene_dists = nearestNeighborDist(t, links, levels, genes)
    samples = sample_nn(t, links, levels, nearest_gene_dists, genesChr,
        args.cluster_size, args.no_of_samples)
   
    out = stdout
    if args.compress:
        from gzip import GzipFile
        out = GzipFile(fileobj=stdout, mode='w')
    print >> out, '\t'.join(map(str, samples))
    out.close()

    LOG.info('done')
