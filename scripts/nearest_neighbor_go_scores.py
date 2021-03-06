#!/usr/bin/env python2

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF
from itertools import chain, combinations, imap
from bisect import bisect
from sys import stdout, stderr, exit, maxsize
from cStringIO import StringIO
from os.path import basename
from random import sample
import logging
import csv

from newick_parser import Tree, Branch, Leaf


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

ROOT_NODE = 'GO:0008150'

def readAssociations(data):
    """ read GO gene association (*.gaf) files """
    res = dict()
    for row in csv.reader(data, delimiter='\t'):
        if len(row) < 11:
            continue
        if row[8] != 'P':
            continue
        goId = row[4]
        for i in row[10].split('|'):
            if i not in res:
                res[i] = list()
            if goId not in res[i]:
                res[i].append(goId)
    return res

def readTerm(data):
    """ helper file of readGO() """
    isTerm = False
    curItem = StringIO()
    for line in data:
        if isTerm:
            if not line or line == '\n':
                isTerm = False
                curItem.reset()
                yield curItem
                curItem.reset()
                curItem.truncate()
            curItem.write(line)
        elif line == '[Term]\n':
            isTerm = True
        # otherwise continue
    
def readGO(data):
    """ read GO data in OBO format """
    res = dict()
    ns_dict = dict()
    for term in readTerm(data):
        goId = None
        isARels = list()
        namespace = None
        for line in term:
            if line.startswith('id:'):
                goId = line[4:-1]
            elif line.startswith('namespace:'):
                namespace = line[11:-1]
            elif line.startswith('is_a:'):
                isARels.append(line[6:line.find(' !')])
        if isARels and namespace == 'biological_process': 
            res[goId] = isARels
        ns_dict[goId] = namespace

    # clean out everything except biological process:
    for k,v in res.items():
        if ns_dict[k] != 'biological_process':
            del res[k]
    return res

def getAllRootPaths(goId, goData):
    """ find all paths to the root for a given GO ID """
    if not goId:
        return list()
    if goId not in goData: 
        return [[goId]]

    res = list()
    for i in goData[goId]:
        paths = getAllRootPaths(i, goData)
        for path in paths:
            res.append([goId] + path)
    return res


def readGenes(data):
    """ read genes from annotation file """
    genes = dict()
    
    for line in csv.reader(data, delimiter='\t'):
        if not genes.has_key(line[0]):
            genes[line[0]] = set()
        genes[line[0]].add(line[3])
        
    return genes

def constructLevelMap(goTree):
    """ top-down assignment of levels for each node in the tree """

    levels = {goTree.subtree: 0}
    queue = [(goTree.subtree, 0)]

    while queue:
        node, l = queue.pop()
        if type(node) == Branch:
            for c in node.subtrees:
                levels[c] = l+1
                queue.append((c, l+1))
    return levels


def constructGOTree(pathDict):
    """ construct GO hierarchy and links from paths of genes that go from their
    leaf representations to the root""" 

    jpaths = list(set(chain(*pathDict.values())))
    max_l = max(map(len, jpaths))
    pNodes = [None] * len(jpaths)
    t = Tree(None)
    l = 0
    while l <= max_l:
        l += 1
        for i in xrange(len(jpaths)):
            p = jpaths[i]
            lastp = len(p) 
            # skip path if too short
            if l > lastp:
                continue 

            # we are at the root
            b = None
            if not (l - 1):
                if t.subtree and t.subtree.label != p[-l]:
                    # ooops
                    import pdb; pdb.set_trace() 
                elif not t.subtree: 
                    t.subtree = Branch(list(), label=p[-l])
                b = t.subtree
            # for all other nodes that have yet not been identified
            elif p[-l] not in [x.label for x in pNodes[i].subtrees]:
                b = Branch(list(), label=p[-l])
                pNodes[i].subtrees.append(b)
                b.setAncestor(pNodes[i])
            # for identified nodes
            else:
                b = [x for x in pNodes[i].subtrees if x.label == p[-l]][0]
            pNodes[i] = b


    leafMapping = dict()
    for i in xrange(len(pNodes)):
        leafMapping[pNodes[i].label] = pNodes[i]

    # construct gene representation link list
    links = dict()
    for gene, reps in pathDict.items():
        gl = list()
        for pLeaf in set(map(lambda x: x[0], reps)):
            l = Leaf(gene)
            parent = leafMapping[pLeaf]
            parent.subtrees.append(l)
            l.setAncestor(parent)
            gl.append(l)
        links[gene] = gl
    return t, links


def nearestNeighborDist(t, links, levels, genes):
    """ compute nearest-neighbor distance in GO hiearchy t for each gene in
    genes """

    res = {gene:float('inf') for gene in genes}


    # sort initial genes into disjoint sets, one for each level; create data
    # structure pushMaps that allows to push elements from bottom to top
    pushMaps = dict()
    levelSets = dict()
    for gene in genes:
        if links.has_key(gene):
            for rep in links[gene]:
                l = levels[rep] 
                if not levelSets.has_key(l):
                    levelSets[l] = set()
                levelSets[l].add(rep)
                pushMaps[rep] = {rep.label:l}

    if not levelSets:
        return levelSets
   
    mx_level = max(levelSets.keys())
    levelSets = [levelSets.get(l, set()) for l in xrange(mx_level+1)]

    # iterate through each level and process elements associated with nodes of
    # the tree, that are successively pushed from bottom to top
    for l in xrange(mx_level, -1, -1):
        for node in levelSets[l]:
            elems = pushMaps[node]
            if len(elems) < 2:
                # if only one element is present, it is by definition the
                # push-up represenative
                ((gene, l1), ) = elems.items()
            else:
                for (el1, l1), (el2, l2) in combinations(elems.items(), 2):
                    d = float(l1+l2-2*l-1)/(l1+l2) 
                    res[el1] = min(res[el1], d)
                    res[el2] = min(res[el2], d)
                # identify push-up represenative
                (gene, l1) = min(elems.items(), key=lambda x: x[1])

            # (gene, l1) will be the representative that is pushed to the
            # ancestral node
            p = node.getAncestor()
            if p != None:
                # if p!= None, then, by definition, l > 0
                levelSets[l-1].add(p)
                if not pushMaps.has_key(p):
                    pushMaps[p] = dict()
                if not pushMaps[p].has_key(gene):
                    pushMaps[p][gene] = l1
                else:
                    pushMaps[p][gene] = min(l1, pushMaps[p][gene])
    return res


def printClusterDistances(t, links, levels, nn_genome, cluster_data, out,
        samples=None):
    """ output GO nearest-neighbor distance offset for each gene cluster """

    isHeader = True
    for line in csv.reader(cluster_data, delimiter='\t'):
        if isHeader:
            isHeader = False
            continue

        genes = set(filter(lambda x: links.has_key(x), line[0].split(';')))
        nn_cluster = nearestNeighborDist(t, links, levels, genes)
        
        score = float('inf')
        p = -1
        k = len(genes)
        if k > 1 and samples and samples.has_key(k):
            score = sum(nn_cluster[g]-nn_genome[g] for g in genes)
            p = float(bisect(samples[k], score))/len(samples[k])
        print >> out, '%s\t%s\t%s\t%s' %(line[0], k, score, p)

if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('-s', '--samples', type=str, 
            help='file containing empirically sampled nearest neighbor scores')
    parser.add_argument('-c', '--is_compressed', action='store_true', 
            help='flag to indicate that file containing empirically ' + \
                    'sampled nn scores is compressed (GZIP)')
    parser.add_argument('go_obo_file', type=str, 
            help='gene Ontology hierarchy in OBO format')
    parser.add_argument('go_associations', type=str,
            help='gene Ontology assocations file that links GO terms ' + \
                    'with genes')
    parser.add_argument('genome_annotations', type=str, 
            help='genome annotation file')
    parser.add_argument('cluster_file', type=str,
            help='gene clusters in proprietary format')

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

    samples = None
    if args.samples:
        LOG.info('reading samples from file %s' %args.samples)
        samples = dict()

        data = None
        if args.is_compressed:
            from gzip import open as gopen
            data = gopen(args.samples)
        else:
            data = open(args.samples)
        
        for line in csv.reader(data, delimiter='\t'):
            samples[int(line[0])] = map(float, line[1:])
            
    LOG.info('nearest neighbor analysis over all annotated genes')
    nearest_gene_dists = nearestNeighborDist(t, links, levels, genes)

    printClusterDistances(t, links, levels, nearest_gene_dists,
            open(args.cluster_file), stdout, samples=samples)
