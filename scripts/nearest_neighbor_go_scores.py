#!/usr/bin/env python2

from sys import stdout, stderr, exit
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF
from os.path import basename
from itertools import chain, combinations
import logging
import csv

from newick_parser import Tree, Branch, Leaf

LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)

ROOT_NODE = 'GO:0008150'

def readAssociations(data):
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

    genes = set()
    for line in csv.reader(data, delimiter='\t'):
        genes.add(line[3])

    return genes

def constructLevelMap(goTree):

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
    # return tree and links from each gene to its leaf representations
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
        for rep in reps:
            l = Leaf(gene)
            parent = leafMapping[rep[0]]
            parent.subtrees.append(l)
            l.setAncestor(parent)
            gl.append(l)
        links[gene] = gl
    return t, links


def nearestNeighborDist(t, links, levels, genes):

    res = {gene:float('inf') for gene in genes}

    levelLists = dict()
    pushLists = dict()

    for gene in genes:
        if links.has_key(gene):
            for rep in links[gene]:
                l = levels[rep] 
                if not levelLists.has_key(l):
                    levelLists[l] = set()
                levelLists[l].add(rep)
                pushLists[rep] = {rep.label:l}
   
    mx_level = max(levelLists.keys())
    levelLists = [levelLists.get(l, set()) for l in xrange(mx_level+1)]

    for l in xrange(mx_level, -1, -1):
        for node in levelLists[l]:
            elems = pushLists[node]
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

            p = node.getAncestor()
            if p != None:
                # if p!= None, then, by definition, l > 0
                levelLists[l-1].add(p)
                if not pushLists.has_key(p):
                    pushLists[p] = dict()
                if not pushLists[p].has_key(gene):
                    pushLists[p][gene] = l1
                else:
                    pushLists[p][gene] = min(l1, pushLists[p][gene])
    return res


def printClusterDistances(t, links, levels, nn_genome, cluster_data):
    for line in csv.reader(cluster_data, delimiter='\t'):
        genes = set(filter(lambda x: links.has_key(x), line[0].split(';')))
        nn_cluster = nearestNeighborDist(t, links, levels, genes)
        print '%s\t%s' %(line[0], sum(nn_cluster[g]-nn_genome[g] for g in
            genes))  
        

if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('go_obo_file', type=str, 
            help='Gene Ontology hierarchy in OBO format')
    parser.add_argument('go_associations', type=str,
            help='Gene Ontology assocations file that links GO terms ' + \
                    'with genes')
    parser.add_argument('genome_annotations', type=str, 
            help='genome annotation file')
    parser.add_argument('cluster_file', type=str,
            help='gene clusters in proprietary format')

    args = parser.parse_args()

    # setup logging
    ch = logging.StreamHandler(stderr)
    ch.setLevel(logging.ERROR)
    ch.setFormatter(logging.Formatter('!! %(message)s'))
    cf = logging.FileHandler('%s.nn_go.log' %basename(args.cluster_file).rsplit('.',
        1)[0], mode='w')
    cf.setLevel(logging.DEBUG)
    cf.setFormatter(logging.Formatter('%(levelname)s\t%(asctime)s\t%(message)s'))
    LOG.addHandler(cf)
    LOG.addHandler(ch)

    LOG.info('loading GO OBO file')
    goData = readGO(open(args.go_obo_file))
    LOG.info('loading GO associations')
    assoData = readAssociations(open(args.go_associations))
    LOG.info('laoding gene annotations')
    genes = readGenes(open(args.genome_annotations))

    pathDict = {gene: list(chain(*[map(tuple, getAllRootPaths(rep, goData)) for
        rep in reps])) for gene, reps in assoData.items() if gene in genes}

    LOG.info('found %s out of %s in GO hierarchy' %(len(pathDict), len(genes)))

    LOG.info('construcing GO hierarchy')
    t, links = constructGOTree(pathDict)
    LOG.info('top-down assignment of levels')
    levels = constructLevelMap(t)

    LOG.info('nearest neighbor analysis over all annotated genes')
    nearest_gene_dists = nearestNeighborDist(t, links, levels, genes)

