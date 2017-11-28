#!/usr/bin/env python2

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF
from sys import stderr, stdout, exit 
import logging

import networkx as nx

LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)


if __name__ == '__main__':

    # setting up the argument parser
    parser = ArgumentParser(
            description='Subtracts the second graph from the first',
            formatter_class=ADHF)
    parser.add_argument('-d', '--delta', type=float, 
            help='distance value for absent edges')
    parser.add_argument('graph_1', type=str, help='First graph')
    parser.add_argument('graph_2', type=str, help='First graph')

    args = parser.parse_args()

    # setup logging
    ch = logging.StreamHandler(stderr)
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(logging.Formatter('%(asctime)s %(levelname)s: %(message)s'))
    LOG.addHandler(ch)

    LOG.info('reading first graph')

    G1 = nx.read_gml(args.graph_1)
    if G1.is_directed():
        G1 = G1.to_undirected()

    LOG.info('reading second graph')

    G2 = nx.read_gml(args.graph_2)
    if G2.is_directed():
        G2 = G2.to_undirected()

    import pdb; pdb.set_trace() 

    LOG.info('producing union set of edges')
    
    edges = set(G1.edges())
    edges.update(G2.edges())

    
    LOG.info('constructing difference graph')

    G = nx.Graph()

    for (u, v) in edges:
        w1= args.delta
        try:
            w1 = G1[u][v]['weight']
        except:
            pass
        w2= args.delta
        try:
            w2 = G2[u][v]['weight']
        except:
            pass
        G.add_edge(u, v, weight=w2-w1)

    LOG.info('output graph')
    G.write_gml(G, stdout) 
    LOG.info('DONE')


