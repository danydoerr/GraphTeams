configfile: 'config.yaml'

from os.path import basename, join, isdir
from itertools import chain
from glob import glob

BIN_DIR = config['bin_dir']

HIC_DATA_DIR = config['hic_data_dir']
ORGANISMS = sorted(basename(x) for x in glob('%s/*' %HIC_DATA_DIR) if isdir(x))
ORG_SHORT = '_'.join(map(lambda x: x[:3].lower(), ORGANISMS))
HIC_MAPS_BASE = sorted(map(lambda y: join(*(y.split('/')[
        len(HIC_DATA_DIR.split('/')):])), chain(*(glob(
        '%s/%s/*%s' %(HIC_DATA_DIR, x, config['hic_file_ending'])) for x in
        ORGANISMS))))
NORM_HIC_DIR = config['normalized_data_dir']
GRAPH_DATA_DIR = config['graph_data_dir']
SEQ_GRAPH_DATA_DIR = config['graph_data_dir'] + '_seq'

GENE_DATA_DIR = config['gene_data_dir']
HOMOLOGY_MAPS = sorted(glob('%s/%s/*%s' %(GENE_DATA_DIR, x,
        config['gene_data_file_ending']))[0] for x in ORGANISMS)

TEAMS_DIR = config['graph_teams_dir']
SEQ_TEAMS_DIR = config['graph_teams_dir'] + '_seq'
DELTA = config['delta']

rule all: 
    input:
        expand('%s/%s_d{delta}.csv' %(TEAMS_DIR, ORG_SHORT), delta=DELTA)

rule all_sequential: 
    input:
        expand('%s/%s_d{delta}.csv' %(SEQ_TEAMS_DIR, ORG_SHORT), delta=DELTA)

rule cutMaps:
    input:
        '%s/{organism}/{hic_map}' %HIC_DATA_DIR
    output:
        temp('%s/{organism}/{hic_map}.truncated' %NORM_HIC_DIR)
    shell:
        'mkdir -p $(dirname {output}); tail -n +3 {input} | cut -f2- > {output}'


rule normalizeMaps:
    input:
        expand('%s/{hic_map}.truncated' %NORM_HIC_DIR, hic_map=HIC_MAPS_BASE)
    output:
        expand('%s/{hic_map}.truncated.dmat' %NORM_HIC_DIR, hic_map=HIC_MAPS_BASE)
    shell:
        '%s/Normalizer.py {input}' %BIN_DIR


rule skipHeader:
    input:
        '%s/{organism}/{gene_data}%s' %(GENE_DATA_DIR,
                config['gene_data_file_ending'])
    output:
        temp('%s/{organism}/{gene_data}.noheader' %GENE_DATA_DIR)
    shell:
        'tail -n +2 {input} > {output}'


rule makeAnnotationFile:
    input:
        '%s/{organism}/{gene_data}.noheader' %GENE_DATA_DIR
    output:
        ann = '%s/{organism}/{gene_data}.annotation' %GENE_DATA_DIR,
	tmp = temp('%s/{organism}/{gene_data}.annotation.tmp' %GENE_DATA_DIR)
    shell:
        '%s/AnFileMaker.py {input} {output.tmp};' %BIN_DIR +
	'sort {output.tmp} | uniq > {output.ann}'


rule makeHomologyTable:
    input:
        map(lambda x: '%s.noheader' %x.rsplit('.', 1)[0], HOMOLOGY_MAPS)
    output:
        table = '%s/homology_%s.csv' %(GENE_DATA_DIR, ORG_SHORT),
        tmp = temp('%s/homology_%s.csv.tmp' %(GENE_DATA_DIR, ORG_SHORT))
    shell:
        '%s/MakeHomList.py {input} {output.tmp};' %BIN_DIR + 
        'sort {output.tmp} | uniq > {output.table};'


rule buildGraphs:
    input:
        hic_dmat = lambda wildcards: expand('%s/{hic_map}.truncated.dmat' %(NORM_HIC_DIR), hic_map=(x for x in HIC_MAPS_BASE if x.split('/', 1)[0] == wildcards.organism)),
        annotation_file = lambda wildcards: ['%s.annotation' %x.rsplit('.', 1)[0] for x in
                HOMOLOGY_MAPS if x.rsplit('/', 2)[-2] == wildcards.organism],
        homology_table = '%s/homology_%s.csv' %(GENE_DATA_DIR, ORG_SHORT)
    params:
        bin_size = config['hic_map_resolution']
    output:
        '%s/{organism}.ml' %(GRAPH_DATA_DIR)
    shell:
        'mkdir -p %s;' %GRAPH_DATA_DIR + 
        '%s/ParseToGraphml.py -s {params.bin_size} ' %BIN_DIR +
        '{input.annotation_file} {input.homology_table} {output} '
        '{input.hic_dmat}'


rule findTeams:
    input:
        graph=expand('%s/{organism}.ml' %GRAPH_DATA_DIR, organism=ORGANISMS)
    output:
        '%s/%s_d{delta}.csv' %(TEAMS_DIR, ORG_SHORT)
    benchmark:
        'benchmarks/teams_spatial_d{delta}.txt'
    shell:
        'mkdir -p %s;' %TEAMS_DIR + 
        '%s/graph_teams.py -d {wildcards.delta} {input} > {output}' %BIN_DIR


rule buildSequentialGraphs:
    input:
        hic_dmat = lambda wildcards: expand('%s/{hic_map}.truncated.dmat' %(NORM_HIC_DIR), hic_map=(x for x in HIC_MAPS_BASE if x.split('/', 1)[0] == wildcards.organism)),
        annotation_file = lambda wildcards: ['%s.annotation' %x.rsplit('.', 1)[0] for x in
                HOMOLOGY_MAPS if x.rsplit('/', 2)[-2] == wildcards.organism],
        homology_table = '%s/homology_%s.csv' %(GENE_DATA_DIR, ORG_SHORT)
    params:
        bin_size = config['hic_map_resolution']
    output:
        '%s/{organism}.ml' %(SEQ_GRAPH_DATA_DIR)
    shell:
        'mkdir -p %s;' %SEQ_GRAPH_DATA_DIR + 
        '%s/SeqGraphMaker.py {params.bin_size} ' %BIN_DIR +
        '{input.annotation_file} {input.homology_table} {output} '
        '{input.hic_dmat}'


rule findStringTeams:
    input:
        graph=expand('%s/{organism}.ml' %SEQ_GRAPH_DATA_DIR, organism=ORGANISMS)
    output:
        '%s/%s_d{delta}.csv' %(SEQ_TEAMS_DIR, ORG_SHORT)
    benchmark:
        'benchmarks/teams_seq_d{delta}.txt'
    shell:
        'mkdir -p %s;' %SEQ_TEAMS_DIR + 
        '%s/graph_teams.py -d {wildcards.delta} {input} > {output}' %BIN_DIR

rule make_stats_file:
    input:
        spatial_teams = expand('%s/%s_d{delta}.csv' %(TEAMS_DIR, ORG_SHORT),
        delta=DELTA),
        seq_teams = expand('%s/%s_d{delta}.csv' %(SEQ_TEAMS_DIR, ORG_SHORT),
        delta=DELTA),
        bench = expand('benchmarks/teams_spatial_d{delta}.txt', delta=DELTA)
    output:
        '%s_stats.csv' %ORG_SHORT
    shell:
        'for i in %s; do ' %(' '.join(map(str, DELTA))) + 
        '   %s/ClusterEvaluator.py -d$i %s/%s_d$i.csv %s/%s_d$i.csv benchmarks/teams_spatial_d$i.txt;' %(BIN_DIR, TEAMS_DIR, ORG_SHORT, SEQ_TEAMS_DIR, ORG_SHORT) + 
        'done > {output}'

rule visualize_stats:
    input:
        '%s_stats.csv' %ORG_SHORT
    output:
        '%s_stats.pdf' %ORG_SHORT
    shell:
        '%s/visualize_cluster_stats.py {input} > {output}' %BIN_DIR
