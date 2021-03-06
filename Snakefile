configfile: 'config.yaml'

from sys import stdout, stderr, exit
from os.path import basename, join, isdir
from itertools import chain
from glob import glob
import re

BIN_DIR = config['bin_dir']
HIC_DATA_DIR = config['hic_data_dir']

ORGANISMS = config.get('compare_only', sorted(basename(x) for x in glob(
        '%s/*' %HIC_DATA_DIR) if isdir(x)))
ORG_SHORT = '_'.join(map(lambda x: x[:3].lower(), ORGANISMS))

ORG_DIFF = [x for x in config.get('diff_pairs', ()) if x[0] in ORGANISMS and
        x[1] in ORGANISMS]

HIC_MAPS_BASE = sorted(map(lambda y: join(*(y.split('/')[
        len(HIC_DATA_DIR.split('/')):])), chain(*(glob(
        '%s/%s/*%s' %(HIC_DATA_DIR, x, config['hic_file_ending'])) for x in
        ORGANISMS))))
HIC_FORMAT = config['hic_map_format']
HIC_ROW_OFFSET = 0
HIC_COL_OFFSET = 0

if HIC_FORMAT == 'DIXON':
    HIC_ROW_OFFSET = config.get('hic_row_offset', 1)
    HIC_COL_OFFSET = config.get('hic_col_offset', 1)
elif HIC_FORMAT == 'HOMER':
    HIC_ROW_OFFSET = config.get('hic_row_offset', 1)
    HIC_COL_OFFSET = config.get('hic_col_offset', 2)

GRAPH_DATA_DIR = config['graph_data_dir']
SEQ_GRAPH_DATA_DIR = config['graph_data_dir'] + '_seq'
DIFF_GRAPH_DATA_DIR = config['graph_data_dir'] + '_diff'

GENE_DATA_DIR = config['gene_data_dir']
HOMOLOGY_MAPS = sorted(glob('%s/%s/*%s' %(GENE_DATA_DIR, x,
        config['gene_data_file_ending']))[0] for x in ORGANISMS)

TEAMS_DIR = config['graph_teams_dir']
SEQ_TEAMS_DIR = config['graph_teams_dir'] + '_seq'
DIFF_TEAMS_DIR = config['graph_teams_dir'] + '_diff'
DELTA = config['delta']
MX_DELTA = config.get('max_delta', max(DELTA))

GO_ANALYSIS_DIR = config['go_analysis_dir']
GO_ASSOC_DATA = next(iter(glob('%s/*.gaf' %config['go_data_dir'])), [])
GO_OBO_FILE = next(iter(glob('%s/*.obo' %config['go_data_dir'])), [])
GO_REF = config['go_reference_species']
GO_SAMPLE_SIZE = config['go_sample_pool_size']

OCCURRENCE_MIN = config.get('occurrence_min', 1)
REAL_GCS = config.get('biological_geneclusters_file', 'biological_geneclusters.txt')

if MX_DELTA < max(DELTA):
    print(('\t!! ERROR: value for \'max_delta\' cannot be set lower than ' + \
            'maximum delta in list \'delta\''))
    exit(1)

if GO_REF and GO_REF not in ORGANISMS:
    print(('\t!! ERROR: GO reference genome %s not found in set ' + \
            'of organisms {%s}. Exiting') %(GO_REF, ', '.join(ORGANISMS)))
    exit(1)

LOG_DIR = config['log_dir']

if not isdir(LOG_DIR):
    from os import makedirs
    makedirs(LOG_DIR)

rule all: 
    input:
        expand('%s/%s_d{delta}.csv' %(TEAMS_DIR, ORG_SHORT), delta=DELTA)


rule all_sequential: 
    input:
        expand('%s/%s_d{delta}.csv' %(SEQ_TEAMS_DIR, ORG_SHORT), delta=DELTA)


rule go_analysis:
    input: 
        expand('%s/{teams_dir}/%s_ref_%s_d{delta}.significant' %(
                GO_ANALYSIS_DIR, ORG_SHORT, GO_REF), teams_dir=(TEAMS_DIR,
                SEQ_TEAMS_DIR), delta=DELTA), 
        '%s_ref_%s_go_score_stats.csv' %(ORG_SHORT, GO_REF)


rule all_diff: 
    input:
        expand('%s/%s_d{delta}.csv' %(DIFF_TEAMS_DIR, '_'.join(map(lambda x:
                '-'.join(x), ORG_DIFF))), delta=DELTA)


rule check_occurrence:
    input:
        'occurrences_%s_m%s.csv' %(ORG_SHORT, OCCURRENCE_MIN)
       

rule check_occurrence_seq:
    input:
        'occurrences_seq_%s_m%s.csv' %(ORG_SHORT, OCCURRENCE_MIN)
        
rule check_occurrence_diff:
    input:
        'occurrences_diff_%s_m%s.csv' %('_'.join(map(lambda x: '-'.join(x),
                ORG_DIFF)), OCCURRENCE_MIN)

rule normalize:
    input:
        '%s/{hic_map}' %HIC_DATA_DIR
    params:
        row_offset = HIC_ROW_OFFSET,
        col_offset = HIC_COL_OFFSET,
        ntype = config.get('normalization_type', 'MAX')
    output:
        '%s/{hic_map}.dmat' %HIC_DATA_DIR 
    log:
        '%s/normalizer_{hic_map}.log' %LOG_DIR
    shell:
        '%s/normalizer.py -t {params.ntype} -x {params.col_offset} ' %BIN_DIR +
        '-y {params.row_offset} {input} 2> {log}'


rule process_annotations:
    input:
        '%s/{organism}/{gene_data}%s' %(GENE_DATA_DIR, config['gene_data_file_ending'])
    params:
        ann_format = config['gene_data_format']
    output:
        '%s/{organism}/{gene_data}.annotation' %GENE_DATA_DIR
    shell:
        '%s/process_annotations.py -f {params.ann_format} ' %BIN_DIR + 
        '{input} > {output}'


if config['gene_data_format'] == 'ENSEMBLE':
    rule makeHomologyTable:
        input:
            HOMOLOGY_MAPS
        output:
            table = '%s/homology_%s.csv' %(GENE_DATA_DIR, ORG_SHORT),
            tmp = temp('%s/homology_%s.csv.tmp' %(GENE_DATA_DIR, ORG_SHORT))
        shell:
            '%s/MakeHomList.py {input} > {output.tmp};' %BIN_DIR + 
            'sort {output.tmp} | uniq > {output.table};'


rule buildGraphs:
    input:
        hic_dmat = lambda wildcards: expand('%s/{hic_map}.dmat' %HIC_DATA_DIR,
                hic_map=(x for x in HIC_MAPS_BASE if x.split('/', 1)[0] ==
                wildcards.organism)),
        annotation_file = lambda wildcards: ['%s.annotation' %x.rsplit('.',
                1)[0] for x in HOMOLOGY_MAPS if x.startswith('%s/%s' %(GENE_DATA_DIR, 
                wildcards.organism))],
        homology_table = '%s/homology_%s.csv' %(GENE_DATA_DIR, ORG_SHORT)
    params:
        mx_delta = MX_DELTA,
        hic_format = HIC_FORMAT,
    output:
        '%s/{organism}_d%s.ml' %(GRAPH_DATA_DIR, MX_DELTA)
    log:
        '%s/build_graphs_{organism}_d%s.log' %(LOG_DIR, MX_DELTA)
    shell:
        'mkdir -p %s;' %GRAPH_DATA_DIR + 
        '%s/hic_to_graph.py -d {params.mx_delta} -f {params.hic_format}' %BIN_DIR +
        ' {input.annotation_file} {input.homology_table} {input.hic_dmat} > '
        '{output} 2> {log}'


rule buildSequentialGraphs:
    input:
        hic_dmat = lambda wildcards: expand('%s/{hic_map}.dmat' %HIC_DATA_DIR,
                hic_map=(x for x in HIC_MAPS_BASE if x.split('/', 1)[0] ==
                wildcards.organism)),
        annotation_file = lambda wildcards: ['%s.annotation' %x.rsplit('.',
                1)[0] for x in HOMOLOGY_MAPS if x.rsplit('/', 2)[-2] ==
                wildcards.organism],
        homology_table = '%s/homology_%s.csv' %(GENE_DATA_DIR, ORG_SHORT)
    params:
        mx_delta = MX_DELTA,
        hic_format = HIC_FORMAT,
    output:
        '%s/{organism}_d%s.ml' %(SEQ_GRAPH_DATA_DIR, MX_DELTA)
    log:
        '%s/build_graphs_seq_{organism}_d%s.log' %(LOG_DIR, MX_DELTA)
    shell:
        'mkdir -p %s;' %SEQ_GRAPH_DATA_DIR + 
        '%s/hic_to_graph.py -d {params.mx_delta} -f {params.hic_format}' %BIN_DIR +
        ' --sequential {input.annotation_file} {input.homology_table} '
        '{input.hic_dmat} > {output} 2> {log}'


rule buildDiffGraphs:
    input:
        o1 = '%s/{o1}_d%s.ml' %(GRAPH_DATA_DIR, MX_DELTA),
        o2 = '%s/{o2}_d%s.ml' %(GRAPH_DATA_DIR, MX_DELTA)
    params:
        mx_delta = MX_DELTA
    output:
        '%s/{o1}-{o2}_d%s.ml' %(DIFF_GRAPH_DATA_DIR, MX_DELTA)
    log:
        '%s/build_diff_graphs_{o1}-{o2}_d%s.log' %(LOG_DIR, MX_DELTA)
    shell:
        'mkdir -p %s;' %DIFF_GRAPH_DATA_DIR + 
        '%s/construct_diff_graph.py -d {params.mx_delta} {input.o1} ' %BIN_DIR + 
        '{input.o2} > {output} 2> {log}'
       

rule findTeams:
    input:
        graph=expand('%s/{organism}_d%s.ml' %(GRAPH_DATA_DIR, MX_DELTA),
                organism=ORGANISMS)
    output:
        '%s/%s_d{delta,[0-9.-]+}.csv' %(TEAMS_DIR, ORG_SHORT)
    benchmark:
        'benchmarks/teams_spatial_%s_d{delta}.txt' %ORG_SHORT
    shell:
        'mkdir -p %s;' %TEAMS_DIR + 
        '%s/graph_teams.py -d {wildcards.delta} {input} > {output}' %BIN_DIR


rule findStringTeams:
    input:
        graph=expand('%s/{organism}_d%s.ml' %(SEQ_GRAPH_DATA_DIR, MX_DELTA),
                organism=ORGANISMS)
    output:
        '%s/%s_d{delta,[0-9.-]+}.csv' %(SEQ_TEAMS_DIR, ORG_SHORT)
    benchmark:
        'benchmarks/teams_seq_%s_d{delta}.txt' %ORG_SHORT
    shell:
        'mkdir -p %s;' %SEQ_TEAMS_DIR + 
        '%s/graph_teams.py -d {wildcards.delta} {input} > {output}' %BIN_DIR


rule findDiffTeams:
    input:
        graphs=expand('%s/{odiff}_d%s.ml' %(DIFF_GRAPH_DATA_DIR, MX_DELTA),
                odiff = map(lambda x: '-'.join(x), ORG_DIFF))
    output:
        '%s/%s_d{delta,[0-9.-]+}.csv' %(DIFF_TEAMS_DIR, '_'.join(map(lambda x:
                '-'.join(x), ORG_DIFF)))
    benchmark:
        'benchmarks/teams_diff_%s_d{delta}.txt' %'_'.join(map(lambda x:
                '-'.join(x), ORG_DIFF))
    shell:
        'mkdir -p %s;' %DIFF_TEAMS_DIR + 
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
        '   %s/ClusterEvaluator.py -d$i %s/%s_d$i.csv %s/%s_d$i.csv benchmarks/teams_spatial_%s_d$i.txt;' %(BIN_DIR, TEAMS_DIR, ORG_SHORT, SEQ_TEAMS_DIR, ORG_SHORT, ORG_SHORT) + 
        'done > {output}'


rule visualize_cluster_stats:
    input:
        '%s_stats.csv' %ORG_SHORT
    output:
        '%s_stats.pdf' %ORG_SHORT
    shell:
        '%s/visualize_cluster_stats.py {input} > {output}' %BIN_DIR

rule sample_go_neighbor_cluster_scores:
    input:
        obo = GO_OBO_FILE,
        assoc = GO_ASSOC_DATA,
        annot = ['%s.annotation' %x.rsplit('.', 1)[0] for x in HOMOLOGY_MAPS if
                x.startswith('%s/%s' %(GENE_DATA_DIR, GO_REF))],
    params:
        pool_size = GO_SAMPLE_SIZE
    output:
        temp('%s/%s_samples_s{cluster_size}_n%s.csv'%(GO_ANALYSIS_DIR, GO_REF,
                GO_SAMPLE_SIZE))
    log:
        '%s/sample_nn_go_ref_%s_n%s.log' %(LOG_DIR, GO_REF, GO_SAMPLE_SIZE)
    shell:
        'mkdir -p %s;' %GO_ANALYSIS_DIR +
        '%s/sample_nn_go_scores.py {input.obo} {input.assoc} ' %BIN_DIR +
        '{input.annot} {wildcards.cluster_size} {params.pool_size} > {output} '
        '2> {log}'

rule all_sample_go_neighbor_cluster_scores:
    input:
        expand('%s/%s_samples_s{cluster_size}_n%s.csv' %(GO_ANALYSIS_DIR,
                GO_REF, GO_SAMPLE_SIZE), cluster_size=range(2,
                config['go_sample_max_cluster_size']+1))
    output:
        '%s/%s_samples_n%s.csv.gz' %(GO_ANALYSIS_DIR, GO_REF, GO_SAMPLE_SIZE)
    shell:
        'OUT="{output}";'
        'for i in `seq 2 %s`; do' %config['go_sample_max_cluster_size'] +
        '   echo -en "$i\\t" >> ${{OUT%%.gz}};' 
        '   cat %s/%s_samples_s${{i}}_n%s.csv >> ${{OUT%%.gz}};' %(GO_ANALYSIS_DIR, GO_REF, GO_SAMPLE_SIZE) + 
        'done; gzip ${{OUT%.gz}};'


rule go_neighbor_cluster_scores:
    input:
        obo = GO_OBO_FILE,
        assoc = GO_ASSOC_DATA,
        annot = ['%s.annotation' %x.rsplit('.', 1)[0] for x in HOMOLOGY_MAPS if
                x.startswith('%s/%s' %(GENE_DATA_DIR, GO_REF))],
        samples = '%s/%s_samples_n%s.csv.gz' %(GO_ANALYSIS_DIR, GO_REF,
                GO_SAMPLE_SIZE),
        teams = '{teams_dir}/%s_d{delta}.csv' %ORG_SHORT
    output:
        '%s/{teams_dir}/%s_ref_%s_d{delta,[0-9.-]+}.csv' %(GO_ANALYSIS_DIR,
                ORG_SHORT, GO_REF)
    log:
        '%s/nearest_neighbor_go_scores_%s_d{delta}_n%s.log' %(LOG_DIR,
                ORG_SHORT, GO_SAMPLE_SIZE)
    shell:
        'mkdir -p %s/{wildcards.teams_dir};' %GO_ANALYSIS_DIR +
        '%s/nearest_neighbor_go_scores.py -c -s {input.samples} ' %BIN_DIR +
        '{input.obo} {input.assoc} {input.annot} {input.teams} > {output} '
        '2> {log}'


rule go_scores_significant:
    input:
        '%s/{teams_dir}/%s_ref_%s_d{delta}.csv' %(GO_ANALYSIS_DIR, ORG_SHORT,
                GO_REF)
    output:
        '%s/{teams_dir}/%s_ref_%s_d{delta,[0-9.-]+}.significant' %(
                GO_ANALYSIS_DIR, ORG_SHORT, GO_REF)
    shell:
        'awk -F "\\t" \'{{if ($4 >= 0 && $4 < 0.05) print $0}}\' {input} | '
        'sort -gk4 > {output}'
        

rule go_score_stats:
    input:
        expand('%s/{teams_dir}/%s_ref_%s_d{delta}.csv' %(GO_ANALYSIS_DIR,
                ORG_SHORT, GO_REF), teams_dir=(TEAMS_DIR, SEQ_TEAMS_DIR),
                delta=DELTA)
    output:
        '%s_ref_%s_go_score_stats.csv' %(ORG_SHORT, GO_REF)
    shell:
        'for d in %s; do' %' '.join(map(str, DELTA)) + 
        '   n=`wc -l %s/%s/%s_ref_%s_d$d.csv| cut -f1 -d\ `; '%(GO_ANALYSIS_DIR, TEAMS_DIR, ORG_SHORT, GO_REF) +
        '   m=`wc -l %s/%s/%s_ref_%s_d$d.csv| cut -f1 -d\ `; '%(GO_ANALYSIS_DIR, SEQ_TEAMS_DIR, ORG_SHORT, GO_REF) +
        '   i=`awk -F "\\t" "{{if (\\\$2 > 0) sum += \\\$3/\\\$2}} END {{print sum/$n}}" %s/%s/%s_ref_%s_d$d.csv`;' %(GO_ANALYSIS_DIR, TEAMS_DIR, ORG_SHORT, GO_REF) +
        '   j=`awk -F "\\t" "{{if (\\\$2 > 0) sum += \\\$3/\\\$2}} END {{print sum/$m}}" %s/%s/%s_ref_%s_d$d.csv`;' %(GO_ANALYSIS_DIR, SEQ_TEAMS_DIR, ORG_SHORT, GO_REF) +
        '   echo -e "$d\t$i\t$j";'
        'done > {output}'

rule check_occurrence_rule:
    input:
        realgc = REAL_GCS,
        teams = expand('%s/%s_d{delta}.csv' %(TEAMS_DIR, ORG_SHORT), delta=DELTA)
    params:
        min_occ = OCCURRENCE_MIN
    output:
        'occurrences_%s_m%s.csv' %(ORG_SHORT, OCCURRENCE_MIN)
    shell:
        '%s/check_occurrence.py -m {params.min_occ} {input.realgc} ' %BIN_DIR +
        '{input.teams} > {output}'

rule check_occurrence_seq_rule:
    input:
        realgc = REAL_GCS,
        teams = expand('%s/%s_d{delta}.csv' %(SEQ_TEAMS_DIR, ORG_SHORT),
                delta=DELTA)
    params:
        min_occ = OCCURRENCE_MIN
    output:
        'occurrences_seq_%s_m%s.csv' %(ORG_SHORT, OCCURRENCE_MIN)
    shell:
        '%s/check_occurrence.py -m {params.min_occ} {input.realgc} ' %BIN_DIR +
        '{input.teams} > {output}'

rule check_occurrence_diff_rule:
    input:
        realgc = REAL_GCS,
        teams = expand('%s/%s_d{delta}.csv' %(DIFF_TEAMS_DIR,
                '_'.join(map(lambda x: '-'.join(x), ORG_DIFF))), delta=DELTA)
    params:
        min_occ = OCCURRENCE_MIN
    output:
        'occurrences_diff_%s_m%s.csv' %('_'.join(map(lambda x: '-'.join(x),
                ORG_DIFF)), OCCURRENCE_MIN)
    shell:
        '%s/check_occurrence.py -m {params.min_occ} {input.realgc} ' %BIN_DIR +
        '{input.teams} > {output}'

