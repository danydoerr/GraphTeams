CONTENTS OF THIS FILE
------------------------------------------------------------------------------

* Requirements
* Installation
* Tutorial


REQUIREMENTS
------------------------------------------------------------------------------

* PYTHON3 and PYTHON2 (https://www.python.org/)
* PYTHON2 libraries:
    - NetworkX
* SNAKEMAKE (http://snakemake.bitbucket.org)


INSTALLATION
------------------------------------------------------------------------------

0. Satisfy above-described dependencies
1. Clone git repository
    $ git clone http://github.com/danydoerr/GraphTeams


TUTORIAL
------------------------------------------------------------------------------

1. For each organism of your dataset:
    * create folders 
        - data/hic_maps/{organism name} and 
        - data/genomic/{organism name}
    * copy all intra-chromosomal Hi-C maps of the organism into
      data/hic_maps/{organism name} (note that the suffix of the file names must
      be as specified in the config.yaml file)
    * copy Ensemble Data table of the organism into /data/genomic/{organism
      name} (data from all chromsomes of the organism should be contained in a
      single file; note that the suffixe of the file name must be as specificed
      in the config.yaml file)

2. Adjust configuration file config.yaml:
    * adjust the resolution of your Hi-C maps (at the moment, all maps must have
      the same resolution);
    * adjust the file endings of the Hi-C maps and the genome data files that
      contain information about the location of genes and their gene family
      membership; 
    * adjust the choices of the delta threshold (for more information on this
      parameter, read our soon-to-be published paper);

3. Run 'snakemake --cores {number of cores you want to allocate for running GraphTeams}
4. Identified graph teams will be written to the graph_teams sub-directory
