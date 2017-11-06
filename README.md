CONTENTS OF THIS FILE
------------------------------------------------------------------------------

* Requirements
* Installation
* Tutorial
* Important Note


REQUIREMENTS
------------------------------------------------------------------------------

* PYTHON3 and PYTHON2 (https://www.python.org/)
* PYTHON2 libraries:
    - NetworkX
* PYTHON3 libraries:
    - psutil
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
    * copy all Hi-C maps of the organism into data/hic_maps/{organism name} 
      (note that all file names must contain the resolution of their 
      corresponding Hi-C map in base pairs directly in front of their suffix 
      as specified in the config.yaml file)
    * file names of intrachromosomal Hi-C maps should contain a substring 
      "chrID." where "ID" is the identifier of the chromosome belonging to the
      map which is also used in the corresponding annotation file
    * file names of interchromosomal Hi-C maps should contain a substring
      "chrID.chrID." where "ID" is the identifier of the chromosomes belonging 
      to the map which is also used in the corresponding annotation file
    * copy Ensemble Data table of the organism into /data/genomic/{organism
      name} (data from all chromsomes of the organism should be contained in a
      single file; note that the suffix of the file name must be as specificed
      in the config.yaml file)

2. Adjust configuration file config.yaml:
    * adjust the file endings of the Hi-C maps and the genome data files that
      contain information about the location of genes and their gene family
      membership; 
    * adjust the choices of the delta threshold (for more information on this
      parameter, we refer to https://link.springer.com/chapter/10.1007/978-3-319-67979-2_11);

3. Run 'snakemake --cores {number of cores you want to allocate for running GraphTeams}
4. Identified graph teams will be written to the graph_teams sub-directory

IMPORTANT NOTE
------------------------------------------------------------------------------

Our workflow now allows the usage of Hi-C maps with different resolutions.
Therefore, the Hi-C maps' resolutions are adjusted with respect to the highest
resolution among all maps. This requires all resolutions to be divisible by the
highest resolution without remainder. Make sure that all Hi-C maps to be used
fulfil this criteria.
