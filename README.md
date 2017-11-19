CONTENTS OF THIS FILE
------------------------------------------------------------------------------

* Requirements
* Installation
* Tutorial
* Important Note
* Generation of Homology Table


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
    * create a folder data/hic_maps/{organism name} where {organism name} is
      the same name for a species used inside the homology table of Ensembl
      and in connection with it.
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

2. Copy a homology table gained by Ensembl into /data/genomic/{organism name}
   where {organism name} is the name of the outgoing species the homology 
   table was generated with. (For detailed instructions on how to generate 
   such a homology table using Ensembl's BioMart see GENERATION OF HOMOLOGY 
   TABLE below.)

3. Adjust configuration file config.yaml:
    * adjust the file endings of the Hi-C maps and the genome data files that
      contain information about the location of genes and their gene family
      membership; 
    * adjust the choices of the delta threshold (for more information on this
      parameter, we refer to https://link.springer.com/chapter/10.1007/978-3-319-67979-2_11);

3. Run 'snakemake --cores {number of cores you want to allocate for running GraphTeams}'
4. Identified graph teams will be written to the graph_teams sub-directory

IMPORTANT NOTE
------------------------------------------------------------------------------

Our workflow now allows the usage of Hi-C maps with different resolutions.
Therefore, the Hi-C maps' resolutions are adjusted with respect to the highest
resolution among all maps. This requires all resolutions to be divisible by the
highest resolution without remainder. Make sure that all Hi-C maps to be used
fulfil this criteria.

GENERATION OF HOMOLOGY TABLE
------------------------------------------------------------------------------

Beside from Hi-C data of the species under investigation, the detection of 
gene clusters additionally requires gene annotations and homology assignments.
We assume that both is queried using Ensembl's BioMart web interface tool 
(available under http://ensembl.org/biomart/martview). Detailed instructions
on how to generate these data is given in the following:

1. Go to http://ensembl.org/biomart/martview

2. Choose "Ensembl Genes 90".

3. Choose one of the species you want to gain information about.

4. Click on “Filters” in the left menu.

5. Unfold the “MULTI SPECIES COMPARISONS” box, tick the “Homologue filters”
   option and choose another species you want to gain information about from 
   the drop-down menu.

6. Click on “Attributes” in the left menu.

7. Click on “Homologues”.

8. Unfold the "GENE" box, make a tick at "Chromosome/scaffold name", 
   "Gene start (bp)", "Gene end (bp)" and "Gene name", and remove the ticks 
   for "Gene stable ID" and "Transcript stable ID".

9. Unfold the "ORTHOLOGUES" box and tick for each species (except for the spe-
   cies chosen in 3.)
    * "{Species name} gene name",
    * "{Species name} chromosome/scaffold name",
    * "{Species name} chromosome/scaffold start (bp)" and
    * "{Species name} chromosome/scaffold end (bp)".

10. Click on the “Results” button (top left).

11. Choose a CSV file as output and click on the "Go" button.

Note: In the above describtion we used the currently available Ensembl data-
      base version "Genes 90". However, available versions may differ in the 
      future.
