# thesis
thesis project: screening of EVEs in spodoptera genomes

This repository contains the snakemake workflows designed for screening assembled genomes and whole genome sequencing data

#CONDA DEPENDENCIES
the following dependencies are obtained via conda during analysis, these do not have to be manually installed beforehand
  base:
      Python2
      Python3
  pipeline 1:
    ete3
    ete_toolchain
    ncbi datasets
    taxonkit 
   
  pipeline 2:
    ngs-bits
    pysam
    mmseqs2
    bwa
    samtools
    
#MANUAL DEPENDENCIES:
  pipeline 1:
    diamond
    Blast_To_Bed.py (EVE_FINDER toolbox)
    Refine_Candidates pipeline
  pipeline 2:
    vyper
    
