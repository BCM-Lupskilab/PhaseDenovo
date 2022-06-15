# PhaseDenovo

This repo includes the R script to phase de novo variant based on long and short-read sequencing data. 


## Input files

1. A bed file includes *de novo* variants:

*  Chrom
*  pos
*  REF
*  ALT

2. vcf file with variant phased using [whatshap(v.1.0)](https://whatshap.readthedocs.io/en/latest/)
run whatshap using short-read and long-read data to get the physical phased variant (0|1, 1|0, etc.)

> [~]$ whatshap phase -o phased.vcf --reference=reference.fasta input.vcf nanopore.bam pacbio.cram illumina.bam 



## Run on Rstudio

[Example](https://cluhaowie.github.io/PhaseDenovo/run_example.html) of running the script in Rstudio.

## Session infor
> 