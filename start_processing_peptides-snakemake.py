#!/bin/bash
#SBATCH --job-name=netCTL  # Job name
#SBATCH --mail-type=ALL           # notifications for job done & fail
#SBATCH --mail-user=eknodel@asu.edu # send-to address
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 96:00:00
#SBATCH --mem=1024

newgrp combinedlab

source activate var_call_env
module load r/latest

#export
PERL5LIB=/packages/6x/vcftools/0.1.12b/lib/per15/site_perl

snakemake --snakefile processing_peptides-snakemake.py -j 15 --keep-target-files --rerun-incomplete --cluster "sbatch -n 1 -c 8 -t 96:00:00"
