#!/bin/bash

#SBATCH --job-name=scQC
#SBATCH --partition=amd
#SBATCH --time=5-00:00:00
#SBATCH --mem=5G

nextflow -log logs/.nextflow.log run main.nf -profile tartu_hpc \
    --samples assets/input_examples/samples.tsv \
    --vcf /path/to/genotypes.vcf.gz \
    --MT-threshold 8 \
    --exonic-regions assets/exonic_regions/hg38exonsUCSC.bed \
    --outdir results \
    -resume
