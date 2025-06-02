# scQC
Nextflow workflow for scRNAseq QC, demultiplexing, and doublet detection

Prerequisites:
* Nextflow
* Singularity
* The Demuxafy singularity image. Install instructions:
    * `wget -O Demuxafy_v2.1.0.sif https://www.dropbox.com/scl/fi/g0cuyjwomdavom6u6kb2v/Demuxafy.sif?rlkey=xfey1agg371jo4lubsljfavkh`
    * Put the downloaded `.sif` file in the `singularity_imgs` directory.
 * cd to the scQC directory

An example of the nextflow run script is in [`run.sh`](run.sh)

Parameters:
* `--samples` TSV containing the sample IDs, nf-core/scrnaseq Cell Ranger output directories, files of individuals in the pools, cell metadata (example: [`assets/input_examples/samples.tsv`](assets/input_examples/samples.tsv))
* `--vcf` VCF containing the genotypes of the individuals in the samples
* `--GT-field` the genotype field in the VCF (default: GT)
* `--MT-threshold` maximum percentage of mitochondrial RNA allowed in a cell
* `--exonic-regions` file containing the exonic regions of genes (example: [`assets/exonic_regions/hg38exonsUCSC.bed`](assets/exonic_regions/hg38exonsUCSC.bed))
* `--outdir` path where to write the outputs

Outputs:
* `preprocess` contains the MT content and MALAT1 expression distribution plots, and the count matrices after MT content and MALAT1 expression filtering
* `demuxafy` contains the outputs of all tools used from the Demuxafy image
* `postprocess` contains the count matrices with donor assignments in `.h5ad` files. Also, the barplots of filtering summaries and comparison to metadata

Note: The cell metadata is used only for comparing filtering and demultiplexing results to metadata, not for assigning donor labels.
