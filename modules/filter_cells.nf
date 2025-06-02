process FILTER_CELLS {
    container "quay.io/peepk/scqc-seurat:scqc-seurat_linux-amd64"
    publishDir "$params.outdir/preprocess/filter_cells/$sample_id", pattern: "*.png", mode: 'copy'
    publishDir "$params.outdir/preprocess/filter_cells/$sample_id", pattern: "$filtered_mtx_outdir", mode: 'copy'
    publishDir "$params.outdir/preprocess/filter_cells/$sample_id", pattern: "$filter_summary_outpath", mode: 'copy'
    tag "$sample_id"
    label "process_medium"

    input:
    tuple val(sample_id), path(cellranger_outdir, name: "outs_symlink")
    val(n_cells_outpath)

    script:
    filtered_mtx_outdir = "filtered_mtx"
    filter_summary_outpath = "filter_summary.tsv"
    """
    echo Hello from filter_cells.nf
    # Make a copy of the cellranger output directory because we will modify it's contents but
    # replace the large bam file with a symlink to save space.
    mkdir -p outs
    cp -aL $cellranger_outdir/* outs/
    ln -sf \$(readlink -f $cellranger_outdir)/possorted_genome_bam.bam outs/possorted_genome_bam.bam

    mkdir -p $filtered_mtx_outdir

    filter_cells.R $params.MTThreshold $filtered_mtx_outdir $filter_summary_outpath $sample_id $n_cells_outpath
    """

    output:
    // Written to results
    tuple val(sample_id), path(filtered_mtx_outdir, type: 'dir'), emit: filtered_mtx_outdir
    path(filter_summary_outpath)
    path("*.png")

    // Used in subsequent processes
    tuple val(sample_id), path("outs", type: "dir"), emit: my_cellranger_outdir  // Contains the filtered mtx in the structure of the original cellranger outdir
    path(n_cells_outpath), emit: n_cells
}
