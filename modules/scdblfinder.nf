process SCDBLFINDER {
    container "$projectDir/singularity_imgs/Demuxafy_v2.1.0.sif"
    publishDir "${params.outdir}/demuxafy/scdblfinder", mode: 'copy'
    tag "$sample_id"
    label "process_medium"

    input:
    tuple val(sample_id), path(countmtx_dir)

    script:
    """
    mkdir -p $sample_id

    scDblFinder.R -o $sample_id -t $countmtx_dir

    n_dbls=\$(awk '\$1=="doublet" {print \$2}' $sample_id/scDblFinder_doublet_summary.tsv)
    """

    output:
    tuple val(sample_id), path(sample_id, type: 'dir'), emit: outdir
    tuple val(sample_id), env(n_dbls), emit: n_dbls
}
