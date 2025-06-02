process ASSIGNMENT {
    container "quay.io/peepk/scqc_py:v1.0.0"
    publishDir "$params.outdir/postprocess/assignment", pattern: "$outdir", mode: 'copy'
    tag "$sample_id"

    input:
    tuple val(sample_id), path(combine_results), path(cell_meta), path(mtx_conversions)
    val(n_singlets)

    script:
    outdir = sample_id
    """
    mkdir -p $outdir
    scqc_assignment.py \
        --sample $sample_id \
        --genotype_tools $params.demuxTools \
        --transcription_tools $params.dblTools \
        --combined_results $combine_results \
        --cell_meta $cell_meta \
        --n_singlets $n_singlets \
        --mtx_conversions $mtx_conversions \
        --outdir $outdir
    """

    output:
    tuple val(sample_id), path(outdir, type: 'dir'), emit: outdir
    path(n_singlets), emit: n_singlets
}
