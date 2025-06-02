process COMBINE_RESULTS {
    container "$projectDir/singularity_imgs/Demuxafy_v2.1.0.sif"
    publishDir "$params.outdir/demuxafy/combine_results", pattern: "$sample_id", mode: 'copy'
    tag "$sample_id"

    input:
    tuple val(sample_id), path(vireo_outdir, stageAs: 'vireo'), path(demuxalot_outdir, stageAs: 'demuxalot'), path(scds_outdir, stageAs: 'scds'), path(scdblfinder_outdir, stageAs: 'scdblfinder'), path(doubletfinder_outdir, stageAs: 'doubletfinder')

    script:
    tsv_out = "$sample_id/combined_results.tsv"
    """
    mkdir -p $sample_id

    Combine_Results.R \
        -o $tsv_out \
        --DoubletFinder doubletfinder \
        --scds scds \
        --scDblFinder scdblfinder \
        --vireo vireo \
        --demuxalot demuxalot \
        --method "MajoritySinglet"
    """

    output:
    tuple val(sample_id), path(sample_id, type:"dir")
    tuple val(sample_id), path(tsv_out), emit: combined_results
    /**/
}
