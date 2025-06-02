process DOUBLETFINDER {
    container "$projectDir/singularity_imgs/Demuxafy_v2.1.0.sif"
    publishDir "${params.outdir}/demuxafy/doubletfinder", mode: 'copy'
    tag "$sample_id"
    //label "process_medium"

    input:
    tuple val(sample_id), path(countmtx), val(scds_n_dbls), val(scdblfinder_n_dbls)

    script:
    n_dbls = [scds_n_dbls, scdblfinder_n_dbls]*.toInteger().average().round()  // The number of expected doublets is the average of the previous methods
    """
    echo "Expected number of doublets for pool $sample_id is $n_dbls"
    mkdir -p $sample_id
    scqc_doubletfinder_wrapper.R $sample_id $countmtx $n_dbls
    """

    output:
    tuple val(sample_id), path(sample_id, type: 'dir'), emit: outdir
}
