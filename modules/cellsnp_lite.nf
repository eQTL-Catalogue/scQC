process CELLSNP_LITE {
    container "$projectDir/singularity_imgs/Demuxafy_v2.1.0.sif"
    publishDir "${params.outdir}/demuxafy/cellsnp", mode: 'copy'
    tag "$sample_id"

    input:
    tuple val(sample_id), path(barcodes), path(bam), path(bam_idx), path(vcf), path(vcf_idx)

    script:
    """
    mkdir -p $sample_id

    cellsnp-lite \
        -s $bam \
        -b $barcodes \
        -O ./$sample_id \
        -R $vcf \
        -p $task.cpus \
        --minMAF 0.1 \
        --minCOUNT 20 \
        --gzip
    """

    output:
    tuple val(sample_id), path(sample_id, type: 'dir')
}
