process VIREO {
    container "$projectDir/singularity_imgs/Demuxafy_v2.1.0.sif"
    publishDir "${params.outdir}/demuxafy/vireo", mode: 'copy'
    tag "$sample_id"
    label "process_low"

    input:
    tuple val(sample_id), path(cellsnp_outdir, name: "cellsnp_outdir"), path(vcf)

    script:
    """
    mkdir -p $sample_id
    
    vireo \
        -c $cellsnp_outdir \
        -d $vcf \
        -t $params.GTField \
        -p $task.cpus \
        -o ./$sample_id \
        --callAmbientRNAs
    """

    output:
    tuple val(sample_id), path(sample_id, type: 'dir')
}
