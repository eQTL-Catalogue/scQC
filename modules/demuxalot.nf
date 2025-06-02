process DEMUXALOT {
    container "$projectDir/singularity_imgs/Demuxafy_v2.1.0.sif"
    publishDir "${params.outdir}/demuxafy/demuxalot", mode: 'copy'
    tag "$sample_id"

    input:
    tuple val(sample_id), path(vcf), path(vcf_idx), path(barcodes), path(bam), path(bam_idx), path(inds)

    script:
    """
    mkdir -p $sample_id

    Demuxalot.py \
        -b $barcodes \
        -a $bam \
        -n $inds \
        -v $vcf \
        -o ./$sample_id \
        -r True
    """
    
    // Optional command to get the summary of the classified cells
    //    bash demuxalot_summary.sh ./$sample/assignments_refined.tsv.gz > ./$sample/demuxalot_summary.tsv

    output:
    tuple val(sample_id), path(sample_id, type: 'dir')
}
