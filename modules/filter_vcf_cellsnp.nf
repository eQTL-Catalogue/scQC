process FILTER_VCF {
    container 'quay.io/biocontainers/bcftools:1.21--h3a4d415_1'
    publishDir "${params.outdir}/demuxafy/filter_vcf/$sample_id", mode: 'copy'
    tag "$sample_id"
    label "filter"

    input:
    tuple val(sample_id), path(cellsnp_outdir), path(vcf), path(vcf_idx)

    script:
    outpath = "common_exonic_SNPs.cellsnp_vars.vcf"
    """
    bcftools view --threads $task.cpus -R $cellsnp_outdir/cellSNP.base.vcf.gz $vcf -Ov -o $outpath
    """

    output:
    tuple val(sample_id), path(outpath)
}
