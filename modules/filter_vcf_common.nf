process FILTER_VCF {
    container "quay.io/biocontainers/bcftools:1.21--h3a4d415_1"
    publishDir "${params.outdir}/preprocess/filter_vcf/$sample_id", mode: 'copy'
    tag "$sample_id"
    label "filter"

    input:
    tuple val(sample_id), path(pool_inds), path(vcf), path(idx)
    each path(regions)

    script:
    outpath = "common_exonic_SNPs.vcf.gz"
    """
    bcftools view --threads $task.cpus --write-index=csi -i 'TYPE="snp" && MAF>0.05' -R $regions -S $pool_inds $vcf -Oz -o $outpath
    """

    output:
    tuple val(sample_id), path(outpath), path("${outpath}.csi"), emit: filtered_vcf
}
