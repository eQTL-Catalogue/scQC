include { VIREO } from "$projectDir/modules/vireo"
include { CELLSNP_LITE } from "$projectDir/modules/cellsnp_lite"
include { FILTER_VCF } from "$projectDir/modules/filter_vcf_cellsnp"

workflow VIREO_WF {
    take:
    vcf_ch
    my_cellranger_outdir_ch

    main:
    //samples_ch.view()
    //vcf_ch.view()

    cellsnp_input_ch = my_cellranger_outdir_ch
        .map{ it -> [
            it[0], 
            "${it[1]}/filtered_feature_bc_matrix/barcodes.tsv.gz", 
            "${it[1]}/possorted_genome_bam.bam", 
            "${it[1]}/possorted_genome_bam.bam.bai"
        ]}
        .join(vcf_ch)
    //cellsnp_input_ch.view()

    // CellSNP
    CELLSNP_LITE(cellsnp_input_ch)
    cellsnp_out_ch = CELLSNP_LITE.out
    //cellsnp_out_ch.view()

    // Filter VCF further to keep only variants that are in the CellSNP output
    filter_vcf_ch = cellsnp_out_ch.join(vcf_ch)
    //filter_vcf_ch.view()

    FILTER_VCF(filter_vcf_ch)
    vcf_cellsnp_ch = FILTER_VCF.out
    //vcf_cellsnp_ch.view()

    vireo_input_ch = cellsnp_out_ch.join(vcf_cellsnp_ch)
    //vireo_input_ch.view()

    VIREO(vireo_input_ch)

    emit:
    vireo_out = VIREO.out
    /**/
}