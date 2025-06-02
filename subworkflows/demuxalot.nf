include { DEMUXALOT } from "$projectDir/modules/demuxalot"

workflow DEMUXALOT_WF {
    take:
    vcf_ch
    my_cellranger_outdir_ch
    inds_ch

    main:
    demuxalot_input_ch = vcf_ch
        .join(
            my_cellranger_outdir_ch.map { it -> [
                it[0], 
                "${it[1]}/filtered_feature_bc_matrix/barcodes.tsv.gz", 
                "${it[1]}/possorted_genome_bam.bam", 
                "${it[1]}/possorted_genome_bam.bam.bai"
            ]}
            .join(inds_ch)
        )
    //demuxalot_input_ch.view()

    DEMUXALOT(demuxalot_input_ch)

    emit:
    demuxalot_out = DEMUXALOT.out
}
