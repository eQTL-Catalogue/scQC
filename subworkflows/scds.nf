include { SCDS } from "$projectDir/modules/scds"

workflow SCDS_WF {
    take:
    my_cellranger_outdir_ch

    main:
    scds_input_ch = my_cellranger_outdir_ch.map{ it -> [
        it[0], 
        "${it[1]}/filtered_feature_bc_matrix"
    ]}

    SCDS(scds_input_ch)

    emit:
    scds_outdir = SCDS.out.outdir
    scds_n_dbls = SCDS.out.n_dbls
}
