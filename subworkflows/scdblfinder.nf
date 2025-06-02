include { SCDBLFINDER } from "$projectDir/modules/scdblfinder"

workflow SCDBLFINDER_WF {
    take:
    my_cellranger_outdir_ch

    main:
    scdblfinder_input_ch = my_cellranger_outdir_ch.map{ it -> [
        it[0], 
        "${it[1]}/filtered_feature_bc_matrix"
    ]}

    SCDBLFINDER(scdblfinder_input_ch)

    emit:
    scdblfinder_outdir = SCDBLFINDER.out.outdir
    scdblfinder_n_dbls = SCDBLFINDER.out.n_dbls
}
