include { DOUBLETFINDER } from "$projectDir/modules/doubletfinder"

workflow DOUBLETFINDER_WF {
    take:
    samples_ch
    scds_n_dbls_ch
    scdblfinder_n_dbls_ch

    main:
    expected_n_dbls_ch = scds_n_dbls_ch.join(scdblfinder_n_dbls_ch)
    //expected_n_dbls_ch.view()

    samples_ch = samples_ch.map{ it -> [
        it[0], 
        "${it[1]}/filtered_cells.norm.rds"
    ]}

    DOUBLETFINDER(samples_ch.join(expected_n_dbls_ch))

    emit:
    doubletfinder_outdir = DOUBLETFINDER.out.outdir
    /**/
}
