include { COMBINE_RESULTS } from "$projectDir/modules/combine_results"
include { ASSIGNMENT } from "$projectDir/modules/assignment"

workflow ASSIGNMENT_WF {
    take:
    vireo_ch
    demuxalot_ch
    scds_ch
    scdblfinder_ch
    doubletfinder_ch
    cell_meta_ch
    mtx_conversions

    main:
    COMBINE_RESULTS(vireo_ch.join(demuxalot_ch).join(scds_ch).join(scdblfinder_ch).join(doubletfinder_ch))
    //COMBINE_RESULTS.out.combined_results.view()

    n_singlets_outpath = "n_singlets.tsv"

    ASSIGNMENT(COMBINE_RESULTS.out.combined_results.join(cell_meta_ch).join(mtx_conversions), n_singlets_outpath)

    n_singlets_ch = ASSIGNMENT.out.n_singlets.collectFile(
        name: n_singlets_outpath, 
        keepHeader: true, 
        skip: 1
    )
    //n_singlets_ch.view()

    emit:
    assignment_outdir = ASSIGNMENT.out.outdir
    n_singlets = n_singlets_ch
}
