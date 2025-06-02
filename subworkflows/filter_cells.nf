include { FILTER_CELLS } from "$projectDir/modules/filter_cells"

workflow FILTER_CELLS_WF {
    take:
    cellranger_outdir_ch

    main:
    n_cells_outpath = "n_cells.tsv"

    FILTER_CELLS(cellranger_outdir_ch, n_cells_outpath)

    n_cells_ch = FILTER_CELLS.out.n_cells.collectFile(
        name: n_cells_outpath, 
        keepHeader: true, 
        skip: 1
    )
    //n_cells_ch.view()

    emit:
    my_cellranger_outdir = FILTER_CELLS.out.my_cellranger_outdir
    filtered_mtx_outdir = FILTER_CELLS.out.filtered_mtx_outdir
    n_cells = n_cells_ch
    /**/
}
