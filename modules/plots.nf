process PLOTS {
    container "quay.io/peepk/scqc_py:v1.0.0"
    publishDir "$params.outdir/postprocess/plots", mode: 'copy'

    input:
    path(n_cells)
    path(n_singlets)

    script:
    """
    summary_plots.py $n_cells $n_singlets
    """

    output:
    path("*.png")
    path("*.tsv")
    /**/
}
