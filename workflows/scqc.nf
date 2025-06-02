include { FILTER_CELLS_WF } from "$projectDir/subworkflows/filter_cells"
include { FILTER_VCF } from "$projectDir/modules/filter_vcf_common"
include { VIREO_WF } from "$projectDir/subworkflows/vireo"
include { DEMUXALOT_WF } from "$projectDir/subworkflows/demuxalot"
include { SCDS_WF } from "$projectDir/subworkflows/scds"
include { SCDBLFINDER_WF } from "$projectDir/subworkflows/scdblfinder"
include { DOUBLETFINDER_WF } from "$projectDir/subworkflows/doubletfinder"
include { ASSIGNMENT_WF } from "$projectDir/subworkflows/assignment"
include { PLOTS } from "$projectDir/modules/plots"

workflow SCQC {
    take:
    vcf_ch  // Tuple: VCF, idx
    exonic_regions_ch  // BED
    samples_ch  // Queue channel of tuples: sample id, cellranger outdir, pool individuals, cell metadata

    main:
    //vcf_ch.view()
    //exonic_regions_ch.view()
    //samples_ch.view()

    // Separate channels for different files
    cellranger_outdir_ch = samples_ch.map { it -> [it[0], it[1]] }
    inds_ch = samples_ch.map { it -> [it[0], it[2]] }
    cell_meta_ch = samples_ch.map { it -> [it[0], it[3]] }

    // Throw out low quality cells and normalize/convert the countmatrices 
    // to various formats needed in subsequent steps.
    FILTER_CELLS_WF(cellranger_outdir_ch)
    qc_passed_cellranger_out_ch = FILTER_CELLS_WF.out.my_cellranger_outdir
    

    // Filter the VCF to keep only each pool's samples 
    // and common (MAF > 0.05) SNPs that overlap genes (exons).
    FILTER_VCF(inds_ch.combine(vcf_ch), exonic_regions_ch)
    filtered_vcf_ch = FILTER_VCF.out.filtered_vcf

    // Run demultiplexing tools
    VIREO_WF(filtered_vcf_ch, qc_passed_cellranger_out_ch)
    DEMUXALOT_WF(filtered_vcf_ch, qc_passed_cellranger_out_ch, inds_ch)

    // Run doublet detection tools
    SCDS_WF(qc_passed_cellranger_out_ch)
    SCDBLFINDER_WF(qc_passed_cellranger_out_ch)
    // SCDS and SCDBLFINDER outputs are used to estimate the number of doublets 
    // which DOUBLETFINDER needs as input.

    DOUBLETFINDER_WF(
        FILTER_CELLS_WF.out.filtered_mtx_outdir, 
        SCDS_WF.out.scds_n_dbls, 
        SCDBLFINDER_WF.out.scdblfinder_n_dbls
    )

    // Add donor labels to HQ cells
    ASSIGNMENT_WF(
        VIREO_WF.out.vireo_out, 
        DEMUXALOT_WF.out.demuxalot_out, 
        SCDS_WF.out.scds_outdir, 
        SCDBLFINDER_WF.out.scdblfinder_outdir, 
        DOUBLETFINDER_WF.out.doubletfinder_outdir, 
        cell_meta_ch, 
        FILTER_CELLS_WF.out.filtered_mtx_outdir
    )

    // Plot summaries of filtering and comparison to metadata
    PLOTS(FILTER_CELLS_WF.out.n_cells, ASSIGNMENT_WF.out.n_singlets)
    /**/
}
