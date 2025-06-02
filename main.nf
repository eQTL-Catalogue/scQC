nextflow.enable.dsl=2

include { SCQC } from "$projectDir/workflows/scqc"

workflow {
    // Even though the command line params are written with a hyphen (e.g. "exonic-regions"),
    // Nextflow implicitly creates a camelCase param ("exonicRegions") as well
    // because hyphenated param names cannot be used in Nextflow code.
    //println params
    
    println("Samplesheet: $params.samples")
    println("VCF: $params.vcf")
    println("Genotype field: $params.GTField")
    println("Exonic regions: $params.exonicRegions")
    println("Mitochondrial content threshold: $params.MTThreshold%")
    println("Demultiplexing tools: $params.demuxTools")
    println("Doublet detection tools: $params.dblTools")
    println("Output directory: $params.outdir\n")

    // Determine the vcf index based on which file (.tbi or .csi) is present in the vcf directory
    vcf_ch = Channel.fromPath(params.vcf, checkIfExists: true)
        .map{ it -> [it, file(it + '.tbi').exists() ? it + '.tbi' : it + '.csi'] }
    //vcf_ch.view()

    exonic_regions_ch = Channel.fromPath(params.exonicRegions, checkIfExists: true)
    
    samples_ch = Channel.fromPath(params.samples, checkIfExists: true)
        .splitCsv(header: true, sep: "\t", strip: true)
        .map( row -> [row.sample, row.cellranger_out, row.inds, row.cell_meta] )

    SCQC(vcf_ch, exonic_regions_ch, samples_ch)
    /**/
}
