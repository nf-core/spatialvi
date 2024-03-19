/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { READ_DATA              } from '../modules/local/read_data'
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { INPUT_CHECK            } from '../subworkflows/local/input_check'
include { SPACERANGER            } from '../subworkflows/local/spaceranger'
include { DOWNSTREAM             } from '../subworkflows/local/downstream'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_spatialtranscriptomics_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SPATIALTRANSCRIPTOMICS {

    take:
    samplesheet // file: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // SUBWORKFLOW: Read and validate samplesheet
    //
    INPUT_CHECK (
        samplesheet
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // MODULE: FastQC
    //
    FASTQC(
        INPUT_CHECK.out.ch_spaceranger_input.map{ it -> [it[0] /* meta */, it[1] /* reads */]}
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})

    //
    // SUBWORKFLOW: Space Ranger raw data processing
    //
    DOWNSTREAM_REQUIRED_SPACERANGER_FILES = [
        "raw_feature_bc_matrix.h5",
        "tissue_positions.csv",
        "scalefactors_json.json",
        "tissue_hires_image.png",
        "tissue_lowres_image.png"
    ]
    SPACERANGER (
        INPUT_CHECK.out.ch_spaceranger_input
    )
    ch_versions = ch_versions.mix(SPACERANGER.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(SPACERANGER.out.sr_dir.collect{it[1]})
    ch_downstream_input = INPUT_CHECK.out.ch_downstream_input.concat(SPACERANGER.out.sr_dir).map{
        meta, outs -> [meta, outs.findAll{ it -> DOWNSTREAM_REQUIRED_SPACERANGER_FILES.contains(it.name) }]
    }

    //
    // MODULE: Read ST data and save as `anndata`
    //
    READ_DATA (
        ch_downstream_input
    )
    ch_versions = ch_versions.mix(READ_DATA.out.versions)

    //
    // SUBWORKFLOW: Downstream analyses of ST data
    //
    DOWNSTREAM (
        READ_DATA.out.sdata_raw
    )
    ch_versions = ch_versions.mix(DOWNSTREAM.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
