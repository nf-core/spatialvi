/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { READ_DATA              } from '../modules/local/read_data'
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { INPUT_CHECK            } from '../subworkflows/local/input_check'
include { SPACERANGER            } from '../subworkflows/local/spaceranger'
include { DOWNSTREAM             } from '../subworkflows/local/downstream'
include { AGGREGATION            } from '../subworkflows/local/aggregation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_spatialvi_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SPATIALVI {

    take:
    samplesheet                    // file   : /path/to/samplesheet
    spaceranger_reference          // dir    : /path/to/reference
    spaceranger_probeset           // file   : /path/to/csv
    hd_bin_size                    // integer: Bin size for Visium HD
    qc_min_counts                  // integer: Minimum UMIs per spot
    qc_min_genes                   // integer: Minimum genes per spot
    qc_min_spots                   // integer: Minimum spots per gene
    qc_mito_threshold              // float  : Maximum mito. content per spot
    qc_ribo_threshold              // float  : Minimum ribo. content per spot
    qc_hb_threshold                // float  : Maximum haem. content per spot
    cluster_n_hvgs                 // integer: Number of HVGs to use
    cluster_resolution             // float  : Spot clustering resolution
    svg_autocorr_method            // string : Autocorrelation method
    n_top_svgs                     // integer: Number of variable genes to plot
    merge_sdata                    // boolean: Whether to merge sdata or not
    integrate_sdata                // boolean: Whether to integrate sdata or not
    integration_cluster_resolution // float  : Integration cluster resolution
    integration_n_hvgs             // integer: Number of HVGs to integrate with
    multiqc_config                 // file   : /path/to/multiqc/config
    multiqc_logo                   // file   : /path/to/multiqc/logo
    multiqc_methods_description    // file   : /path/to/multiqc/description
    outdir                         // dir    : /path/to/output/directory

    main:

    ch_multiqc_files = channel.empty()

    //
    // SUBWORKFLOW: Read and validate samplesheet
    //
    INPUT_CHECK (
        samplesheet,
        hd_bin_size
    )

    //
    // MODULE: FastQC
    //
    FASTQC(
        INPUT_CHECK.out.ch_spaceranger_input.map{ it -> [it[0] /* meta */, it[1] /* reads */]}
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{ it -> it[1] })

    //
    // SUBWORKFLOW: Space Ranger raw data processing
    //
    SPACERANGER (
        INPUT_CHECK.out.ch_spaceranger_input,
        spaceranger_reference,
        spaceranger_probeset,
    )
    ch_multiqc_files = ch_multiqc_files.mix(SPACERANGER.out.sr_dir.collect{ it -> it[1] })
    ch_downstream_input = INPUT_CHECK.out.ch_downstream_input
        .mix(SPACERANGER.out.sr_dir)

    //
    // MODULE: Read ST data and save as `SpatialData`
    //
    READ_DATA (
        ch_downstream_input,
        hd_bin_size
    )

    //
    // SUBWORKFLOW: Downstream analyses of ST data
    //
    DOWNSTREAM (
        READ_DATA.out.sdata_raw,
        qc_min_counts,
        qc_min_genes,
        qc_min_spots,
        qc_mito_threshold,
        qc_ribo_threshold,
        qc_hb_threshold,
        cluster_n_hvgs,
        cluster_resolution,
        svg_autocorr_method,
        n_top_svgs,
    )

    //
    // SUBWORKFLOW: Sample aggregation (optional)
    //
    if (merge_sdata || integrate_sdata) {
        AGGREGATION (
            DOWNSTREAM.out.svg_sdata,
            merge_sdata,
            integrate_sdata,
            integration_cluster_resolution,
            integration_n_hvgs,
        )
    }

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(topic_versions.versions_file)
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_' + 'spatialvi_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = multiqc_config ?
        channel.fromPath(multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = multiqc_logo ?
        channel.fromPath(multiqc_logo, checkIfExists: true) :
        channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = multiqc_methods_description ?
        file(multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(
        DOWNSTREAM.out.qc_mqc.map{ it -> it[1] }.collect()
    )
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: [ multiqc_report.html ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
