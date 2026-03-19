/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SDATA_READ_VISIUM      } from '../modules/local/sdata/read_visium/main'
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
    samplesheet                    //    file: /path/to/samplesheet
    spaceranger_reference          //     dir: /path/to/reference
    spaceranger_probeset           //    file: /path/to/csv
    hd_bin_size                    // integer: Bin size for Visium HD
    qc_min_counts                  // integer: Minimum UMIs per spot
    qc_min_genes                   // integer: Minimum genes per spot
    qc_min_spots                   // integer: Minimum spots per gene
    qc_mito_threshold              //   float: Maximum mito. content per spot
    qc_ribo_threshold              //   float: Minimum ribo. content per spot
    qc_hb_threshold                //   float: Maximum haem. content per spot
    normalize_target_sum           //  string: Target sum of total count normalization
    n_highly_variable_genes        // integer: Number of HVGs to use
    hvg_flavor                     //  string: Flavor for HVG calculations
    n_principal_components         // integer: Number of principal components to compute
    pca_use_highly_variable        // boolean: Whether to only use highly variable genes for PCA
    n_neighbours                   // integer: Number of nearest neighbours to compute
    neighbours_n_pcs               // integer: Number of PCs to use for nearest neighbours
    umap_min_dist                  //   float: Minimum distance between embedded points
    umap_spread                    //   float: Scale of embedded points
    cluster_resolution             //   float: Spot clustering resolution
    cluster_key_added              //  string: Obs key where cluster labels are added
    rank_genes_group_by            //  string: Column name to group by for differential expression testing
    rank_genes_method              //  string: Method to use for differential expression testing
    spatial_coord_type             //  string: Type of spatial coordinate system
    spatial_n_neighbours           // integer: Number of spatial neighbourhoods
    spatial_cluster_key            //  string: Obs key where spatial cluster labels are added
    svg_autocorr_method            //  string: Autocorrelation method
    n_top_svgs                     // integer: Number of variable genes to plot
    merge_sdata                    // boolean: Whether to merge sdata or not
    integrate_sdata                // boolean: Whether to integrate sdata or not
    integration_cluster_resolution //   float: Integration cluster resolution
    integration_n_hvgs             // integer: Number of HVGs to integrate with
    multiqc_config                 //    file: /path/to/multiqc/config
    multiqc_logo                   //    file: /path/to/multiqc/logo
    multiqc_methods_description    //    file: /path/to/multiqc/description
    outdir                         //     dir: /path/to/output/directory

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
        INPUT_CHECK.out.ch_spaceranger_input
            .map { it -> [it[0], it[1]] } // [ meta, reads ]
    )
    ch_multiqc_files = ch_multiqc_files
        .mix(FASTQC.out.zip.collect { it -> it[1] })

    //
    // SUBWORKFLOW: Space Ranger raw data processing
    //
    SPACERANGER (
        INPUT_CHECK.out.ch_spaceranger_input,
        spaceranger_reference,
        spaceranger_probeset,
    )
    ch_multiqc_files = ch_multiqc_files
        .mix(SPACERANGER.out.sr_dir.collect { it -> it[1] })

    //
    // Combine pre-existing spaceranger outputs with newly processed ones
    //
    ch_spaceranger_dir = INPUT_CHECK.out.ch_downstream_input
        .mix(SPACERANGER.out.sr_dir)

    //
    // MODULE: Read ST data and save as SpatialData
    //
    SDATA_READ_VISIUM (
        ch_spaceranger_dir,
        hd_bin_size
    )

    //
    // SUBWORKFLOW: Downstream analyses of ST data
    // This includes: QC, filtering, normalization,
    // dimensionality reduction, clustering, and spatial analysis
    //
    DOWNSTREAM (
        SDATA_READ_VISIUM.out.sdata,
        qc_min_counts,
        qc_min_genes,
        qc_min_spots,
        qc_mito_threshold,
        qc_ribo_threshold,
        qc_hb_threshold,
        normalize_target_sum,
        n_highly_variable_genes,
        hvg_flavor,
        n_principal_components,
        pca_use_highly_variable,
        n_neighbours,
        neighbours_n_pcs,
        umap_min_dist,
        umap_spread,
        cluster_resolution,
        cluster_key_added,
        rank_genes_group_by,
        rank_genes_method,
        spatial_coord_type,
        spatial_n_neighbours,
        spatial_cluster_key,
        svg_autocorr_method,
        n_top_svgs
    )

    //
    // SUBWORKFLOW: Sample aggregation (optional)
    //
    // TODO: add back integration
    /*if (merge_sdata || integrate_sdata) {
        AGGREGATION (
            DOWNSTREAM.out.sdata_svg,
            merge_sdata,
            integrate_sdata,
            integration_cluster_resolution,
            integration_n_hvgs
        )
    }*/

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
            name: 'nf_core_spatialvi_software_mqc_versions.yml',
            sort: true,
            newLine: true
        )
        .set { ch_collated_versions }

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
    ch_multiqc_files    = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))

    ch_multiqc_custom_methods_description = multiqc_methods_description ?
        file(multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml",
             checkIfExists: true)
    ch_methods_description = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    //
    // Add filter statistics to MultiQC (JSON format needs custom config)
    //
    ch_multiqc_files = ch_multiqc_files.mix(
        DOWNSTREAM.out.filter_stats.map { _meta, json -> json }.collect()
    )

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
    // Final SpatialData outputs
    sdata_raw       = SDATA_READ_VISIUM.out.sdata      // channel: [ meta, zarr ]
    sdata_svg       = DOWNSTREAM.out.sdata_svg         // channel: [ meta, zarr ]

    // Reports TODO
    //qc_html         = DOWNSTREAM.out.qc_html           // channel: [ meta, html ]
    //clustering_html = DOWNSTREAM.out.clustering_html   // channel: [ meta, html ]
    //svg_html        = DOWNSTREAM.out.svg_html          // channel: [ meta, html ]

    // SVG results
    svg_csv         = DOWNSTREAM.out.svg_csv           // channel: [ meta, csv ]

    // MultiQC
    multiqc_report  = MULTIQC.out.report.toList()      // channel: [ multiqc_report.html ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
