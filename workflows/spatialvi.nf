/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SDATA_READ_VISIUM         } from "../modules/local/sdata/read_visium"
include { FASTQC                    } from "../modules/nf-core/fastqc"
include { SCANPY_RANK_GENES_GROUPS  } from "../modules/local/scanpy/rank_genes_groups"
include { SDATA_MERGE               } from "../modules/local/sdata/merge"
include { SDATA_TO_LEGACY_ANNDATA   } from "../modules/local/sdata/to_legacy_anndata"
include { MULTIQC                   } from "../modules/nf-core/multiqc"
include { QUARTO_NOTEBOOK as REPORT } from "../modules/nf-core/quarto/notebook"
include { SDATA_UPDATE_TABLE        } from "../modules/local/sdata/update_table"

include { INPUT_CHECK               } from "../subworkflows/local/input_check"
include { SPACERANGER               } from "../subworkflows/local/spaceranger"
include { PREPROCESSING             } from "../subworkflows/local/preprocessing"
include { CLUSTERING                } from "../subworkflows/local/clustering"
include { SPATIAL                   } from "../subworkflows/local/spatial"
include { INTEGRATION               } from "../subworkflows/local/integration"
include { paramsSummaryMultiqc      } from "../subworkflows/nf-core/utils_nfcore_pipeline"
include { paramsSummaryMap          } from "plugin/nf-schema"
include { softwareVersionsToYAML    } from "../subworkflows/nf-core/utils_nfcore_pipeline"
include { methodsDescriptionText    } from "../subworkflows/local/utils_nfcore_spatialvi_pipeline"

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
    skip_downstream                // boolean: Whether to skip downstream steps or not
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
    n_neighbors                    // integer: Number of nearest neighbors to compute
    neighbors_n_pcs                // integer: Number of PCs to use for nearest neighbors
    umap_min_dist                  //   float: Minimum distance between embedded points
    umap_spread                    //   float: Scale of embedded points
    cluster_resolution             //   float: Spot clustering resolution
    rank_genes_method              //  string: Method to use for differential expression testing
    spatial_coord_type             //  string: Type of spatial coordinate system
    spatial_n_neighbors            // integer: Number of spatial neighborhoods
    svg_autocorr_method            //  string: Autocorrelation method
    n_top_svgs                     // integer: Number of variable genes to plot
    skip_integration               // boolean: Whether to integrate sdata or not
    integration_method             //  string: Integration method to use
    integration_cluster_resolution //   float: Integration cluster resolution
    multiqc_config                 //    file: /path/to/multiqc/config
    multiqc_logo                   //    file: /path/to/multiqc/logo
    multiqc_methods_description    //    file: /path/to/multiqc/description
    outdir                         //     dir: /path/to/output/directory

    main:

    def ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    // =========================================================================
    // VALIDATION AND SPACE RANGER PROCESSING
    // =========================================================================

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
    // MODULE: Read ST data and save as SpatialData
    //
    ch_spaceranger_dir = INPUT_CHECK.out.ch_downstream_input
        .mix(SPACERANGER.out.sr_dir)
    SDATA_READ_VISIUM (
        ch_spaceranger_dir,
        hd_bin_size
    )
    ch_sdata_raw = SDATA_READ_VISIUM.out.sdata

    // =========================================================================
    // PER-SAMPLE ANALYSES
    // =========================================================================

    ch_adata              = channel.empty()
    ch_svg_csv            = channel.empty()
    ch_sdata_output       = channel.empty()
    ch_report_html        = channel.empty()
    ch_report_notebook    = channel.empty()
    ch_report_params_yaml = channel.empty()
    ch_report_artifacts   = channel.empty()

    if (!skip_downstream) {

        //
        // MODULE: Extract legacy AnnData for scanpy processing
        //
        SDATA_TO_LEGACY_ANNDATA (
            ch_sdata_raw
        )

        //
        // SUBWORKFLOW: Pre-processing
        //
        PREPROCESSING (
            SDATA_TO_LEGACY_ANNDATA.out.adata,
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
            pca_use_highly_variable
        )
        ch_multiqc_files = ch_multiqc_files
            .mix(PREPROCESSING.out.filter_stats.collect { it -> it[1] })

        //
        // SUBWORKFLOW: Clustering
        //
        use_rep = ''
        umap_key_added = 'X_umap'
        leiden_key_added = 'clusters'
        CLUSTERING (
            PREPROCESSING.out.adata,
            n_neighbors,
            neighbors_n_pcs,
            use_rep,
            umap_min_dist,
            umap_spread,
            umap_key_added,
            cluster_resolution,
            leiden_key_added
        )

        //
        // MODULE: Differential expression analysis
        //
        rank_genes_group_by = 'clusters'
        SCANPY_RANK_GENES_GROUPS (
            CLUSTERING.out.adata,
            rank_genes_group_by,
            rank_genes_method
        )

        //
        // SUBWORKFLOW: Spatial analyses
        //
        SPATIAL (
            SCANPY_RANK_GENES_GROUPS.out.adata,
            spatial_coord_type,
            spatial_n_neighbors,
            svg_autocorr_method
        )
        ch_adata   = SPATIAL.out.adata
        ch_svg_csv = SPATIAL.out.svg_csv

        //
        // MODULE: Update SpatialData with AnnData results
        //
        SDATA_UPDATE_TABLE (
            ch_sdata_raw.join(ch_adata),
            ''
        )
        ch_sdata_output = SDATA_UPDATE_TABLE.out.sdata

        //
        // MODULE: Per-sample Quarto reports
        //
        report_notebook = file("${projectDir}/bin/report.qmd", checkIfExists: true)
        extensions = channel.fromPath("${projectDir}/assets/_extensions").collect()
        ch_report_input_data = ch_sdata_output
            .map { it -> it[1] }
        ch_report_notebook = ch_sdata_output
            .map { it -> it[0] }
            .combine(channel.value(report_notebook))
            .map { meta, notebook -> tuple(meta, notebook) }
        ch_report_params = ch_sdata_output
            .map { _meta, sdata ->
                [
                    input_sdata  : sdata.name,
                    n_top_svgs   : n_top_svgs,
                    artifact_dir : "artifacts"
                ]
            }
        REPORT (
            ch_report_notebook,
            ch_report_params,
            ch_report_input_data,
            extensions
        )
        ch_report_html        = REPORT.out.html
        ch_report_notebook    = REPORT.out.notebook
        ch_report_params_yaml = REPORT.out.params_yaml
        ch_report_artifacts   = REPORT.out.artifacts

    } else {
        ch_sdata_output = ch_sdata_raw
    }

    // =========================================================================
    // MERGING AND INTEGRATION
    // =========================================================================

    //
    // MODULE: Merge per-sample SpatialData objects into one
    //
    ch_sdata_output_sorted = ch_sdata_output
        .toSortedList { a, b -> a[0].id <=> b[0].id }
        .flatMap()
        .map { _meta, zarr -> zarr }
        .collect()
    SDATA_MERGE (
        ch_sdata_output_sorted
    )
    ch_sdata_merged = SDATA_MERGE.out.sdata

    ch_integration_html        = channel.empty()
    ch_integration_notebook    = channel.empty()
    ch_integration_params_yaml = channel.empty()
    ch_integration_artifacts   = channel.empty()
    ch_sdata_integrated        = channel.empty()
    ch_adata_integrated        = channel.empty()

    //
    // SUBWORKFLOW: Sample aggregation (optional)
    //
    if (!skip_downstream && !skip_integration) {
        INTEGRATION (
            ch_sdata_merged,
            ch_adata,
            integration_method,
            n_neighbors,
            neighbors_n_pcs,
            umap_min_dist,
            umap_spread,
            integration_cluster_resolution
        )
        ch_integration_html        = INTEGRATION.out.html
        ch_integration_notebook    = INTEGRATION.out.notebook
        ch_integration_params_yaml = INTEGRATION.out.params_yaml
        ch_integration_artifacts   = INTEGRATION.out.artifacts
        ch_sdata_integrated        = INTEGRATION.out.sdata
        ch_adata_integrated        = INTEGRATION.out.adata
    }

    // =========================================================================
    // FINALISATION AND MULTIQC
    // =========================================================================

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

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_'  +  'spatialvi_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'spatialvi'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )

    emit:
    // SpatialData outputs
    sdata_raw               = ch_sdata_raw               // channel: [ meta, zarr ]
    sdata                   = ch_sdata_output            // channel: [ meta, zarr ]
    sdata_merged            = ch_sdata_merged            // channel: [ zarr ]
    sdata_integrated        = ch_sdata_integrated        // channel: [ zarr ]

    // AnnData outputs
    adata                   = ch_adata                   // channel: [ meta, h5ad ]
    adata_integrated        = ch_adata_integrated        // channel: [ h5ad ]

    // Per-sample report outputs
    report_html             = ch_report_html             // channel: [ meta, html ]
    report_notebook         = ch_report_notebook         // channel: [ meta, qmd ]
    report_params_yaml      = ch_report_params_yaml      // channel: [ meta, yml ]
    report_artifacts        = ch_report_artifacts        // channel: [ meta, dir ]

    // SVG results
    svg_csv                 = ch_svg_csv                 // channel: [ meta, csv ]

    // Integration report outputs
    integration_html        = ch_integration_html        // channel: [ meta, html ]
    integration_notebook    = ch_integration_notebook    // channel: [ meta, qmd ]
    integration_params_yaml = ch_integration_params_yaml // channel: [ meta, yml ]
    integration_artifacts   = ch_integration_artifacts   // channel: [ meta, dir ]

    // MultiQC
    multiqc_report          = MULTIQC.out.report.toList() // channel: [ html ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
