//
// Subworkflow for downstream analyses of ST data
//

include { QUARTONOTEBOOK as ST_QUALITY_CONTROLS } from '../../modules/nf-core/quartonotebook/main'
include { QUARTONOTEBOOK as ST_SPATIAL_DE       } from '../../modules/nf-core/quartonotebook/main'
include { QUARTONOTEBOOK as ST_CLUSTERING       } from '../../modules/nf-core/quartonotebook/main'

workflow ST_DOWNSTREAM {

    take:
    st_adata_raw

    main:

    ch_versions = Channel.empty()

    //
    // Quarto reports and extension files
    //
    quality_controls_file = file("${projectDir}/bin/st_quality_controls.qmd", checkIfExists: true)
    clustering_file = file("${projectDir}/bin/st_clustering.qmd", checkIfExists: true)
    spatial_de_file = file("${projectDir}/bin/st_spatial_de.qmd", checkIfExists: true)
    extensions = Channel.fromPath("${projectDir}/assets/_extensions").collect()

    // input:
    // tuple val(meta), path(notebook)
    // val parameters
    // path input_files
    // path extensions

    // output:
    // tuple val(meta), path("*.html")     , emit: html
    // tuple val(meta), path("${notebook}"), emit: notebook
    // tuple val(meta), path("artifacts/*"), emit: artifacts, optional: true
    // tuple val(meta), path("params.yml") , emit: params_yaml, optional: true
    // tuple val(meta), path("_extensions"), emit: extensions, optional: true
    // path "versions.yml"                 , emit: versions

    //
    // Quality controls and filtering
    //
    ch_quality_controls_input_data = st_adata_raw
        .map { it -> it[1] }
    ch_quality_controls_notebook = st_adata_raw
        .map { tuple(it[0], quality_controls_file) }
    quality_controls_params = [
        input_raw_data: "st_adata_raw.h5ad",
        min_counts: params.st_qc_min_counts,
        min_genes: params.st_qc_min_genes,
        min_spots: params.st_qc_min_spots,
        mito_threshold: params.st_qc_mito_threshold,
        ribo_threshold: params.st_qc_ribo_threshold,
        hb_threshold: params.st_qc_hb_threshold,
        output_adata_filtered: "st_adata_filtered.h5ad"
    ]
    ST_QUALITY_CONTROLS (
        ch_quality_controls_notebook,
        quality_controls_params,
        ch_quality_controls_input_data,
        extensions
    )
    ch_versions = ch_versions.mix(ST_QUALITY_CONTROLS.out.versions)

    //
    // Normalisation, dimensionality reduction and clustering
    //
    ch_clustering_input_data = ST_QUALITY_CONTROLS.out.artifacts
        .map { it -> it[1] }
    ch_clustering_notebook = ST_QUALITY_CONTROLS.out.artifacts
        .map { tuple(it[0], clustering_file) }
    clustering_params = [
        input_adata_filtered: "st_adata_filtered.h5ad",
        cluster_resolution: params.st_cluster_resolution,
        n_hvgs: params.st_cluster_n_hvgs,
        output_adata_processed: "st_adata_processed.h5ad"
    ]
    ST_CLUSTERING (
        ch_clustering_notebook,
        clustering_params,
        ch_clustering_input_data,
        extensions
    )
    ch_versions = ch_versions.mix(ST_CLUSTERING.out.versions)

    //
    // Spatial differential expression
    //
    ch_spatial_de_input_data = ST_CLUSTERING.out.artifacts
        .map { it -> it[1] }
    ch_spatial_de_notebook = ST_CLUSTERING.out.artifacts
        .map { tuple(it[0], spatial_de_file) }
    spatial_de_params = [
        input_adata_processed: "st_adata_processed.h5ad",
        n_top_spatial_degs: params.st_n_top_spatial_degs,
        output_spatial_degs: "st_spatial_de.csv"
    ]
    ST_SPATIAL_DE (
        ch_spatial_de_notebook,
        spatial_de_params,
        ch_spatial_de_input_data,
        extensions
    )
    ch_versions = ch_versions.mix(ST_SPATIAL_DE.out.versions)

    emit:
    st_qc_html             = ST_QUALITY_CONTROLS.out.html        // channel: [ meta, html ]
    st_adata_filtered      = ST_QUALITY_CONTROLS.out.artifacts   // channel: [ meta, h5ad ]
    st_qc_notebook         = ST_QUALITY_CONTROLS.out.notebook    // channel: [ meta, qmd ]
    st_qc_params           = ST_QUALITY_CONTROLS.out.params_yaml // channel: [ meta, yml ]

    st_clustering_html     = ST_CLUSTERING.out.html              // channel: [ html ]
    st_adata_processed     = ST_CLUSTERING.out.artifacts         // channel: [ meta, h5ad]
    st_clustering_notebook = ST_CLUSTERING.out.notebook          // channel: [ meta, qmd ]
    st_clustering_params   = ST_CLUSTERING.out.params_yaml       // channel: [ meta, yml ]

    st_spatial_html        = ST_SPATIAL_DE.out.html              // channel: [ meta, html ]
    st_degs                = ST_SPATIAL_DE.out.artifacts         // channel: [ meta, csv ]
    st_spatial_notebook    = ST_SPATIAL_DE.out.notebook          // channel: [ meta, qmd ]
    st_spatial_params      = ST_SPATIAL_DE.out.params_yaml       // channel: [ meta, yml ]

    versions               = ch_versions                         // channel: [ versions.yml ]
}
