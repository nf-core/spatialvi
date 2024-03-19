//
// Subworkflow for downstream analyses of ST data
//

include { QUARTONOTEBOOK as ST_QUALITY_CONTROLS         } from '../../modules/nf-core/quartonotebook/main'
include { QUARTONOTEBOOK as ST_SPATIALLY_VARIABLE_GENES } from '../../modules/nf-core/quartonotebook/main'
include { QUARTONOTEBOOK as ST_CLUSTERING               } from '../../modules/nf-core/quartonotebook/main'

workflow ST_DOWNSTREAM {

    take:
    st_sdata_raw

    main:

    ch_versions = Channel.empty()

    //
    // Quarto reports and extension files
    //
    quality_controls_notebook = file("${projectDir}/bin/st_quality_controls.qmd", checkIfExists: true)
    clustering_notebook = file("${projectDir}/bin/st_clustering.qmd", checkIfExists: true)
    spatially_variable_genes_notebook = file("${projectDir}/bin/st_spatially_variable_genes.qmd", checkIfExists: true)
    extensions = Channel.fromPath("${projectDir}/assets/_extensions").collect()

    //
    // Quality controls and filtering
    //
    ch_quality_controls_input_data = st_sdata_raw
        .map { it -> it[1] }
    ch_quality_controls_notebook = st_sdata_raw
        .map { tuple(it[0], quality_controls_notebook) }
    quality_controls_params = [
        input_sdata: "st_sdata_raw.zarr",
        min_counts: params.st_qc_min_counts,
        min_genes: params.st_qc_min_genes,
        min_spots: params.st_qc_min_spots,
        mito_threshold: params.st_qc_mito_threshold,
        ribo_threshold: params.st_qc_ribo_threshold,
        hb_threshold: params.st_qc_hb_threshold,
        artifact_dir: "artifacts",
        output_adata: "st_adata_filtered.h5ad",
        output_sdata: "st_sdata_filtered.zarr",
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
        .map { tuple(it[0], clustering_notebook) }
    clustering_params = [
        input_sdata: "st_sdata_filtered.zarr",
        cluster_resolution: params.st_cluster_resolution,
        n_hvgs: params.st_cluster_n_hvgs,
        artifact_dir: "artifacts",
        output_adata: "st_adata_processed.h5ad",
        output_sdata: "st_sdata_processed.zarr",
    ]
    ST_CLUSTERING (
        ch_clustering_notebook,
        clustering_params,
        ch_clustering_input_data,
        extensions
    )
    ch_versions = ch_versions.mix(ST_CLUSTERING.out.versions)

    //
    // Spatially variable genes
    //
    ch_spatially_variable_genes_input_data = ST_CLUSTERING.out.artifacts
        .map { it -> it[1] }
    ch_spatially_variable_genes_notebook = ST_CLUSTERING.out.artifacts
        .map { tuple(it[0], spatially_variable_genes_notebook) }
    spatially_variable_genes_params = [
        input_sdata: "st_sdata_processed.zarr",
        n_top_spatial_degs: params.st_n_top_spatial_degs,
        artifact_dir: "artifacts",
        output_csv: "st_spatially_variable_genes.csv",
        output_adata: "st_adata_spatially_variable_genes.h5ad",
        output_sdata: "st_sdata.zarr",
    ]
    ST_SPATIALLY_VARIABLE_GENES (
        ch_spatially_variable_genes_notebook,
        spatially_variable_genes_params,
        ch_spatially_variable_genes_input_data,
        extensions
    )
    ch_versions = ch_versions.mix(ST_SPATIALLY_VARIABLE_GENES.out.versions)

    emit:
    st_qc_html             = ST_QUALITY_CONTROLS.out.html                // channel: [ meta, html ]
    st_sdata_filtered      = ST_QUALITY_CONTROLS.out.artifacts           // channel: [ meta, h5ad ]
    st_qc_notebook         = ST_QUALITY_CONTROLS.out.notebook            // channel: [ meta, qmd ]
    st_qc_params           = ST_QUALITY_CONTROLS.out.params_yaml         // channel: [ meta, yml ]

    st_clustering_html     = ST_CLUSTERING.out.html                      // channel: [ html ]
    st_sdata_processed     = ST_CLUSTERING.out.artifacts                 // channel: [ meta, h5ad]
    st_clustering_notebook = ST_CLUSTERING.out.notebook                  // channel: [ meta, qmd ]
    st_clustering_params   = ST_CLUSTERING.out.params_yaml               // channel: [ meta, yml ]

    st_svg_html            = ST_SPATIALLY_VARIABLE_GENES.out.html        // channel: [ meta, html ]
    st_output              = ST_SPATIALLY_VARIABLE_GENES.out.artifacts   // channel: [ meta, csv ]
    st_svg_notebook        = ST_SPATIALLY_VARIABLE_GENES.out.notebook    // channel: [ meta, qmd ]
    st_svg_params          = ST_SPATIALLY_VARIABLE_GENES.out.params_yaml // channel: [ meta, yml ]

    versions               = ch_versions                                 // channel: [ versions.yml ]
}
