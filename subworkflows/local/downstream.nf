//
// Subworkflow for downstream analyses of ST data
//

include { QUARTONOTEBOOK as QUALITY_CONTROLS         } from '../../modules/nf-core/quartonotebook/main'
include { QUARTONOTEBOOK as SPATIALLY_VARIABLE_GENES } from '../../modules/nf-core/quartonotebook/main'
include { QUARTONOTEBOOK as CLUSTERING               } from '../../modules/nf-core/quartonotebook/main'

workflow DOWNSTREAM {

    take:
    sdata_raw

    main:

    ch_versions = Channel.empty()

    //
    // Quarto reports and extension files
    //
    quality_controls_notebook = file("${projectDir}/bin/quality_controls.qmd", checkIfExists: true)
    clustering_notebook = file("${projectDir}/bin/clustering.qmd", checkIfExists: true)
    spatially_variable_genes_notebook = file("${projectDir}/bin/spatially_variable_genes.qmd", checkIfExists: true)
    extensions = Channel.fromPath("${projectDir}/assets/_extensions").collect()

    //
    // Quality controls and filtering
    //
    ch_quality_controls_input_data = sdata_raw
        .map { it -> it[1] }
    ch_quality_controls_notebook = sdata_raw
        .map { tuple(it[0], quality_controls_notebook) }
    quality_controls_params = [
        input_sdata: "sdata_raw.zarr",
        min_counts: params.qc_min_counts,
        min_genes: params.qc_min_genes,
        min_spots: params.qc_min_spots,
        mito_threshold: params.qc_mito_threshold,
        ribo_threshold: params.qc_ribo_threshold,
        hb_threshold: params.qc_hb_threshold,
        artifact_dir: "artifacts",
        output_adata: "adata_filtered.h5ad",
        output_sdata: "sdata_filtered.zarr",
    ]
    QUALITY_CONTROLS (
        ch_quality_controls_notebook,
        quality_controls_params,
        ch_quality_controls_input_data,
        extensions
    )
    ch_versions = ch_versions.mix(QUALITY_CONTROLS.out.versions)
    ch_qc = QUALITY_CONTROLS.out.artifacts
        | map { meta, artifacts -> [meta, artifacts[0], meta, artifacts[1]] }
        | flatten
        | collate ( 2 )
        | branch {
            sdata: it[1].name.endsWith('.zarr')
            mqc: it[1].name.endsWith('.csv')
        }
    ch_qc_sdata = ch_qc.sdata
    ch_qc_mqc   = ch_qc.mqc
    ch_qc_html  = QUALITY_CONTROLS.out.html
    ch_qc_nb    = QUALITY_CONTROLS.out.notebook
    ch_qc_yml   = QUALITY_CONTROLS.out.params_yaml

    //
    // Normalisation, dimensionality reduction and clustering
    //
    ch_clustering_input_data = QUALITY_CONTROLS.out.artifacts
        .map { it -> it[1] }
    ch_clustering_notebook = QUALITY_CONTROLS.out.artifacts
        .map { tuple(it[0], clustering_notebook) }
    clustering_params = [
        input_sdata: "sdata_filtered.zarr",
        cluster_resolution: params.cluster_resolution,
        n_hvgs: params.cluster_n_hvgs,
        artifact_dir: "artifacts",
        output_adata: "adata_processed.h5ad",
        output_sdata: "sdata_processed.zarr",
    ]
    CLUSTERING (
        ch_clustering_notebook,
        clustering_params,
        ch_clustering_input_data,
        extensions
    )
    ch_versions = ch_versions.mix(CLUSTERING.out.versions)

    //
    // Spatially variable genes
    //
    ch_spatially_variable_genes_input_data = CLUSTERING.out.artifacts
        .map { it -> it[1] }
    ch_spatially_variable_genes_notebook = CLUSTERING.out.artifacts
        .map { tuple(it[0], spatially_variable_genes_notebook) }
    spatially_variable_genes_params = [
        input_sdata: "sdata_processed.zarr",
        svg_autocorr_method: params.svg_autocorr_method,
        n_top_svgs: params.n_top_svgs,
        artifact_dir: "artifacts",
        output_csv: "spatially_variable_genes.csv",
        output_adata: "adata_spatially_variable_genes.h5ad",
        output_sdata: "sdata.zarr",
    ]
    SPATIALLY_VARIABLE_GENES (
        ch_spatially_variable_genes_notebook,
        spatially_variable_genes_params,
        ch_spatially_variable_genes_input_data,
        extensions
    )
    ch_versions = ch_versions.mix(SPATIALLY_VARIABLE_GENES.out.versions)

    emit:
    qc_html           = ch_qc_html  // channel: [ meta, html ]
    qc_sdata          = ch_qc_sdata // channel: [ meta, zarr ]
    qc_mqc            = ch_qc_mqc   // channel: [ meta, csv ]
    qc_nb             = ch_qc_nb    // channel: [ meta, qmd ]
    qc_params         = ch_qc_yml   // channel: [ meta, yml ]

    clustering_html   = CLUSTERING.out.html                      // channel: [ html ]
    clustering_sdata  = CLUSTERING.out.artifacts                 // channel: [ meta, h5ad]
    clustering_nb     = CLUSTERING.out.notebook                  // channel: [ meta, qmd ]
    clustering_params = CLUSTERING.out.params_yaml               // channel: [ meta, yml ]

    svg_html          = SPATIALLY_VARIABLE_GENES.out.html        // channel: [ meta, html ]
    svg_csv           = SPATIALLY_VARIABLE_GENES.out.artifacts   // channel: [ meta, csv ]
    svg_nb            = SPATIALLY_VARIABLE_GENES.out.notebook    // channel: [ meta, qmd ]
    svg_params        = SPATIALLY_VARIABLE_GENES.out.params_yaml // channel: [ meta, yml ]

    versions          = ch_versions                              // channel: [ versions.yml ]
}
