
// Subworkflow for downstream analyses of ST data
//

include { QUARTONOTEBOOK as REPORT     } from '../../../modules/nf-core/quartonotebook/main'
include { SCANPY_CALCULATE_QC_METRICS  } from '../../../modules/local/scanpy/calculate_qc_metrics/main'
include { SCANPY_FILTER                } from '../../../modules/local/scanpy/filter/main'
include { SCANPY_HIGHLY_VARIABLE_GENES } from '../../../modules/local/scanpy/highly_variable_genes/main'
include { SCANPY_LEIDEN                } from '../../../modules/local/scanpy/leiden/main'
include { SCANPY_LOG1P                 } from '../../../modules/local/scanpy/log1p/main'
include { SCANPY_NEIGHBORS             } from '../../../modules/local/scanpy/neighbors/main'
include { SCANPY_NORMALIZE_TOTAL       } from '../../../modules/local/scanpy/normalize_total/main'
include { SCANPY_PCA                   } from '../../../modules/local/scanpy/pca/main'
include { SCANPY_RANK_GENES_GROUPS     } from '../../../modules/local/scanpy/rank_genes_groups/main'
include { SCANPY_UMAP                  } from '../../../modules/local/scanpy/umap/main'
include { SDATA_UPDATE_TABLE           } from '../../../modules/local/sdata/update_table/main'
include { SQUIDPY_INTERACTION_MATRIX   } from '../../../modules/local/squidpy/interaction_matrix/main'
include { SQUIDPY_NHOOD_ENRICHMENT     } from '../../../modules/local/squidpy/nhood_enrichment/main'
include { SQUIDPY_SPATIAL_AUTOCORR     } from '../../../modules/local/squidpy/spatial_autocorr/main'
include { SQUIDPY_SPATIAL_NEIGHBORS    } from '../../../modules/local/squidpy/spatial_neighbors/main'

workflow DOWNSTREAM {

    take:
    ch_sdata_input          // channel: [ meta, zarr ]
    ch_adata_input          // channel: [ meta, h5ad ]
    qc_min_counts           // integer: Minimum UMIs per spot
    qc_min_genes            // integer: Minimum genes per spot
    qc_min_spots            // integer: Minimum spots per gene
    qc_mito_threshold       //   float: Maximum mitochondrial content per spot
    qc_ribo_threshold       //   float: Minimum ribosomal content per spot
    qc_hb_threshold         //   float: Maximum haemoglobin content per spot
    normalize_target_sum    //  string: Target sum of total count normalization
    n_highly_variable_genes // integer: Number of highly variable genes to use
    hvg_flavor              //  string: Flavor for HVG calculations
    n_principal_components  // integer: Number of principal components to compute
    pca_use_highly_variable // boolean: Whether to only use highly variable genes for PCA
    n_neighbours            // integer: Number of nearest neighbours to compute
    neighbours_n_pcs        // integer: Number of PCs to use for nearest neighbours
    umap_min_dist           //   float: Minimum distance between embedded points
    umap_spread             //   float: Scale of embedded points
    cluster_resolution      //   float: Spot clustering resolution
    rank_genes_method       //  string: Method to use for differential expression testing
    spatial_coord_type      //  string: Type of spatial coordinate system
    spatial_n_neighbours    // integer: Number of spatial neighbours to use
    svg_autocorr_method     //  string: Spatial autocorrelation method ('moran' or 'geary')
    n_top_svgs              // integer: Number of top spatially variable genes to report

    main:

    // =========================================================================
    // QUALITY CONTROL AND FILTERING
    // =========================================================================

    //
    // Quality control metrics
    //
    SCANPY_CALCULATE_QC_METRICS (
        ch_adata_input
    )

    //
    // Filtering
    //
    SCANPY_FILTER (
        SCANPY_CALCULATE_QC_METRICS.out.adata,
        qc_min_counts,
        qc_min_genes,
        qc_min_spots,
        qc_mito_threshold,
        qc_ribo_threshold,
        qc_hb_threshold
    )

    // =========================================================================
    // NORMALIZATION AND FEATURE SELECTION
    // =========================================================================

    //
    // Normalization
    //
    SCANPY_NORMALIZE_TOTAL (
        SCANPY_FILTER.out.adata,
        normalize_target_sum
    )

    //
    // Log transformation
    //
    SCANPY_LOG1P (
        SCANPY_NORMALIZE_TOTAL.out.adata
    )

    //
    // Highly variable genes
    //
    SCANPY_HIGHLY_VARIABLE_GENES (
        SCANPY_LOG1P.out.adata,
        n_highly_variable_genes,
        hvg_flavor
    )

    // =========================================================================
    // DIMENSIONALITY REDUCTION AND CLUSTERING
    // =========================================================================

    //
    // PCA
    //
    SCANPY_PCA (
        SCANPY_HIGHLY_VARIABLE_GENES.out.adata,
        n_principal_components,
        pca_use_highly_variable
    )

    //
    // Neighbors graph (for UMAP and Leiden)
    //
    use_rep = ''
    SCANPY_NEIGHBORS (
        SCANPY_PCA.out.adata,
        n_neighbours,
        neighbours_n_pcs,
        use_rep
    )

    //
    // UMAP
    //
    umap_key_added = 'X_umap'
    SCANPY_UMAP (
        SCANPY_NEIGHBORS.out.adata,
        umap_min_dist,
        umap_spread,
        umap_key_added
    )

    //
    // Leiden clustering
    //
    leiden_key_added = 'clusters'
    SCANPY_LEIDEN (
        SCANPY_UMAP.out.adata,
        cluster_resolution,
        leiden_key_added
    )

    // =========================================================================
    // DIFFERENTIAL EXPRESSION AND SPATIAL ANALYSIS
    // =========================================================================

    //
    // Spatial differential expression analysis (rank genes by cluster)
    //
    rank_genes_group_by = 'clusters'
    SCANPY_RANK_GENES_GROUPS (
        SCANPY_LEIDEN.out.adata,
        rank_genes_group_by,
        rank_genes_method
    )

    //
    // Spatial neighbors
    //
    SQUIDPY_SPATIAL_NEIGHBORS (
        SCANPY_RANK_GENES_GROUPS.out.adata,
        spatial_coord_type,
        spatial_n_neighbours
    )

    //
    // Spatial neighbourhood enrichment analysis
    //
    cluster_key = 'clusters'
    SQUIDPY_NHOOD_ENRICHMENT (
        SQUIDPY_SPATIAL_NEIGHBORS.out.adata,
        cluster_key
    )

    //
    // Spatial interaction matrix
    //
    SQUIDPY_INTERACTION_MATRIX (
        SQUIDPY_NHOOD_ENRICHMENT.out.adata,
        cluster_key
    )

    //
    // Spatial autocorrelation (spatially variable genes)
    //
    SQUIDPY_SPATIAL_AUTOCORR (
        SQUIDPY_INTERACTION_MATRIX.out.adata,
        svg_autocorr_method
    )

    //
    // Update SpatialData with SVG results (final checkpoint)
    //
    SDATA_UPDATE_TABLE (
        ch_sdata_input.join(SQUIDPY_SPATIAL_AUTOCORR.out.adata),
        ''
    )
    ch_sdata_output = SDATA_UPDATE_TABLE.out.sdata

    // =========================================================================
    // REPORT
    // =========================================================================

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

    emit:
    // Final data outputs
    sdata              = ch_sdata_output                    // channel: [ meta, zarr ]
    adata              = SQUIDPY_SPATIAL_AUTOCORR.out.adata // channel: [ meta, h5ad ]

    // Intermediate AnnData outputs (useful for debugging)
    adata_qc           = SCANPY_FILTER.out.adata            // channel: [ meta, h5ad ]
    adata_clustered    = SCANPY_LEIDEN.out.adata            // channel: [ meta, h5ad ]

    // Filter statistics (for MultiQC)
    filter_stats       = SCANPY_FILTER.out.stats            // channel: [ meta, json ]

    // SVG results
    svg_csv            = SQUIDPY_SPATIAL_AUTOCORR.out.csv   // channel: [ meta, csv ]

    // Report outputs
    report_html        = REPORT.out.html                    // channel: [ meta, html ]
    report_notebook    = REPORT.out.notebook                // channel: [ meta, qmd ]
    report_params_yaml = REPORT.out.params_yaml             // channel: [ meta, yml ]
    report_artifacts   = REPORT.out.artifacts               // channel: [ meta, artifacts ]
}
