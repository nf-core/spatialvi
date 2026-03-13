//
// Subworkflow for downstream analyses of ST data
//

include { SDATA_TO_LEGACY_ANNDATA                        } from '../../../modules/local/sdata/to_legacy_anndata/main'
include { SCANPY_CALCULATE_QC_METRICS                    } from '../../../modules/local/scanpy/calculate_qc_metrics/main'
include { SCANPY_FILTER                                  } from '../../../modules/local/scanpy/filter/main'
include { SCANPY_NORMALIZE_TOTAL                         } from '../../../modules/local/scanpy/normalize_total/main'
include { SCANPY_LOG1P                                   } from '../../../modules/local/scanpy/log1p/main'
include { SCANPY_HIGHLY_VARIABLE_GENES                   } from '../../../modules/local/scanpy/highly_variable_genes/main'
include { SCANPY_PCA                                     } from '../../../modules/local/scanpy/pca/main'
include { SCANPY_NEIGHBORS                               } from '../../../modules/local/scanpy/neighbors/main'
include { SCANPY_UMAP                                    } from '../../../modules/local/scanpy/umap/main'
include { SCANPY_LEIDEN                                  } from '../../../modules/local/scanpy/leiden/main'
include { SCANPY_RANK_GENES_GROUPS                       } from '../../../modules/local/scanpy/rank_genes_groups/main'
include { SQUIDPY_SPATIAL_NEIGHBORS                      } from '../../../modules/local/squidpy/spatial_neighbors/main'
include { SQUIDPY_NHOOD_ENRICHMENT                       } from '../../../modules/local/squidpy/nhood_enrichment/main'
include { SQUIDPY_INTERACTION_MATRIX                     } from '../../../modules/local/squidpy/interaction_matrix/main'
include { SQUIDPY_SPATIAL_AUTOCORR                       } from '../../../modules/local/squidpy/spatial_autocorr/main'
include { SDATA_UPDATE_TABLE as SDATA_UPDATE_TABLE_QC         } from '../../../modules/local/sdata/update_table/main'
include { SDATA_UPDATE_TABLE as SDATA_UPDATE_TABLE_CLUSTERING } from '../../../modules/local/sdata/update_table/main'
include { SDATA_UPDATE_TABLE as SDATA_UPDATE_TABLE_SVG        } from '../../../modules/local/sdata/update_table/main'

include { QUARTONOTEBOOK as REPORT                       } from '../../../modules/nf-core/quartonotebook/main'

workflow DOWNSTREAM {

    take:
    ch_sdata_raw        // channel: [ meta, sdata.zarr ]
    qc_min_counts       // integer: Minimum UMIs per spot
    qc_min_genes        // integer: Minimum genes per spot
    qc_min_spots        // integer: Minimum spots per gene
    qc_mito_threshold   // float  : Maximum mitochondrial content per spot
    qc_ribo_threshold   // float  : Minimum ribosomal content per spot
    qc_hb_threshold     // float  : Maximum haemoglobin content per spot
    cluster_n_hvgs      // integer: Number of highly variable genes to use
    cluster_resolution  // float  : Spot clustering resolution
    svg_autocorr_method // string : Spatial autocorrelation method ('moran' or 'geary')
    n_top_svgs          // integer: Number of top spatially variable genes to report

    main:

    //
    // Quarto report template and extensions
    //
    report_notebook = file("${projectDir}/bin/report.qmd", checkIfExists: true)
    extensions = channel.fromPath("${projectDir}/assets/_extensions").collect()

    // =========================================================================
    // DATA LOADING
    // =========================================================================

    //
    // Extract legacy AnnData for scanpy processing
    //
    SDATA_TO_LEGACY_ANNDATA(ch_sdata_raw)

    // =========================================================================
    // QUALITY CONTROL AND FILTERING
    // =========================================================================

    //
    // Quality control metrics
    //
    SCANPY_CALCULATE_QC_METRICS(SDATA_TO_LEGACY_ANNDATA.out.adata)

    //
    // Filtering
    //
    SCANPY_FILTER(
        SCANPY_CALCULATE_QC_METRICS.out.adata,
        qc_min_counts,
        qc_min_genes,
        qc_min_spots,
        qc_mito_threshold,
        qc_ribo_threshold,
        qc_hb_threshold
    )

    //
    // Update SpatialData with filtered AnnData (checkpoint for QC report)
    //
    ch_sdata_for_qc_update = ch_sdata_raw
        .join(SCANPY_FILTER.out.adata)

    SDATA_UPDATE_TABLE_QC(ch_sdata_for_qc_update)

    ch_qc_sdata = SDATA_UPDATE_TABLE_QC.out.sdata

    // =========================================================================
    // NORMALIZATION AND FEATURE SELECTION
    // =========================================================================

    //
    // Normalization
    //
    SCANPY_NORMALIZE_TOTAL(SCANPY_FILTER.out.adata)

    //
    // Log transformation
    //
    SCANPY_LOG1P(SCANPY_NORMALIZE_TOTAL.out.adata)

    //
    // Highly variable genes
    //
    SCANPY_HIGHLY_VARIABLE_GENES(
        SCANPY_LOG1P.out.adata,
        cluster_n_hvgs
    )

    // =========================================================================
    // DIMENSIONALITY REDUCTION AND CLUSTERING
    // =========================================================================

    //
    // PCA
    //
    SCANPY_PCA(SCANPY_HIGHLY_VARIABLE_GENES.out.adata)

    //
    // Neighbors graph (for UMAP and Leiden)
    //
    SCANPY_NEIGHBORS(SCANPY_PCA.out.adata)

    //
    // UMAP
    //
    SCANPY_UMAP(SCANPY_NEIGHBORS.out.adata)

    //
    // Leiden clustering
    //
    SCANPY_LEIDEN(
        SCANPY_UMAP.out.adata,
        cluster_resolution
    )

    //
    // Update SpatialData with clustered AnnData
    //
    ch_sdata_for_clustering_update = ch_qc_sdata
        .join(SCANPY_LEIDEN.out.adata)

    SDATA_UPDATE_TABLE_CLUSTERING(ch_sdata_for_clustering_update)

    ch_clustering_sdata = SDATA_UPDATE_TABLE_CLUSTERING.out.sdata

    // =========================================================================
    // DIFFERENTIAL EXPRESSION AND SPATIAL ANALYSIS
    // =========================================================================

    //
    // Differential expression analysis (rank genes by cluster)
    //
    SCANPY_RANK_GENES_GROUPS(SCANPY_LEIDEN.out.adata)

    //
    // Spatial neighbors (for spatial analyses)
    //
    SQUIDPY_SPATIAL_NEIGHBORS(SCANPY_RANK_GENES_GROUPS.out.adata)

    //
    // Neighborhood enrichment analysis
    //
    SQUIDPY_NHOOD_ENRICHMENT(SQUIDPY_SPATIAL_NEIGHBORS.out.adata)

    //
    // Interaction matrix
    //
    SQUIDPY_INTERACTION_MATRIX(SQUIDPY_NHOOD_ENRICHMENT.out.adata)

    //
    // Spatial autocorrelation (spatially variable genes)
    //
    SQUIDPY_SPATIAL_AUTOCORR(
        SQUIDPY_INTERACTION_MATRIX.out.adata,
        svg_autocorr_method
    )

    //
    // Update SpatialData with SVG results (final checkpoint)
    //
    ch_sdata_for_svg_update = ch_clustering_sdata
        .join(SQUIDPY_SPATIAL_AUTOCORR.out.adata)

    SDATA_UPDATE_TABLE_SVG(ch_sdata_for_svg_update)

    ch_svg_sdata = SDATA_UPDATE_TABLE_SVG.out.sdata

    // =========================================================================
    // REPORT
    // =========================================================================

    ch_report_input_data = ch_svg_sdata
        .map { it -> it[1] }
    ch_report_input_data.view { "DEBUG ch_report_input_data: $it" }
    
    ch_report_notebook = ch_svg_sdata
        .map { it -> it[0] }
        .combine(channel.value(report_notebook))
        .map { meta, notebook -> tuple(meta, notebook) }
    
    ch_report_notebook.view { "DEBUG ch_report_notebook: $it" }
    
    ch_report_params = ch_svg_sdata
        .map { meta, sdata -> 
            [
                input_sdata  : sdata.name,
                artifact_dir : "artifacts"
            ]
        }

    REPORT(
        ch_report_notebook,
        ch_report_params,
        ch_report_input_data,
        extensions
    )

    ch_report_html      = REPORT.out.html
    ch_report_nb        = REPORT.out.notebook
    ch_report_yml       = REPORT.out.params_yaml
    ch_report_artifacts = REPORT.out.artifacts

    emit:
    // SpatialData outputs
    sdata_qc              = ch_qc_sdata                              // channel: [ meta, zarr ]
    sdata_clustered       = ch_clustering_sdata                      // channel: [ meta, zarr ]
    sdata_svg             = ch_svg_sdata                             // channel: [ meta, zarr ]

    // AnnData outputs (intermediate, useful for debugging)
    adata_qc              = SCANPY_FILTER.out.adata                  // channel: [ meta, h5ad ]
    adata_clustered       = SCANPY_LEIDEN.out.adata                  // channel: [ meta, h5ad ]
    adata_svg             = SQUIDPY_SPATIAL_AUTOCORR.out.adata       // channel: [ meta, h5ad ]

    // Filter statistics (for MultiQC)
    filter_stats          = SCANPY_FILTER.out.stats                  // channel: [ meta, json ]

    // SVG results
    svg_csv               = SQUIDPY_SPATIAL_AUTOCORR.out.csv         // channel: [ meta, csv ]

    // Report outputs
    report_html           = ch_report_html                           // channel: [ meta, html ]
    report_notebook       = ch_report_nb                             // channel: [ meta, qmd ]
    report_params_yaml    = ch_report_yml                            // channel: [ meta, yml ]
    report_artifacts      = ch_report_artifacts                      // channel: [ meta, artifacts ]
}