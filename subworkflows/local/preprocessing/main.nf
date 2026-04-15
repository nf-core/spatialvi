include { SCANPY_CALCULATE_QC_METRICS  } from "../../../modules/local/scanpy/calculate_qc_metrics"
include { SCANPY_FILTER                } from "../../../modules/local/scanpy/filter"
include { SCANPY_HIGHLY_VARIABLE_GENES } from "../../../modules/local/scanpy/highly_variable_genes"
include { SCANPY_LOG1P                 } from "../../../modules/local/scanpy/log1p"
include { SCANPY_NORMALIZE_TOTAL       } from "../../../modules/local/scanpy/normalize_total"
include { SCANPY_PCA                   } from "../../../modules/local/scanpy/pca"

workflow PREPROCESSING {

    take:
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

    main:

    //
    // MODULE: Calculate quality control metrics
    //
    SCANPY_CALCULATE_QC_METRICS (
        ch_adata_input
    )

    //
    // MODULE: Filtering
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

    //
    // MODULE: Normalization
    //
    SCANPY_NORMALIZE_TOTAL (
        SCANPY_FILTER.out.adata,
        normalize_target_sum
    )

    //
    // MODULE: Log-transformation
    //
    SCANPY_LOG1P (
        SCANPY_NORMALIZE_TOTAL.out.adata
    )

    //
    // MODULE: Highly variable gene selection
    //
    SCANPY_HIGHLY_VARIABLE_GENES (
        SCANPY_LOG1P.out.adata,
        n_highly_variable_genes,
        hvg_flavor
    )

    //
    // MODULE: Principal Component Analysis
    //
    SCANPY_PCA (
        SCANPY_HIGHLY_VARIABLE_GENES.out.adata,
        n_principal_components,
        pca_use_highly_variable
    )


    emit:
    adata          = SCANPY_PCA.out.adata          // channel: [ meta, h5ad ]
    filter_stats   = SCANPY_FILTER.out.stats       // channel: [ meta, csv ]
}
