include { SCANPY_HIGHLY_VARIABLE_GENES } from "../../../modules/local/scanpy/highly_variable_genes"
include { SCANPY_LOG1P                 } from "../../../modules/local/scanpy/log1p"
include { SCANPY_NORMALIZE_TOTAL       } from "../../../modules/local/scanpy/normalize_total"
include { SCANPY_PCA                   } from "../../../modules/local/scanpy/pca"

workflow PREPROCESSING {

    take:
    ch_adata_input          // channel: [ meta, h5ad ]
    normalize_target_sum    //  string: Target sum of total count normalization
    n_highly_variable_genes // integer: Number of highly variable genes to use
    hvg_flavor              //  string: Flavor for HVG calculations
    hvg_batch_key           //  string: Column in `obs` to select HVGs per batch, or ''
    n_principal_components  // integer: Number of principal components to compute
    pca_use_highly_variable // boolean: Whether to only use highly variable genes for PCA
    pca_key_added           //  string: Key for the PCA results

    main:

    //
    // MODULE: Normalization
    //
    SCANPY_NORMALIZE_TOTAL (
        ch_adata_input,
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
        hvg_flavor,
        hvg_batch_key,
        true // write_adata
    )

    //
    // MODULE: Principal Component Analysis
    //
    SCANPY_PCA (
        SCANPY_HIGHLY_VARIABLE_GENES.out.adata,
        n_principal_components,
        pca_use_highly_variable,
        pca_key_added,
        true // write_adata
    )

    emit:
    adata = SCANPY_PCA.out.adata // channel: [ meta, h5ad ]
}
