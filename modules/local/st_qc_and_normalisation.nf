//
// Spatial data pre-processing
//
process ST_QC_AND_NORMALISATION {

    // TODO: Add final Conda/container directive
    // TODO: Export versions

    tag "${meta.id}"
    label "process_low"

    container "cavenel/spatialtranscriptomics"

    input:
    path(report)
    tuple val(meta), path(st_raw, stageAs: "adata_raw.h5ad")
    path(mito_data)

    output:
    tuple val(meta), path("*/st_adata_norm.h5ad")           , emit: st_data_norm
    tuple val(meta), path("*/st_adata_plain.h5ad")          , emit: st_data_plain
    tuple val(meta), path("*/st_qc_and_normalisation.html") , emit: html
    // path("versions.yml")                                    , emit: versions

    script:
    """
    quarto render ${report} \
        --output st_qc_and_normalisation.html \
        -P rawAdata:${st_raw} \
        -P mitoFile:${mito_data} \
        -P pltFigSize:${params.st_preprocess_fig_size} \
        -P minCounts:${params.st_preprocess_min_counts} \
        -P minGenes:${params.st_preprocess_min_genes} \
        -P minCells:${params.st_preprocess_min_cells} \
        -P histplotQCmaxTotalCounts:${params.st_preprocess_hist_qc_max_total_counts} \
        -P histplotQCminGeneCounts:${params.st_preprocess_hist_qc_min_gene_counts} \
        -P histplotQCbins:${params.st_preprocess_hist_qc_bins} \
        -P nameDataPlain:st_adata_plain.h5ad \
        -P nameDataNorm:st_adata_norm.h5ad

    mkdir -p ${meta.id}
    mv st_adata_plain.h5ad ${meta.id}/st_adata_plain.h5ad
    mv st_adata_norm.h5ad ${meta.id}/st_adata_norm.h5ad
    mv st_qc_and_normalisation.html ${meta.id}/st_qc_and_normalisation.html
    """
}
