//
// Spatial data pre-processing
//
process ST_QC_AND_NORMALISATION {

    // TODO: Add a better description
    // TODO: Find solution for Quarto with Conda

    tag "${meta.id}"
    label "process_low"

    container "docker.io/erikfas/spatialtranscriptomics"

    input:
    path(report)
    tuple val(meta), path(st_raw, stageAs: "adata_raw.h5ad")

    output:
    tuple val(meta), path("st_adata_norm.h5ad")           , emit: st_data_norm
    tuple val(meta), path("st_adata_plain.h5ad")          , emit: st_data_plain
    tuple val(meta), path("st_qc_and_normalisation.html") , emit: html
    tuple val(meta), path("st_qc_and_normalisation_files"), emit: html_files
    path("versions.yml")                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    quarto render ${report} \
        --output st_qc_and_normalisation.html \
        -P rawAdata:${st_raw} \
        -P pltFigSize:${params.st_preprocess_fig_size} \
        -P minCounts:${params.st_preprocess_min_counts} \
        -P minGenes:${params.st_preprocess_min_genes} \
        -P minCells:${params.st_preprocess_min_cells} \
        -P histplotQCmaxTotalCounts:${params.st_preprocess_hist_qc_max_total_counts} \
        -P histplotQCminGeneCounts:${params.st_preprocess_hist_qc_min_gene_counts} \
        -P histplotQCbins:${params.st_preprocess_hist_qc_bins} \
        -P nameDataPlain:st_adata_plain.h5ad \
        -P nameDataNorm:st_adata_norm.h5ad

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quarto: \$(quarto -v)
        scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
    END_VERSIONS
    """
}
