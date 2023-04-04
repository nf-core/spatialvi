//
// Spatial data pre-processing
//
process ST_QC_AND_NORMALISATION {

    // TODO: Add final Conda/container directive
    // TODO: Export versions

    tag "${sample_id}"
    label "process_low"

    container "cavenel/spatialtranscriptomics"

    input:
    path(report)
    tuple val(sample_id), path(st_raw, stageAs: "adata_raw.h5ad")
    path(mito_data)

    output:
    tuple val(sample_id), path("*.st_adata_norm.h5ad")           , emit: st_data_norm
    tuple val(sample_id), path("*.st_adata_plain.h5ad")          , emit: st_data_plain
    tuple val(sample_id), path("*.st_qc_and_normalisation.html") , emit: html
    tuple val(sample_id), path("st_qc_and_normalisation_files/*"), emit: html_files
    // path("versions.yml")                                         , emit: versions

    script:
    """
    quarto render ${report} \
        --output ${sample_id}.st_qc_and_normalisation.html \
        -P rawAdata:${st_raw} \
        -P mitoFile:${mito_data} \
        -P pltFigSize:${params.STpreprocess_pltFigSize} \
        -P minCounts:${params.STpreprocess_minCounts} \
        -P minGenes:${params.STpreprocess_minGenes} \
        -P minCells:${params.STpreprocess_minCells} \
        -P histplotQCmaxTotalCounts:${params.STpreprocess_histplotQCmaxTotalCounts} \
        -P histplotQCminGeneCounts:${params.STpreprocess_histplotQCminGeneCounts} \
        -P histplotQCbins:${params.STpreprocess_histplotQCbins} \
        -P nameDataPlain:st_adata_plain.h5ad \
        -P nameDataNorm:st_adata_norm.h5ad

    mv st_adata_plain.h5ad ${sample_id}.st_adata_plain.h5ad
    mv st_adata_norm.h5ad ${sample_id}.st_adata_norm.h5ad
    """
}
