//
// Spatial data pre-processing
//
process ST_PREPROCESS {

    // TODO: Add final Conda/container directive
    // TODO: Export versions

    tag "${sample_id}"
    label "process_low"

    container "erikfas/spatialtranscriptomics"

    input:
    tuple val(sample_id), path(st_raw), path(st_factors)
    path(mito_data)

    output:
    tuple val(sample_id), path("*_norm.h5ad")       , emit: st_data_norm
    tuple val(sample_id), path("*_plain.h5ad")      , emit: st_data_plain
    tuple val(sample_id), path("*.st_adata_x.npz")  , emit: st_adata_x
    tuple val(sample_id), path("*.st_adata_var.npz"), emit: st_adata_var
    tuple val(sample_id), path("*.st_adata_obs.npz"), emit: st_adata_obs
    tuple val(sample_id), path("*.png")             , emit: figures
    // path("versions.yml")                            , emit: versions

    script:
    """
    stPreprocess.py \
        --npFactorsOutputName ${st_factors} \
        --rawAdata ${st_raw} \
        --mitoFile ${mito_data} \
        --pltFigSize ${params.STpreprocess_pltFigSize} \
        --minCounts ${params.STpreprocess_minCounts} \
        --minGenes ${params.STpreprocess_minGenes} \
        --minCells ${params.STpreprocess_minCells} \
        --histplotQCmaxTotalCounts ${params.STpreprocess_histplotQCmaxTotalCounts} \
        --histplotQCminGeneCounts ${params.STpreprocess_histplotQCminGeneCounts} \
        --histplotQCbins ${params.STpreprocess_histplotQCbins} \
        --nameDataPlain ${sample_id}.st_adata_plain.h5ad \
        --nameDataNorm ${sample_id}.st_adata_norm.h5ad \
        --nameX ${sample_id}.st_adata_x.npz \
        --nameVar ${sample_id}.st_adata_var.npz \
        --nameObs ${sample_id}.st_adata_obs.npz
    """
}
