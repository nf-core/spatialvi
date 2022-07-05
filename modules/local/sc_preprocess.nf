//
// Single cell data pre-processing
//
process SC_PREPROCESS {

    tag "${meta}"
    label "process_low"

    // TODO: Add Conda/container directive

    input:
    tuple val(meta), path(sc_raw), path(sc_factors)
    path(mito_data)

    output:
    tuple val(meta), path("*_norm.h5ad")       , emit: sc_data_norm
    tuple val(meta), path("*_plain.h5ad")      , emit: sc_data_plain
    tuple val(meta), path("*.sc_adata_x.npz")  , emit: sc_adata_x
    tuple val(meta), path("*.sc_adata_var.npz"), emit: sc_adata_var
    tuple val(meta), path("*.sc_adata_obs.npz"), emit: sc_adata_obs
    tuple val(meta), path("*.png")             , emit: figures

    script:
    """
    scPreprocess.py \
        --npFactorsOutputName=${sc_factors} \
        --rawAdata=${sc_raw} \
        --mitoFile=${mito_data} \
        --pltFigSize=${params.SCpreprocess_pltFigSize} \
        --minCounts=${params.SCpreprocess_minCounts} \
        --minGenes=${params.SCpreprocess_minGenes} \
        --minCells=${params.SCpreprocess_minCells} \
        --histplotQCmaxTotalCounts=${params.SCpreprocess_histplotQCmaxTotalCounts} \
        --histplotQCminGeneCounts=${params.SCpreprocess_histplotQCminGeneCounts} \
        --histplotQCbins=${params.SCpreprocess_histplotQCbins} \
        --nameDataPlain=${meta.id}.sc_adata_plain.h5ad \
        --nameDataNorm=${meta.id}.sc_adata_norm.h5ad \
        --nameX ${meta.id}.sc_adata_x.npz \
        --nameVar ${meta.id}.sc_adata_var.npz \
        --nameObs ${meta.id}.sc_adata_obs.npz
    """
}
