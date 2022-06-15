//
// Single cell data pre-processing
//
process SC_PREPROCESS {

    label "python_process"

    input:
    tuple val(sample_id), path(sc_raw), path(sc_factors)
    path(mito_data)

    output:
    tuple val(sample_id), path("*_norm.h5ad")       , emit: sc_data_norm
    tuple val(sample_id), path("*_plain.h5ad")      , emit: sc_data_plain
    tuple val(sample_id), path("*.sc_adata_x.npz")  , emit: sc_adata_x
    tuple val(sample_id), path("*.sc_adata_var.npz"), emit: sc_adata_var
    tuple val(sample_id), path("*.sc_adata_obs.npz"), emit: sc_adata_obs
    tuple val(sample_id), path("*.png")             , emit: figures

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
        --nameDataPlain=${sample_id}.sc_adata_plain.h5ad \
        --nameDataNorm=${sample_id}.sc_adata_norm.h5ad \
        --nameX ${sample_id}.sc_adata_x.npz \
        --nameVar ${sample_id}.sc_adata_var.npz \
        --nameObs ${sample_id}.sc_adata_obs.npz
    """
}
