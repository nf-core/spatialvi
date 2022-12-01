//
// Clustering etc. TODO: better description
//
process ST_CLUSTERING {

    label "python_process"

    // TODO: Add Conda/container directive
    container "erikfas/spatialtranscriptomics"

    input:
    tuple val(sample_id), path(st_adata_norm), path(sc_adata_norm)

    output:
    tuple val(sample), path("*.st_*.h5ad"), emit: st_adata_processed
    tuple val(sample), path("*.sc_*.h5ad"), emit: sc_adata_processed
    path("*.png")                         , emit: figures, optional: true

    script:
    """
    stClusteringWorkflow.py \
        --fileNameST ${st_adata_norm} \
        --fileNameSC ${sc_adata_norm} \
        --resolution=$params.Clustering_resolution \
        --saveFileST ${sample_id}.st_adata_processed.h5ad \
        --saveFileSC ${sample_id}.sc_adata_processed.h5ad
    """
}
