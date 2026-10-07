process SCANPY_PCA {
    tag "${meta.id}"
    label 'process_medium'

    container "community.wave.seqera.io/library/harmonypy_scanorama_gcc_gxx_pruned:95f731fde0b9ddef"

    input:
    tuple val(meta), path(adata, stageAs: "input.h5ad", arity: '1')
    val n_pcs
    val use_highly_variable
    val key_added
    val write_adata

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: adata, optional: true
    tuple val(meta), path("obsm/*.pkl")    , emit: obsm
    tuple val(meta), path("varm/*.pkl")    , emit: varm
    tuple val(meta), path("uns/*.pkl")     , emit: uns
    path "versions.yml"                    , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'pca.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def touch_adata = write_adata.toString() == 'true' ? "touch ${prefix}.h5ad" : ''
    def varm_key = key_added == 'X_pca' ? 'PCs'   : key_added
    def uns_key  = key_added == 'X_pca' ? 'pca'   : key_added
    """
    ${touch_adata}
    mkdir -p obsm varm uns
    touch obsm/${key_added}.pkl
    touch varm/${varm_key}.pkl
    touch uns/${uns_key}.pkl
    touch versions.yml
    """
}
