process SCANPY_HARMONY {
    tag "${meta.id}"
    label 'process_medium'

    container "community.wave.seqera.io/library/harmonypy_scanorama_gcc_gxx_pruned:95f731fde0b9ddef"

    input:
    tuple val(meta), path(h5ad, stageAs: "input.h5ad", arity: '1')
    val key
    val embedding_added
    val write_adata

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: adata, optional: true
    tuple val(meta), path("obsm/*.pkl")    , emit: obsm
    path "versions.yml"                    , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'harmony.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def touch_adata = write_adata.toString() == 'true' ? "touch ${prefix}.h5ad" : ''
    """
    ${touch_adata}
    mkdir -p obsm
    touch obsm/${embedding_added}.pkl
    touch versions.yml
    """
}
