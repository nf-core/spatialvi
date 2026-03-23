process SCANPY_LEIDEN {
    tag "${meta.id}"
    label 'process_medium'

    container "community.wave.seqera.io/library/harmonypy_scanorama_gcc_gxx_pruned:95f731fde0b9ddef"

    input:
    tuple val(meta), path(adata, stageAs: "input.h5ad")
    val resolution
    val key_added

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: adata
    path "versions.yml"                    , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'leiden.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad
    touch versions.yml
    """
}
