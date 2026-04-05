process SCANPY_SCANORAMA {
    tag "${prefix}"
    label 'process_medium'

    container "community.wave.seqera.io/library/harmonypy_scanorama_gcc_gxx_pruned:95f731fde0b9ddef"

    input:
    path h5ad
    val key
    val embedding_added

    output:
    path "${prefix}.h5ad", emit: adata
    path "versions.yml"  , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "merged"
    template 'scanorama.py'

    stub:
    prefix = task.ext.prefix ?: "merged"
    """
    touch ${prefix}.h5ad
    touch versions.yml
    """
}
