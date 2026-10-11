process ADATA_MERGE {
    label 'process_low'

    container "community.wave.seqera.io/library/bin2cell_harmonypy_scanorama_imagecodecs_pruned:d457f4c6ebb7b325"

    input:
    path(h5ad, arity: '1..*')
    val join
    val label
    val preserve_spatial
    val layer

    output:
    path "${prefix}.h5ad", emit: adata
    path "versions.yml"  , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "merged"
    template 'merge.py'

    stub:
    prefix = task.ext.prefix ?: "merged"
    """
    touch ${prefix}.h5ad
    touch versions.yml
    """
}
