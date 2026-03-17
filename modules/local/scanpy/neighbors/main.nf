process SCANPY_NEIGHBORS {
    tag "${meta.id}"
    label 'process_medium'

    container "community.wave.seqera.io/library/harmonypy_scanorama_gcc_gxx_pruned:95f731fde0b9ddef"

    input:
    tuple val(meta), path(adata)

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: adata
    path "versions.yml",                     emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}_neighbors"
    n_neighbors = task.ext.n_neighbors ?: 15
    n_pcs = task.ext.n_pcs ?: "null"
    use_rep = task.ext.use_rep ?: "null"
    template 'neighbors.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_neighbors"
    """
    touch ${prefix}.h5ad
    touch versions.yml
    """
}
