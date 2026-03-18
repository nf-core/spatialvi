process SCANPY_NEIGHBORS {
    tag "${meta.id}"
    label 'process_medium'

    container "community.wave.seqera.io/library/harmonypy_scanorama_gcc_gxx_pruned:95f731fde0b9ddef"

    input:
    tuple val(meta), path(adata)
    val n_neighbors
    val neighbours_n_pcs
    val neighbours_use_rep

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: adata
    path "versions.yml",                     emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}_neighbors"
    template 'neighbors.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_neighbors"
    """
    touch ${prefix}.h5ad
    touch versions.yml
    """
}
