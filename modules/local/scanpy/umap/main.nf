process SCANPY_UMAP {
    tag "${meta.id}"
    label 'process_medium'

    container "community.wave.seqera.io/library/harmonypy_scanorama_gcc_gxx_pruned:95f731fde0b9ddef"

    input:
    tuple val(meta), path(adata)

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: adata
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}_umap"
    min_dist = task.ext.min_dist ?: 0.5
    spread = task.ext.spread ?: 1.0
    n_components = task.ext.n_components ?: 2
    template 'umap.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_umap"
    """
    touch ${prefix}.h5ad
    touch versions.yml
    """
}