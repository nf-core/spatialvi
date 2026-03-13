process SCANPY_PCA {
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
    prefix = task.ext.prefix ?: "${meta.id}_pca"
    n_comps = task.ext.n_comps ?: 50
    use_highly_variable = task.ext.use_highly_variable ?: true
    template 'pca.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_pca"
    """
    touch ${prefix}.h5ad
    touch versions.yml
    """
}