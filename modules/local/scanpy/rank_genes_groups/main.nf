process SCANPY_RANK_GENES_GROUPS {
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
    prefix = task.ext.prefix ?: "${meta.id}_deg"
    groupby = task.ext.groupby ?: "clusters"
    method = task.ext.method ?: "t-test"
    template 'rank_genes_groups.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_deg"
    """
    touch ${prefix}.h5ad
    touch versions.yml
    """
}