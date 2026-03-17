process SQUIDPY_NHOOD_ENRICHMENT {
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
    prefix = task.ext.prefix ?: "${meta.id}_nhood_enrichment"
    cluster_key = task.ext.cluster_key ?: "clusters"
    template 'nhood_enrichment.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_nhood_enrichment"
    """
    touch ${prefix}.h5ad
    touch versions.yml
    """
}
