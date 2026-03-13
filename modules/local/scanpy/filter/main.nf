process SCANPY_FILTER {
    tag "${meta.id}"
    label 'process_single'

    container "community.wave.seqera.io/library/harmonypy_scanorama_gcc_gxx_pruned:95f731fde0b9ddef"

    input:
    tuple val(meta), path(adata)
    val min_counts
    val min_genes
    val min_spots
    val mito_threshold
    val ribo_threshold
    val hb_threshold

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: adata
    tuple val(meta), path("${prefix}_stats.json"), emit: stats
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}_filtered"
    template 'filter.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_filtered"
    """
    touch ${prefix}.h5ad
    echo '{}' > ${prefix}_stats.json
    touch versions.yml
    """
}