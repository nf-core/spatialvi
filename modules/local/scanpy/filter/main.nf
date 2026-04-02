process SCANPY_FILTER {
    tag "${meta.id}"
    label 'process_single'

    container "community.wave.seqera.io/library/harmonypy_scanorama_gcc_gxx_pruned:95f731fde0b9ddef"

    input:
    tuple val(meta), path(adata, stageAs: "input.h5ad")
    val min_counts
    val min_genes
    val min_obs
    val mito_threshold
    val ribo_threshold
    val hb_threshold

    output:
    tuple val(meta), path("${prefix}.h5ad")      , emit: adata
    tuple val(meta), path("${prefix}_stats.json"), emit: stats
    path "versions.yml"                          , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'filter.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad
    echo '{}' > ${prefix}_stats.json
    touch versions.yml
    """
}
