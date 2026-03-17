process SQUIDPY_SPATIAL_AUTOCORR {
    tag "${meta.id}"
    label 'process_high'

    container "community.wave.seqera.io/library/harmonypy_scanorama_gcc_gxx_pruned:95f731fde0b9ddef"

    input:
    tuple val(meta), path(adata)
    val mode

    output:
    tuple val(meta), path("${prefix}.h5ad"),    emit: adata
    tuple val(meta), path("${prefix}_svg.csv"), emit: csv
    path "versions.yml",                        emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}_spatial_autocorr"
    n_top_genes = task.ext.n_top_genes ?: 50
    template 'spatial_autocorr.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_spatial_autocorr"
    """
    touch ${prefix}.h5ad
    touch ${prefix}_svg.csv
    touch versions.yml
    """
}
