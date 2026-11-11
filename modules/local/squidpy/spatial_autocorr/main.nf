process SQUIDPY_SPATIAL_AUTOCORR {
    tag "${meta.id}"
    label 'process_high'

    container "community.wave.seqera.io/library/harmonypy_scanorama_gcc_gxx_pruned:95f731fde0b9ddef"

    input:
    tuple val(meta), path(adata, stageAs: "input.h5ad", arity: '1')
    val mode
    val write_adata

    output:
    tuple val(meta), path("${prefix}.h5ad")   , emit: adata, optional: true
    tuple val(meta), path("${prefix}_svg.csv"), emit: csv
    tuple val(meta), path("uns/*.pkl")        , emit: uns
    path "versions.yml"                       , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'spatial_autocorr.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def touch_adata = write_adata.toString() == 'true' ? "touch ${prefix}.h5ad" : ''
    def uns_key = mode == 'moran' ? 'moranI' : 'gearyC'
    """
    ${touch_adata}
    touch ${prefix}_svg.csv
    mkdir -p uns
    touch uns/${uns_key}.pkl
    touch versions.yml
    """
}
