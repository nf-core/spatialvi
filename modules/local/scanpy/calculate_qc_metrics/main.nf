process SCANPY_CALCULATE_QC_METRICS {
    tag "${meta.id}"
    label 'process_single'

    container "community.wave.seqera.io/library/harmonypy_scanorama_gcc_gxx_pruned:95f731fde0b9ddef"

    input:
    tuple val(meta), path(adata, stageAs: "input.h5ad", arity: '1')
    val write_adata

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: adata, optional: true
    tuple val(meta), path("obs/*.pkl")     , emit: obs
    tuple val(meta), path("var/*.pkl")     , emit: var
    path "versions.yml"                    , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'calculate_qc_metrics.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def touch_adata = write_adata.toString() == 'true' ? "touch ${prefix}.h5ad" : ''
    """
    ${touch_adata}
    mkdir -p obs var
    touch obs/qc_metrics.pkl
    touch var/qc_metrics.pkl
    touch versions.yml
    """
}
